(* This file is part of Biocaml.
   @author Christophe Troestler *)

open Bigarray
open Printf

type proba = float

type vec = (float, float64_elt, fortran_layout) Array1.t
type mat = (float, float64_elt, fortran_layout) Array2.t
type int_vec = (int, int_elt, fortran_layout) Array1.t

let vec_create n = Array1.create float64 fortran_layout n
let mat_create n m = Array2.create float64 fortran_layout n m

let vec_copy v =
  let v' = vec_create (Array1.dim v) in
  Array1.blit v v';
  v'

let mat_copy a =
  let a' = mat_create (Array2.dim1 a) (Array2.dim2 a) in
  Array2.blit a a';
  a'

type viterbi_work = { delta: mat;
                      phi: (int, int_elt, fortran_layout) Array2.t }

type baum_welch_work = { alpha: mat;  (* forward quantities *)
                         beta: mat;   (* backward quantities *)
                       }

module type T =
  sig
    type t
    val a : t -> mat
    val b : t -> mat
    val init : t -> vec
    val copy : t -> t

    type state_seq
    type obs_seq
    val length_obs : obs_seq -> int

    val proba_obs : ?work:vec -> t -> obs_seq -> proba
    val viterbi_work : n_states: int -> n_obs: int -> viterbi_work
    val viterbi : ?work:viterbi_work -> ?states: state_seq ->
                  t -> obs_seq -> state_seq * proba
    val baum_welch_work : n_states: int -> n_obs: int -> baum_welch_work
    val baum_welch : ?work:baum_welch_work -> ?max_iter: int ->
                     t -> obs_seq -> unit
    val forward : ?alpha:mat -> t -> obs_seq -> mat
    val backward : ?beta:mat -> t -> obs_seq -> mat
  end

type t = { a: mat;          (* transition proba, square matrix *)
           b: mat;
           init: vec }

let a hmm = hmm.a
let b hmm = hmm.b
let init hmm = hmm.init
let copy hmm = { a = mat_copy hmm.a;
                 b = mat_copy hmm.b;
                 init = vec_copy hmm.init }

type state_seq = int_vec
let length_state = Array1.dim
type obs_seq = int_vec
let length_obs = Array1.dim;;

let make ~n_states ~n_obs =
  let a = mat_create n_states n_states
  and b = mat_create n_states n_obs
  and init = vec_create n_states in
  let equi_state = 1. /. float n_states
  and equi_obs = 1. /. float n_obs in
  Array2.fill a equi_state;
  Array2.fill b equi_obs;
  Array1.fill init equi_state;
  { a; b; init }

(* Forward
 ***********************************************************************)

DEFINE FORWARD(n_obs, get, obs) =
  let alpha = match alpha with
    | None -> mat_create (Array2.dim2 hmm.a) n_obs
    | Some m ->
       if Array2.dim1 m < Array2.dim2 hmm.a then
         invalid_arg "Biocaml_HMM.forward: dim1(alpha) too small";
       if Array2.dim2 m < n_obs then
         invalid_arg "Biocaml_HMM.forward: dim2(alpha) too small";
       m in
  (* α.{x,1} = b.{x, o₁} π_x(1)  *)
  let obs_1 = get obs 1 in
  for x = 1 to Array2.dim2 hmm.a do
    alpha.{x,1} <- hmm.init.{x} *. hmm.b.{x, obs_1}
  done;
  (* α.{x, n} = ∑ α.{y, n-1} a.{x, y} b.{x, obs.{n}} on all y ∈ E
     where a.{x, y} = Prob(X_n = x | X_(n-1) = y) *)
  for n = 2 to n_obs do
    let n_1 = n - 1 in
    let obs_n = get obs n in
    for x = 1 to Array2.dim1 hmm.a do
      let sum = ref 0. in
      for y = 1 to Array2.dim2 hmm.a do
        sum := !sum +. alpha.{y, n_1} *. hmm.a.{x,y}
      done;
      alpha.{x, n} <- !sum *. hmm.b.{x, obs_n}
    done
  done;
  alpha;;

let forward ?alpha hmm (obs: int_vec) =
  let n_obs = Array1.dim obs in
  FORWARD(n_obs, Array1.get, obs)
;;

(* Proba of obs
 ***********************************************************************)

(* Same as forward but we only need the last column so we do not
   allocate the full matrix [alpha]. *)
DEFINE PROBA_OBS(n_obs, get, obs) =
  let n_states = Array2.dim2 hmm.a in
  let alpha = match work with
    | None -> vec_create (2 * n_states)
    | Some v ->
       if Array1.dim v < 2 * n_states then
         invalid_arg "Biocaml_HMM.proba_obs: dim(work) too small";
       v in
  for x = 1 to Array2.dim2 hmm.a do
    alpha.{x} <- hmm.b.{x, get obs 1} *. hmm.init.{x}
  done;
  let curr = ref 0 and next = ref n_states in (* offset to avoid copying *)
  for n = 2 to n_obs do
    let obs_n = get obs n in
    for x = 1 to n_states do
      let sum = ref 0. in
      for y = 1 to n_states do
        sum := !sum +. alpha.{y + !curr} *. hmm.a.{x,y}
      done;
      alpha.{x + !next} <- !sum *. hmm.b.{x, obs_n}
    done;
    let will_be_next = !curr in
    curr := !next;
    next := will_be_next
  done;
  (* ∑ alpha *)
  let sum = ref 0. in
  for x = 1 + !curr to Array1.dim alpha + !curr do
    sum := !sum +. alpha.{x}
  done;
  !sum *. 1. (* force unboxing *) ;;

let proba_obs ?work hmm (obs:int_vec) =
  let n_obs = Array1.dim obs in
  PROBA_OBS(n_obs, Array1.get, obs)
;;

(* Backward
 ***********************************************************************)

(* β.{x, n} = ∑_{y ∈ E}  b.{y, obs.{n+1}} β.{y, n+1} a.{y, x} *)
DEFINE BACKWARD(n_obs, get, obs) =
  let n_states = Array2.dim2 hmm.a in
  let beta = match beta with
    | None -> mat_create n_states n_obs
    | Some m ->
      if Array2.dim1 m < n_states then
        invalid_arg "Biocaml_HMM.backward: dim1(beta) too small";
      if Array2.dim2 m < n_obs then
        invalid_arg "Biocaml_HMM.backward: dim2(beta) too small";
      m in
  Array2.fill beta 0.;
  for x = 1 to n_states do beta.{x, n_obs} <- 1. done;
  for n = n_obs - 1 downto 1 do
    let o = get obs (n+1) in
    for x = 1 to n_states do
      for y = 1 to n_states do
        beta.{x, n} <- beta.{x, n} +. hmm.b.{y, o} *. beta.{y, n+1} *. hmm.a.{y, x}
      done;
    done;
  done;
  beta;;

let backward ?beta hmm (obs:int_vec) =
  let n_obs = Array1.dim obs in
  BACKWARD(n_obs, Array1.get, obs)
;;

(* Viterbi
 ***********************************************************************)

let viterbi_work ~n_states ~n_obs =
  { delta = mat_create n_states n_obs;
    phi = Array2.create int fortran_layout n_states n_obs }
;;

DEFINE VITERBI(n_obs, get, obs, create_state, length_state, set_state) =
  let n_states = Array2.dim1 hmm.a in
  let seq_states = match states with
    | None -> create_state n_obs
    | Some w ->
      if length_state w < n_obs then
        invalid_arg(sprintf "Biocaml_HMM.viterbi: dim states = %i, need >= %i"
                            (length_state w) n_obs);
      w in
  let { delta; phi } = match work with
    | None -> viterbi_work n_states n_obs
    | Some w ->
      if Array2.dim1 w.delta < n_states then
        invalid_arg(sprintf "Biocaml_HMM.viterbi: n_state(work) = %i, \
                             need >= %i" (Array2.dim1 w.delta) n_states);
      if Array2.dim2 w.delta < n_obs then
        invalid_arg(sprintf "Biocaml_HMM.viterbi: n_obs(work) = %i, need >= %i"
                            (Array2.dim2 w.delta) n_obs);
      w in
  (* δ.{x, 1} = P(Y₁ = o₁, X₁ = x) = b.{x, o₁} π_x(1) *)
  for x = 1 to n_states do
    delta.{x,1} <- hmm.b.{x,1} *. hmm.init.{x}
  done;
  (* FIXME: only 2 consecutive cols of δ are needed *)
  for n = 1 to n_obs - 1 do
    let o = get obs (n+1) in
    for x = 1 to n_states do
      (* δ.{x, n+1} = max_{y ∈ E}  δ.{y, n} a.{x, y} b.{x, o_{n+1}} *)
      let m = ref neg_infinity
      and arg_m = ref 0 in
      for y = 1 to n_states do
        let m' = delta.{y, n} *. hmm.a.{x,y} in
        if m' > !m then (m := m';  arg_m := y)
      done;
      delta.{x, n+1} <- !m *. hmm.b.{x, o};
      phi.{x, n+1} <- !arg_m
    done
  done;
  (* Determine the most probable sequence by backtracking. *)
  let state_max = ref 1
  and max = ref delta.{1, n_obs} in
  for x = 2 to n_states do
    let d = delta.{x, n_obs} in
    if d > !max then (max := d; state_max := x)
  done;
  set_state (seq_states: state_seq) n_obs !state_max;
  let proba = !max *. 1. in
  for n = n_obs - 1 downto 1 do
    state_max := phi.{!state_max, n+1}; (* s(n) ← φ(s(n+1), n+1) *)
    set_state (seq_states: state_seq) n !state_max;
  done;
  seq_states, proba


let viterbi ?work ?states hmm (obs: int_vec) =
  let n_obs = Array1.dim obs in
  (* To trigger the optimization, [Array1.set] must be used
     explicitely and its first argument must be type annotated (done
     in the macro). *)
  VITERBI(n_obs, Array1.get, obs,
          Array1.create int fortran_layout, Array1.dim, Array1.set)
;;

(* Baum Welch
 ***********************************************************************)

let baum_welch_work ~n_states ~n_obs =
  if n_states <= 0 then
    invalid_arg "Biocaml_HMM.baum_welch_work: number of states <= 0";
  if n_obs <= 0 then
    invalid_arg "Biocaml_HMM.baum_welch_work: number of observations <= 0";
  { alpha = mat_create n_states n_obs;
    beta = mat_create n_states n_obs;
  }

let sum_col (m: mat) col =
  let sum = ref 0. in
  for x = 1 to Array2.dim1 m do sum := !sum +. m.{x, col} done;
  !sum *. 1. (* force unboxing *)

(* Implementation of Baum-Welch algorithm based on:

   Hidden Markov Models and the Baum–Welch Algorithm, IEEE Information
   Theory Society Newsletter, Dec. 2003. [REDO the computs!]
   http://www-rcf.usc.edu/~lototsky/MATH508/Baum-Welch.pdf
   http://courses.cs.tamu.edu/rgutier/cpsc689_s07/welch2003baumWelch.pdf
 *)
DEFINE BAUM_WELCH(n_obs, get, obs) =
  let n_states = Array2.dim1 hmm.a in
  let { alpha; beta } = match work with
    | None -> baum_welch_work n_states n_obs
    | Some w ->
       (* Worspace can only be created by special fun. Dims correlated. *)
       if Array2.dim1 w.alpha < n_states then
         invalid_arg "Biocaml_HMM.baum_welch: number of states of \
                      \"work\" too small";
       if Array2.dim2 w.alpha < n_obs then
         invalid_arg "Biocaml_HMM.baum_welch: number of observations of \
                      \"work\" too small";
       w in
  Printf.printf "alpha... %!";
  ignore(forward ~alpha hmm obs);       (* modify [alpha] *)
  Printf.printf "beta... %!";
  ignore(backward ~beta hmm obs);       (* modify [beta] *)
  let a = hmm.a and b = hmm.b and init = hmm.init in
  let better = ref true in
  let n_iter = ref 0 in
  let p_obs = ref(sum_col alpha n_obs) in
  while !better && !n_iter < max_iter && !p_obs > 0. do
    (* Baum-Welch re-estimation (computations done for our case), N=n_obs:

       init.{x} ← α₁(x) β₁(x) / P(obs)

                         ∑_{n=1,...,N-1} α.{x,n} b.{y, o_n+1} β.{y, n+1}
       a.{y,x} ← a.{y,x} ------------------------------------------------
                         ∑_{n=1,...,N-1} α.{x,n} β.{x,n}

                 ∑_{n=1,...,N | o_n = k} α.{x,n} β.{x,n}
       b.{x,k} ← ---------------------------------------
                 ∑_{n=1,...,N}           α.{x,n} β.{x,n}
     *)
    Printf.printf "re-estimation %i (P(obs)=%g)\n%!" !n_iter !p_obs;
    Array2.fill b 0.;
    for x = 1 to n_states do
      init.{x} <- alpha.{x,1} *. beta.{x,1} /. !p_obs;
      let den_a = ref 0. in (* ∑_{n=1,...,N-1} α.{x,n} β.{x,n} *)
      for n = 1 to n_obs - 1 do
        let ab = alpha.{x,n} *. beta.{x,n} in
        den_a := !den_a +. ab;
        let k = get obs n in
        b.{x, k} <- b.{x, k} +. ab;
      done;
      let ab = alpha.{x, n_obs} *. beta.{x, n_obs} in
      let den_b = !den_a +. ab in (* ∑_{n=1,...,N} α.{x,n} β.{x,n} *)
      let k = get obs n_obs in
      b.{x, k} <- b.{x, k} +. ab;
      (* Re-estimate [b.{x,_}]. *)
      if den_b > 0. then (
        for k = 1 to Array2.dim2 b do b.{x,k} <- b.{x,k} /. den_b done
      )
      else (
        (* ∀n P(X_n = x, obs) = 0, i.e. the observations say nothing
           about the emissions from [x].  Equi-distribution is then "best". *)
        let equi = 1. /. float(Array2.dim2 b) in
        for k = 1 to Array2.dim2 b do b.{x,k} <- equi done
      );
      (* Re-estimate [a.{_,x}]. *)
      if !den_a > 0. then (
        for y = 1 to n_states do
          if a.{y,x} > 0. then (
            let sum = ref 0. in
            for n = 2 to n_obs do
              sum := !sum +. alpha.{x,n-1} *. b.{y, get obs n} *. beta.{y,n}
            done;
            a.{y,x} <- a.{y,x} *. !sum /. !den_a
          )
        done
      )
    done;

    Printf.printf "alpha %i... %!" !n_iter;
    ignore(forward ~alpha hmm obs); (* overwrite [alpha] *)
    Printf.printf "beta... %!";
    ignore(backward ~beta hmm obs); (* overwrite [beta] *)
    let p_obs_new = sum_col alpha n_obs in
    Printf.printf "P(obs_new) = %g\n%!" p_obs_new;
    better := p_obs_new > 1.001 *. !p_obs;
    p_obs := p_obs_new;
    incr n_iter;
  done
  ; Printf.printf "P(obs) = %g \n%!" !p_obs
;;

let baum_welch_max_iter = 20

let baum_welch ?work ?(max_iter=baum_welch_max_iter) hmm (obs:int_vec) =
  let n_obs = Array1.dim obs in
  BAUM_WELCH(n_obs, Array1.get, obs)

(* Functorial interface
 ***********************************************************************)

module type STATE_SEQ =
  sig
    type t
    val create : int -> t
    val length : t -> int
    val set : t -> int -> state:int -> unit
  end

module type OBS_SEQ =
  sig
    type t
    val n_obs : int
    val length : t -> int
    val get : t -> int -> int
  end

module Make(S: STATE_SEQ)(O: OBS_SEQ) = struct
  type hmm = t
  type t = hmm

  let a = a and b = a and init = init and copy = copy (* re-export *)

  type state_seq = S.t
  type obs_seq = O.t
  let length_obs = O.length

  let make ~n_states = make ~n_states ~n_obs:O.n_obs

  let forward ?alpha hmm obs =
    let n_obs = O.length obs in
    FORWARD(n_obs, O.get, obs)

  let proba_obs ?work hmm obs =
    let n_obs = O.length obs in
    PROBA_OBS(n_obs, O.get, obs)

  let backward ?beta hmm obs =
    let n_obs = O.length obs in
    BACKWARD(n_obs, O.get, obs)

  let viterbi_work = viterbi_work       (* re-export *)

  let viterbi ?work ?states hmm obs =
    let n_obs = O.length obs in
    VITERBI(n_obs, O.get, obs,  S.create, S.length, S.set)

  let baum_welch_work = baum_welch_work (* reexport *)

  let baum_welch ?work ?(max_iter=baum_welch_max_iter) hmm obs =
    let n_obs = O.length obs in
    BAUM_WELCH(n_obs, O.get, obs)
end

(* Some helpers to create SEQ *)

module Int_state = struct
  type t = int_vec
  let create n = Array1.create int fortran_layout n
  let length = Array1.dim
  let set (v: int_vec) i ~state = v.{i} <- state
end

module type STRING = sig val chars : string end

module Char_state(C: STRING) =
  struct
    type t = string

    let to_char = String.copy C.chars
    let create = String.create
    let length = String.length
    let set v i ~state = v.[i] <- to_char.[state]
  end

module Char_obs(C: STRING) =
  struct
    let n_obs = String.length C.chars
    let to_int = Array.make 256 0 (* waste a bit of memory but be fast *)
    let () =
      for k = 0 to n_obs - 1 do
        to_int.(Char.code C.chars.[k]) <- k + 1
      done;

    type t = string
    let n_obs = n_obs
    let length = String.length
    let get v i = to_int.(Char.code v.[i])
  end
