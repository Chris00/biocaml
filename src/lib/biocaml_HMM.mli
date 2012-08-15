
(** Hidden Markov models.

  A Hidden Markov Model is given by two sequences of random variables,
  the states (X_n) and the observables (Y_n), n = 1, 2,...

  References:
  - {{:http://courses.cs.tamu.edu/rgutier/cpsc689_s07/welch2003baumWelch.pdf}
    Hidden Markov Models and the Baum–Welch Algorithm}, IEEE Information
    Theory Society Newsletter, Dec. 2003.
  - {{:http://courses.media.mit.edu/2010fall/mas622j/ProblemSets/ps4/tutorial.pdf}
    Some Mathematics for HMM}, Dawei Shen, 2008.
 *)

open Bigarray

type proba = float

type vec = (float, float64_elt, fortran_layout) Array1.t
type mat = (float, float64_elt, fortran_layout) Array2.t
type int_vec = (int, int_elt, fortran_layout) Array1.t

type viterbi_work
(** A workspace for {!T.viterbi}.  Can be used with HMMs instanced
    with various data if desired. *)

type baum_welch_work
(** A workspace for {!T.baum_welch}.  Can be used with HMMs instanced
    with various data if desired. *)

module type T =
  sig
    type t
    (** A hidden Markov model (mutable).  Such a model is given by three
        matrices [a], [b] and [init] (see below). *)
    val a : t -> mat
    (** [a] is the transition matrix, i.e.
        {[
          a.{x, y} = P(X_n = x | X_(n-1) = y)    (* it is independent of n *).
        ]}
        where [x] and [y] run through all states (converted to integers). *)
    val b : t -> mat
    (** [b] contains the emission probabilities:
        {[
          b.{x, k} = P(Y_n = o_k | X_n = x)      (* it is independent of n *)
        ]}
        where [x] runs though all states and [k] through all possible
        observations (both converted to integers). *)
    val init : t -> vec
    (** [init] the initial probabilities π_x(1): [init.{x}] = P(X₁ = [x]). *)

    val copy : t -> t
    (** [copy hmm] returns a copy of [hmm]. *)

    type state_seq
    (** A sequence of states. *)

    type obs_seq
    (** A sequence of observations. *)

    val length_obs : obs_seq -> int
    (** Length of the sequence of observations. *)


    (** {3 Standard HMM functions} *)

    val proba_obs : ?work:vec -> t -> obs_seq -> proba
    (** [proba_obs hmm obs] returns the probability of the sequence of
        observations [obs] for the model [hmm].
        @param work a vector with length greater or equal to twice the
        number of states. *)

    val viterbi_work : n_states: int -> length_obs: int -> viterbi_work
    (** Allocates enough memory for {!T.viterbi} to be run with a HMM of
        [n_states] states and sequences of at most [length_obs]
        observations. *)

    val viterbi : ?work:viterbi_work -> ?states: state_seq ->
                  t -> obs_seq -> state_seq * proba
    (** [viterbi hmm obs] returns [(s, p)] where [s] the most probable
        sequence of states given the observations [obs] for the hidden
        Markov model [hmm].  More precisely, it returns the sequence
        x₁x₂...xN, N = [length_obs obs], that maximizes the probability:
        {[
          P(X₁ = x₁,..., X_N = x_N | Y₁ = o₁,..., Y_N = o_N)
        ]}
        where [obs] = o₁...o_N.  The probability [p] is
        {[
          p = P(X₁ = x₁,..., X_N = x_N and Y₁ = o₁,..., Y_N = o_N)
        ]}
        Beware that it is the {i joint} distribution, not the
        conditional one.

        @param states preallocated sequence.  Will be populated with
        the most probable sequence and returned.
        @param viterbi_work the workspace needed for [viterbi] to
        work.  Can be allocated with {!viterbi_work}.  *)

    val baum_welch_work : n_states: int -> length_obs: int -> baum_welch_work
    (** Allocates enough memory for {!baum_welch} to be run with a HMM
        of [n_states] states and sequences of at most [length_obs]
        observations. *)

    val baum_welch : ?work:baum_welch_work -> ?max_iter: int ->
                     t -> obs_seq -> unit
    (** [baum_welch hmm seq] reestimate the [hmm] parameters
        (modifying [hmm] in place) to maximize the probability that
        the sequence of observations [seq] occurs.

        @param max_iter the maximum number of iterations of the
        Baum-Welch algorithm.  Default: [20].  *)

    val forward : ?alpha:mat -> t -> obs_seq -> mat
    (** [forward hmm obs] efficiently computes the matrix [α.{x,n}],
        where [x] varies among the states (1,..., dim1 a) and [n =
        1,..., length_obs obs], defined by
        {[
          α.{x, n} = P(Y₁ = o₁,..., Y_n = o_n, X_n = x)
        ]}
        where [obs] = o₁...o_N and N = [length_obs obs].  In
        particular, note that (using [Lacaml]):
        P(Y₁ = o₁,..., Y_n = o_n) = [Vec.sum(Mat.col alpha n)]. *)

    val backward : ?beta:mat -> t -> obs_seq -> mat
    (** [backward a b obs] efficiently computes the matrix [β.{x,n}]
        where [x] varies among the states (1,..., dim1 a) and [n =
        1,..., N] defined by
        {[
          β.{x,n} = P(Y_{n+1} = o_{n+1},..., Y_N = o_N | X_n = x)
        ]}
        where N = [length_obs obs] and [obs] = o₁,...,o_N.  *)
    ;;
  end

(** {2 HMM with integer states and observations} *)

include T with type state_seq = int_vec and type obs_seq = int_vec

val make : n_states: int -> n_obs: int -> t
(** Create an untrained HMM with [n_states] (numbered [1, 2,...,
    n_states]) that can handle [n_obs] (numbered [1,..., n_obs]).  *)

val of_mat : ?check: bool -> a: mat -> b: mat -> vec -> t
(** [of_mat a b init] returns a HMM with transition probabilities [a],
    emission probabilities [b] and initial probabilities [init].  See
    {!T.a}, {!T.b} and {!T.init} for the content of these matrices.
    @raise Invalid_argument if the dimensions of [a], [b] and [init]
    are not coherent.

    @param check verifies that the matrices satisfy the probability
    constraints.  Default: [true]. *)


(** {2 HMM with arbitrary states and observations} *)

(** Specification of sequence of states. *)
module type STATE_SEQ =
  sig
    type t  (** A sequence of states *)

    val create : int -> t
    (** [create n] create a sequence of states of length [n].  The
        sequence does not need to be initialized. *)

    val length : t -> int
    (** [length seq] returns the number of states in the sequence [seq]. *)

    val set : t -> int -> state:int -> unit
    (** [set seq i state] set the element [i] of the sequence [seq] to
        the state numbered [state].  States are numbered [1,...,N]. *)
  end

(** Specification of sequence of observations. *)
module type OBS_SEQ =
  sig
    type t  (** A sequence of observations *)

    val n_obs : int
    (** The number of different values observations can take. *)

    val length : t -> int
    (** [length seq] returns the number of observations in the
        sequence [seq]. *)

    val get : t -> int -> int
    (** [get seq i] returns the [i]th ([i=1,..., length seq])
        observation in the sequence [seq].  This function must take
        care of converting the observation to an integer in the range
        [1,..., n_obs] (otherwise the HMM functions will raise
        exceptions, for example [Invalid_argument]).  *)
  end

(** Create a HMM module for the sequences of states specified by [S]
    and the sequences of observations specified by [O]. *)
module Make(S: STATE_SEQ)(O: OBS_SEQ) : sig
  include T with type state_seq = S.t and type obs_seq = O.t

  val make : n_states: int -> t

  val of_mat : ?check: bool -> a: mat -> b: mat -> vec -> t
  (** Same as {!Biocaml_HMM.of_mat} with the additional check that [b]
      has the right size fo the number of possible observations
      {!OBS_SEQ.n_obs}. *)
end


module Int_state : STATE_SEQ with type t = int_vec

module type STRING = sig val chars : string end

module Char_state : functor (C:STRING) -> STATE_SEQ with type t = string
(** [char_obs chars] returns a module where states are strings
    (seen as sequences of char elements) taking values in [chars]. *)

module Char_obs : functor (C: STRING) -> OBS_SEQ with type t = string
(** [char_obs chars] returns a module where observations are strings
    (seen as sequences of char elements) taking values in [chars]. *)

