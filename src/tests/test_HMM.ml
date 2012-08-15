open Bigarray
open Printf
open OUnit

let mat a = Array2.of_array float64 fortran_layout a
let vec v = Array1.of_array float64 fortran_layout v
let int_vec v = Array1.of_array int fortran_layout v

let assert_almostequal ?(epsilon=1e-13) name x y =
  assert_bool name (abs_float(x -. y) < epsilon)

let assert_equalmat ?(epsilon=1e-13) a_name a a_exact =
  for x = 1 to Array2.dim1 a do
    for k = 1 to Array2.dim2 a do
      assert_bool (sprintf "%s.{%i,%i}" a_name x k)
                  (abs_float(a.{x,k} -. a_exact.{x,k}) < epsilon)
    done
  done


module Obs = struct let chars = "RCS" end (* Rainy, Cloudy, Sunny *)
module HMM = Biocaml_HMM.Make (Biocaml_HMM.Int_state)
                              (Biocaml_HMM.Char_obs(Obs))

let test_Pobs () =
  let a = mat [| [| 0.4; 0.2; 0.1 |];
                 [| 0.3; 0.6; 0.1 |];
                 [| 0.3; 0.2; 0.8 |] |]
  (* State and observations are in 1-1 correspondance (thus not really a
     hidden model). *)
  and b = mat [| [| 1.; 0.; 0. |];
                 [| 0.; 1.; 0. |];
                 [| 0.; 0.; 1. |] |]
  and init = vec [| 0.; 0.; 1. (* S *) |] in
  let hmm = HMM.of_mat ~a ~b init in
  assert_almostequal "S" (HMM.proba_obs hmm "S") 1.;
  assert_almostequal "SSSRRSCS" (HMM.proba_obs hmm "SSSRRSCS") 1.536e-4

let test_forward_backward () =
  (* http://en.wikipedia.org/wiki/Forward%E2%80%93backward_algorithm#Example *)
  let a = mat [| [| 0.7; 0.3 |]; [| 0.3; 0.7 |] |]
  and b = mat [| [| 0.9; 0.1 |]; [| 0.2; 0.8 |] |]
  and init = vec [| 0.5; 0.5 |] in
  let hmm = Biocaml_HMM.of_mat a b init in
  let alpha = Biocaml_HMM.forward hmm (int_vec [|1; 1; 2; 1; 1 |]) in
  assert_equalmat "alpha" alpha ~epsilon:1e-4
                  (mat [|[| 0.45; 0.310475; 0.0229; 0.0407; 0.0297 |];
                         [| 0.1;  0.040975; 0.0974; 0.0150; 0.0045 |] |]);
  let beta = Biocaml_HMM.backward hmm (int_vec [|1; 1; 2; 1; 1 |]) in
  assert_equalmat "beta" beta ~epsilon:1e-4
                  (mat [| [| 0.0661; 0.09060; 0.45925; 0.69; 1. |];
                          [| 0.0455; 0.15028; 0.24365; 0.41; 1. |] |])

let test_viterbi () =
  let module Obs = struct let chars = "HT" end in
  let module HMM = Biocaml_HMM.Make (Biocaml_HMM.Int_state)
                                    (Biocaml_HMM.Char_obs(Obs)) in
  let hmm = HMM.make ~n_states:3 in
  let obs1 = "HHHHTHTTTT" and obs2 = "HTTHTHHTTH" in
  let st1, p1 = HMM.viterbi hmm obs1 in
  printf ""


let test_baum_welch () =
  ()

let tests = "HMM" >::: [
            "Proba(obs)" >:: test_Pobs;
            "Forward/Backward" >:: test_forward_backward;
            "Viterbi" >:: test_viterbi;
            "Baum-Welch" >:: test_baum_welch;
          ]
