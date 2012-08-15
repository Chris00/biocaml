open Bigarray
open Printf
open OUnit

let assert_almostequal ?epsilon x y =
  assert_equal ~cmp:(cmp_float ?epsilon) x y

module Obs = struct let chars = "RCS" end (* Rainy, Cloudy, Sunny *)
module HMM = Biocaml_HMM.Make (Biocaml_HMM.Int_state)
                              (Biocaml_HMM.Char_obs(Obs))

let a = Array2.of_array float64 fortran_layout [| [| 0.4; 0.2; 0.1 |];
                                                  [| 0.3; 0.6; 0.1 |];
                                                  [| 0.3; 0.2; 0.8 |] |]
(* State and observations are in 1-1 correspondance (thus not really a
   hidden model). *)
let b = Array2.of_array float64 fortran_layout [| [| 1.; 0.; 0. |];
                                                  [| 0.; 1.; 0. |];
                                                  [| 0.; 0.; 1. |] |]

let init = Array1.of_array float64 fortran_layout [| 0.; 0.; 1. (* S *) |]

let test_Pobs () =
  let hmm = HMM.of_mat ~a ~b init in
  assert_almostequal (HMM.proba_obs hmm "S") 1.;
  assert_almostequal (HMM.proba_obs hmm "SSSRRSCS") 1.536e-4

let tests = "HMM" >::: [
            "Proba(obs)" >:: test_Pobs;
          ]
