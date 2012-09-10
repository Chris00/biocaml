(** BED data. A BED file is in the format shown below, where columns
    must be separted by a tab character.

    {v
    chrA   lo1   hi1
    chrA   lo2   hi2
    .      .     .
    .      .     .
    .      .     .
    chrB   lo1   hi1
    chrB   lo2   hi2
    .      .     .
    .      .     .
    .      .     .
    v}

    The definition is that intervals are zero based and half-open. So by
    default the line "chrA lo hi" is parsed to the interval [\[lo + 1,
    hi\]] on chromosome [chrA]. Similarly, when printing, the default
    is to print [\[lo - 1, hi\]]. The optional argument
    [increment_lo_hi] allows changing this behavior for non-conformant
    files. In addition, the optional argument [chr_map] is a [string
    -> string] function that allows changing of the chromosome name to
    a specified format, and defaults to [identity].

    Some tools require that the set of intervals do not overlap within
    each chromosome. This is not enforced, but you can use
    [any_overlap] to verify this property when needed.
*)


type t = string * int * int * [`Float of float| `Int of int | `String of string] list
with sexp
(** The type of BED data stream items. *)
  
type parse_error =
[ `not_a_float of Biocaml_pos.t * string
| `not_an_int of Biocaml_pos.t * string
| `wrong_number_of_columns of Biocaml_pos.t * string list
| `incomplete_input of Biocaml_pos.t * string list * string option
]
with sexp
(** The possible parsing errors. *)

type parsing_spec = [
| `enforce of [ `float | `int | `string ] list
| `strings
| `best_effort
]
with sexp
(** The specification of how to parse the remaining columns. *)

module Transform: sig
  val string_to_t:
    ?filename:string ->
    ?more_columns:parsing_spec ->
    unit ->
    (string, (t, parse_error) Core.Result.t) Biocaml_transform.t
(** Create a [Biocaml_transform.t] parser, while providing the format of the
    additional columns (default [`best_effort]). *)

  val t_to_string: unit ->
    (t, string) Biocaml_transform.t
end 
