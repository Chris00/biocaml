
type raw_alignment = {
  qname : string;
  flag : int;
  ref_id: int;
  pos : int; (** 0-based, -1 if undefined*)
  mapq : int;
  bin: int;
  cigar : string;
  next_ref_id : int;
  pnext : int;
  tlen : int;
  seq : string;
  qual : int array;
  optional : string;
} with sexp

type raw_item =
[ `alignment of raw_alignment
| `header of string
| `reference_information of (string * int) array ]
with sexp
  
module Transform: sig
  type raw_bam_error = [
  | `read_name_not_null_terminated of string
  | `reference_information_name_not_null_terminated of string
  | `reference_information_overflow of int * string
  | `wrong_magic_number of string
  | `wrong_int32 of string
  ]
  with sexp

  val string_to_raw:
    ?zlib_buffer_size:int ->
    unit ->
    (string,
     (raw_item, [> `unzip of Biocaml_zip.Transform.unzip_error
                | `bam of raw_bam_error ] )
       Core.Result.t)
      Biocaml_transform.t

  type parse_optional_error = [
  | `wrong_auxiliary_data of
      [ `array_size of int
      | `null_terminated_hexarray
      | `null_terminated_string
      | `out_of_bounds
      | `wrong_int32 of string
      | `unknown_type of char ] * string
  ]
  with sexp

  val parse_optional: ?pos:int -> ?len:int -> string ->
    (Biocaml_sam.optional_content, parse_optional_error) Core.Result.t
    
  type parse_cigar_error = [
  | `wrong_cigar of string
  | `wrong_cigar_length of int ]
  with sexp

  val parse_cigar: ?pos:int -> ?len:int -> string ->
    (Biocaml_sam.cigar_op array, parse_cigar_error) Core.Result.t
      
  type raw_to_item_error = [
  | `header_line_not_first of int
  | `header_line_without_version of (string * string) list
  | `header_line_wrong_sorting of string
  | `invalid_header_tag of int * string
  | `invalid_tag_value_list of int * string list
  | `reference_sequence_not_found of raw_alignment
  | parse_optional_error
  | parse_cigar_error
  | `wrong_flag of raw_alignment
  | `wrong_mapq of raw_alignment
  | `wrong_pnext of raw_alignment
  | `wrong_pos of raw_alignment
  | `wrong_qname of raw_alignment
  | `wrong_tlen of raw_alignment ]
  with sexp
    
  val raw_to_item :
    unit ->
    (raw_item, (Biocaml_sam.item, [> raw_to_item_error]) Core.Result.t)
      Biocaml_transform.t

  type item_to_raw_error =
  [ `cannot_get_sequence of Biocaml_sam.alignment
  | `header_line_not_first of string
  | `reference_name_not_found of Biocaml_sam.alignment * string ]
  with sexp

  val item_to_raw :
    unit ->
    (Biocaml_sam.item,
     (raw_item, item_to_raw_error) Core.Result.t) Biocaml_transform.t
       
      
  val raw_to_string :
    ?zlib_buffer_size:int ->
    unit ->
    (raw_item, string) Biocaml_transform.t
end
    
