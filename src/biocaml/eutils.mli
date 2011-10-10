(** Entrez Utilities API

    This modules provides a partial access to Entrez databases
    such as Pubmed, Gene or Protein. The API proposed by the 
    NCBI is based on HTTP requests, and this modules contains
    a couple of functions to ease the construction of appropriate
    URLs. It is thus of rather low-level; in particular, there is
    no support for parsing the answers, which are simple strings.

    Databases in Entrez can be seen as collections of records, each
    record representing an object of the database. The basic usage
    of the API is first to search a database with the esearch utility.
    Given a query string, esearch will return a collection of 
    identifiers. These identifiers are then used to fetch the actual
    records with the efetch utility.
*)

type database = [
| `gene
| `genome
| `geodatasets
| `geoprofiles
| `protein
| `pubmed
| `pubmedcentral
| `sra
| `unigene
| `taxonomy
]
(** Represents available databases *)


(** 4 Low level access 
    
    For a documentation of the parameters, see http://www.ncbi.nlm.nih.gov/books/NBK25499/
*)



val esearch_url : 
  ?retstart:int -> ?retmax:int -> 
  ?rettype:[`uilist | `count] -> 
  ?field:string ->
  ?datetype:[`pdat | `mdat | `edat] ->
  ?reldate:int ->
  ?mindate:string -> ?maxdate:string ->
  database -> string -> string
(** Construction of esearch URLs. *)

val efetch_url : 
  ?rettype:string -> ?retmode:string ->
  database -> string list -> string
(** Construction of efetch URLs. Note that this access method 
    does not support more than 200 ids *)
