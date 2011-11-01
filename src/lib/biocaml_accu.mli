(** A datastructure (based on Hashtbl) to accumulate values *)
open Batteries

type ('instance,'bin,'increment,'accu) t

val create : ?n:int -> 'd -> ('a -> 'b) -> ('c -> 'd -> 'd) -> ('a,'b,'c,'d) t
(** [n] is the approximate size of the domain *)

val add : ('a,'b,'c,'d) t -> 'a -> 'c -> unit

val enum : ('a,'b,'c,'d) t -> ('b * 'd) Enum.t

val get  : ('a,'b,'c,'d) t -> 'b -> 'd




type 'instance counter = ('instance, 'instance, int, int) t

val counter : ?n:int -> unit -> 'a counter

val count : 'a counter -> 'a -> int -> unit

