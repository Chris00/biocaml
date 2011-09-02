type 'a location = 'a * Range.t

module type Domain = sig
  type 'a t

  val of_locations : 'a location Enum.t -> 'a t

  val inter : 'a t -> 'a t -> 'a t
  val diff : 'a t -> 'a t -> 'a t
  val size : 'a t -> int

  val intersects : 'a location -> 'a t -> bool
  val intersection_size : 'a location -> 'a t -> int

  val enum : 'a t -> 'a location Enum.t
end


module type Signal = sig
  type ('a,'b) t
  val make : ('b list -> 'c) -> ('a location * 'c) Enum.t -> ('a,'c) t

  val eval : 'a -> int -> ('a,'b) t -> 'b

  val fold : ('a -> Range.t -> 'b -> 'c -> 'c) -> ('a,'b) t -> 'c -> 'c

  val enum : ('a,'b) t -> ('a location * 'b) Enum.t
end

module type LSet = sig
  type 'a t

  val make : 'a location Enum.t -> 'a t


  val fold : ('a -> Range.t -> 'b -> 'b) -> 'a t -> 'b -> 'b

  val intersects : 'a location -> 'a t -> bool

  val enum : 'a t -> 'a location Enum.t

  val union : 'a t -> 'a t -> 'a t
  val add : 'a location -> 'a t -> 'a t
end

module type LMap = sig
  type ('a,'b) t

  val make : ('a location * 'b) Enum.t -> ('a,'b) t

  val fold : ('a -> Range.t -> 'b -> 'c -> 'c) -> ('a,'b) t -> 'c -> 'c

  val pwfold : ('a -> Range.t -> 'b list -> 'c -> 'c) -> ('a, 'b) t -> 'c -> 'c

  val intersects : 'a location -> ('a,'b) t -> bool

  val enum : ('a,'b) t -> ('a location * 'b) Enum.t

  val union : ('a,'b) t -> ('a,'b) t -> ('a,'b) t
  val add : 'a location -> 'b -> ('a,'b) t -> ('a,'b) t
end

