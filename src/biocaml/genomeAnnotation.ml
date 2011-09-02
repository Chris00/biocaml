open Batteries

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

let rec iset_intersects_range i j s = ISet.(BatAvlTree.(
  if i > j then raise (Invalid_argument "iset_intersects_range") ;
  if is_empty s then false
  else
    let v1, v2 = root s in 
    if j < v1 then iset_intersects_range i j (left_branch s)
    else if v2 < i then iset_intersects_range i j (right_branch s)
    else true
))

module ISetDomain = struct
  type 'a t = ('a, ISet.t) PMap.t

  let of_locations e = 
    let accu = Accu.create ISet.empty fst (fun (_,r) -> Range.(ISet.add_range r.lo r.hi)) in
    Enum.iter (fun loc -> Accu.add accu loc loc ) e ;
    PMap.of_enum (Accu.enum accu)

  let inter u v =
    PMap.foldi
      (fun k set_u accu ->
	 try 
	   let set_v = PMap.find k v in
	   PMap.add k (ISet.inter set_u set_v) accu
	 with Not_found -> accu)
      u PMap.empty

  let diff u v =
    PMap.foldi
      (fun k set_u accu ->
	 let set_u' = 
	   try 
	     let set_v = PMap.find k v in
	     ISet.diff set_u set_v
	   with Not_found -> set_u
	 in PMap.add k set_u' accu)
      u PMap.empty

  let size x = 
    PMap.foldi (fun _ set accu -> ISet.cardinal set + accu) x 0

  let intersection_size (k,r) dom = ISet.(
    try 
      inter
	Range.(add_range r.lo r.hi empty)
	(PMap.find k dom)
      |> cardinal
    with Not_found -> 0
  )

  let intersects (k,r) dom = 
    try Range.(iset_intersects_range r.lo r.hi (PMap.find k dom))
    with Not_found -> false
    
  let enum dom = 
    (PMap.enum dom) /@ (fun (k,s) -> Enum.map (fun (lo,hi) -> k, Range.make lo hi) (ISet.enum s))
    |> Enum.concat
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

