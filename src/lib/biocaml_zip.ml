open Biocaml_internal_pervasives

let dbg fmt = Debug.make "ZIP" fmt

module Transform = struct
  type unzip_error =
  [ `garbage_at_end_of_compressed_data of string
  | `zlib of string
  | `wrong_gzip_header of
      [ `compression_method | `flags | `magic_number ] * int ]
  with sexp

  open Result
  let try_skip_gzip_header_exn buffer =
    let bytes_read = ref 0 in
    let ignore_bytes n = bytes_read := !bytes_read + n in
    let assert_byte f error =
      let byte = (Char.to_int buffer.[!bytes_read]) in
      if f byte
      then (incr bytes_read; return byte)
      else fail (`wrong_gzip_header (error, byte)) in
    let assert_byte_eq b error = assert_byte ((=) b) error in
    let read_byte () =
      let byte = (Char.to_int buffer.[!bytes_read]) in
      (incr bytes_read; return byte) in
    let rec skip_null_terminated () =
      read_byte () >>= fun x ->
      if x = 0 then return () else skip_null_terminated () in
    assert_byte_eq 0x1f `magic_number >>= fun _ ->
    assert_byte_eq 0x8b `magic_number >>= fun _ ->
    assert_byte_eq 8 `compression_method >>= fun _ ->
    assert_byte (fun b -> b land 0xe0 = 0) `flags >>= fun flags ->
    ignore_bytes 6;
    (if flags land 0x04 <> 0 then (
    (* Skip extra data *)
      read_byte () >>= fun len1 ->
      read_byte () >>= fun len2 ->
      ignore_bytes (len1 + len2 lsl 8);
      return ()
     ) else return ())
    >>= fun () ->
    (if flags land 0x08 <> 0 then skip_null_terminated () else return ())
    >>= fun () ->
    (if flags land 0x10 <> 0 then skip_null_terminated () else return ())
    >>= fun () ->
    if flags land 0x02 <> 0 then begin
    (* Skip header CRC *)
      ignore_bytes 2;
    end;
    if !bytes_read > String.length buffer then
      failwith "NOT-YET"
    else
      return !bytes_read

  let inflate_as_much_as_possible in_buffer buffer
      zstream zlib_write_buffer zlib_buffer_size format =
    let len = String.length buffer in
    let try_to_output  used_out =
      if used_out > 0
      then (`output (Ok (String.sub zlib_write_buffer 0 used_out)))
      else (`not_ready) in
    let (finished, used_in, used_out) =
      Zlib.inflate zstream buffer 0 len
        zlib_write_buffer 0 zlib_buffer_size
        Zlib.Z_SYNC_FLUSH in
    dbg "unzip len: %d, finished: %b, used_in: %d, used_out: %d" len
      finished used_in used_out;
    if used_in < len then (
      if finished then (
        match format with
        | `gzip when len - used_in >= 8 ->
        (* The 8-bytes CRC must be skipped then another Gzip header +
           stream can follow *)
          let continue_after_gzip = len - used_in > 8 in
          if continue_after_gzip then
            Buffer.add_string in_buffer
              String.(sub buffer (used_in + 8) (len - used_in - 8));
          `finished_gzip (try_to_output used_out)
        | `gzip ->
          Buffer.add_string in_buffer
            String.(sub buffer used_in (len - used_in));
          try_to_output used_out
        | _ ->
          `error (`garbage_at_end_of_compressed_data
                     String.(sub buffer used_in (len - used_in)))
      ) else (
        Buffer.add_string in_buffer
          String.(sub buffer used_in (len - used_in));
        try_to_output used_out
      )
    ) else try_to_output used_out

  let unzip ?(format=`raw) ?(zlib_buffer_size=4096) () =
    let zstream =  ref (Zlib.inflate_init false) in
    let in_buffer = Buffer.create 42 in
    let zlib_write_buffer = String.create zlib_buffer_size  in
    dbg "unzip init: zlib_buffer_size: %d" zlib_buffer_size;
    let current_state =
      ref (match format with `gzip -> `gzip_header | `raw -> `inflate) in
    let rec next stopped =
      let buffered = Buffer.contents in_buffer in
      let len = String.length buffered in
      Buffer.clear in_buffer;
      dbg "unzip: len: %d" len;
      begin match len with
      | 0 -> if stopped then `end_of_stream else `not_ready
      | _ ->
        begin match !current_state with
        | `inflate ->
          begin
            try 
              let inflation =
                inflate_as_much_as_possible in_buffer buffered
                  !zstream zlib_write_buffer zlib_buffer_size  format in
              begin match inflation with
              | `output o -> `output o
              | `error e -> `output (Error e)
              | `not_ready -> `not_ready
              | `finished_gzip out ->
                current_state := `gzip_header;
                zstream :=  Zlib.inflate_init false;
                out
              end
            with
            | e ->
              dbg "E: %s; len: %d" (Exn.to_string e) len;
              `output (Error (`zlib (Exn.to_string e)))
          end
        | `gzip_header ->
          begin
            try
              match try_skip_gzip_header_exn buffered with
              | Ok bytes_read ->
                current_state := `inflate;
                dbg "header bytes_read: %d" bytes_read;
                Buffer.add_string in_buffer
                  String.(sub buffered bytes_read (len - bytes_read));
                next stopped
              | Error e -> `output (Error e)
            with
              e -> Buffer.add_string in_buffer buffered; `not_ready
          end
        end
      end
    in
    Biocaml_transform.make_stoppable_with_error ()
      ~feed:(fun string -> Buffer.add_string in_buffer string;) ~next

  let gzip_default_header =
  (* ID1, Id2, compr meth, flags,
     mtime (4 bytes),
     xflags, OS (unknown) *)
    "\x1F\x8B\x08\x00\
   \x00\x00\x00\x00\
   \x00\xff"

  let zip ?(format=`raw) ?(level=8) ?(zlib_buffer_size=4096) () =
    let zstream =  ref (Zlib.deflate_init level false) in
    let in_buffer = Buffer.create 42 in
    let zlib_write_buffer = String.create zlib_buffer_size  in
    dbg "zip init: level: %d, zlib_buffer_size: %d" level zlib_buffer_size;
    let state =
      ref (match format with `raw -> `deflating | `gzip -> `gzip_header) in
    let this_is_the_end = ref false in
    let update_crc, output_crc =
      match format with
      | `raw -> ((fun _ _ -> ()), (fun () -> ""))
      | `gzip ->
        let gzip_crc = ref 0l in
        let gzip_size = ref 0l in
        ((fun buf used_in  ->
          gzip_crc := Zlib.update_crc !gzip_crc buf 0 used_in;
          gzip_size := Int32.(!gzip_size + (of_int_exn used_in));
          ()),
         (fun () ->
           let buf = String.create 8 in
           dbg "crc: %ld, size: %ld" !gzip_crc !gzip_size;
           Binary_packing.pack_signed_32
             ~byte_order:`Little_endian ~pos:0 ~buf !gzip_crc;
           Binary_packing.pack_signed_32
             ~byte_order:`Little_endian ~pos:4 ~buf !gzip_size;
           buf))
    in
    let next stopped =
      match state.contents with
      | `gzip_header -> state := `deflating; `output gzip_default_header
      | `deflating ->
        let buffered = Buffer.contents in_buffer in
        begin match String.length buffered with
        | 0 ->
          if stopped
          then begin
            if !this_is_the_end then `end_of_stream
            else begin
              let (_, _, used_out) =
                Zlib.deflate !zstream "" 0 0
                  zlib_write_buffer 0 zlib_buffer_size Zlib.Z_FINISH in
              let to_output =
                if used_out < zlib_buffer_size then (
                  this_is_the_end := true;
                  Zlib.deflate_end !zstream;
                  String.(sub zlib_write_buffer 0 used_out ^ output_crc ())
                ) else (
                  String.(sub zlib_write_buffer 0 used_out)
                ) in
              dbg "at the pseudo-end... used_out: %d, outputting %d bytes"
                used_out String.(length to_output);
              `output to_output
            end
          end
          else `not_ready
        | len ->
          Buffer.clear in_buffer;
          dbg "zip: len: %d" len;
          let (_, used_in, used_out) =
            Zlib.deflate !zstream buffered 0 len
              zlib_write_buffer 0 zlib_buffer_size
            (* Zlib.Z_NO_FLUSH *)
              Zlib.Z_SYNC_FLUSH
          in
          dbg "used_in: %d -- used_out: %d" used_in used_out;
          update_crc buffered used_in;
          if used_in < len
          then (Buffer.add_substring in_buffer
                  buffered used_in (String.length buffered - used_in));
          if used_out > 0 then
            `output String.(sub zlib_write_buffer 0 used_out)
          else
            `not_ready
        end
    in
    Biocaml_transform.make_stoppable ()
      ~feed:(fun string -> Buffer.add_string in_buffer string;) ~next 
end
