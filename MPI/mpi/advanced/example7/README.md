Example 7: Composite Messages
=============================

Sending many messages between processes can significantly
slow your programs.
Especially if the communications network
has a high latency.

If you need to send many pieces of information, one approach
to improving efficiency is to combine several pieces of
information into a single message.
MPI offers us two ways to
agglomerate our variables.
The first is to create a user
derived type which can be a compound of intrinsic types.
The
second is to pack a number of variables into a message buffer
prior to transmission and unpack them again at the receiving
end.

Both of these approaches can be very effective.
However, care
must be taken.
Pitfalls include:

1. If we try to agglomerate too much information on a process
   before transmission, we could end up delaying the
   receiving process which could have been kept busy if only
   information had been sent to them earlier.

2. We should also be mindful of the computation cost of the
   agglomeration on the sending process.  If the communication
   network has a very low latency, we may have done better
   by sending the information piece-meal.

pack-unpack
-----------

Here we pack variables into a single message buffer on the
sending process and unpack them on the receiving process.
This is a direct and no-nonsense approach.
However, this
approach requires that the programmer explicitly carry
out the packing and unpacking operations, which can be
error prone.
Also this approach is computationally more
expensive than its alternative should a similar message
be sent several times in a program.

derived-type
------------

This approach requires more setup but is computationally
more efficient if used several times and is less error
prone as the compound message is automatically broken into
its constituent parts.

The derived-type is based around the offsets in memory
between the storage locations of the variables concerned.

The function `MPI_Address()` is provided so that the offsets
can be reliably computed—pointer arithmetic is not strictly
guaranteed outside of the same array.

### Exercise

Devise programs which employ the other methods of constructing user-derived types:
- `MPI_Type_contiguous()`
- `MPI_Type_indexed()`
