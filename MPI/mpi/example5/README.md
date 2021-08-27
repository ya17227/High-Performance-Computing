Example 5: Communication Pattern Skeletons
===========================================

In this directory you will find a number of example programs showing different
methods of communication: simple send & receive, combined send-receive, asynchronous send
& receive and, finally, one-sided communications.

You may find these communication 'skeletons' useful for your assignment.
You may copy
the skeleton and adapt it to the requirements of the assignment.

I've also included an example where you can attach a debugger to a running process to
find a deadlock.


skeleton1
---------

I have used `MPI_Ssend()` to emphasise the importance of the communication pattern.
- See what happened if you, for example, change the even processes to send to the left
  before the right.

skeleton2
---------

In this case I have used `MPI_Sendrecv()`.
This pattern is harder to break.
You could, of course, set `left` or `right` to an incorrect rank, and that will break it.

skeleton2-simple2d
------------------

Here, we see `MPI_Sendrecv()` used for a halo-exchange on a 2d-grid.
This is getting more
like what is required for the assignment.
This is a very simple program that does only
one halo exchange and then examines the result.

skeleton2-heated-plate
----------------------

Even more like the assignment, we now see `MPI_Sendrecv()` used to calculate values on a
2d-grid following a time-stepping algorithm.
Each process uses a stencil to compute
new values for the heat distribution across a rectangular plate and the time-stepping
proceeds until some form of convergence is reached (in this case, a fixed number of
steps is performed).

skeleton3
---------

This example uses a similar pattern to skeleton2, but now uses asynchronous comms.
**Beware, this pattern is very fragile!**
For example, try commenting out the call to `MPI_Waitall()`, recompile and re-run.
Errors due to missing required synchronisations will be compounded if your comms are in a loop.

skeleton4
---------

This skeleton uses RMA, i.e. _one-sided_, MPI calls to achieve a similar pattern to the other skeletons.

debugging a deadlock
--------------------

An example gdb session:

```
> cd debugging/mpi
> mpirun -np 2 ./deadlock.exe 1
> gdb
```

With gdb started, you can begin to interrogate the program
(see [gdb_cheat_sheet.pdf](gdb_cheat_sheet.pdf) in this dir for more on gdb commands):

```
(gdb) attach <process-id>
(gdb) print rank
(gdb) set debugWait = 0
(gdb) next
(gdb) next
(gdb) next
(gdb) next
...
```

You'll see the line on which the program hangs.
(Note that if a function does not return,
then gdb itself will be hung too.)

To detach the debugger from the running process:

```
(gdb) detach
```

