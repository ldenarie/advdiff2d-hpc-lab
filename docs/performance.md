# Observations

The first [reference run](../results/20261109T141104Z_gha-ubuntu-24.04) reports
the following data:

| kernel | threads | grid | steps | reps | median MLUP/s | min–max | spread | eff. GB/s |
|---|---|---|---|---|---|---|---|---|
| seq | 1 | 128×128 | 12208 | 7 | 173.8 | 166.5–174.7 | 4.7 % | 4.17 |
| seq | 1 | 512×512 | 763 | 7 | 173.2 | 166.9–173.9 | 4.0 % | 4.16 |
| seq | 1 | 2048×2048 | 48 | 7 | 173.5 | 170.9–174.6 | 2.1 % | 4.16 |
| seq | 1 | 4096×4096 | 12 | 7 | 173.5 | 163.4–179.5 | 9.3 % | 4.16 |

We can observe that the throughput is approximately constant at ~173 MLUP/s,
meaning that the increase in the grid size does not slow down the processing.

**Hypothesis**: This suggests that the limit is in the processor operation, and
not in the memory bandwidth (otherwise, we would see an decrease when the size
of the grid increases, as the larger, but slower memory layers (caches) are
involved.).

# Analysis

## On the validity of performance comparison

The benchmarks run on GitHub Actions. As such, it is impossible to ensure that
two differents job run on the same machine. Even on a signle machine, it is
impossible to know the amount of ressources that is allocated to our job. In
order to make fair comparisons between the different version of our kernel, one
should only compare runs that were performed on the same job. All comparison
made further in this document are made this way. In particular, this explains
why, depending on the table, the same kernel can have different measures.

## Count instructions read per cell update

In order to analyse where the time is spent, I used the callgrind tool from
valgrind.

First, the driver needs to be compiled with the debugging information. The usual
usage is to compile with `-DCMAKE_BUILD_TYPE=RelWithDebInfo`, but by defaut,
this option compiles with the optimization `O2`. In order to profile the same
program that we want to run at the end, it is necessary to compile with `O3`.
This can be done by just adding the option `-g`.

```
cmake -S . -B build-profile -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS=-g
cmake --build build-profile/
```

Then, I can run callgrind with a short execution of the program, and get the
output.

```
valgrind --tool=callgrind ./build-profile/advdiff2d nx=256 ny=256 steps=20 reps=1 warmup=0
callgrind_annotate
```

Here is a sample of the output showcasing the function call with the most
instruction reads:

```
--------------------------------------------------------------------------------
Ir                   file:function
--------------------------------------------------------------------------------
 45,875,220 (58.07%)  src/include/advdiff/details/kernel.hpp:advdiff::step_seq(advdiff::Field const&, advdiff::Field&, advdiff::Grid const&, advdiff::Physics const&, double)
 7,926,280 (10.03%)  src/kernel_seq.cpp:advdiff::step_seq(advdiff::Field const&, advdiff::Field&, advdiff::Grid const&, advdiff::Physics const&, double) [/workspaces/advdiff2d-hpc-lab/build-profile/advdiff2d]
 7,884,800 ( 9.98%)  include/advdiff/problem.hpp:advdiff::step_seq(advdiff::Field const&, advdiff::Field&, advdiff::Grid const&, advdiff::Physics const&, double)
 6,553,640 ( 8.30%)  include/advdiff/field.hpp:advdiff::step_seq(advdiff::Field const&, advdiff::Field&, advdiff::Grid const&, advdiff::Physics const&, double)
 2,162,688 ( 2.74%)  ./math/../sysdeps/ieee754/dbl-64/e_exp.c:__ieee754_exp_fma [/usr/lib/x86_64-linux-gnu/libm.so.6]
 1,310,860 ( 1.66%)  /usr/include/c++/13/bits/stl_vector.h:advdiff::step_seq(advdiff::Field const&, advdiff::Field&, advdiff::Grid const&, advdiff::Physics const&, double)
 1,310,720 ( 1.66%)  ./math/./w_exp_template.c:exp@@GLIBC_2.29 [/usr/lib/x86_64-linux-gnu/libm.so.6]
   788,220 ( 1.00%)  src/main.cpp:main [/workspaces/advdiff2d-hpc-lab/build-profile/advdiff2d]
   592,128 ( 0.75%)  include/advdiff/problem.hpp:main
```

We can see that most of the time is spent in the kernel, for a total of
68,239,940 instruction reads. Relative to the number of cell update (nx * ny *
steps), it makes 52 instruction read per cell update.

# What did the compiler optimized

I temporarily added the `-fopt-info-vec-optimized -fopt-info-vec-missed` flags
to the compilation of `kernel_seq.cpp`.

The compilation of `kernel_seq.cpp` now outputs:

```
/src_dir/src/kernel_seq.cpp:16:23: missed: couldn't vectorize loop
/src_dir/src/kernel_seq.cpp:17:21: missed: not vectorized: no vectype for stmt: _3 = MEM[(double &)_57];
 scalar_
/src_dir/src/kernel_seq.cpp:12:23: optimized: loop vectorized using 16 byte vectors
/src_dir/src/kernel_seq.cpp:12:23: optimized:  loop versioned for vectorization because of possible aliasing
/src_dir/src/kernel_seq.cpp:23:23: missed: couldn't vectorize loop
/src_dir/src/kernel_seq.cpp:23:23: missed: not vectorized: control flow in loop.
/src_dir/src/kernel_seq.cpp:24:27: missed: couldn't vectorize loop
/src_dir/src/kernel_seq.cpp:24:27: missed: not vectorized: control flow in loop.
/src_dir/src/kernel_seq.cpp:32:23: missed: couldn't vectorize loop
/src_dir/src/kernel_seq.cpp:32:23: missed: not vectorized: control flow in loop.
/src_dir/src/kernel_seq.cpp:33:27: missed: couldn't vectorize loop
/src_dir/src/kernel_seq.cpp:33:27: missed: not vectorized: control flow in loop.
/src_dir/src/kernel_seq.cpp:49:23: missed: couldn't vectorize loop
/src_dir/src/kernel_seq.cpp:49:23: missed: not vectorized: control flow in loop.
/usr/inctl_numeric.h:140:22: optimized: loop vectorized using 16 byte vectors
/src_dir/src/kernel_seq.cpp:57:23: missed: couldn't vectorize loop
/src_dir/src/kernel_seq.cpp:57:23: missed: not vectorized: loop nest containing two or more consecutive inner loops cannot be vectorized
/usr/inctl_algo.h:5650:24: missed: couldn't vectorize loop
/src_dir/src/kernel_seq.cpp:55:8: missed: not vectorized: unsupported use in stmt.
/src_dir/src/kernel_seq.cpp:66:23: missed: couldn't vectorize loop
/src_dir/src/kernel_seq.cpp:66:23: missed: not vectorized: loop nest containing two or more consecutive inner loops cannot be vectorized
/usr/inctl_algo.h:5714:24: missed: couldn't vectorize loop
/src_dir/src/kernel_seq.cpp:64:8: missed: not vectorized: unsupported use in stmt.
```

The loops involved in the kernel functions are the one in lines 32 and 33. Both
were not vectorized because of "control flow in loop". One possible reason is
that the content of the loop involves a conditional execution (the advection
term depends on the signs of the velocity components u and v).

The fact that the compiler does not vectorize the loop only explais part of why
the thoughput is so low. Another possible explanation is the nubmer of operation
in each cell update. For every cell, multiple terms are computed, and these
computation involves double floating point divisions, which are very costly
compared to other operation.

# Optimization

## Remove conditional in for loop

Let's try to remove the conditional in the for loop so that the compiler can
vectorize it.  There are several possibilities. Probably, the one with the most
performance is isolate each possible branch in its own function, and call a
different for loop in each case. But this will be my second choice. My first
choice will be to add both term instead of using the conditional, and add
factors 1. and 0. that make the kernel only add the right value. This is my
first choice because that will allow me to check that the for loop is indeed
factorized.

After this change, the compiler reports that the inner loop (the one in `j`) was
vectorized. A new call to valgrind callgrind tool shows that the total number of
instruction read is 28'683'600. That makes around 22 instruction read per cell
update. Which is much better than the 52 that we had in the beginning.

There are still a lot of computation, including divisions, performed in each
cell update. This is particularly true because the kernel now computes useless
terms which are then canceled after a multiplication by 0.

## Remove conditional in the for loop

Now that we showed that the for loop can be vectorized and that it makes a big
difference, we can implement a version that does not compute the terms that will
be canceled by the multiplication by a factor 0. What we definitely do not want
is to write 4 different for loops for each combination of u and v signs. This
solution would work, but it would never scale for an equation more complex.

The main thing to notice is that the condition is not dependent on the cell. So
the result will be exactly the same for each cell in a single field update. And
the formula to update a cell is just a linear combination of the value of the
cell and the values of its neighbors. The factors associated to each cell in the
linear combination is the same for each cell. As such, they can be computed
outside of the loop. Actually, this is the formulation of the update written in
[the equation at §3 in the model
description](model.md#stability-and-maximum-principle). Seen individually, this
equation seems to involve more computation than the original formulation, with
the conditionals, but as the weight can be computed only one time for all the
cells, it is much more efficient to do it this way.

After implementing this version, the compiler still vectorizes the loop. The
number of instruction read per cell update now drops to 9.

But with this implementation the test that checks that a uniform field is
unchanged now fails. The sum of the factors is indeed not perfectly 1 (because
of numerical errors), and as such, the uniform field is not perfectly unchanged.
It is in theory possible to make the computation enforce the fact that a
uniform field is unchanged. By writing the update in the following manner:

```
c_new(i,j) = c(i,j) + (Cx+Dx) (c(i-1,j) - c(i,j)) + Dx (c(i+1,j) - c(i,j)) + (Cy+Dy) (c(i,j-1) - c(i,j)) + Dy (c(i,j+1) - c(i,j))
```

This incremental formulation ensure that a uniform field is unchanged as each
difference is 0. The price for this added precision is that there is now 4
additional substraction in place of one substraction and a multiplication. This
is a compromise that must be chosen using a proper measure. Running valgrind
callgrind tool, it shows that the number of instruction read per cell update is
now 10, which is a bit higher than 9. Running the benchmark on the CI, it show
the following numbers:

| kernel | threads | grid | steps | reps | median MLUP/s | min–max | spread | eff. GB/s |
|---|---|---|---|---|---|---|---|---|
| incremental | 1 | 128×128 | 12208 | 7 | 1142.5 | 1140.0–1143.8 | 0.3 % | 27.42 |
| incremental | 1 | 512×512 | 763 | 7 | 1145.4 | 1131.5–1147.4 | 1.4 % | 27.49 |
| incremental | 1 | 2048×2048 | 48 | 7 | 1124.2 | 1105.2–1134.9 | 2.6 % | 26.98 |
| incremental | 1 | 4096×4096 | 12 | 7 | 1127.9 | 1119.2–1138.6 | 1.7 % | 27.07 |
| naive | 1 | 128×128 | 12208 | 7 | 117.7 | 117.2–117.8 | 0.5 % | 2.83 |
| naive | 1 | 512×512 | 763 | 7 | 117.6 | 117.3–117.6 | 0.3 % | 2.82 |
| naive | 1 | 2048×2048 | 48 | 7 | 117.4 | 117.3–117.5 | 0.1 % | 2.82 |
| naive | 1 | 4096×4096 | 12 | 7 | 117.4 | 116.1–117.5 | 1.1 % | 2.82 |
| weighted | 1 | 128×128 | 12208 | 7 | 1377.3 | 1375.1–1379.8 | 0.3 % | 33.06 |
| weighted | 1 | 512×512 | 763 | 7 | 1437.8 | 1431.9–1439.3 | 0.5 % | 34.51 |
| weighted | 1 | 2048×2048 | 48 | 7 | 1303.6 | 1240.8–1314.1 | 5.6 % | 31.29 |
| weighted | 1 | 4096×4096 | 12 | 7 | 1258.1 | 1248.0–1264.2 | 1.3 % | 30.19 |

Going from the naive to the optimized version goes from 117 to 1437 MLUP/s (for
a 512 times 512 grid). While going from the optimized version to the one that is
optimized but is more precise only drops down to 1145 MLUP/s. For sake of
precision (which is important on long run simulation), I would argue that the
version that keeps a uniform field unchanged is worth it. This is espacially
true because for larger grids, the limitation is the memory bandwidth, so the
difference between the weighted and incremental version is even lower. Also, one
thing that must be noted is that the kernel I decided to call more precise seems
to be indeed more precise. The total mass relative error for the incremental
kernel shows that in the data. On a coarser grid (128 x 128), the total mass
relative error stays in the order of magnitude of 1e-16. For the weighted
kernel, this error goes up to 1e-12 (in order of magnitude). For the sake of
this experiment, I will just keep the two kernel versions available.

