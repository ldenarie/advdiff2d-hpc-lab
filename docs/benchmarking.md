# Benchmarking protocol

How performance is measured, and what the numbers can and cannot say.

## 1. Driver contract

`bench/run_bench.sh` calls the driver once per configuration and expects:

```sh
advdiff2d kernel=seq nx=N ny=N steps=S reps=R warmup=W output=csv
```

| Option | Meaning |
|---|---|
| `kernel` | kernel to run: `seq` for now; unknown values are an error |
| `steps` | time steps in one repetition |
| `reps` | number of timed repetitions |
| `warmup` | number of untimed repetitions run before the timed ones |
| `output` | `text` (default, human-readable) or `csv` |

With `output=csv`, standard output holds exactly one header line and one line
per timed repetition, and nothing else (diagnostics go to standard error):

```
kernel,nx,ny,steps,rep,time_s,mlups,mass_rel_err
seq,512,512,763,0,0.4421,449.8,4.4e-16
...
```

- `time_s`: wall-clock time of the `steps` time steps of one repetition,
  measured with `std::chrono::steady_clock`. It includes the periodic boundary
  update, which is part of a step; it excludes allocation and initialization.
- `mlups`: `nx * ny * steps / time_s / 1e6`, million lattice (cell) updates
  per second.
- `mass_rel_err`: relative change of the total mass during the repetition,
  `|M_end - M_start| / M_start`, where `M_start` is the mass at the start of
  that repetition. A sanity check that the run computed something meaningful.
- The field is allocated and initialized once; repetitions continue the same
  simulation. Warm-up repetitions absorb first-touch page faults and cache
  effects of the first steps.

`run_bench.sh` prepends `label,commit,threads` to each line and checks the
header and the number of lines.

## 2. Protocol

- **Release build, no sanitizers.** The script refuses other builds unless
  `-f` is given.
- **Repetitions.** Each configuration runs 1 warm-up and 7 timed repetitions by
  default. Results report the **median**, with min and max to show the spread.
- **Duration.** The number of steps is chosen so that each repetition performs
  about 2·10⁸ cell updates (tenths of a second to a few seconds), well above
  the clock resolution and the timer overhead.
- **Grid sizes** cover the memory hierarchy: 128² (two fields ≈ 0.3 MB,
  in L2), 512² (≈ 4 MB, around L2/L3), 2048² (≈ 67 MB) and 4096² (≈ 270 MB),
  far above the last-level cache, in main memory.
- **Test case.** Driver defaults: Gaussian of width 0.05 L centred in the
  domain, u = 1, v = 0.5, kappa = 1e-4, S = 0.9. The kernel is memory-bound
  and has no data-dependent branch, so its speed hardly depends on these
  values; they are fixed so that every run measures the same computation.
- **Traceability.** Each run writes `results/<UTC time>_<label>/` with
  `raw.csv` (every repetition) and `env.txt` (commit, compiler, flags, CPU
  model, `lscpu`, memory, `OMP_*` variables).

## 3. Derived quantities

- **Effective bandwidth** = MLUP/s × bytes per update. The default of 24 bytes
  assumes the minimum traffic per cell: read `c`, write `c_new`, and the
  write-allocate read of `c_new`; neighbours are assumed to come from cache.
  It is a model, not a measurement, and is labelled as such.
- **Speedup** S(p) = T(1) / T(p), and **efficiency** E(p) = S(p) / p, from
  median times, against the same kernel on 1 thread (not against `seq`, unless
  stated).

## 4. Limits of these measurements

- **Shared virtual machines.** GitHub-hosted runners are VMs on shared hosts:
  the CPU model changes between jobs and neighbours add noise. Only results
  from **the same job** are compared; `env.txt` records the machine.
- **Few cores.** Public-repository runners have 4 vCPUs, which may be 2
  physical cores with 2 hardware threads each (see `Thread(s) per core` in
  `env.txt`). Scaling beyond the number of physical cores is not expected to
  be linear, and nothing is extrapolated beyond 4.
- **One node.** MPI runs on one machine communicate through shared memory: no
  network latency or bandwidth is measured.
- **Frequency.** Turbo and frequency scaling are not controlled on these VMs.
