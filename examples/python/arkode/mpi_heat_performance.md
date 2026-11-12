# MPI heat performance on a shared GPU

For these examples, CPU wins on small vectors and CUDA wins on the largest batch.
Two MPI ranks sharing one GPU are slower than one GPU rank. These results do not
measure multi-GPU scaling or CUDA-aware MPI.

## Configuration

Measured on November 12, 2026 with a Xeon W-2255 (10 physical cores), one RTX A4000
(16 GiB), NVIDIA driver 620.70.06, and CUDA toolkit 13.5. The GPU also drives the
desktop. The build uses RelWithDebInfo, CUDA architecture 86, and double precision.
Versions are Open MPI 5.0.5, NumPy 2.5.1, JAX 0.12.0.1, and Torch 2.15.0+cu132.

All configurations use 201 interior grid rows, independent heat fields, two sine
modes, `k=0.5`, final time `0.1`, relative tolerance `1e-6`, and absolute tolerance
`1e-10`. The solver is LSRKStep RKC2 with an analytic spectral bound and the I
controller. The default unpreconditioned implicit example is not this benchmark.
Torch uses `torch.compile`; JAX uses `jax.jit`. Each process warms one full solve
before measurement. Framework compilation and imports are outside integration time.

CPU math-library thread counts are one. Open MPI uses `--bind-to none`. JAX GPU
preallocation is disabled. GPU ranks share visible device 0 and stage halo rows
through host buffers. Backend order rotates over five samples per configuration.

## Synchronized integration times

Times are seconds. Each cell is the median followed by the minimum and maximum
of five runs. The timer includes complete `ARKodeEvolve` execution and final GPU
synchronization, using the slowest rank. Setup, warmup, and validation are outside.

| Ranks | Unknowns | NumPy CPU | JAX CUDA | Torch CUDA |
|---|---:|---:|---:|---:|
| 1 | 201 | 0.0180 (0.0180–0.0181) | 0.734 (0.729–0.741) | 0.335 (0.329–0.339) |
| 1 | 51,456 | 0.354 (0.346–0.361) | 0.766 (0.765–0.778) | 0.354 (0.348–0.356) |
| 1 | 823,296 | 12.070 (11.974–12.131) | 1.191 (1.185–1.192) | 0.803 (0.801–0.805) |
| 2 | 201 | 0.0211 (0.0209–0.0216) | 1.315 (1.297–1.323) | 0.931 (0.924–0.931) |
| 2 | 51,456 | 0.179 (0.178–0.183) | 1.382 (1.378–1.425) | 0.993 (0.984–1.006) |
| 2 | 823,296 | 5.362 (5.329–5.393) | 1.854 (1.850–1.862) | 1.425 (1.409–1.434) |

At 823,296 unknowns, CUDA reduces warmed integration time by about 10.1 times
with JAX and 15.0 times with Torch at one rank. At two ranks, the factors are 2.9
and 3.8. CPU and Torch have no measurable difference at the one-rank middle size
because their ranges overlap. GPU overhead dominates the small cases.

All 90 samples completed without errors, convergence failures, or stderr output.
The largest analytic error was `3.73e-6`. Backends have identical work counts for
each configuration: 95 steps and 1,140 RHS calls for one field, and 103 steps and
1,188 RHS calls for the larger batches. The timed work is not a cached solution.

## Startup and warmup change the comparison

Median process times include launcher startup, imports, preparation, warmup,
integration, validation, printing, and teardown:

| Ranks | Unknowns | NumPy CPU | JAX CUDA | Torch CUDA |
|---|---:|---:|---:|---:|
| 1 | 823,296 | 24.200 s | 3.630 s | 8.131 s |
| 2 | 823,296 | 10.982 s | 5.153 s | 9.375 s |

Torch's shorter warmed solve does not mean its fresh process finishes first.
Compiler disk caches were warm, so these are not first-ever installation timings.

## What limits the solve

Separate Python profiles at the largest one-rank size attribute about 6.55 s of
the CPU run to NumPy RHS callbacks and 5.52 s to solver execution outside those
callbacks. The large CPU workload is vector and stencil work, not halo exchange.
This profile does not distinguish CPU memory bandwidth from native arithmetic.

For GPU variants, callbacks dominate the host profile. JAX spends time in DLPack
conversion, reshape dispatch, host/device halo transfers, and immutable-result
publication. Torch spends time in compiled-kernel dispatch, DLPack conversion,
halo transfers, and synchronization. Each solve makes 1,188 synchronous RHS
callbacks. Adding ranks increases this coordination while sharing the same GPU.
These profiles locate expensive code paths but do not prove a kernel-level
hardware bandwidth limit.

Nsight Systems tracing failed inside CUPTI due to conflicting injected/runtime
libraries. No CUDA trace timings are used here. Python profiling runs are separate
from the 90 reported samples.

## Reproduce

Follow [Run MPI heat examples](README_mpi.md), then run:

```sh
python scripts/benchmark_sundials4py_mpi.py \
	--ranks 1 2 --fields 1 256 4096 --repeats 5 \
	--output /tmp/mpi-heat-results.jsonl
```

The original raw samples are at
`/tmp/opencode/mpi-heat-performance-synchronized.jsonl` on the test workstation.
Profiles and detailed notes are under `/tmp/opencode/mpi-heat-*`. The comparison
applies to these batched heat equations and their current callback interfaces,
not to arbitrary coupled PDEs or framework workloads.
