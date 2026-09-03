# Prototype

The first, flatter pass at the same method, kept for reference. Each operator
lives in its own file and the experiments are single scripts with the settings
inlined:

| File | Role |
|------|------|
| `S_operator.m` | sample-and-hold |
| `P_operator.m` | modular band-limiter |
| `ideal_lowpass_filter.m` | the LPF used inside **P** |
| `classical_hybrid_interpolator.m` | relaxation update |
| `chebyshev_hybrid_interpolator.m` | Chebyshev-accelerated update |
| `main_test.m`, `test_final.m`, `validation.m` | early experiments |
| `search_AB.m` | first version of the frame-bound grid search |

The code in the repository root supersedes all of this — it is the version to
read and to run.
