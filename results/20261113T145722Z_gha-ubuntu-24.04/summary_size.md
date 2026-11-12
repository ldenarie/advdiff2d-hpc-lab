Results for gha-ubuntu-24.04 @ b597fd5e08bf. Effective bandwidth assumes 24 bytes per cell update.
Ceiling: element rate of the STREAM scale kernel in the same job (main memory). It bounds grids larger than the last-level cache; cache-resident grids can exceed it.

| kernel | threads | grid | steps | reps | median MLUP/s | min–max | spread | eff. GB/s | % of ceiling |
|---|---|---|---|---|---|---|---|---|---|
| incremental | 1 | 128×128 | 12208 | 7 | 1409.2 | 1406.2–1418.0 | 0.8 % | 33.82 | 70 % |
| incremental | 1 | 512×512 | 763 | 7 | 1357.3 | 1353.7–1358.7 | 0.4 % | 32.57 | 68 % |
| incremental | 1 | 2048×2048 | 48 | 7 | 946.5 | 937.6–968.1 | 3.2 % | 22.72 | 47 % |
| incremental | 1 | 4096×4096 | 12 | 7 | 897.5 | 885.5–900.1 | 1.6 % | 21.54 | 45 % |
| naive | 1 | 128×128 | 12208 | 7 | 135.5 | 135.3–135.5 | 0.2 % | 3.25 | 7 % |
| naive | 1 | 512×512 | 763 | 7 | 135.1 | 135.0–135.1 | 0.1 % | 3.24 | 7 % |
| naive | 1 | 2048×2048 | 48 | 7 | 135.1 | 135.1–135.2 | 0.1 % | 3.24 | 7 % |
| naive | 1 | 4096×4096 | 12 | 7 | 135.2 | 131.6–135.3 | 2.7 % | 3.24 | 7 % |
| weighted | 1 | 128×128 | 12208 | 7 | 1738.2 | 1697.6–1768.5 | 4.1 % | 41.72 | 87 % |
| weighted | 1 | 512×512 | 763 | 7 | 1763.9 | 1755.0–1769.3 | 0.8 % | 42.33 | 88 % |
| weighted | 1 | 2048×2048 | 48 | 7 | 979.3 | 968.2–987.8 | 2.0 % | 23.50 | 49 % |
| weighted | 1 | 4096×4096 | 12 | 7 | 944.0 | 940.8–952.8 | 1.3 % | 22.66 | 47 % |
