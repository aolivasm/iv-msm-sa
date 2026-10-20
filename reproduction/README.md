# Reproducing the paper

This directory contains paper-specific analysis code, not a second R package.
For the user-facing estimators and examples, use [`ivmsm`](../packages/ivmsm/README.md).

From the repository root, run `cd reproduction` in your terminal before using
the commands in the [reproduction guide](../docs/reproduction.md).

- `experiments/`: continuous- and binary-outcome simulation drivers.
- `scripts/`: analysis, collection, and plotting commands.
- `R/`: supporting standalone implementations used by those scripts.
- `tests/`: small checks that do not launch the full studies.
- `data/`: instructions for obtaining OHIE data separately.

Results are generated locally and are not distributed. The analysis source
layout is preserved here so moving the directory does not change the methods,
seeds, or experiment settings.
