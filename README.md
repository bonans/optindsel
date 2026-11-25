# Selecting 2 × 2 submatrices with optimally bounded inverses

[![DOI][doi-badge]][doi-link]
[![arXiv][arxiv-badge]][arxiv-link]
[![MATLAB][code-badge]][code-link]
[![License: MIT][license-badge]](LICENSE.txt)

[doi-badge]: https://img.shields.io/badge/DOI-10.48550/arXiv.XXXX.XXXXX-blue
[doi-link]: https://doi.org/10.48550/arXiv.XXXX.XXXXX

[arxiv-badge]: https://img.shields.io/badge/arXiv-XXXX.XXXXX-red
[arxiv-link]: https://arxiv.org/abs/XXXX.XXXXX

[code-badge]: https://img.shields.io/badge/MATLAB-≥R2022a-blue.svg
[code-link]: https://www.mathworks.com/products/matlab.html

[license-badge]: https://img.shields.io/badge/License-MIT-green.svg

This repository contains the MATLAB implementations and reproducibility scripts for the following paper:

Bonan Sun<br>
*Selecting* $2 \times 2$ *submatrices with optimally bounded inverses*<br>
Preprint, 2026

## Requirements

- **MATLAB** version R2022a or later (verified in MATLAB R2025a; no external toolboxes required)

## Quick Start

### Using `opt2row`

Given a matrix $U \in \mathbb{R}^{n \times 2}$ ($n \ge 2$) with orthonormal columns ($U^\top U = I_2$), `opt2row` selects a two-element row index set $I = \lbrace i, j \rbrace$ ($i < j$) in a single streaming pass ($\mathcal{O}(n)$ operations, $\mathcal{O}(1)$ auxiliary memory):

```matlab
setpath();

[U, ~] = qr(randn(1000, 2), 'econ');
[I, wnorm, invnorm] = opt2row(U);          % guarantees ||U(I,:)^(-1)||_2 <= sqrt(n)

% General positive row weights w in (0, 1)^n with sum(w) <= 1:
w = rand(1000, 1); w = w / sum(w);
[I_w, wnorm_w] = opt2row(U, w);            % guarantees ||U(I_w,:)^(-1) diag(sqrt(w(I_w)))||_2 <= 1
```

### Reproducing the Paper Results

Run `reproduce_all` in MATLAB (~1 minute on a modern machine):

```matlab
>> reproduce_all                      % Compute experiments and save figure
>> reproduce_all('mode', 'compute')   % Headless compute only (writes data/opt2_experiments.mat)
>> reproduce_all('mode', 'plot')      % Plot figures/opt2_experiments.pdf from saved data
```

From the command line (headless batch mode):

```bash
matlab -batch "reproduce_all"
```

## Execution Modes

| Mode | Compute & save `.mat` | Plot & save `.pdf` |
|------|-----------------------|--------------------|
| `'compute&plot'` (default) | Yes | Yes |
| `'compute'` | Yes | No |
| `'plot'` | No (falls back if missing) | Yes |

## Output

| File | Description |
|------|-------------|
| [`data/opt2_experiments.mat`](data/opt2_experiments.mat) | Saved numerical results for both experiments in Section 5 |
| [`figures/opt2_experiments.pdf`](figures/opt2_experiments.pdf) | Fig. 1 (left & right): Perturbed 3-block matrices ($8 \le n \le 1024$) and Proposition 4.2 noise transition ($n = 512$) |

**Note on figures:** Figure 1 in the paper is rendered with TikZ/pgfplots directly from the values stored in `data/opt2_experiments.mat`. Running `figures_opt2.m` (or `reproduce_all`) generates the equivalent MATLAB vector figure `figures/opt2_experiments.pdf` from the same data.

## Code Structure

```
optindsel/
├── reproduce_all.m              # Main driver to run experiments and generate figures
├── setpath.m                    # Adds src/ and experiments/ to MATLAB path
├── src/                         # Core algorithms and baseline implementations
│   ├── opt2row.m                # Algorithm 2: single-pass O(n) streaming 2x2 row selection
│   ├── inv2norm.m               # Proposition 3.2: closed-form ||B^(-1)||_2 for 2x2 matrices
│   ├── threeblock2row.m         # Section 4: three-block orthonormal matrix generator
│   ├── exactmaxvol2row.m        # Exhaustive global maximum-volume 2x2 row selection
│   ├── exact2row.m              # Exhaustive minimum ||U_I^(-1) diag(sqrt(w_I))||_2 search
│   ├── maxvol2.m                # Local maximum-volume swapping (adapted from TT-Toolbox)
│   └── bsqr.m                   # Bischof-Stewart pivoted QR (adapted from Damle 2026)
├── experiments/                 # Reproducibility scripts for Section 5
│   ├── run_opt2.m               # Runs both experiments and saves data/opt2_experiments.mat
│   └── figures_opt2.m           # Plots figures/opt2_experiments.pdf from saved data
├── data/                        # Saved numerical records
│   └── opt2_experiments.mat     # Precomputed R2025a results for Section 5
└── figures/                     # Output figures
    └── opt2_experiments.pdf     # MATLAB reproduction of Figure 1
```

## Third-Party Attributions

- `src/maxvol2.m` is adapted from [`oseledets/TT-Toolbox`](https://github.com/oseledets/TT-Toolbox) (MIT License, Copyright (c) 2009–2012 Ivan Oseledets et al.).
- `src/bsqr.m` is adapted in MATLAB from [`asdamle/BischofStewartQR`](https://github.com/asdamle/BischofStewartQR) (MIT License, Copyright (c) 2026 Anil Damle).

## License

MIT License — see [LICENSE.txt](LICENSE.txt).

## Citation

If you use this code, please cite:

```bibtex
@misc{Sun2026opt2,
  title        = {Selecting {$2 \times 2$} submatrices with optimally bounded inverses},
  author       = {Sun, Bonan},
  howpublished = {arXiv preprint arXiv:XXXX.XXXXX},
  year         = {2026}
}
```
