# Universal Scaling in One-Dimensional Non-Reciprocal Matter

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

> Code associated with the paper "Universal scaling in one-dimensional non-reciprocal matter" by Shuoguang Liu, Ryo Hanai, and Peter B. Littlewood ([arXiv:2503.14384](https://arxiv.org/abs/2503.14384)).


## 📖 Overview

This repository contains the code needed to reproduce the main results of the paper:

> **Universal scaling in one-dimensional non-reciprocal matter**
>
> Shuoguang Liu, Ryo Hanai, and Peter B. Littlewood, 2025.

The scripts are written in **MATLAB 2024** and generate both simulation data and figures used in the main text and supplementary materials.


## 🚀 Features

- Full codebase to reproduce figures from the main text and supplementary materials.
- Simulations of critical fluctuations and dynamical pattern formation.
- Calculation of correlation functions in both time and frequency domains.


## 🛠️ Requirements

- MATLAB 2024 or later


## 💡 Usage

1. Clone the repository:

   ```bash
   git clone https://github.com/<your-username>/<your-repo>.git
   cd <your-repo>
   ```

2. Open MATLAB and set the path to the cloned folder.

3. To generate data:

   Run:
   ```matlab
   main_computation
   ```

4. To compute correlation functions:

   Use:
   ```matlab
   avg_corr_psi
   FFT
   ```


## 📊 Figures Correspondence

### Main Text Figures

- **Fig. 1:** (no code association)
- **Fig. 2:**
  - (a) `xover_jp.m`
  - (b) `collapse_CEP.m` and `collapse_EW.m`
- **Fig. 3:**
  - (a) `alpha_CEP.m`
  - (b) `fit_peak.m` and `z_CEP.m`
  - (c-d) same as Fig. 2(b)
- **Fig. 4:**
  - (a) `log_scaling.m` and `power_scaling.m`
  - (b) `xover_noise.m`
  - (c) `beta_pattern.m` and `phase_diagram.m`
- **Fig. 5:**
  - (a) `pattern_formation.m`
  - (b) `topo_vortex_evolve.m`


### Supplementary Information Figures

- **Fig. S1:** `collapse_CEP.m`
- **Fig. S2:** `waiting_time_CEP.m` and `waiting_time_EW.m`
- **Fig. S3-S5:** `collapse_width_CEP.m` (reused)
- **Fig. S6:** `freq_noise.m` and `freq_size.m`
- **Fig. S7:** `amp_pattern.m`
- **Fig. S8:** `psi_beta_jp.m` and `beta_pattern.m`
- **Fig. S9:** same as main text Fig. 5(b)
- **Fig. S10:** `main_computation.m` and `pattern_fit_peak.m`
- **Fig. S11:** `diffusion_freq.m`
- **Fig. S12:** `compare_correlator.m`
- **Fig. S13:** `FFT_scaling_width.m`


## 📚 Citation

If you use this code in your work, please cite:

> Shuoguang Liu, Ryo Hanai, Peter B. Littlewood, "Universal scaling in one-dimensional non-reciprocal matter," 2025, [arXiv:2503.14384](https://arxiv.org/abs/2503.14384).

BibTeX:

```bibtex
@article{liu2025universal,
  title={Universal scaling in one-dimensional non-reciprocal matter},
  author={Liu, Shuoguang and Hanai, Ryo and Littlewood, Peter B.},
  year={2025},
  eprint={2503.14384},
  archivePrefix={arXiv},
  primaryClass={cond-mat.stat-mech}
}
```


## ⚖️ License

This project is licensed under the [MIT License](LICENSE).


## 👥 Authors & Acknowledgments

- **Shuoguang Liu** – Lead Developer
- **Ryo Hanai** – Coauthor and Advisor
- **Peter B. Littlewood** – Principal Investigator

Special thanks to our collaborators and the University of Chicago for providing computational resources.

---

<div align="center">
  _"Exploring the universal dynamics beyond equilibrium."_
</div>
