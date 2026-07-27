# Universal Scaling in One-Dimensional Non-Reciprocal Matter

[![License](https://img.shields.io/badge/license-BSD%202--Clause-blue.svg)](LICENSE)

> Code associated with the paper **“Universal scaling in one-dimensional non-reciprocal matter”** by Shuoguang Liu, Peter B. Littlewood and Ryo Hanai ([arXiv:2503.14384](https://arxiv.org/abs/2503.14384)).


## 📖 Overview

This repository contains the code needed to reproduce the main results of the paper:

> **Universal scaling in one-dimensional non-reciprocal matter**  
> Shuoguang Liu, Ryo Hanai, and Peter B. Littlewood

The scripts are written for **MATLAB R2024**.

All raw simulation data presented in this work can be generated using the scripts in the `/computation` directory. The scripts for data analysis and figure generation are provided in the `/analysis` directory. The processed data required to reproduce the figures are available in the `/data` directory.

Below is a list of the files corresponding to each figure.


## 📊 Figures Correspondence

### Main Text

- **Fig. 1**
  - (a) `main1.m` — for insets
  - (b–c) no code associated

- **Fig. 2**
  - (a) `main1.m`, `xover_chi_all.m`
  - (b) `main1.m`, `size_CEP.m`, `collapse_CEP.m`, `size_EW.m`, and `collapse_EW.m`

- **Fig. 3**
  - (a) `width_CEP`, `alpha_CEP.m`
  - (b) `Fit_peak.m`, `FFT_scaling.m`, `z_CEP.m`
  - (c–d) same as Fig. 2(b)

- **Fig. 4**
  - (a) `phaseA_evolve.m`
  - (b) `size_pattern.m`, `log_scaling.m`, and `power_scaling.m`
  - (c) `xover_freq_all.m`


### Supplementary Information

- **Fig. S1**
  - (a) `main2.m`, `CEP_approx.m`
  - (b) `main2.m`, `EW_approx.m`

- **Fig. S2**
  - `compare_correlator_EW.m`

- **Fig. S3**
  - (a–b) `compare_correlator_CEP.m`
  - (c) `main3.m`, `ensemble_CEP.m`

- **Fig. S4**
  - (a–b) `main2.m`, `compare_correlator_pattern.m`
  - (c) `main3.m`, `ensemble_pattern.m`

- **Fig. S5**
  - (a–b) `main4.m`, `collapse_CEP_CBB.m`

- **Fig. S6**
  - (a) `waiting_time_CEP.m`
  - (b) `waiting_time_EW.m`

- **Fig. S7**
  - (a–b) `amp.m`

- **Fig. S8**
  - `alpha_CEP_robust.m`

- **Fig. S9**
  - (a–b) `collapse_width_CEP.m`

- **Fig. S10**
  - `Fit_peak.m`, `FFT_scaling_width.m`

- **Fig. S11**
  - (a) `collapse_width_CEP.m`
  - (b) `alpha_CEP.m`

- **Fig. S12**
  - (a) `collapse_width_EW.m`
  - (b) `size_EW.m`

- **Fig. S13**
  - (a–b) `amp_pattern.m`

- **Fig. S14**
  - (a) `freq_size.m`
  - (b) `freq_noise.m`

- **Fig. S15**
  - (a) `main1.m`
  - (b–c) no code associated

- **Fig. S16**
  - (a–c) `log_scaling.m`


## 📚 Citation

If you use this code in your work, please cite:

> Shuoguang Liu, Ryo Hanai, and Peter B. Littlewood,  
> **“Universal scaling in one-dimensional non-reciprocal matter,”**  
> [arXiv:2503.14384](https://arxiv.org/abs/2503.14384).

```bibtex
@article{liu2025universal,
  title         = {Universal scaling in one-dimensional non-reciprocal matter},
  author        = {Liu, Shuoguang and Hanai, Ryo and Littlewood, Peter B.},
  year          = {2025},
  eprint        = {2503.14384},
  archivePrefix = {arXiv},
  primaryClass  = {cond-mat.stat-mech}
}
```


## ⚖️ License

This project is licensed under the [BSD 2-Clause License](LICENSE).


## 👥 Authors

- **Shuoguang Liu**  
  [shuoguang@uchicago.edu](mailto:shuoguang@uchicago.edu)  
  James Franck Institute and Department of Physics, University of Chicago, Chicago IL 60637, USA

- **Ryo Hanai**  
  [hanai.r.7e4b@m.isct.ac.jp](mailto:hanai.r.7e4b@m.isct.ac.jp)  
  Department of Physics, Institute of Science Tokyo, 2-12-1 Ookayama Meguro-ku, Tokyo, 152-8551, Japan

- **Peter B. Littlewood**  
  [littlewood@uchicago.edu](mailto:littlewood@uchicago.edu)  
  James Franck Institute and Department of Physics, University of Chicago, Chicago IL 60637, USA  
  School of Physics and Astronomy, The University of St Andrews, St Andrews, KY16 9AJ, United Kingdom


## 🙏 Acknowledgments

This research benefited from Physics Frontier Center for Living Systems funded by the National Science Foundation (PHY-2317138).

RH was supported by a Grant in Aid for Transformative Research Areas (No. 25H01364), for Scientific Research (B) (General) (No. 25K00935), and for Research Activity Start-up from JSPS in Japan (No. 23K19034), and the National Research Foundation (NRF) funded by the Ministry of Science of Korea (Grant No. RS-2023-00249900).

The computation benefited from Research Computing Center at the University of Chicago.
