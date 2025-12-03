# Universal Scaling in One-Dimensional Non-Reciprocal Matter

[![License](https://img.shields.io/badge/license-BSD%202--Clause-blue.svg)](LICENSE)

> Code associated with the paper "Universal scaling in one-dimensional non-reciprocal matter" by Shuoguang Liu, Ryo Hanai, and Peter B. Littlewood ([arXiv:2503.14384](https://arxiv.org/abs/2503.14384)).


## 📖 Overview

This repository contains the code needed to reproduce the main results of the paper:

> **Universal scaling in one-dimensional non-reciprocal matter**
>
> Shuoguang Liu, Ryo Hanai, and Peter B. Littlewood, 2025.

The scripts are written in **MATLAB** and generate both simulation data and figures used in the main text and supplementary information.

To generate data, run `main_computation.m`.
To compute correlation functions, run `avg_corr_psi.m` (time domain) or `FFT.m` (frequency domain).



## 📊 Figures Correspondence

### Main Text Figures

- **Fig. 1:**
  - (a) `main_computation.m` (for insets)
  - (b-c) no code associated
- **Fig. 2:**
  - (a) `xover_chi_all.m`
  - (b) `collapse_CEP.m` and `collapse_EW.m`
- **Fig. 3:**
  - (a) `alpha_CEP.m`
  - (b) `fit_peak.m` and `z_CEP.m`
  - (c-d) same as Fig. 2(b)
- **Fig. 4:**
  - (a) `phaseA_evolve.m`
  - (b) `log_scaling.m` and `power_scaling.m`
  - (c) `xover_freq_all.m`


### Supplementary Information Figures

- **Fig. S1:** `collapse_CEP.m`
- **Fig. S2:** `waiting_time_CEP.m` and `waiting_time_EW.m`
- **Fig. S3:** `amp_pattern.m`
- **Fig. S4:** `collapse_CEP.m`
- **Fig. S5, S6:** `collapse_width_CEP.m`
- **Fig. S7:** `collapse_width_EW.m`
- **Fig. S8:** `amp_pattern.m`
- **Fig. S9:** `freq_size.m` and `freq_noise.m`
- **Fig. S10:** `main_computation.m` (for a); no code for (c, d)
- **Fig. S11:** `log_scaling.m`
- **Fig. S12:** `CEP_approx.m` and `EW_approx.m`
- **Fig. S13:** `compare_correlator_EW.m`
- **Fig. S14:** `compare_correlator_CEP.m` and `ensemble_CEP.m`
- **Fig. S15:** `compare_correlator_pattern.m` and `ensemble_pattern.m`
- **Fig. S16:** `FFT_scaling_width.m`


##  📚 Citation

If you use this code in your work, please cite:

> Shuoguang Liu, Peter B. Littlewood, Ryo Hanai, "Universal scaling in one-dimensional non-reciprocal matter," 2025, [arXiv:2503.14384](https://arxiv.org/abs/2503.14384).



## 👥 Authors

- **Shuoguang Liu**  
  [shuoguang@uchicago.edu](mailto:shuoguang@uchicago.edu)  
  James Franck Institute and Department of Physics, University of Chicago, Chicago IL 60637, USA

- **Peter B. Littlewood**  
  [littlewood@uchicago.edu](mailto:littlewood@uchicago.edu)  
  James Franck Institute and Department of Physics, University of Chicago, Chicago IL 60637, USA  
  School of Physics and Astronomy, The University of St Andrews, St Andrews, KY16 9AJ, United Kingdom
  
- **Ryo Hanai**  
  [hanai.r.7e4b@m.isct.ac.jp](mailto:hanai.r.7e4b@m.isct.ac.jp)  
  Department of Physics, Institute of Science Tokyo, 2-12-1 Ookayama Meguro-ku, Tokyo, 152-8551, Japan


##  Acknowledgments

This research benefited from Physics Frontier Center for Living Systems funded by the National Science Foundation (PHY-2317138). RH was supported by a Grant in Aid for Transformative Research Areas (No. 25H01364), for Scientific Research (B) (General) (No. 25K00935), and for Research Activity Start-up from JSPS in Japan (No. 23K19034) and the National Research Foundation (NRF) funded by the Ministry of Science of Korea (Grant No. RS-2023-00249900). The computation benefited from Research Computing Center at the University of Chicago.
