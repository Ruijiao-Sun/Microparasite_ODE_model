# A continuous-time microparasite model incorporating infection intensity and parasite aggregation
Ruijiao Sun<sup>1</sup>, Jason Cosens Walsman<sup>1</sup>, Mark Wilber<sup>2</sup>, Cheryl J. Briggs<sup>1</sup>
1. Marine Science Institute, University of California Santa Barbara, Santa Barbara, CA, USA
2. School of Natural Resources, University of Tennessee Institute of Agriculture, Knoxville, Tennessee 37996

Please contact Ruijiao Sun (ruijiaos@ucsb.edu) for suggestions, comments, improvements, etc.

## Overview
This repository contains MATLAB code for our microparasite model that incorporates infection intensity and pathogen aggregation, with application to the amphibian chytrid fungus disease system.

## Folder
A brief description of the contents within each folder is provided below.


### PDE model
Contains MATLAB code for numerical simulations of the full partial differential equation (PDE) model.

### 1st order approximation
Includes MATLAB code for the first-order approximation, which tracks changes in mean pathogen load while assuming a fixed aggregation level.
- **equilibrium.m** performs analyses to evaluate how aggregation, load-dependent mortality, and within-host pathogen growth influence host–parasite dynamics. This code generates Figure 3 in the main text.
- **trade_off_analysis.m** analyzes virulence-transmission trade-offs under different model scenarios. This code generates Figure 4.
- **Recruit.m** - Density-dependent recruitment function.
- **load_death.m** — Load-dependent mortality function.
- **Growth.m** — Within-host pathogen growth function.
- **Bd_system.m** — Wrapper function that combines model components and is called by equilibrium.m and trade_off_analysis.m.
- **Bd_stability.m** — Stability analysis of equilibrium points.
- **Bd_R0.m** — Calculates the basic reproduction number R0.
- **Bd_model_fixedsigma.m** - Main model ODE equations
- **Bd_jacobian.m** — Jacobian matrix calculation for stability analysis.
- **Bd_equlibrium.m** — Computes equilibrium conditions.

### 2nd order approximation
Includes MATLAB code for the second-order approximation, which tracks changes in both mean pathogen load and aggregation (variance).
- **equilibrium.m** 
- **trade_off_analysis.m**
- **Recruit.m**
- **load_death.m**
- **Growth.m**
- **dist_plot.m**
- **Bd_system.m**
- **Bd_stability.m**
- **Bd_R0.m**
- **Bd_model_varysigma.m**
- **Bd_jacobian.m**
- **Bd_equlibrium.m**
