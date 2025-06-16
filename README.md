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

### 2nd order approximation
Includes MATLAB code for the second-order approximation, which tracks changes in both mean pathogen load and aggregation (variance).
