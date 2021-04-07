# HRP
Code for Hierarchical Riemannian Pursuit

## 1.0 Prerequisites
+ **Matlab**
+ **SPAMS Matlab toolbox v2.6**
Download SPAMS from  http://spams-devel.gforge.inria.fr/downloads.html .
https://github.com/xhm1014/spams-matlab-install-on-win10 shows how to install it.
+ **KSVD Matlab toolbox**
Download KSVD v13 from https://www.cs.technion.ac.il/~ronrubin/software.html
(OMP-Box v10 is required)

## 2.0 Generate results for the convergence curves
Run   `TSP_gap_sim.m` in the folder `curve_convergence`.

## 3.0 Generate results for the sample complexity curves
Run  `TSP_sample_sim.m` in the folder `curve_samplecomplexity`.

## 4.0 Generate results for the heatmap with synthetic data
+ Unzip the .zip files in the folder `heatmap_synthetic`
+ Run `TSP_rev1_syndata_main.m` in the folder `heatmap_synthetic`

## 5.0 Generate results for the table with real-world sensor data
+ Copy all the toolboxes and methods from the folder `heatmap_synthetic` to the folder `method` under `table_sensor`
+ Run `TSP_Rev_Real.m` in the folder `table_sensor`
