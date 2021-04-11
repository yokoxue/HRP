# HRP
Code for paper "Efficient Sparse Coding using Hierarchical Riemannian Pursuit" Ye Xue, Vincent Lau, and Songfu Cai (submited)

## 1.0 Prerequisites
+ **Matlab**
+ **SPAMS Matlab toolbox v2.6**
Download SPAMS from  http://spams-devel.gforge.inria.fr/downloads.html .
https://github.com/xhm1014/spams-matlab-install-on-win10 shows how to install it.
+ **KSVD Matlab toolbox**
Download KSVD v13 from https://www.cs.technion.ac.il/~ronrubin/software.html
(OMP-Box v10 is required).
+ **CVX Matlab toolbox**
Download CVX toobox from http://cvxr.com/cvx/ and install.

## 2.0 Generate the results for the convergence curves
Run   `Converge_sim.m` in the folder `curve_convergence`.

## 3.0 Generate the results for the sample complexity curves
Run  `Sample_sim.m` in the folder `curve_samplecomplexity`.

## 4.0 Generate the results for the RMSE heatmap with synthetic data
+ Unzip the .zip files in the folder `heatmap_synthetic`.
+ Run `Syndata_main.m` in the folder `heatmap_synthetic`.

## 5.0 Generate the results for the table with real-world sensor data
+ Unzip all the .zip files in the folder `table_sensor`.
+ Run `Sensor_Data_main.m` in the folder `table_sensor`.
+ Raw data of the sensor readings of the Airly network can be downloaded from https://www.kaggle.com/datascienceairly/air-quality-data-from-extensive-network-of-sensors.
