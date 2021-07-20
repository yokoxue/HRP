# HRP
Code for paper "Efficient Sparse Coding using Hierarchical Riemannian Pursuit," in IEEE Transactions on Signal Processing, Y. Xue, V. K. N. Lau and S. Cai, doi: 10.1109/TSP.2021.3093769.[[paper]](https://arxiv.org/pdf/2104.10314.pdf)

If you find they are useful, please cite:
```
@ARTICLE{9470930,
  author={Xue, Ye and Lau, Vincent K. N. and Cai, Songfu},
  journal={IEEE Transactions on Signal Processing}, 
  title={Efficient Sparse Coding using Hierarchical Riemannian Pursuit}, 
  year={2021},
  volume={},
  number={},
  pages={1-1},
  doi={10.1109/TSP.2021.3093769}}
```
## 1.0 Prerequisites
+ **Matlab**

+ **KSVD Matlab toolbox （for Baseline 1)**

Download KSVD v13 from https://www.cs.technion.ac.il/~ronrubin/software.html and install
(OMP-Box v10 is required).

+ **SPAMS Matlab toolbox v2.6 （for Baseline 2)**

Download SPAMS from  http://spams-devel.gforge.inria.fr/downloads.html.
Follow the steps in https://github.com/xhm1014/spams-matlab-install-on-win10 to install.


+ **CVX Matlab toolbox （for Baseline 4)**

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
