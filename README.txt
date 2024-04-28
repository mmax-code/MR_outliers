This repository contains material to reproduce the results of the manuscript "Outlier detection in Mendelian Randomisaton" by Maximilian Mandl, Anne-Laure Boulesteix, Stephen Burgess, and Verena Zuber. 
 
Session Information:

R version 4.0.4 (2021-02-15)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Debian GNU/Linux 11 (bullseye)

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0

locale:
 [1] LC_CTYPE=de_DE.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=de_DE.UTF-8        LC_COLLATE=de_DE.UTF-8    
 [5] LC_MONETARY=de_DE.UTF-8    LC_MESSAGES=de_DE.UTF-8   
 [7] LC_PAPER=de_DE.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=de_DE.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets 
[7] methods   base     

other attached packages:
 [1] MendelianRandomization_0.9.0 RadialMR_1.1                
 [3] scales_1.1.1                 RcppEigen_0.3.3.9.1         
 [5] doParallel_1.0.16            iterators_1.0.13            
 [7] foreach_1.5.1                dplyr_1.0.8                 
 [9] caret_6.0-86                 ggplot2_3.3.3               
[11] lattice_0.20-41              mvtnorm_1.1-1               

     

Structure:

- Univariable:
    - Real Data
        - Univariable_Data.R
        - VitD.RData 
    - Simulation
        - Univariable_Simulation.R 

- Multivariable:
    - Real Data
        - Multivariable_Data.R
        - chd_amd_alz.csv
    - Simulation
        - Multivariable_Simulation.R

        
Description:

- Univariable_Data.R: File for running the univariable MR data part. 
- VitD.RData: Summary-level data for MS-VitD association.
- Univariable_Simulation.R: Runs the simulation analysis for the univariable settings. 
- Multivariable_Data.R: File for running the multivariable MR data part. 
- chd_amd_alz.csv: Summary-level data for chd,amd,alz and ldl,hdl,triglycerids associations. 
- Multivariable_Simulation.R:  Runs the simulation analysis for the multivariable settings. 
