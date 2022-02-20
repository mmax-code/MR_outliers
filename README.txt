 This repository contains material to reproduce the results of the manuscript "Outlier detection in Mendelian Randomisaton" by Maximilian Mandl, Anne-Laure Boulesteix, Stephen Burgess, and Verena Zuber. 
 
Session Information:

R version 4.0.4 (2021-02-15)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Debian GNU/Linux 11 (bullseye)

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0

Random number generation:
 RNG:     Mersenne-Twister 
 Normal:  Inversion 
 Sample:  Rounding 
 
locale:
 [1] LC_CTYPE=de_DE.UTF-8       LC_NUMERIC=C               LC_TIME=de_DE.UTF-8       
 [4] LC_COLLATE=de_DE.UTF-8     LC_MONETARY=de_DE.UTF-8    LC_MESSAGES=de_DE.UTF-8   
 [7] LC_PAPER=de_DE.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
[10] LC_TELEPHONE=C             LC_MEASUREMENT=de_DE.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] scales_1.1.1        RcppEigen_0.3.3.9.1 doParallel_1.0.16   iterators_1.0.13   
 [5] foreach_1.5.1       MRPRESSO_1.0        dplyr_1.0.4         caret_6.0-86       
 [9] ggplot2_3.3.3       lattice_0.20-41     mvtnorm_1.1-1      

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.6           lubridate_1.7.9.2    prettyunits_1.1.1    class_7.3-18        
 [5] ps_1.5.0             assertthat_0.2.1     rprojroot_2.0.2      ipred_0.9-9         
 [9] R6_2.5.0             plyr_1.8.6           stats4_4.0.4         pillar_1.4.7        
[13] rlang_0.4.10         data.table_1.13.6    callr_3.5.1          rpart_4.1-15        
[17] Matrix_1.3-2         desc_1.2.0           devtools_2.3.2       splines_4.0.4       
[21] gower_0.2.2          stringr_1.4.0        munsell_0.5.0        compiler_4.0.4      
[25] pkgconfig_2.0.3      pkgbuild_1.2.0       Metrics_0.1.4        nnet_7.3-15         
[29] tidyselect_1.1.0     tibble_3.0.6         prodlim_2019.11.13   codetools_0.2-18    
[33] crayon_1.4.0         withr_2.4.1          MASS_7.3-53.1        recipes_0.1.15      
[37] ModelMetrics_1.2.2.2 grid_4.0.4           nlme_3.1-152         gtable_0.3.0        
[41] lifecycle_0.2.0      DBI_1.1.1            magrittr_2.0.1       pROC_1.17.0.1       
[45] cli_2.3.0            stringi_1.5.3        cachem_1.0.3         reshape2_1.4.4      
[49] fs_1.5.0             remotes_2.2.0        testthat_3.0.1       timeDate_3043.102   
[53] ellipsis_0.3.1       generics_0.1.0       vctrs_0.3.6          lava_1.6.8.1        
[57] tools_4.0.4          glue_1.4.2           purrr_0.3.4          processx_3.4.5      
[61] pkgload_1.1.0        fastmap_1.1.0        survival_3.2-7       colorspace_2.0-0    
[65] sessioninfo_1.1.1    memoise_2.0.0        usethis_2.0.0    


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
- Univariable_Simulation.R: Runs the simulation analysis for the univariable setting. 
- Multivariable_Data.R: File for running the multivariable MR data part. 
- chd_amd_alz.csv: Summary-level data for chd,amd,alz and ldl,hdl,triglycerids associations. 
- Multivariable_Simulation.R:  Runs the simulation analysis for the multivariable setting. 


For the reproduction of the univariable results, please adjust the working directory and run Univariable_Data.R (for the real data application) and Univariable_Simulation.R (for the simulation results). For the reproduction of the multivariable results, please adjust the working directory and run Multivariable_Data.R (for the real data application) and Multivariable_Simulation.R (for the simulation results).


