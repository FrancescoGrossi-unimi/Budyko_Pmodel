# AET and WHC comparison

The vignette [AET_WHC_comparison.Rmd](https://github.com/FrancescoGrossi-unimi/Budyko_Pmodel/blob/main/vignettes/AET_WHC_comparison.Rmd) highighits the difference between the p_model run with the water holding capacity (WHC) estimated with the water root zone kept at 2m and the new estimates of WHC developed by [Stocker, 2023](https://www.nature.com/articles/s41561-023-01125-2).

Is possible to run directly the vignette (driver data not included at the moment).

To run all the analysis from scratch is necessary to download and insert in the folder [data](https://github.com/FrancescoGrossi-unimi/Budyko_Pmodel/tree/main/data) the data from [Zenodo](https://zenodo.org/records/8403081) and insert in two folder named CSV and LSM according to the data type, and the file [cwdx80.nc](https://zenodo.org/records/5515246) to update the WHC. It is also need the whc at 2m (ask Beni)

After downoladed the file, is possible to run the script [create_driver.R](https://github.com/FrancescoGrossi-unimi/Budyko_Pmodel/blob/main/R/create_driver.R). This script will create a driver data that will be used to run the p_model in the further script.

Is possible to run the p_model using the driver data created before, the script [create_file_for_vignette.R](https://github.com/FrancescoGrossi-unimi/Budyko_Pmodel/blob/main/R/create_file_for_vignette.R) will create two dataframes that will be compared in the vignette.
