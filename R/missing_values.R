
# run all  and jump to line 75

#-----------------------------------

library(dplyr)
library(tidyr)
library(rsofun)
library(here)
library(ncdf4)

# insert correct whc and run p model (not interesting)

# paramter for p model
params_modl <- list(
  kphio              = 0.04998,
  kphio_par_a        = 0.0,
  kphio_par_b        = 1.0,
  soilm_thetastar    = 0.6 * 240,
  soilm_betao        = 0.0,
  beta_unitcostratio = 146.0,
  rd_to_vcmax        = 0.014,
  tau_acclim         = 30.0,
  kc_jmax            = 0.41
)

driver <- readRDS(here("data","rsofun_driver_data_clean.rds"))

# updating whc, insert the file whc_2m.nc in your data folder or directory of interest

nc_file <- nc_open(here("data","whc_2m.nc"))

whc = ncvar_get(nc_file, "whc_2m")
lons = ncvar_get(nc_file, "lon")
lats = ncvar_get(nc_file, "lat")

geo <- driver |>
  unnest(site_info) |>
  select(lon  , lat)

geo$sitename <- driver$sitename

n <- 1 # parameter to select size of slice to average

whc <- lapply(geo$sitename, function(x){
  tmp <- geo[geo$sitename == x,]
  lonid <- which(lons > tmp$lon)[1]
  latid <- which(lats > tmp$lat)[1]
  whc_grid <- whc[(lonid-n):(lonid+n), (latid-n):(latid+n)]
  whc_site <- mean(as.numeric(whc_grid, na.rm=T))
  return(whc_site)
})

whc = unlist(whc)

for(i in 1:241){
  driver$site_info[i][[1]][4] <- whc[i]
}
rm(whc) # for menory

# run p model

output <- rsofun::runread_pmodel_f(
  driver,
  par = params_modl
)

# some simulations failed due to missing values in the forcing. If that happens,
# the number of rows in the output data is 1.
output <- output |>
  mutate(len = purrr::map_int(data, ~nrow(.))) |>
  filter(len != 1) |>
  select(-len)

#-----------------------------------
# look at missing value

output_missing = output |>
  unnest(data)  |>
  select(aet,pet)

bool = as.data.frame(is.na(output_missing))
# line missing == 2 =  both value are missing
line_missing = apply(bool, 1, sum)
table(line_missing)

col_missing = apply(bool, 2, sum)
print(col_missing)
# 36865 estimations have both aet and pet missing
# since there are 36865 missing pet, when pet is missing also aet is missing

# this will reflect in aet missing also in the dataframe (38 sites have missing values)
# see AET_WHC comparison, in the world map the eastern region present missing data

whc_missing = output |>
  unnest(site_info) |>
  select(whc)
table(is.na(whc_missing))

# 7 whc values are missing

forcing_missing = driver |>
  unnest(forcing) |>
  # remove the nested dataframe and non numeric variable
  select(-sitename, -params_siml, -site_info, -date)

table(is.na(forcing_missing))
# 2% of the forcing data are missing

bool_forcing = as.data.frame(is.na(forcing_missing))
forcing_line_missing = apply(bool_forcing, 1, sum)
# 13 == all the value are missing in the row
table(forcing_line_missing)
# no row as all the value missing, maximum 3 values per row

# scan every column in forcing_missing
forcing_lcolum_missing = apply(bool_forcing, 2, sum)
print(forcing_lcolum_missing)

# the net radiation seems to be the problem, having 146730 values missing
# gpp, which is the second-most column whit missing values, has 6205 of them
