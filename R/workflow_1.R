# this script is the merge of create_driver and create_file for vignette
# the final output will be the two plots that compare the observed and predicted ET with different whc

library(tidyverse)
library(reshape2)
library(lubridate)
library(ncdf4)
library(here)
library(cwd)

# driver creation, skip if already present
#-----------

# acquire metadata, function taken from jaideep repo rsofun, "read_meta_fdk.R"

files_csv = list.files(here("data-raw/CSV"))
files_lsm = list.files(here("data-raw/LSM"))

files_lsm <- files_lsm[grep("Flux.nc", files_lsm)]
files_csv = files_csv[grep("HH", files_csv)]
sites <- lapply(file_csv,function(x) {substr(x,5,10)})
check_order <- lapply(files_lsm,function(x) {substr(x,1,6)}) # check if order is the same

identical(sites,check_order)

params_siml = list( tibble(spinup = TRUE, spinupyears = 10,recycle=1,outdt= 1,
                           ltre=FALSE, ltne = FALSE, ltrd=FALSE, ltnd= FALSE, lgr3= TRUE,
                           lgn3= FALSE,lgr4= FALSE))

nc = nc_open(here("data","cwdx80.nc"))
lons = ncvar_get(nc, "lon")
lats = ncvar_get(nc, "lat")
S80 = ncvar_get(nc, "cwdx80")

# for cycle, select the column of interest and take the daily average length(sites)
df <- NULL

for(i in 1:3){

  site <- sites[[i]]
  csv <- files_csv[i]
  lsm <- files_lsm[i]

  # acquire metadata

  meta <- nc_open(here("data-raw/LSM",lsm))
  lat = as.numeric(ncvar_get(meta,"latitude"))
  long = as.numeric(ncvar_get(meta,"longitude"))
  elv = as.numeric(ncvar_get(meta,"elevation"))

  hhdf <- readr::read_csv(paste0(here("data-raw/CSV"),"/",csv))

  # if LE_CORR not present skip the site

  if (!("LE_CORR" %in% colnames(hhdf))) next

  # Add date and time columns to hhdf for easier further processing.

  hhdf <- hhdf |>
    mutate(time = lubridate::as_datetime(as.character(TIMESTAMP_START), tz = "GMT", format="%Y%m%d%H%M")) |>
    mutate(date = lubridate::as_date(time))

  #add NERTRAD and other columns if absent, otherwise p model won't work
  # fiter out also this (?)

  if (!("NETRAD" %in% colnames(hhdf))) {
    hhdf$NETRAD = NA
  }

  if (!("FPAR" %in% colnames(hhdf))) {
    hhdf$FPAR = NA
  }

  if (!("GPP_DT_VUT_REF" %in% colnames(hhdf))) {
    hhdf$GPP_DT_VUT_REF = NA
  }

  # Aggregate to daily 24-hr means  -----
  ddf_24hr_mean <-
    try(
      hhdf |>
        group_by(date) |>
        select(time, date, TA_F_MDS, VPD_F_MDS, SW_IN_F_MDS, NETRAD, PA_F, P_F,
               CO2_F_MDS, GPP_DT_VUT_REF,FPAR,LE_CORR) |>
        summarize_all(.funs = mean, na.rm = TRUE)
    ) |>
    mutate(LE_CORR = convert_et(LE_CORR,TA_F_MDS,PA_F*1000)) #pressure is in kPa

  # tmax and tmin
  tmaxmin <-
    hhdf |>
    group_by(date) |>
    summarize(
      tmax = max(TA_F_MDS),
      tmin = min(TA_F_MDS)
    )

  # Creating driver object

  lonid = which(lons > long)[1]
  latid = which(lats > lat)[1]
  n = 1 # parameter to select the slice to average
  S80_slice = S80[(lonid-n):(lonid+n), (latid-n):(latid+n)]
  whc_site = mean(as.numeric(S80_slice, na.rm=T))

  site_info = list(tibble(
    lon= long,
    lat= lat,
    elv =  elv,
    whc = whc_site))

  forcing = ddf_24hr_mean |>
    dplyr::filter(!(lubridate::mday(date) == 29 & lubridate::month(date) == 2)) |>
    left_join(tmaxmin) |>
    group_by(date) |>
    summarize(
      date = date,
      temp = TA_F_MDS,
      vpd = VPD_F_MDS * 100,
      ppfd = SW_IN_F_MDS * 2.04 * 1e-06, # 2.04 is kfFEC
      netrad = NETRAD,
      patm = PA_F * 1000,
      snow = 0,
      rain = P_F * 48 /(60 * 60 * 24), # P_F [mm timestep-1] * 48 [timesteps day-1] / 86400 [secs day-1 ]
      tmin = tmin, # TMIN_F_MDS,
      tmax = tmax, # TMAX_F_MDS,
      fapar = FPAR,
      co2 = CO2_F_MDS,
      ccov = 0,
      gpp = GPP_DT_VUT_REF,
      ET = LE_CORR
    ) |>
    mutate(gpp = gpp*86400/1e6*12)  |>    # convert [umol m-2 s-1] to [gC m-2 day-1]
    list()

  tmp <- data.frame(sitename = site)

  tmp$params_siml = params_siml
  tmp$site_info = site_info
  tmp$forcing = forcing

  df <- rbind(df,tmp)
}

saveRDS(df, here("data","new_driver_data_obs.rds"), compress="xz")
#-----------

# run p model with different whc, skip if already present

#-----------


# function used to create dataframes
get_annual_aet_pet <- function(df){
  df |>
    mutate(year = lubridate::year(date)) |>
    group_by(year) |>
    summarise(aet = sum(aet),
              pet = sum(pet)) |>
    ungroup() |>
    summarise(aet = mean(aet),
              pet = mean(pet))
}

get_annual_prec_cond <- function(df){
  df |>
    mutate(year = lubridate::year(date)) |>
    group_by(year) |>
    summarise(prec_cond = sum(prec_cond)) |>
    ungroup() |>
    summarise(prec_cond = mean(prec_cond))
}

# some sites have missing value within the year
# wheighted average to ensure correct ET
get_annual_ET <-  function(df){
  df |>

    mutate(year = lubridate::year(date)) |>
    group_by(year) |>
    summarise(valid_days   = sum(!is.na(ET)),
              ET   = sum(ET,na.rm=TRUE)) |>
    ungroup() |>
    summarise(ET = sum(ET*valid_days/365)/(sum(valid_days/365)))
}

# paramter for p model
params_modl <- list(
  kphio              = 0.04998,    # setup ORG in Stocker et al. 2020 GMD
  kphio_par_a        = 0.0,        # set to zero to disable temperature-dependence of kphio
  kphio_par_b        = 1.0,
  soilm_thetastar    = 0.6 * 240,  # to recover old setup with soil moisture stress
  soilm_betao        = 0.0,
  beta_unitcostratio = 146.0,
  rd_to_vcmax        = 0.014,      # value from Atkin et al. 2015 for C3 herbaceous
  tau_acclim         = 30.0,
  kc_jmax            = 0.41
)

old_driver <- new_driver <- readRDS(here("data","new_driver_data_obs.rds"))

# change whc to previous result

nc_file <- nc_open(here("data","whc_2m.nc"))

whc = ncvar_get(nc_file, "whc_2m")
lons = ncvar_get(nc_file, "lon")
lats = ncvar_get(nc_file, "lat")

geo <- old_driver |>
  unnest(site_info) |>
  select(lon  , lat)

geo$sitename <- old_driver$sitename

n <- 1 # parameter to select size of slice to average

old_whc <- lapply(geo$sitename, function(x){
  tmp <- geo[geo$sitename == x,]
  lonid <- which(lons > tmp$lon)[1]
  latid <- which(lats > tmp$lat)[1]
  whc_grid <- whc[(lonid-n):(lonid+n), (latid-n):(latid+n)]
  whc_site <- mean(as.numeric(whc_grid, na.rm=T))
  return(whc_site)
})

old_whc = unlist(old_whc)

for(i in 1:dim(old_driver)[1]){
  old_driver$site_info[i][[1]][4] <- old_whc[i]
}
rm(whc) # for memory

# run p model

old_output <- rsofun::runread_pmodel_f(
  old_driver,
  par = params_modl
)

# some simulations failed due to missing values in the forcing. If that happens,
# the number of rows in the output data is 1.
old_output <- old_output |>
  mutate(len = purrr::map_int(data, ~nrow(.))) |>
  filter(len != 1) |>
  select(-len)

# create dataframe
old_adf <- old_output |>
  mutate(old_adf = purrr::map(data, ~get_annual_aet_pet(.))) |>
  unnest(old_adf) |>
  select(sitename, aet,pet)

old_adf <- old_driver |>
  unnest(forcing) |>
  left_join(
    old_output |>
      unnest(data) |>
      select(-snow, -netrad, -fapar, gpp_pmodel = gpp),
    by = c("sitename", "date")
  ) |>
  mutate(prec = (rain + snow) * 60 * 60 * 24) |>
  mutate(prec_cond = prec + cond) |>
  group_by(sitename) |>
  nest() |>
  mutate(old_adf = purrr::map(data, ~get_annual_prec_cond(.))) |>
  unnest(old_adf) |>
  select(sitename, prec_cond) |>
  right_join(
    old_adf,
    by = "sitename"
  )

whc <-  old_output |>
  unnest(site_info) |>
  select(whc)

old_adf$whc <- whc[[1]]

old_adf <- old_driver |>
  mutate(ET = purrr::map(forcing, ~get_annual_ET(.))) |>
  unnest(ET) |>
  select(sitename, ET) |>
  right_join(
    old_adf,
    by = "sitename"
  )

# run p model
new_output <- rsofun::runread_pmodel_f(
  new_driver,
  par = params_modl
)

# some simulations failed due to missing values in the forcing. If that happens,
# the number of rows in the output data is 1.
new_output <- new_output |>
  mutate(len = purrr::map_int(data, ~nrow(.))) |>
  filter(len != 1) |>
  select(-len)

# create dataframe
new_adf <- new_output |>
  mutate(new_adf = purrr::map(data, ~get_annual_aet_pet(.))) |>
  unnest(new_adf) |>
  select(sitename, aet, pet)

new_adf <- new_driver |>
  unnest(forcing) |>
  left_join(
    new_output |>
      unnest(data) |>
      select(-snow, -netrad, -fapar, gpp_pmodel = gpp),
    by = c("sitename", "date")
  ) |>
  mutate(prec = (rain + snow) * 60 * 60 * 24) |>
  mutate(prec_cond = prec + cond) |>
  group_by(sitename) |>
  nest() |>
  mutate(new_adf = purrr::map(data, ~get_annual_prec_cond(.))) |>
  unnest(new_adf) |>
  select(sitename, prec_cond) |>
  right_join(
    new_adf,
    by = "sitename"
  )

whc <-  new_output |>
  unnest(site_info) |>
  select(whc)

new_adf$whc <- whc[[1]]

new_adf <- new_driver |>
  mutate(ET = purrr::map(forcing, ~get_annual_ET(.))) |>
  unnest(ET) |>
  select(sitename, ET) |>
  right_join(
    new_adf,
    by = "sitename"
  )

saveRDS(old_adf,here("data","budyko_old_whc_ET.rds"))
saveRDS(new_adf,here("data","budyko_new_whc_ET.rds"))
saveRDS(geo,here("data","site_coordinates.rds"))
#-----------

# data visualization

old_adf <- readRDS(here("data","budyko_old_whc_ET.rds"))
new_adf <- readRDS(here("data","budyko_new_whc_ET.rds"))
geo <- readRDS(here("data","site_coordinates.rds"))

ggplot(old_adf,aes(x= ET,y=aet)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  labs(title = "Evapotranspiration vs predicted AET with new WHC",
       subtitle = paste0("R^2 = ",cor(old_adf$ET,old_adf$aet, use = "complete.obs")^2))

ggplot(new_adf,aes(x= ET,y=aet)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  labs(title = "Evapotranspiration vs predicted AET with old WHC",
       subtitle = paste0("R^2 = ",cor(new_adf$ET,new_adf$aet, use = "complete.obs")^2))
