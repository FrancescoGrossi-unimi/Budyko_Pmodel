library(tidyverse)
library(reshape2)
library(lubridate)
library(ncdf4)
library(here)
library(cwd)
source(here("R","read_meta_fdk.R"))

files_csv = list.files(here("data-raw/CSV"))
files_lsm = list.files(here("data-raw/LSM"))

file_csv = files_csv[grep("HH", files_csv)]
sites <- lapply(file_csv,function(x) {substr(x,5,10)})


params_siml = list( tibble(spinup = TRUE, spinupyears = 10,recycle=1,outdt= 1,
                           ltre=FALSE, ltne = FALSE, ltrd=FALSE, ltnd= FALSE, lgr3= TRUE,
                           lgn3= FALSE,lgr4= FALSE))

nc = nc_open(here("data","cwdx80.nc"))
lons = ncvar_get(nc, "lon")
lats = ncvar_get(nc, "lat")
S80 = ncvar_get(nc, "cwdx80")

# for cycle
df <- NULL

for(i in 1:length(sites)){

  site <- sites[[i]]
  csv <- file_csv[i]

  # acquire metadata

  meta <- suppressWarnings(
    try(
      read_meta_fdk(
        site = site,
        path = here("data-raw/LSM"),
        meta_data = T
      )
    )
  )

  longitude  <- meta[[1]]$longitude
  latitude  <- meta[[1]]$latitude
  elevation  <- meta[[1]]$elevation

  hhdf <- readr::read_csv(paste0(here("data-raw/CSV"),"/",csv))

  # if LE_CORR not present skip the site

  if (!("LE_CORR" %in% colnames(hhdf))) next

  # Add date and time columns to hhdf for easier further processing.
  # ---------------------------------------------------------
  hhdf <- hhdf |>
    mutate(time = lubridate::as_datetime(as.character(TIMESTAMP_START), tz = "GMT", format="%Y%m%d%H%M")) |>
    mutate(date = lubridate::as_date(time))

  #add NERTRAD and other columns if absent

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

  # Creating driver object  ---------------

  lonid = which(lons > longitude)[1]
  latid = which(lats > latitude)[1]
  n = 1 # parameter to select the slice to average
  S80_slice = S80[(lonid-n):(lonid+n), (latid-n):(latid+n)]
  whc_site = mean(as.numeric(S80_slice, na.rm=T))

  site_info = list(tibble(
    lon= longitude,
    lat= latitude,
    elv =  elevation,
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
