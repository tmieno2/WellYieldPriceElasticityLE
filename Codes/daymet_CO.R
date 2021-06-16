######################################
# Download and assemble Daymet data
######################################

#===================================
# Preparation
#===================================
library(daymetr)
library(readr)
library(here)

setwd(here())

#--------------------------
# get well location data
#--------------------------
wells_loc <- readRDS('./Data/wells_loc_sf.rds') %>%
  cbind(.,st_coordinates(.)) %>%
  setnames(c('X','Y'),c('long','lat')) %>%
  data.table() 

#===================================
# Download the data by well
#===================================
wells_len <- nrow(wells_loc)
# i <- 1

daymet_download <- function(i){
	print(i)
	well_temp <- wells_loc[i,]
	temp_wdid <- well_temp[,wdid]
	temp_data <- download_daymet(
		  site = temp_wdid,
      lat = well_temp[,lat],
      lon = well_temp[,long],
      start = 2011,
      end = 2017,
      internal = TRUE
    )
  temp_daymet_data <- temp_data$data %>%
    data.table() %>%
	  .[,wdid:=temp_wdid] %>%
    .[,altitude:=temp_data$altitude]

	return(temp_daymet_data)
}

# daymet_all <- mclapply(1:wells_len,daymet_download,mc.cores=8) %>%
# 	rbindlist()

daymet_all <- lapply(1:wells_len,daymet_download) %>%
  rbindlist()

saveRDS(daymet_all,'./Data/dyamet_raw.rds')

#===================================
# Put together data
#===================================
daymet_all <- readRDS('./Data/dyamet_raw.rds')

#--- source functions ---#
source('./Codes/ASCE_PenmanMonteith.R')
MeanWspd <- 3.25

daymet_sum <- daymet_all %>%
	setnames(
    c("year","yday","dayl..s.","prcp..mm.day.","srad..W.m.2.","swe..kg.m.2.","tmax..deg.c.","tmin..deg.c.","vp..Pa."),
    c("year","year.day","daylight","precip","solar.rad","snow.water.equiv","max.temp","min.temp","vapor.pressure")
  ) %>%
  .[,`:=`(
    daily_gdd = (max.temp + min.temp)/2 - 10,
    precip = precip*0.0393701,
    solar.rad = (solar.rad*daylight)/1000000,  # Converted to match Tim's ET process form matlab
    vapor.pressure = vapor.pressure/1000,  # Converted to match Tim's ET process form matlab
    Wspd = MeanWspd,
    date = as.Date(year.day - 1, origin = paste(year,"-01-01",sep=""))
    )] %>%
  .[,`:=`(
    day = as.numeric(substring(date, 9, 10)),
    month = as.numeric(substring(date, 6, 7))
	)] %>%
	wells_loc[.,on='wdid'] %>%
  .[,ET0 := ASCE_PenmanMonteith(
		max.temp,min.temp,solar.rad,vapor.pressure,
		Wspd,precip,year,month,day,lat,altitude
	)*0.0393701]

#--------------------------
# summarize the data
#--------------------------
daymet_season <- daymet_sum[month %in% c(5,6,7,8,9),.(precip = sum(precip), et = sum(ET0), gdd = sum(daily_gdd)),by=.(wdid,year)]

saveRDS(daymet_season,'./Data/daymet_season_CO.rds')

