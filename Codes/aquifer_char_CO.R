#' ---
#' title: "Identify aquifer characteristics from USGS"
#' author: "Taro Mieno"
#' output:
#'   html_document:
#'     number_sections: yes
#'     theme: flatly
#'     highlight: zenburn
#'     toc_float: yes
#'     toc: yes
#'     toc_depth: 3
#' geometry: margin=1in
#' ---


#/*=================================================*/
#' # Preparation
#/*=================================================*/
setwd('~/Dropbox/GroundwaterWithMani/WellYieldPriceElasticity')

#=== library ===#
source('~/Box Sync/R_libraries_load/library.R')
library(lfe)
library(raster)
library(sf)
library(dplyr)
library(stringr)
library(parallel)

#/*----------------------------------*/
#' ## Wells location
#/*----------------------------------*/
wells_location <- readRDS('./Data/CO/wells_loc_sf.rds') 

#/*=================================================*/
#' # Gepgraphic matching
#/*=================================================*/
# --------------------------
# Identify saturated thickness at each well
#--------------------------
tif_file_names <- list.files(path="~/Dropbox/GroundwaterWithMani/Data/HPA_sat_thickness_inputs" , pattern="*.tif", full.names=TRUE) %>%
	.[!str_detect(.,'tif.')]
tif_len <- length(tif_file_names)

crs_to_use <- raster(tif_file_names[1])@crs@projargs

wells_for_sat <- wells_location %>%
	st_transform(crs_to_use) %>%
	as('Spatial')

get_sat <- function (i){
	print(i)
	#--- file name ---#
	temp_filename <- tif_file_names[i]

	#--- import the data ---#
	temp_raster <- raster(temp_filename)

	#--- extract and form a data.table ---#
	temp_data <- data.table(
		wdid = wells_for_sat@data$wdid,
		sat_thickness = raster::extract(temp_raster,wells_for_sat),
		year = str_extract(temp_filename,"[0-9]{4}")
		)
	return(temp_data)
}

all_sat <- mclapply(1:tif_len,get_sat,mc.cores=6) %>%
	rbindlist() %>%
	.[,sat_thickness:=sat_thickness*3.28084] %>%  # meter to feet
	.[,unit:='feet']

#--- save the data ---#
saveRDS(all_sat,'./Data/CO/all_sat_CO.rds')

# all_sat[sat_thickness<0,sat_thickness] %>% hist
# nrow(all_sat)
