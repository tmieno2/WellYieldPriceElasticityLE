#' ---
#' title: "Data management for CO"
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
#=== wd ===#
setwd('~/Dropbox/GroundwaterWithMani/WellYieldPriceElasticity/LandEconVersion')
#=== packages ===#
library(data.table)
library(magrittr)
library(sf)
library(stargazer)
library(lfe)

#/*=================================================*/
#' # Initial Data Preparation
#/*=================================================*/

#=== water use, well capacity, pumping cost efficiency, price ===#
#' PricePanelData.csv was provided by Aaron Hrozencik.
data <- fread('./Data/PriceData.csv') %>%
  .[,.(wdid,Year,PCC,TrueMP,TrueAP,WellCap,pumpingAF,REAname,numPCC,p1,p2,p3,pump_hp)] %>%
  setnames('Year','year')

#=== well location ===#
wells_location <- st_read("./Data/FinalPermit.shp") %>%
  setnames(names(.),tolower(names(.))) %>%
  data.table() %>%
  .[,.(wdid,latdecdeg,longdecdeg)] %>%
  .[,wdid:=as.numeric(as.character(wdid))] %>%
  .[!is.na(wdid),]

#=== combine ===#
data_comp <- wells_location[data,on='wdid']

#=== sf version ===#
wells_loc_sf <- data_comp %>%
  .[,.(wdid,latdecdeg,longdecdeg)] %>%
  unique(by='wdid') %>%
  .[!is.na(latdecdeg),] %>%
  .[!is.na(longdecdeg),] %>%
  st_as_sf(coords=c('longdecdeg','latdecdeg')) %>%
  st_set_crs(4269)

#===  save ===#
saveRDS(data_comp,'./Data/data_CO.rds')
saveRDS(wells_loc_sf,'./Data/wells_loc_sf.rds')

#/*=================================================*/
#' # Combine datasets
#/*=================================================*/
#/*----------------------------------*/
#' ## Import datasets 
#/*----------------------------------*/
#=== water use, well capacity, pumping cost efficiency, price ===#
data <- readRDS('./Data/data_CO.rds')

#=== daymet ===#
daymet_season <- readRDS('./Data/daymet_season_CO.rds')

#=== saturated thickness ===#
all_sat <- readRDS('./Data/all_sat_CO.rds') %>%
  .[,year:=as.numeric(year)] %>%
  .[year>=2011,]

#/*----------------------------------*/
#' ## Combine all
#/*----------------------------------*/
data_reg <- daymet_season[data,on=c('wdid','year')] %>%
  #=== at least three pumping tests ===#
  .[numPCC>=3,] %>%
  #=== sat thickness available only until 2016 ===#
  .[year<=2016,] %>%
  #=== merge with sat ===#
  all_sat[.,on=c('wdid','year')] 

#/*----------------------------------*/
#' ## Merge with USGS dwt data
#/*----------------------------------*/
#--- import CO dwt data ---#
CO_gwl_sf <- readRDS("./Data/CO_gwl_sf.rds") 
CO_gwl_dt <- data.table(CO_gwl_sf) %>% 
  .[, geometry := NULL]

#--- transform ir wells to have the same CRS as the dwt data ---#
wells_loc_sf <- st_transform(wells_loc_sf, st_crs(CO_gwl_sf))

#--- keep only the ones that are inside the boundary of ir wells ---#
CO_gwl_wells <- distinct(CO_gwl_sf, geometry, .keep_all = TRUE) %>% 
  st_crop(st_bbox(wells_loc_sf))

# ggplot() +
#   geom_sf(data = wells_loc_sf, col = "red") +
#   geom_sf(data = CO_gwl_sf, col = "blue") 

#--- find ir wells that have at least one USGS well within its 2 mile radius ---#
ir_wells_in_buffer_id <- st_buffer(CO_gwl_wells, dist = 3200) %>% 
  st_intersection(wells_loc_sf, .) %>% 
  pull(wdid)

ir_wells_in_buffer <- filter(wells_loc_sf, wdid %in% ir_wells_in_buffer_id)

#--- find the nearest USGS wells for each of the ir wells ---#
nearest_index <- st_nearest_feature(ir_wells_in_buffer, CO_gwl_wells)

#--- ir-usge correspondence ---#
ir_usgs <- ir_wells_in_buffer %>% 
  data.table() %>% 
  .[, site_no := CO_gwl_wells[nearest_index, ] %>%  pull(site_no)] %>% 
  .[, geometry := NULL]

saveRDS(ir_usgs, "./Data/ir_usgs.rds")

#--- merge the data ---#
data_reg <- ir_usgs[data_reg, on = "wdid"] %>% 
  .[,usgs_near := TRUE] %>% 
  .[is.na(site_no), usgs_near := FALSE] %>% 
  CO_gwl_dt[., on = c("site_no", "year")]

# data_reg[!is.na(dwt), ]
# data_reg[!is.na(site_no), ]

#/*=================================================*/
#' # Other operations
#/*=================================================*/

data_reg <- setnames(data_reg, c('WellCap','PCC'),c('gpm','etw')) %>%
    #=== drop years prior 2013 for Highline as it has a different rate structure ===#
  .[,REA:='Y-W'] %>%
  .[REAname=='Highline' & year<=2012,REA:='Highline-1'] %>%
  .[REAname=='Highline' & year > 2012,REA:='Highline-2'] %>%
  # .[!(year < 2013 & REAname=='Highline'),] %>%
  .[etw < 1200,] %>%
  .[,pc:=etw/12*TrueMP] %>%
  .[,tier:=1] %>%
  .[TrueMP==p2,tier:=2] %>%
  .[TrueMP==p3,tier:=3] %>%
  .[,wdid_tier_hp:=paste0(wdid,tier,pump_hp)] %>%
  .[,wdid_tier:=paste0(wdid,tier)] %>%
  .[,year_REA:= paste0(year,REAname)] %>%
  #=== within-transformed gpm and etw ===#
  .[,`:=`(
    w_gpm=gpm-mean(gpm),
    w_etw=etw-mean(etw),
    w_sat=sat_thickness-mean(sat_thickness),
    w_pc=pc-mean(pc),
    w_et=et-mean(et),
    w_precip=precip-mean(precip)
    ),
    by=wdid_tier_hp] %>%
  #=== get rid of outliers ===#
  .[!(pumpingAF>300 & gpm<500),] %>%
  .[!(REA=='Y-W' & abs(w_pc)>2.5),] %>%
  .[!(REA=='Highline' & abs(w_pc)>1.5),] %>%
  # .[,any_gpm_over_2000:=sum(gpm>2000),by=wdid] %>%
  # .[any_gpm_over_2000==0,] %>%
  .[,any_pump_0:=sum(pumpingAF<10),by=wdid] %>%
  .[any_pump_0==0,]
  # %>%
  # .[!(pumpingAF>760 & year==2012),]

saveRDS(data_reg,'./Data/final_data_CO.rds')
