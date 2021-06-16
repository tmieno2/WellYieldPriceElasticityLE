#' ---
#' title: "Get depth to water table data from USGS"
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

#+ setup, include=FALSE
knitr::opts_chunk(
  echo = TRUE,
  cache = FALSE,
  comment = NA,
  message = FALSE,
  warning = FALSE,
  tidy = FALSE,
  cache.lazy = FALSE,
  #--- figure ---#
  dpi=400,
  fig.width=7.5,
  fig.height=5,
  out.width="750px",
  out.height="500px"
)

#/*=================================================*/
#' # Preparation
#/*=================================================*/
library(here)
library(dataRetrieval)
library(magrittr)
library(sf)
library(dplyr)
library(data.table)

CO_gwl_raw <- readNWISdata(
    stateCd="Colorado", 
    startDate = "2006-01-01", 
    endDate = "2016-12-31", 
    service = "gwlevels"
  ) %>% 
	data.table() %>% 
	.[, .(site_no, lev_dt, lev_va)] %>% 
  setnames(c("lev_dt", "lev_va"), c("date", "dwt")) 

CO_gwl_March <- CO_gwl_raw %>% 
  .[, date := as.Date(date)] %>% 
  .[, `:=`(month = month(date), year = year(date))] %>% 
  .[month %in% c(1, 2, 3),] %>% 
  .[, .(dwt  = mean(dwt)), by = .(site_no, year)]

CO_site_ls <- CO_gwl_March$site_no %>% unique()

sites_info <- readNWISsite(siteNumbers = CO_site_ls) %>% 
  dplyr::select(site_no, dec_lat_va, dec_long_va) %>% 
  #--- turn the data into an sf object ---#
  st_as_sf(coords = c("dec_long_va", "dec_lat_va")) %>% 
  #--- NAD 83 ---#
  st_set_crs(4269) %>% 
  #--- project to WGS UTM 14 ---#
  st_transform(32614) 

CO_gwl_sf <- left_join(sites_info, CO_gwl_March, by = "site_no") %>% 
	na.omit()

saveRDS(CO_gwl_sf, "./Data/CO_gwl_sf.rds")
