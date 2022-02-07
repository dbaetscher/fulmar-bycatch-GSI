#Load R libraries
library(tidyverse)
library(tidyr)
library(DBI)
library(odbc)
library(lubridate)

#  Connect to AKFIN database
# con <- dbConnect(odbc::odbc(), "akfin", 
#                  UID=rstudioapi::askForPassword("Enter AKFIN Username"), 
#                  PWD= rstudioapi::askForPassword("Enter AKFIN Password"))
# 
# #  Query the comprehensive_norpac table in the Council scheme in the AKFIN database.
# #  gear_type=8 is longline
# data <- dbFetch(dbSendQuery(con,"select cruise, 
#                                   permit, 
#                                   haul_seq, 
#                                   haul_date,
#                                   duration_in_min as duration, 
#                                   official_total_catch as otc, 
#                                   mm_percent_monitored, 
#                                   bird_deterrence, 
#                                   percent_retained,
#                                   extrapolated_number,
#                                   extrapolated_weight,
#                                   trip_target_code,
#                                   haul_target_code,
#                                   ves_akr_name,
#                                   londd_end,
#                                   londd_start,
#                                   latdd_end,
#                                   latdd_start,
#                                   total_hooks_pots,
#                                   akfin_year,
#                                   avg_sst_celsius
#                                   from council.comprehensive_norpac
#                                   where year>2005
#                                   and year<2018
#                                   and gear_type=8")) %>% 
#   rename_all(tolower)
# 
# saveRDS(data,file="diana_longline_ancillary.RDS")


#  Round spatial grids to half degrees.
rounding=0.5
#  Sum the retained groundfish catch for each haul
#  For those sets with both deploy and retrieve coordinates, average the two. 
#  If there is only one (typically the retrieval location), use that location.
sumdat <- readRDS("diana_longline_ancillary.RDS") %>% 
  mutate(month=month(haul_date),
         retained_wt=(percent_retained/100)*extrapolated_weight) %>% 
  group_by(cruise,permit,haul_seq,akfin_year,month) %>% 
  summarise(retained_gf_weight=sum(retained_wt,na.rm=TRUE),
            duration=duration[1],
            vessel=ves_akr_name[1],
            londd_end=londd_end[1],
            londd_start=londd_start[1],
            latdd_end=latdd_end[1],
            latdd_start=latdd_start[1],
            hooks=total_hooks_pots[1],
            sst=mean(avg_sst_celsius,na.rm=TRUE)) %>% 
  ungroup %>% 
  rowwise() %>% 
  mutate(lon=sum(c(londd_end,londd_start),na.rm=TRUE)/sum(!is.na(c(londd_end,londd_start))), #if deploy and retrieval locations exist, average them
         lat=sum(c(latdd_end,latdd_start),na.rm=TRUE)/sum(!is.na(c(latdd_end,latdd_start))),
         lonr=round(lon/rounding)*rounding, #create our binned version for half degree grids
         latr=round(lat/rounding)*rounding)

# The above summarized haul level metadata and created half degree bins.
# Now summarize the data by bin, month, and year.
sumdat %>% 
  group_by(akfin_year,month,lonr,latr) %>% 
  summarise(retained_gf_weight=sum(retained_gf_weight,na.rm=TRUE),
            duration=sum(duration,na.rm=TRUE),
            vessels=length(unique(vessel)),
            sst=mean(sst),
            sst_sd=sd(sst),
            hooks=sum(hooks,na.rm=TRUE)) %>% 
  filter(vessels>2) %>% 
  saveRDS("longline_effort_halfdegree_bins_diana.RDS")
