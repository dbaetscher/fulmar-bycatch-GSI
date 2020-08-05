# nofu functions
library(dplyr)
library(tidyr)
library(rubias)

# @param base_ids: colony and sample id
# @param n_dsample: numeric, number to downsample to, e.g., 58
# @param baseline: baseline genotypes dataframe, 2-col format

downsample_baseline <- function(base_ids, n_dsample, base_genos, no_hi_missers){
  
    downsampled_three_colonies <- base_ids %>%
    group_by(repunit) %>%
    filter(repunit != "Chagulak") %>%
    sample_n(., n_dsample, replace = FALSE) #%>%
  #ungroup()
  
  # grab the Chagulak ids from the df
  chag_ids <- base_ids %>%
    filter(repunit == "Chagulak")
  
  # put the df back together
  baseline_dsampled_ids <- bind_rows(downsampled_three_colonies, chag_ids)
  
  # Use a semi_join to keep just those samples in the reference baseline genotypes
  dbaseline <- base_genos %>%
    semi_join(., baseline_dsampled_ids, by = "NMFS_DNA_ID")
  
  # Now, make a combined data frame with the bycatch and baseline samples to keep the integer values for haplotypes consistent.
  # here I want to combine the genotypes
  all_genos <- bind_rows(dbaseline, no_hi_missers)
  
  
  # first make integers of the alleles
  alle_idxs <- all_genos %>% 
    dplyr::select(NMFS_DNA_ID, locus, gene_copy, allele) %>%
    group_by(locus) %>%
    mutate(alleidx = as.integer(factor(allele, levels = unique(allele)))) %>%
    ungroup() %>%
    arrange(NMFS_DNA_ID, locus, alleidx) 
  
  # select just the columns to retain 
  alle_idx2 <- alle_idxs[,-4]
  
  # spread the alleles
  two_col <- alle_idx2 %>%
    unite(loc, locus, gene_copy, sep = ".") %>%
    spread(loc, alleidx)
  
  
  # I can use joins to split the combined data frame
  # first the mixture
  mix_df <- mixture_ids %>%
    left_join(two_col) %>%
    rename(indiv = NMFS_DNA_ID) %>%
    ungroup() # apparently the indiv column was still grouped and causing problems with rubias 
  
  mixture <- as_tibble(mix_df)
  mixture$repunit <- as.character(mixture$repunit)
  
  # and then the baseline
  base_df <- baseline_dsampled_ids %>%
    inner_join(., two_col) %>%
    rename(indiv = NMFS_DNA_ID) %>%
    ungroup()
  
  # perform mixture-assignment 
  mix_assign <- infer_mixture(reference = base_df, mixture = mixture, gen_start_col = 5, method = "MCMC", reps = 10000, burn_in = 1000)
  
  # individual data
  top_assign <- mix_assign$indiv_posteriors %>%
    group_by(indiv) %>%
    top_n(., 1, PofZ) %>%
    ungroup()
  
  # 90% assignment threshold
  # high-likelihood assignments with outlier z-scores?
  assign90 <- top_assign %>%
    filter(z_score < 2.5 & z_score > -2.5) %>%
    filter(PofZ > 0.90) %>%
    select(indiv, collection, PofZ, z_score)

  return(list(
    assign90
    
    
  ))
  
}


## BA spatial overlap from Abram
# https://github.com/abfleishman/trakR/blob/master/R/kernalOverlapBA_p.R

kernalOverlapBA_p<-function (tracks, tripid, groupid,
                             lon,lat,
                             colonyLon, colonyLat,
                             its, h, ud.grid, Plot){
  
  # get unique trips
  UniTripID<-unique(tracks[[tripid]]) #unique trips
  
  # get the unique levels of the group var
  GroupID_levels<-unique(tracks[[groupid]])
  
  # get n trips for each group and the total n trips
  GroupA_length<-length(unique(tracks[[ tripid ]][ tracks[[groupid]]==GroupID_levels[1] ]))
  GroupB_length<-length(unique(tracks[[tripid]][ tracks[[groupid]]==GroupID_levels[2] ]))
  n_UniTripID<-GroupA_length+GroupB_length
  
  # make a trip and group ID column for later use in the SPDF
  tracks$groupid<-tracks[[groupid]]
  tracks$trip.id<-tracks[[tripid]]
  
  # Warn and remove NAs
  if(nrow(tracks[is.na(tracks[[lat]]),])>0){
    lenny<-nrow(tracks)
    tracks<- tracks[!is.na(tracks[[lat]]),]
    warning(paste("There are",lenny-nrow(tracks)," NAs in  your position column. Removing!"))
  }
  
  # Calculate BA for the real groups ------------------
  
  # Make a SPDF
  tracks.spdf <- SpatialPointsDataFrame(coords=cbind(tracks[[lon]],
                                                     tracks[[lat]]),
                                        data=data.frame(id=tracks$groupid),
                                        proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  
  # Project data into laea
  tracks.spdf.t <- spTransform(tracks.spdf,
                               # CRS(paste("+proj=laea +units=km +lon_0=",
                               #           colonyLon, " +lat_0=", colonyLat	, sep="")))
                               CRS("+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs +towgs84=0,0,0"))
  
  # calculate kernelUD
  ud <- kernelUD(tracks.spdf.t, h = h,grid=ud.grid)
  
  # Plot the actual overlap
  # make empty df to solve cmd check issue unbound var
  # uds<-data.frame(long = numeric(0),
  #                 lat = numeric(0),
  #                 order = integer(0),
  #                 hole = logical(0),
  #                 piece = character(0),
  #                 id = character(0),
  #                 group = character(0),
  #                 ud = character(0))
  
  uds<-suppressWarnings(bind_rows(
    mutate(suppressMessages(fortify(getverticeshr(ud, percent=95, standardize=T))), ud="95"),
    mutate(suppressMessages(fortify(getverticeshr(ud, percent=50, standardize=T))), ud="50")))
  
  print(ggplot(data=uds)+
          geom_polygon(aes_string(x="long",y="lat",group="group",fill="id"),alpha=.5)+
          facet_wrap(~ud)+
          labs(fill=groupid))
  
  
  # 50% overlap
  BA_o50<-kerneloverlaphr(ud , method="BA", percent=50, conditional=TRUE)
  BA_o50<-BA_o50[1,2]
  
  # 95% overlap
  BA_o95<-kerneloverlaphr(ud , method="BA", percent=95, conditional=TRUE)
  BA_o95<-BA_o95[1,2] # saving these values from this matrix
  
  # Calculate BA for the Randomized groups ------------------
  
  # create data out structures
  RandomIndices_50<-numeric(length = its)
  RandomIndices_95<-numeric(length = its)
  
  for (i in 1:its){
    
    print(paste("iteration:" ,i,"of",its))
    
    # Shuffle the order of the UniTripIDs
    Idx<-sample(UniTripID, n_UniTripID, replace = FALSE, prob = NULL)
    
    # take the 1st x UniTripIDs where x is the sample size for group A
    UniTripID_A<-Idx[1:GroupA_length]
    
    # take the last x UniTripIDs where x is the sample size for group b. this
    # assumes that all UniTripID are classified into a group
    UniTripID_B<-Idx[(GroupA_length+1):GroupB_length]
    
    # add the random grouping to the data
    tracks$groupRan<-NA
    tracks$groupRan[tracks$trip.id%in%UniTripID_A]<-"A"
    tracks$groupRan[tracks$trip.id%in%UniTripID_B]<-"B"
    
    # Make SPDF with the random grouping as the ID column
    tracks.spdf <- SpatialPointsDataFrame(coords=cbind(tracks[[lon]],tracks[[lat]]),
                                          data=data.frame(id=tracks$groupRan),
                                          proj4string = CRS("+proj=longlat +ellps=WGS84
                                                            +datum=WGS84 +no_defs"))
    # Project data into laea
    tracks.spdf.t <- spTransform(tracks.spdf,
                                 # CRS(paste("+proj=laea +units=km +lon_0=", colonyLon,
                                 #           " +lat_0=", colonyLat	, sep="")))
                                 CRS("+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs +towgs84=0,0,0"))
    
    
    
    # make the kernelUD
    ud1 <- kernelUD(tracks.spdf.t, h = h,grid=ud.grid)
    
    # Plot each if desired (only for diagnostics?)
    if(Plot==T){
      # make empty df to solve cmd check issue unbound var
      # uds<-data.frame(long = numeric(0),
      #                 lat = numeric(0),
      #                 order = integer(0),
      #                 hole = logical(0),
      #                 piece = character(0),
      #                 id = character(0),
      #                 group = character(0),
      #                 ud = character(0))
      uds<-bind_rows(
        mutate(fortify(getverticeshr(ud1, percent=95, standardize=T)), ud="95"),
        mutate(fortify(getverticeshr(ud1, percent=50, standardize=T)), ud="50"))
      head(uds)
      print(ggplot(data=uds)+
              geom_polygon(aes_string(x="long",y="lat",group="group",fill="id"),alpha=.5)+
              facet_wrap(~ud))
    }
    
    # Calculate BA
    BA_50<-kerneloverlaphr(ud1 , method="BA", percent=50, conditional=TRUE)
    BA_95<-kerneloverlaphr(ud1 , method="BA", percent=95, conditional=TRUE)
    
    # add each iterations result
    RandomIndices_50[i]<-BA_50[1,2] # saving the overlap between group 1 & 2 (?)
    RandomIndices_95[i]<-BA_95[1,2]
  }
  
  # Calculate the P
  pval_50<-length(RandomIndices_50[RandomIndices_50<BA_o50])/its # this would need to be a combinatorial - currently % of times index observed is < than random?
  pval_95<-length(RandomIndices_95[RandomIndices_95<BA_o95])/its
  
  # make a results table
  results<-rbind(data.frame(ud=50,p=pval_50,BA=BA_o50,BA_rand_mean=mean(RandomIndices_50),
                            BA_rand_min=min(RandomIndices_50),BA_rand_max=max(RandomIndices_50),BA_rand_sd=sd(RandomIndices_50)),
                 data.frame(ud=95,p=pval_95,BA=BA_o95,BA_rand_mean=mean(RandomIndices_95),
                            BA_rand_min=min(RandomIndices_95),BA_rand_max=max(RandomIndices_95),BA_rand_sd=sd(RandomIndices_95)))
  
  return(results)
  # could have it return all values - for more information for indexes
  # make a list of the results with RandomIndices_95 and RandomIndices_50
  
}