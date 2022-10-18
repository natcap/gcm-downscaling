## GCMClimTool Library
## Angarita H., Yates D., Depsky N. 2014-2021
## This library contains some functions to pertorm the following tasks:
## - Read data of IPCCC GCM models from NetCDF files, in a user defined bounding box
## - Compare GCm results of a given pixel to  local regionalized climate 
## - Perform a Downcaling based on kNN-Boostrapping 
##
# R packages to install

library(chron)
library(fields)
library(ncdf4)
library(evir)
library(maptools)
library(stringr) 
library(dplyr)   # provides more flexible summarise functions for aggregating daily to monthly to yearly
library(plotrix) # multhist function for plotting super-imposed temperature histograms
library(MASS)    # kde2d 2-D Kernel Density Function
library(rgdal)   # readOGR function for updated shapefile importing
library(raster)  # not sure
#library(vgam) #generalized pareto distribution
#library(RNetCDF)

filter <- stats::filter # to unmask the masked filter function by dplyr

#----------------------------------------------------------------------
# *** this library contains all the functions of the GCM Climtool ***
#----------------------------------------------------------------------

#-----------------------------------------------
# I. GENERAL PURPOSE FUNCTIONS * LOAD FIRST
#-----------------------------------------------

val2pctl <- function(val, series) {
  length(series[series < val])/length(series)
}

# Linear Interpolation Function using two vectors of values. It works in cases where you have 2 vectors, one of say Flow = c(0, 5, 10, 15) and one of corresponding Depths = c(0, 2, 3, 3.5). Given a Flow value of say 7.2, to find what the corresponding Depth would be, this function interpolates between the appropriate Depth Values 2 and 3, which correspond to Flows of 5 and 10 to find the correct Depth value #
linapprox <- function(val, valvec, newvalvec)
{
  nearval <- valvec[which.min(abs(valvec - val))]
  valdiff <- val - nearval
  
  if(val < min(valvec) | val > max(valvec))
  {
    newvalest <- NA
    newvalest
  } else {
    if(valdiff > 0)
    {
      if((valvec[which.min(abs(valvec - val)) + 1] - valvec[which.min(abs(valvec - val))]) == 0)
      {
        factor <- 0
      } else {
        factor <- valdiff/(valvec[which.min(abs(valvec - val)) + 1] - valvec[which.min(abs(valvec - val))])
      }
      
      nearnewval <- newvalvec[which.min(abs(valvec - val))]
      newvaldif <- newvalvec[which.min(abs(valvec - val)) + 1] - newvalvec[which.min(abs(valvec - val))]
      newvalest <- nearnewval + factor*newvaldif
    } else if (valdiff < 0)
    {
      if((valvec[which.min(abs(valvec - val))] - valvec[which.min(abs(valvec - val)) - 1]) == 0)
      {
        factor <- 0
      } else {
        factor <- valdiff/(valvec[which.min(abs(valvec - val))] - newvalvec[which.min(abs(valvec - val)) - 1])
      }
      
      nearnewval <- newvalvec[which.min(abs(valvec - val))]
      newvaldif <- newvalvec[which.min(abs(valvec - val))] - newvalvec[which.min(abs(valvec - val)) - 1]
      newvalest <- nearnewval + factor*newvaldif
    } else {
      factor <- 0 
      newvalest <- newvalvec[which.min(abs(valvec - val))]
    }
    
    round(newvalest,3) 
  }
}

# Find Mode of Dataset #
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#number of wet days
wet_fraction = function(tSeries,threshold=0){
  a = which(tSeries > threshold)
  wetP = dim(t(a))[2]
  totalP = dim(t(tSeries[which(!is.na(tSeries))]))[2]
  wet_frac = wetP / totalP
  return(wet_frac)
}

# average length of wet/dry periods
wet_dry_spell_length = function(tSeries,threshold=0){
  
    nDat = dim(t(tSeries))[2]
    
    a = which(tSeries <= threshold)  #dry periods
    b <- a # to be filled with dry period IDs
    b[1] <- 1
    
    dLength = dim(t(a))[2]
    c = which(tSeries > threshold)  #wet periods
    wLength = dim(t(c))[2]
    d <- rep(NA, length(c)) # to be filled with wet period IDs
    d[1] <- 1
    d[which(c < a[1])] <- 1
    
    if (min(a,c) <= threshold) {
      numDrySpell = 1
      numWetSpell = 0
    } else {
      numWetSpell = 1
      numDrySpell = 0
    }
    
    for (i in 2:dLength){
      
      if (a[i-1] != a[i] - 1){
        numWetSpell = numWetSpell + 1
        numDrySpell = numDrySpell + 1 
        b[i] <- b[i-1] + 1
        d[which(c > a[i-1] & c < a[i])] <- max(d[!is.na(d)]) + 1
      }	else {
        b[i] <- b[i-1]
        d[which(c > a[i-1] & c < a[i])] <- max(d[!is.na(d)]) + 1
      }
    }

    maxDryLength <- length(b[which(b == getmode(b))])
    sdDryLength <- sd(sapply(unique(b), function(x) length(b[which(b == x)])))
    p90DryLength <- percentile(sapply(unique(b), function(x) length(b[which(b == x)])), 0.1)
    
    maxWetLength <- length(d[which(d == getmode(d))])
    sdWetLength <- sd(sapply(unique(d), function(x) length(d[which(d == x)])))
    p90WetLength <- percentile(sapply(unique(d), function(x) length(d[which(d == x)])), 0.1)
    
    avgD = dLength/numDrySpell 
    avgW = wLength/numWetSpell
    
    return (c(wetSpells = numWetSpell, drySpells = numDrySpell, wetLength = wLength, dryLength = dLength,avgDry = avgD, avgWet = avgW, maxDry = maxDryLength, maxWet = maxWetLength, StdDevDry = sdDryLength, StdDevWet = sdWetLength, p90Dry = p90DryLength, p90Wet = p90WetLength))
}

percentile = function(tSeries,prob_excedence,plt = FALSE){
  
  a = sort(tSeries,decreasing = TRUE)
  l = dim(t(a))[2]
  p = (1:l)/l
  
  i = tail(which(p<prob_excedence),1)
  
  if (plt == TRUE) {
    plot(p,a,type="l",col = "red",xlab="% Time exceeded", ylab="Value") 
    grid()
  }
  
  return(a[i])
}

bias = function(refSeries,newSeries){
  mref = mean(refSeries,na.rm=TRUE)
  mnew = mean(newSeries,na.rm=TRUE)
  
  return((mnew/mref - 1))
}


tri_state_joint_probability = function(timeSeries,lowerBound,higherBound){
  
  
  nData = dim(t(timeSeries))[2]
  
  cSeries = timeSeries[1:nData-1]
  fSeries = timeSeries[2:nData]
  
  jp_count = mat.or.vec(3, 3)
  
  # current state  = 1: lower or equal than lower bound
  sel = which(cSeries <= lowerBound)
  nSel = dim(t(sel))[2] 
  
  b = which(fSeries[sel] <= lowerBound)
  jp_count[1,1] = max(0,dim(t(b))[2])
  
  b = which(fSeries[sel] > higherBound)
  jp_count[1,3] = max(0,dim(t(b))[2]) 
  
  jp_count[1,2] = nSel - jp_count[1,1] - jp_count[1,3]
  
  
  # current state  = 3: higher than higher bound
  sel = which(cSeries > higherBound)
  nSel = dim(t(sel))[2] 
  
  b = which(fSeries[sel] <= lowerBound)
  jp_count[3,1] = max(0,dim(t(b))[2]) 
  
  b = which(fSeries[sel] > higherBound)
  jp_count[3,3] = max(0,dim(t(b))[2]) 
  
  jp_count[3,2] = nSel - jp_count[3,1] - jp_count[3,3]
  
  
  # current state  = 2: between lower and higher bound
  sel = which(cSeries <= higherBound & cSeries > lowerBound)
  nSel = dim(t(sel))[2] 
  
  b = which(fSeries[sel] <= lowerBound)
  jp_count[2,1] = max(0,dim(t(b))[2]) 
  
  b = which(fSeries[sel] > higherBound)
  jp_count[2,3] = max(0,dim(t(b))[2]) 
  
  jp_count[2,2] = nSel - jp_count[2,1] - jp_count[2,3]
  
  
  jp = jp_count / sum(jp_count)
  
  return(jp)
  
}


tri_state_joint_prob_monthly = function(timeSeriesdf, lowerBounds, higherBounds, monthly_thresholds = TRUE){ # requires a time series object input df with 3 columns --> year, month, value, and vectors of monthly lower and upper bounds, monthly_thresholds dictates whether or not to use monthly changing lower and upper bound thresholds, if FALSE lowerBounds and higherBounds are a single constant value 
  
  colnames(timeSeriesdf) <- c("Year", "Month", "Value")
  
  #in knn_bootstrap function, timeSeriesdf <- data_d
  
  jp_count = mat.or.vec(3, 3*12)
  jp = jp_count
  
  
  if(monthly_thresholds == FALSE)
  {
    lowerBounds <- rep(lowerBounds, 12)
    higherBounds <- rep(higherBounds, 12)
  }
  
  for(m in 1:12)
  {
    TS <- timeSeriesdf[-nrow(timeSeriesdf),]
    TS$Shift <- timeSeriesdf[2:nrow(timeSeriesdf), "Value"]
    monTS <- TS[TS$Month == m,]
    nData <- dim(t(monTS))[2]
    
    # current state  = 1: lower or equal than lower bound
    sel <- which(monTS$Value <= lowerBounds[m]) 
    
    nSel <- dim(t(sel))[2] 
    
    b <- which(monTS[sel, "Shift"] <= lowerBounds[m])
    jp_count[1,(m*3 - 2)] <- max(0,dim(t(b))[2])
    
    b <- which(monTS[sel, "Shift"] > higherBounds[m])
    jp_count[1,(m*3)] <- max(0,dim(t(b))[2]) 
    
    jp_count[1,(m*3 - 1)] <- nSel - jp_count[1,(m*3 - 2)] - jp_count[1,(m*3)]
    
    
    # current state  = 3: higher than higher bound
    sel <- which(monTS$Value > higherBounds[m])
    nSel <- dim(t(sel))[2] 
    
    b <- which(monTS[sel, "Shift"] <= lowerBounds[m])
    jp_count[3,(m*3 - 2)] <- max(0,dim(t(b))[2]) 
    
    b <- which(monTS[sel, "Shift"] > higherBounds[m])
    jp_count[3,(m*3)] <- max(0,dim(t(b))[2]) 
    
    jp_count[3,(m*3 - 1)] <- nSel - jp_count[3,(m*3 - 2)] - jp_count[3,(m*3)]
    
    # current state  = 2: between lower and higher bound
    sel <- which(monTS$Value <= higherBounds[m] & monTS$Value > lowerBounds[m])
    nSel <- dim(t(sel))[2] 
    
    b = which(monTS[sel, "Shift"] <= lowerBounds[m])
    jp_count[2,(m*3 - 2)] <- max(0,dim(t(b))[2]) 
    
    b <- which(monTS[sel, "Shift"] > higherBounds[m])
    jp_count[2,(m*3)] <- max(0,dim(t(b))[2]) 
    
    jp_count[2,(m*3 - 1)] <- nSel - jp_count[2,(m*3 - 2)] - jp_count[2,(m*3)]
    
    jp[,(m*3 - 2):(m*3)] <- round(jp_count[,(m*3 - 2):(m*3)] / sum(jp_count[,(m*3 - 2):(m*3)]),5)
    
    jp <- as.data.frame(jp)
    row.names(jp) <- c("from_Low", "from_Mid", "from_High")
    colnames(jp) <- as.vector(sapply(month.abb, function(x) paste0(c("to_Low", "to_Mid", "to_High"), x)))
  }
  
  return(jp)
  
  #pdf(file = paste0(getwd(), "/JP_Bootstrap_Monthly_plots.pdf"))
  par(mfrow = c(3, 3))
  for(i in 1:3)
  {
    plot(1:12,jp[i,seq(1,34,3)], type = 'l', xlab = "Mes", ylab = "Joint Probability", main = paste(row.names(jp)[i], "toLow"))
    plot(1:12,jp[i,seq(2,35,3)], type = 'l', xlab = "Mes", ylab = "Joint Probability", main = paste(row.names(jp)[i], "toMid")) 
    plot(1:12,jp[i,seq(3,36,3)], type = 'l', xlab = "Mes", ylab = "Joint Probability", main = paste(row.names(jp)[i], "toHigh"))
  }
  mtext(paste("Low Threshold =", round(lowerBounds,2),    "   Upper Threshold =", round(higherBounds,2)), outer = TRUE, cex = 0.7, line = -1)
  dev.off()
}


#----------------------------------------------------------------------------
# II. FUNCTION TO READ GCM AND OBSERVATIONS AND PRODUCE CORRESPONDING INTEGRATED DATA SETS
#----------------------------------------------------------------------------
read_GCM_NCDF = function(boundingBox = c(latMin = 1,latMax = 3,lonMin = -76,lonMax = -74),    # refDate = c(month = 1, day = 1, year = 2006)    # Start date of timestamp in GCM .nc file
                         movingDateStep = 0,                                                  # IMPORTANT: 1 in case reference date changes in each file, (such as in CMCC), otherwise 0
                         main_folder = "F:/GCM_ClimateTool",       
                         results_folder = "/Results",                                         # this path is relative to the main_directory
                         Experiment = "hist-rcp85",                                           # Default values are just for illustration purpose
                         Ensemble = "r1ip1",
                         variableName = "tas",
                         variableLabel = "Temperature [k]",
                         shp_file = "/GIS/WEAP_CATCHMENTS_V14_61_WGS84",                      # these paths are relative to the main_directory
                         obs_files = c(daily = "/IDEAM_OBS_DATA/R/SERIES_tas_DIARIO.csv",
                                       monthly = "/IDEAM_OBS_DATA/R/SERIES_tas.csv",
                                       stationsCatalog = "/IDEAM_OBS_DATA/R/CNE_IDEAM_V5_MAYO_6.csv"),
                         sta_bbox = "region", # determines where to search for observation stations, either in the user-defined bbox region, or in the GCMpixel cell (CHOICES: "region" | "pixel")
                         minObsYear = 1970,
                         minGCMYear = 2006) {

  #---------------------------------------------------------
  # PART 1. READS A GCM IN A BOUNDING BOX
  #----------------------------------------------------------  
  
  # Each folder contains ONE realization of the GCM
  fFolders = paste(main_folder,"GCMs",sep="")
  
  # IMPORTANT: SHP FILE with the boundary of the analysis area:
  gor=readOGR(paste(main_folder,shp_file,sep="")) 
  
  if(boundingBox[1] != -999) {
   latMin =  boundingBox[1]
   latMax =  boundingBox[2]  
   lonMin = boundingBox[3]
   lonMax = boundingBox[4]
  } else{ #get from shapefile
   latMin <- min(gor$Latitude); latMax <- max(gor$Latitude); lonMin <- min(gor$Longitude); lonMax <- max(gor$Longitude)
  }
  
  lonCentr <- mean(c(lonMin,lonMax)) # lon coordinates of centroid of bounding box
  latCentr <- mean(c(latMin,latMax)) # lat coordinates of centroid of bounding box
  
  Experiments <- Experiment
  
  for (Experiment in Experiments) {
    
    ensamble_day <-NULL
    ensamble_year <- NULL
    esm_count = 0;
  
    fList = list.files(path=paste(mFolder,"/required_files/GCMs/",sep=""),
                     pattern= paste("*",variableName,"_day_",model,"_",Experiment,"_",Ensemble,"_*",sep=""),full.names=TRUE)

    fName <- fList[1]
  
# initializes the data frame 
    esm_count = esm_count + 1;
    results_daily <- NULL   
    
    control = nc_open(fName, write=FALSE, readunlim=FALSE)
    modelName = model #att.get.nc(control, "NC_GLOBAL", "model_id")
    fConc= Experiment
        
    calendar <- ncatt_get(control, "time", "calendar")
    baseDate <- as.character(ncatt_get(control, "time", "units"))
    
    baseDate <- gsub("days since ","",baseDate[2])
    baseDate <- gsub(" 00:00:00","",baseDate)
    if(nchar(baseDate) >= 10) {  # format is 1985-01-01 or 1985-01-01 UTC
      baseYear <- as.numeric(substr(baseDate, 1, 4))
      baseMonth <- as.numeric(substr(baseDate, 6, 7))
      baseDay <- as.numeric(substr(baseDate, 9, 10))
    } else  {                    # format is 1985-1-1
      baseYear <- as.numeric(substr(baseDate, 1, 4))
      baseMonth <- as.numeric(substr(baseDate, 6,6))
      baseDay <- as.numeric(substr(baseDate, 8,8))
    }
    
    df$BaseDay <- baseDay
    df$BaseMonth <- baseMonth
    df$BaseYear<- baseYear

    df3$BaseDay <- baseDay
    df3$BaseMonth <- baseMonth
    df3$BaseYear <- baseYear
        
    #  if(experiment == "hist-rcp85"){
    #      rDate = c(month = baseMonth, day = baseDay, year = baseYear)
          #timeIdxs = var.get.nc(control, "time") #RNetCDF
    #      timeIdxs = ncvar_get(control, "time")  #ncdf4
    #      tCount = dim(timeIdxs)
    #  } else {
    #      rDate = c(month = 1, day = 1, year = minGCMYear[1])
    #      #timeIdxs = var.get.nc(control, "time")
    #      timeIdxs = ncvar_get(control, "time")
    #      timeIdxs = 0:(length(timeIdxs) - 1)
    #      tCount = length(timeIdxs)
    #  }
    
    rDate = c(month = baseMonth, day = baseDay, year = baseYear)
    timeIdxs = ncvar_get(control, "time")  #ncdf4
    tCount = dim(timeIdxs)

        # Reads coordinates and creates Lat-Lon grids, for plotting purposes:
        latmat = ncvar_get(control, "lat")
        lonmat = ncvar_get(control, "lon")
        
        # lonmat is required to be in the range -180,180. corrects if the netCDF comes in the 0-360 format,
        neglon = which(lonmat>180)
        lonmat[neglon] = lonmat[neglon]-360
        
        #creates Lat-Lon grids, for plotting purposes:
        latgrd <- outer(lonmat, latmat, FUN=function(xx,yy) { yy })
        longrd <- outer(lonmat, latmat, FUN=function(xx,yy) { xx }) 
        
        # finds the indexes of the grid around the bounding box: This is for One GCM pixel!!
        lonMin_idx = which.min(abs(lonCentr - lonmat)) # now finds the GCM pixel closest to the center of the bounding box
        lonMax_idx = which.min(abs(lonCentr - lonmat)) 
        latMin_idx = which.min(abs(latCentr - latmat)) 
        latMax_idx = which.min(abs(latCentr - latmat))
        
         #This is for all the GCM Pixels inside the bounding box
         Lat_idx_Range <- c(which.min(abs(latMax-latmat)):which.min(abs(latMin-latmat)))
         Lon_idx_Range <- c(which.min(abs(lonMax-lonmat)):which.min(abs(lonMin-lonmat)))
                                               
   
        if(length(timeIdxs)%%365 == 0) # checks if is using 365 day or not, assumes if not 365 day is using leap days
        {
          LeapFix <- "Yes"
          df$LeapFix <- "Yes"
          
        } else {
          LeapFix <- "No"
          df$LeapFix <- "No"
        }
        
        startTime = 1

        # chooses all the pixels in the bounding box
        st = c(min(Lon_idx_Range),min(Lat_idx_Range),startTime)        # lon, lat, time
        ct = c(length(Lon_idx_Range),length(Lat_idx_Range),tCount)
        tmp = ncvar_get(control, as.character(variableName), st, ct) #ncdf4
        data <- as.matrix(colMeans(colMeans(tmp))) #compute the mean of all the pixels

        # fixes funny behavior of var.get.nc function: trims the dimensions of the results if the count of the dimension = 1
        # should only need if reading single pixel.. 
        px1 = 1
        py1 = 1
        var_val = data[,1]
        
        # does some unit conversion:
        if (as.character(variableName) == "pr"){
          var_val = var_val * 86400                                # conversion from kg/ m^2 / s to mm/day 
        } else if (as.character(variableName) == "tas") {
          var_val = var_val #- 273.15                              # conversion from k to c                     
        }
        
        if(LeapFix == "Yes")
        {
          leap_ind_add <- 0
          time_year = month.day.year(timeIdxs + leap_ind_add,origin=rDate)$year
          time_month = month.day.year(timeIdxs + leap_ind_add,origin=rDate)$month
          time_day = month.day.year(timeIdxs + leap_ind_add,origin=rDate)$day
          
          leapday_ind <- which(time_month == 2 & time_day == 29)
          for(i in leapday_ind)
          {
            timeIdxs[i] <- timeIdxs[i] + 1
            timeIdxs[(i+1):length(timeIdxs)] <- timeIdxs[(i+1):length(timeIdxs)] + 1
          }
          
          time_year = month.day.year(timeIdxs + leap_ind_add,origin=rDate)$year
          time_month = month.day.year(timeIdxs + leap_ind_add,origin=rDate)$month
          time_day = month.day.year(timeIdxs + leap_ind_add,origin=rDate)$day
          
          leap_ind_add <- leap_ind_add + length(leapday_ind)
        } else {
          
          time_year = month.day.year(timeIdxs,origin=rDate)$year
          time_month = month.day.year(timeIdxs,origin=rDate)$month
          time_day = month.day.year(timeIdxs,origin=rDate)$day
        }
        
        results_daily_l <-NULL
        results_daily_l = data.frame(timeIdxs,time_year,time_month,time_day,var_val) 
        results_daily_l$timeIdxs <- julian(results_daily_l$time_month,results_daily_l$time_day,results_daily_l$time_year,rDate)  
        
        #quick fix of any strange negative values of PRECIP in GCMs:
        if (as.character(variableName) == "pr"){
          a = which(results_daily_l$var_val < 0)
          results_daily_l$var_val[a] = 0 #NaN
        }
        
        results_daily = rbind(results_daily,results_daily_l) #this will now only be called once, so no really needed
        print(fName)
        
        nc_close(control)
        
        # this moves the reference date one year  
        if (movingDateStep > 0) {
          rDate = c(rDate[1],rDate[2],rDate[3]+movingDateStep)
        }  
        
    seriesName <- paste(as.character(variableName),"_E",esm_count,sep="")
    colnames(results_daily)[which(colnames(results_daily) == "var_val")] <- seriesName
    
    if (as.character(variableName) == "pr") { 
      results_year <-aggregate(results_daily, by=list(results_daily$time_year), FUN=sum)
    } else {
      results_year <-aggregate(results_daily, by=list(results_daily$time_year), FUN=mean)      
    }
    
    results_year$time_year <- results_year$Group.1
    results_year <- results_year[,c("time_year", colnames(results_year)[length(colnames(results_year))])]
    
    if (esm_count == 1) {
      ensemble_day = results_daily
      ensemble_year = results_year
    } else {
      
      ensemble_day[seriesName] <- NULL
      ensemble_day[seriesName] <- results_daily[seriesName]
      
      ensemble_year[seriesName] <- NULL
      ensemble_year[seriesName] <- results_year[seriesName]
    }  
    
    # generates the df of the ensemble at monthly timestep
  
  ensemble_day$month_idx = (ensemble_day$time_year - rDate[3]) * 12 + ensemble_day$time_month
  
  if (as.character(variableName) == "pr") {
    
    ensemble_month <- aggregate(ensemble_day, by = list(ensemble_day$month_idx), FUN=sum)
    ensemble_month2 <- aggregate(ensemble_day, by = list(ensemble_day$month_idx), FUN=mean)   #dummy dataframe
    
    ensemble_month$time_year <- ensemble_month2$time_year
    ensemble_month$time_month <- ensemble_month2$time_month
    ensemble_month <- ensemble_month[,c("time_year","time_month",colnames(ensemble_month)[length(colnames(ensemble_month)) - 1])]
    
  } else {
    ensemble_month <- aggregate(ensemble_day, by = list(ensemble_day$month_idx), FUN=mean)   #dummy dataframe   
    ensemble_month <- ensemble_month[,c("time_year","time_month",colnames(ensemble_month)[length(colnames(ensemble_month)) - 1])]
  }
  
  ## Bounding box around the GCM pixel. So far the 1,1 pixel.
  lon_pixel_space <- lonmat[2] - lonmat[1]
  lat_pixel_space <- latmat[2] - latmat[1]
  
  pix_lonMin = lonmat[lonMin_idx] - lon_pixel_space / 2 #- 180
  pix_lonMax = lonmat[lonMin_idx] + lon_pixel_space / 2 #- 180
  
  pix_latMin = latmat[latMin_idx] - lat_pixel_space / 2
  pix_latMax = latmat[latMin_idx] + lat_pixel_space / 2
  
  
  # -----------------------------------
  # saves dataframe with the GCM query, of each experiment
  
  save(ensemble_year,file=paste(main_folder,results_folder,"/",paste(modelName,fConc,variableName,"ensemble_year.Rda",sep="_"),sep=""))
  save(ensemble_day,file=paste(main_folder,results_folder,"/",paste(modelName,fConc,variableName,"ensemble_day.Rda",sep="_"),sep=""))
  save(ensemble_month,file=paste(main_folder,results_folder,"/",paste(modelName,fConc,variableName,"ensemble_month.Rda",sep="_"),sep=""))
  
  }
  
  
  #---------------------------------------------------------
  # PART 2. READS OBSERVATIONS WITHIN THE BOUNDING BOX OF THE ANALYSIS AREA
  #----------------------------------------------------------
  oPath = paste(main_folder,obs_files,sep="")
  # Reads all historic records from a CSV file:
  VAR_day = read.csv(oPath[1],check.names=F,dec=".") ###EDIT MANON TO ADD IN CHECK.NAMES
  VAR_month = read.csv(oPath[2],check.names=F,dec=".")     # csv file with the variable pr, ts, etc ###EDIT MANON TO ADD IN CHECK.NAMES
  
  if (as.character(variableName) == "pr") {            # careful. TS variables are averaged instead of summed
    VAR_year <-aggregate(VAR_month, by=list(VAR_month$Year), FUN=sum) 
    VAR_year$Year <- VAR_year$Group.1
    VAR_year$Group.1 <- NULL
    VAR_year$timeStamp <- NULL
    VAR_year$Month <- NULL
  } else {
    VAR_year <-aggregate(VAR_month, by=list(VAR_month$Year), FUN=mean)
    VAR_year$Year <- VAR_year$Group.1
    VAR_year$Group.1 <- NULL
    VAR_year$timeStamp <- NULL
    VAR_year$Month <- NULL
  }
  
  availableStages = names(VAR_year)  # codes of stages with data
  
  # Reads station data (id,lat,lon, etc):
  catalog = read.csv(oPath[3],header=TRUE, sep=",")
  
    # CORRECTS THE BOUNDING BOX FORMAT:

    if (pix_lonMin > 180) {
      pix_lonMin = pix_lonMin - 360
    }
    
    if (pix_lonMax > 180) {
      pix_lonMax = pix_lonMax - 360  
    }
    
    if (lonMin > 180) {
      lonMin = lonMin - 360
    } 
    
    if (lonMax > 180) {
      lonMax = lonMax - 360  
    } 

  # Sub-catalog in the bounding box
  if(sta_bbox == "region")
  {
    catalog_bb <- subset(catalog, catalog$longitud >= lonMin & catalog$longitud <= lonMax & catalog$latitud >= latMin & catalog$latitud <= latMax)
  } else {
    
    catalog_bb <- subset(catalog, catalog$longitud >= pix_lonMin & catalog$longitud <= pix_lonMax & catalog$latitud >= pix_latMin & catalog$latitud <= pix_latMax)
  }

  #catalog_bb  <- catalog ## override, just take all stations
  stage_codes = catalog_bb$CODIGO_CAT
  
  # generates the catalog of required AND available data:
  # monthly
  stages = NULL
  stages_x = NULL
  stages_y = NULL
  
  for (st in stage_codes){
    if(paste(as.character(variableName),'_',st,sep = "") %in% availableStages){
      
      id = which(catalog_bb$CODIGO_CAT == st)
      stages_x = c(stages_x,catalog_bb$longitud[id])
      stages_y = c(stages_y,catalog_bb$latitud[id])
      stages = c(stages,paste(as.character(variableName),'_',st,sep = ""))
    }  
  }
  
  # daily:
  availableStages_d = names(VAR_day)  # codes of stages with data
  
  stages_d = NULL
  stages_d_x = NULL
  stages_d_y = NULL
  
  for (st in stage_codes){
    if(paste(as.character(variableName),'_',st,sep = "") %in% availableStages_d){
      id = which(catalog_bb$CODIGO_CAT == st)
      stages_d_x = c(stages_d_x,catalog_bb$longitud[id])
      stages_d_y = c(stages_d_y,catalog_bb$latitud[id])
      stages_d = c(stages_d,paste(as.character(variableName),'_',st,sep = ""))
    }	
  }
  
  mc = c("Year")
  mcm = c("Year","Month")
  mcd = c("timeStamp","Year","Month","Day")
  
  #computes average of stations. missing values are not included
  if(variableName == "pr")
  {
    data_d <- subset(VAR_day, select = c(mcd,stages_d))
    data_d$avg_value <- rowMeans(data.matrix(subset(data_d, select = stages_d)), na.rm = TRUE)
    data_m <- subset(VAR_month, select = c(mcm,stages))
    data_m$avg_value <- rowMeans(subset(data_m, select = stages), na.rm = TRUE)
    
    data_y <- subset(VAR_year, select = c(mc,stages))
    data_y$avg_value <- rowMeans(subset(data_y, select = stages), na.rm = TRUE)
    
  } else {
    data_d <- subset(VAR_day, select = c(mcd,stages_d))
    data_d$avg_value <- rowMeans(data.matrix(subset(data_d, select = stages_d)), na.rm = TRUE)
    data_m <- subset(VAR_month, select = c(mcm,stages))
    data_m$avg_value <- rowMeans(subset(data_m, select = stages), na.rm = TRUE)
    data_y <- subset(VAR_year, select = c(mc,stages))
    data_y$avg_value <- rowMeans(subset(data_y, select = stages), na.rm = TRUE)
  }

  Rcd <- as.data.frame(matrix(NA, nrow = length(data_d$timeStamp), ncol = length(grep(variableName,colnames(data_d))))) # creates blank matrix for plotting non-NAs in period of record plot, will consist of binary 1/0 values indicating present/missing data
  
  date_vec <- as.Date(paste0(data_d$Month,"/",data_d$Day,"/",data_d$Year), format = "%m/%d/%Y")
  
  colnames(Rcd) <- c(colnames(data_d)[grep(variableName,colnames(data_d))]) # assigns the station names to the newly-created Rcd matrix
  
  Rcd[,] <- data_d[,grep(variableName,colnames(data_d))]
  Rcd[!is.na(Rcd)] <- 1
  
  # -----------------------------------
  # saves dataframes with the results of historic observations
  
  save(data_y,file=paste(main_folder,results_folder,"/",paste(modelName,"observed",variableName,"year.Rda",sep="_"),sep=""))
  save(data_m,file=paste(main_folder,results_folder,"/",paste(modelName,"observed",variableName,"month.Rda",sep="_"),sep=""))
  save(data_d,file=paste(main_folder,results_folder,"/",paste(modelName,"observed",variableName,"day.Rda",sep="_"),sep=""))
  
  # --------------------------------------
  # PART 3 MAKES A FEW PLOTS, EXPORTED IN PDF FORMAT:
  # --------------------------------------
  
  pdf(file=paste(main_folder,results_folder,"/",paste(modelName,fConc,variableName,'PLOT.pdf',sep="_"),sep=""), height=8.5, width=11, onefile=TRUE, family='Helvetica', paper="a4r", pointsize=12)   
  par(mfrow = c(3,1))
  par(mar = c(4,4,3,2))
  ## a. plots the context map:
  plot(catalog$longitud,catalog$latitud,col = "grey" ,pch = 16, cex = .75,asp=1,xlab="Lon", ylab="Lat", main = modelName)
  world(add=T)
  plot(gor,border = "red",add=T)
  grid()
  points(stages_d_x,stages_d_y,col = "blue" ,pch = 15, cex = .75)
#  rect(pix_lonMin, pix_latMin, pix_lonMax, pix_latMax, border = "green", lwd = 1.5) # bbox of GCM pixel
  rect(lonMin, latMin, lonMax, latMax, lwd = 1.5, border = "turquoise1") # bbox of region
  
   legend( x="topright", 
#          legend=c("All stations","bbox GCM pixel", paste(variableName,"Stations in GCM bbox"), paste("bbox",region), "Shapefile Polygons"), 
          legend=c("All stations",paste(variableName,"Stations in GCM bbox"), paste("bbox",region), "Shapefile Polygons"), 
#          col=c("grey","green", "blue", "turquoise1", "red"), lwd=1, lty=c(NA,1,NA,1,1), 
          col=c("grey","blue", "turquoise1", "red"), lwd=1, lty=c(NA,NA,1,1), 
#          pch=c(16,NA,15,NA,NA),  merge=FALSE, bg = "white")
          pch=c(16,15,NA,NA),  merge=FALSE, bg = "white")
  
  
  ## b. plots the map with the GCM grid and the local stations:
  plot(stages_d_x,stages_d_y,col = "blue" ,pch = 16, cex = 0.75,asp=1, 
       xlim=c(lonMin, lonMax), ylim=c(min(latMin, boundingBox["latMin"]), max(latMax, boundingBox["latMax"])*1.01),
       xlab="Lon", ylab="Lat")
  
  if (longrd[lonMin_idx:lonMax_idx,latMin_idx:latMax_idx] > 180) {
    lngrd <- (longrd[lonMin_idx:lonMax_idx,latMin_idx:latMax_idx]-360)
  } else {
    lngrd <- longrd[lonMin_idx:lonMax_idx,latMin_idx:latMax_idx]
  }
  world( add=T )
  plot(gor,add=T)
  grid()
  #rect(pix_lonMin, pix_latMin, pix_lonMax, pix_latMax, border = "green", lwd = 3)
  rect(lonMin, latMin, lonMax, latMax, lwd = 3, border = "turquoise1") # bbox of region
  points(lngrd,latgrd[lonMin_idx:lonMax_idx,latMin_idx:latMax_idx],pch = 17,cex = 1.5)
  longrd[longrd > 180] <- longrd[longrd > 180] -360
  points(longrd, latgrd, pch = 17, cex = 1.5)
  
  
  legend( x="topright", 
#          legend=c(paste("Pixel GCM - ", modelName),"bbox GCM pixel", paste(variableName,"stations in bbox"), paste("bbox",region), "Shapefile Polygons"), 
          legend=c(paste("GCM Points - ", modelName),paste(variableName,"stations in bbox"), 
                   paste("bbox",region), "Shapefile Polygons"), 
#         col=c("black","green", "blue", "turquoise1","red"), lwd=1, lty=c(NA,1,NA,1,1), 
          col=c("black","blue", "turquoise1","red"), lwd=1, lty=c(NA,NA,1,1), 
#          pch=c(17,NA,15,NA,NA), merge=FALSE, bg = "white")
           pch=c(17,15,NA,NA), merge=FALSE, bg = "white")

  ## c. plots the historical and GCM time series
  par(xpd = TRUE)
  
  pallete = terrain.colors(dim(t(fFolders))[2]+1)
  ay = which(data_y$Year >= minObsYear)   # filters out years previous to tmin
  
  nplot_exp = 0
  
  for (Experiment in Experiments) {

    fConc = Experiment
    load(file=paste(main_folder,results_folder,"/",paste(modelName,fConc,variableName,"ensemble_year.Rda",sep="_"),sep=""))
    
    nplot_exp = nplot_exp +  1
    
    if (nplot_exp == 1)
      {
        plot(ensemble_year$time_year,ensemble_year[,2],type="l", xlim = c(1850, 2100), ylim=c(min(ensemble_year[,2][is.finite(ensemble_year[,2])], data_y$avg_value[ay][is.finite(data_y$avg_value[ay])]),1.01*max(ensemble_year[,2][is.finite(ensemble_year[,2])], data_y$avg_value[ay][is.finite(data_y$avg_value[ay])])),col = 'green',xlab="Year", ylab=variableLabel)
      } else {
      lines(ensemble_year$time_year,ensemble_year[,2])
      }
    
  }
  
  
  lines(data_y$Year[ay], data_y$avg_value[ay], col = "blue",pin=c(12,8))
  
  legend( x="top", 
          legend=c("Observation Average", paste("GCM", "-", modelName)), 
          col=c("blue",pallete[1]), lwd=1, lty=c(1,1), 
          pch=c(NA,NA), merge=FALSE, bg = "white", horiz = TRUE,inset = c(0, -0.155))
  grid()
  par(xpd = FALSE)
  
  
  # d. plots n days 
  n = 1095   # i.e. 3 years
  if (fConc == "hist-rcp85"){
    tmin = max(minObsYear,min(data_d$Year))
  } else {
    tmin = 1990  # 2008
  }
  tmax = n / 365 + tmin
  
  if (variableName == "pr") {
    pmin = 0
    pmax = 100
  } else if (variableName == "tas") {
    pmin = 0
    pmax = 50   
  }
  
  ayo = which(data_d$Year == tmin)[1]   # filters out years previous to tmin
  aym = which(ensemble_day$time_year==tmin)[1]
  
  #Plots Observed Station Period of Record #
  
  par(mfrow=c(1,1))
  
  plot(date_vec, Rcd[,1], lwd = 3,  col = 4, ylab = '', xlab = '', ylim = c(0, dim(Rcd)[2]), yaxt = 'n', type = 'l', main = paste0("Period of Record Observed ", variableLabel," Stations")) # plots base graph
  
  for (l in 1:dim(Rcd)[2]) # overlays original record lines in blue
  {
    lines(date_vec, l*Rcd[,l], lwd = 3,  col = "dodgerblue")
  }
  
  par(las = 1, cex.axis = 0.50) # sets the y labels to be horizontal
  axis(2, 1:dim(Rcd)[2], labels = colnames(Rcd)) # places station names as the y-axis labels
  
  dev.off()
  
  save(df, file = paste0(mFolder,rFolder,"df.Rda"))
}

#----------------------------------------------------------------------------------
# III. FUNCTION TO GENERATE COMPARATIVE STATS OF A GIVEN GCM TO LOCAL OBSERVATIONS 
#----------------------------------------------------------------------------------

compare_GCM = function(refDate = c(month = 1, day = 1, year = 1850),
                       main_folder = "d:/ncar/dropbox/GCM_ClimateTool",       
                       results_folder = "/Results/chinchina",
                       enso_signal_file = "/ENSO_SIGNAL/ONI_SIGNAL.csv",
                       svd_analysis_file = "/SVD/SVD_Table_Comp1_Pac_Atl.csv",
                       modelName = "CCSM4",
                       historical = "historical",
                       futures = c("rcp85"),
                       varName = "pr",
                       varLabel = "Precipitation (mm)",
                       minObsYear = 1970,
                       maxObsYear = 2000,
                       minGCMYear = 2040,
                       maxGCMYear = 2100,
                       alignHistYears = TRUE) { 
  
  # reads future
  fPath = paste(main_folder,results_folder,"/",sep="")
  fConc = futures
  
  fName = paste(fPath,paste(modelName,fConc,varName,"ensemble_day.Rda",sep="_"),sep="")
  load(fName,verbose=TRUE) #ensemble_day
  
## model Future
  ensemble_day$month_idx = (ensemble_day$time_year - refDate[3]) * 12 + ensemble_day$time_month
  ensemble_day_f <- ensemble_day[ensemble_day$time_year >= minGCMYear & ensemble_day$time_year <= maxGCMYear,]
  #ensemble_day_h <- ensemble_day[ensemble_day$time_year >= minObsYear & ensemble_day$time_year <= maxObsYear,]
  
  fName = paste(fPath,paste(modelName,fConc,varName,"ensemble_year.Rda",sep="_"),sep="")
  load(fName,verbose=TRUE)   
  ensemble_year_f <- ensemble_year[ensemble_year$time_year >= minGCMYear & ensemble_year$time_year <= maxGCMYear,]
  #ensemble_year_h <- ensemble_year[ensemble_year$time_year >= minObsYear & ensemble_year$time_year <= maxObsYear,]
  
  fName = paste(fPath,paste(modelName,fConc,varName,"ensemble_month.Rda",sep="_"),sep="")
  load(fName,verbose=TRUE)   
  ensemble_month_f <- ensemble_month[ensemble_month$time_year >= minGCMYear & ensemble_month$time_year <= maxGCMYear,]
  #ensemble_month_h <- ensemble_month[ensemble_month$time_year >= minObsYear & ensemble_month$time_year <= maxObsYear,]

  
  # model historic
  fName = paste(fPath,paste(modelName,"historical",varName,"ensemble_day.Rda",sep="_"),sep="")
  load(fName,verbose=TRUE) 
  ensemble_day_h <- ensemble_day[ensemble_day$time_year >= minObsYear & ensemble_day$time_year <= maxObsYear,]
  
  fName = paste(fPath,paste(modelName,"historical",varName,"ensemble_year.Rda",sep="_"),sep="")             
  load(fName,verbose=TRUE)
  ensemble_year_h <- ensemble_year[ensemble_year$time_year >= minObsYear & ensemble_year$time_year <= maxObsYear,]
  
  fName = paste(fPath,paste(modelName,"historical",varName,"ensemble_month.Rda",sep="_"),sep="")             
  load(fName,verbose=TRUE)  
  ensemble_month_h <- ensemble_month[ensemble_month$time_year >= minObsYear & ensemble_month$time_year <= maxObsYear,]
  
    
# observed historic
  fName = paste(fPath,paste(modelName, "observed", varName,"day.Rda",sep="_"),sep="")
  load(fName,verbose=TRUE)
  
  fName = paste(fPath,paste(modelName, "observed", varName,"year.Rda",sep="_"),sep="")
  load(fName,verbose=TRUE)
  
  fName = paste(fPath,paste(modelName, "observed", varName,"month.Rda",sep="_"),sep="")
  load(fName,verbose=TRUE)
  
  if(alignHistYears == TRUE)
  {
    data_d <- data_d[data_d$Year >= minObsYear & data_d$Year <= max(ensemble_year_h$time_year),]
    data_m <- data_m[data_m$Year >= minObsYear & data_m$Year <= max(ensemble_year_h$time_year),]
    data_y <- data_y[data_y$Year >= minObsYear & data_y$Year <= max(ensemble_year_h$time_year),]
    
    ensemble_day_h <- ensemble_day_h[ensemble_day_h$time_year >= minObsYear,]
    ensemble_month_h <- ensemble_month_h[ensemble_month_h$time_year >= minObsYear,]
    ensemble_year_h <- ensemble_year_h[ensemble_year_h$time_year >= minObsYear,]
  }
  
  #----------------------------------------------------------------------
  # Reads NINO: Oceanic nino oscilation (ONI) series:
  #----------------------------------------------------------------------
  
  df_oni_m = read.csv(paste(main_folder,enso_signal_file,sep=""))
  df_oni_y = aggregate(df_oni_m$ONI,by=list(df_oni_m$Year),FUN=mean)
  
  #----------------------------------------------------------------------
  # Reads SVD analysis of 1856-2015 Indo-Pacific & Atlantic SSTa series:
  #----------------------------------------------------------------------
  
  svd_m = read.csv(paste(main_folder,svd_analysis_file,sep=""))
  
  if(alignHistYears == TRUE)
  {
    svd_m <- svd_m[svd_m$Year >= minObsYear,]
  }
  
  svd_y = aggregate(svd_m,by=list(svd_m$Year),FUN=mean)
  svd_y$Group.1 <- NULL
  svd_y$Month <- NULL
  
  #----------------------------------------------------------------------
  # Calculates anomalies with respect to refPeriod starting in Ref Year and ending now
  #----------------------------------------------------------------------
  
  aym = which((ensemble_year_h$time_year >= minObsYear)) 
  gcm_avg = mean(ensemble_year_h$pr_E1[aym],na.rm=TRUE)
  
  anm = which(is.na(ensemble_year_h$pr_E1))                    # completes missing data with the mean
  ensemble_year_h$pr_E1[anm] = gcm_avg
  
  ensemble_year_h$anomaly <- ensemble_year_h$pr_E1 - gcm_avg     # historic anomaly
  ensemble_year_f$anomaly <- ensemble_year_f$pr_E1 - gcm_avg 
  
  ayo = which((data_y$Year >= minObsYear)) 
  obs_avg = mean(data_y$avg_value[ayo])
  data_y$anomaly <- data_y$avg_value - obs_avg            # observed anomaly
  
  #----------------------------------------------------------------------
  # Plots to show the seasonal patterns in the observation and in the GCM
  #----------------------------------------------------------------------
  
  pdf(file=paste(main_folder,results_folder,"/",modelName,'_',fConc,"_",varName,'_Intra_Intro_Anomaly.pdf',sep=""), height=8, width=12, onefile=TRUE, family='Helvetica', paper="a4r", pointsize=12) 
  
  par(mfrow = c(3,1), mar = c(2,4,5,1))
  maxy <- 1.03*max(ensemble_month_h$pr_E1[is.finite(ensemble_month_h$pr_E1)], 
                   data_m$avg_value[is.finite(data_m$avg_value)], 
                   ensemble_month_f$pr_E1[is.finite(ensemble_month_f$pr_E1)])
  
  # a. box plot of stationality
  ay = which(data_m$Year >= minObsYear) #1970)
  nam = c("","","","","","","","","","","","")
  
  boxplot(pr_E1~time_month,data=ensemble_month_h, main=modelName, 
          xlab="Month", ylab="Precipitation (mm)", ylim = c(0,maxy), col="cyan", boxwex  = 0.2, at = 1:12 - 0.0,yaxs    = "i" ) # model 1980 -2015  
  
  boxplot(avg_value[ay]~Month[ay],data=data_m, 
          xlab="Month",names = nam, ylab="Precipitation (mm)",col="green3",  boxwex  = 0.2, at = 1:12 - 0.3,add=TRUE)  # historic model
  
  boxplot(pr_E1~time_month,data=ensemble_month_f, 
          xlab="Month",names = nam, ylab="Precipitation (mm)",col="orange", boxwex  = 0.2, at = 1:12 + 0.3,add=TRUE)    # model 2005-2100
  
  AvgMonth_h=aggregate(ensemble_month_h[,3], by=list(Month=ensemble_month_h$time_month),mean,na.rm=FALSE)    # historic model
  AvgMonth_O=aggregate(data_m[ay,ncol(data_m)], by=list(Month=data_m$Month[ay]),mean,na.rm=TRUE)             # observed
  AvgMonth_f=aggregate(ensemble_month_f[,3], by=list(Month=ensemble_month_f$time_month),mean,na.rm=TRUE)     # prospective model
  
  lines(AvgMonth_O$Month,AvgMonth_O$x,type="l",col = "green3",lwd = 2)
  lines(AvgMonth_h$Month,AvgMonth_h$x,type="l",col = "cyan",lwd = 2)
  lines(AvgMonth_f$Month,AvgMonth_f$x,type="l",col = "orange",lwd = 2)
  
  
  par(xpd = TRUE)
  legend(x="top", legend = c(paste0("Observed ", minObsYear, "-", max(data_m$Year)),
                             paste0("GCM ", min(ensemble_month_h$time_year), "-",max(ensemble_year_h$time_year)),
                             paste0("GCM ",min(ensemble_year_f$time_year), "-",maxGCMYear),
                             paste0("Mean Observed"),
                             paste0("Mean GCM Hist"),
                             paste0("Mean GCM Fut")),
  fill= c("green3","cyan","orange","green3","cyan","orange"), horiz = T, inset = c(0,-0.12), bty = 'n', text.width=rep(1.4,6))
  
  # b. Plots  periodogram
  par(xpd = FALSE)
  par(mar = c(2,4,4,1))
  ay = which(data_y$Year>= minObsYear) #1970)
  sp_obs = spectrum(data_y$avg_value[ay],plot=F)
  sp_gcm = spectrum(ensemble_year_h$pr_E1,plot=F)
  sp_oni = spectrum(df_oni_y$x,plot=F)
  sp_svd_pac <- spectrum(svd_y$PACIFIC_comp_1, plot = F)
  sp_svd_atl <- spectrum(svd_y$ATLANTIC_comp_1, plot = F)
  smo = max(sp_obs$spec)
  
  plot(1/sp_obs$freq,sp_obs$spec/smo,
       type="l",log="x",xlim=c(1, 150),ylim=c(0, 1),col = "blue",xlab="Period (years)", ylab="Power Spectrum std.")
  par(new=T)
  
  smg = max(sp_gcm$spec)
  
  plot(1/sp_gcm$freq,sp_gcm$spec/smg,
       type="l",log="x",xlim=c(1, 150),ylim=c(0, 1),col = "green3",xlab="", ylab="")
  par(new=T)
  
  smo = max(sp_oni$spec)
  
  plot(1/sp_oni$freq,sp_oni$spec/smo,
       type="l",log="x",xlim=c(1, 150),ylim=c(0, 1),col = "red",xlab="", ylab="")
  par(new=T)
  
  sm_svd_p = max(sp_svd_pac$spec)
  
  plot(1/sp_svd_pac$freq,sp_svd_pac$spec/sm_svd_p,
       type="l",log="x",xlim=c(1, 150),ylim=c(0, 1),col = "orange",xlab="", ylab="", lty = "dashed", lwd = 1.5)
  par(new=T)
  
  sm_svd_a = max(sp_svd_atl$spec)
  
  plot(1/sp_svd_atl$freq,sp_svd_atl$spec/sm_svd_a,
       type="l",log="x",xlim=c(1, 150),ylim=c(0, 1),col = "magenta",xlab="", ylab="", lty = "dashed", lwd = 1.5)
  par(new=T)
  
  par(xpd = TRUE)
  legend( x="topright", 
          legend=c("Observed Periodogram      ","Periodogram GCM Hist", "ENSO:ONI", "SVD Pacific", "SVD Atlantic"), 
          col=c("blue","green3","red","orange","magenta"), lwd=1.5, lty=c(1,1,1,2,2), 
          pch=c(NA,NA,NA,NA,NA), merge=FALSE, bty = 'n')
  par(xpd = FALSE)
  
  grid()
  
  # c.Anomalies 
  tmin = minObsYear
  tmax = maxGCMYear
  pmin = min(ensemble_year_h$anomaly,ensemble_year_f$anomaly,data_y$anomaly[ayo], na.rm = T)
  pmax = max(ensemble_year_h$anomaly,ensemble_year_f$anomaly,data_y$anomaly[ayo], na.rm = T)
  
  ayo = which(data_y$Year >= minObsYear)
  
  plot(ensemble_year_h$time_year,ensemble_year_h$anomaly,type="l",col="cyan",
       xlim=c(tmin, tmax), ylim=c(pmin, pmax),xlab="Year", ylab="Precipitation - Avg. [mm]")
  par(new=T)
  plot(ensemble_year_f$time_year,ensemble_year_f$anomaly,type="l",col="orange",
       xlim=c(tmin, tmax), ylim=c(pmin, pmax),xlab="", ylab="")
  par(new=T)
  plot(data_y$Year[ayo],data_y$anomaly[ayo],type="l",col="green3",
       xlim=c(tmin, tmax), ylim=c(pmin, pmax),xlab="", ylab="")
  CurveFit <- lowess(ensemble_year_f$time_year,ensemble_year_f$anomaly,f=1) ## linear fit
  lines(CurveFit$x,CurveFit$y,lty=2,col="grey50")
  
  par(xpd = TRUE)
  legend( x="top", 
          legend=c("Observed      ","GCM Historical    ", "GCM Future"), 
          col=c("green3","cyan","orange"), lwd=1.5, lty=c(1,1,1), 
          pch=c(NA,NA,NA), merge=FALSE, inset = c(0,-0.12), horiz = T, bty = 'n')
  par(xpd = FALSE)
  
  grid()
  
  dev.off()
  
  #----------------------------------------------------------------------
  # Calculates other metrics  
  #----------------------------------------------------------------------
  
  wf_model <-NULL
  wet_spell_model <-NULL
  dry_spell_model <-NULL
  p_05_model <-NULL
  b_model <- NULL
  mean_model <-NULL
  
  k = 1
  idx_results <- data.frame(model = character(),
                            parameter=character(),
                            start_year=as.numeric(character()), 
                            end_year=as.numeric(character()),
                            value=as.numeric(character()),
                            stringsAsFactors=FALSE) 
  
  ayo = which(data_d$Year >= minObsYear)
  
  wf_obs = wet_fraction(data_d$avg_value[ayo],1)
  
  spellLengths_obs = wet_dry_spell_length(data_d$avg_value[ayo],1)
  wet_spell_obs <-spellLengths_obs[6]
  dry_spell_obs <-spellLengths_obs[5]
  p_05_obs = percentile(data_d$avg_value[ayo],0.05)
  mean_obs = mean(data_d$avg_value[ayo],na.rm=TRUE)
  
  
  idx_results[k,1] = "obs"; idx_results[k,2] = "wet_fraction"; 
  idx_results[k,3] = minObsYear; 
  idx_results[k,4] = min(minObsYear + 40, max(ensemble_year_f$time_year));
  idx_results[k,5] = wf_obs; k = k + 1
  
  idx_results[k,1] = "obs"; idx_results[k,2] = "wet_spell_length"
  idx_results[k,3] = minObsYear; idx_results[k,4] = min(minObsYear + 40, max(ensemble_year_f$time_year))
  idx_results[k,5] = wet_spell_obs; k = k + 1
  
  idx_results[k,1] = "obs"; idx_results[k,2] = "dry_spell_length"
  idx_results[k,3] = minObsYear; idx_results[k,4] = min(minObsYear + 40, max(ensemble_year_f$time_year))
  idx_results[k,5] = dry_spell_obs; k = k + 1
  
  idx_results[k,1] = "obs"; idx_results[k,2] = "Percentile_5"
  idx_results[k,3] = minObsYear; idx_results[k,4] = min(minObsYear + 40, max(ensemble_year_f$time_year)) 
  idx_results[k,5] = p_05_obs; k = k + 1
  
  idx_results[k,1] = "obs"; idx_results[k,2] = "mean"
  idx_results[k,3] = minObsYear; idx_results[k,4] = min(minObsYear + 40, max(ensemble_year_f$time_year)) 
  idx_results[k,5] = mean_obs ; k = k + 1
  
  
  # HISTORIC GCM:
  
  for (i in 1:ceiling((max(ensemble_year_h$time_year) - minObsYear)/30)){
    bYear = minObsYear + (i - 1) * 30
    eYear = minObsYear  + (i) * 30
    
    aym_b = which(ensemble_day_h$time_year >= bYear)
    aym_b = aym_b[1]
    aym_e = which(ensemble_day_h$time_year < eYear)
    aym_e = tail(aym_e, n=1)
    
    aym = aym_b:aym_e
    
    wf_model[i] = wet_fraction(ensemble_day_h$pr_E1[aym],1)
    
    spellLengths_model = wet_dry_spell_length(ensemble_day_h$pr_E1[aym],1)
    wet_spell_model[i] = spellLengths_model[6]
    dry_spell_model[i] = spellLengths_model[5]
    
    p_05_model[i] = percentile(ensemble_day_h$pr_E1[aym],0.05)
    b_model[i] = bias(data_d$avg_value[ayo],ensemble_day_h$pr_E1[aym])
    mean_model[i] = mean(ensemble_day_h$pr_E1[aym],na.rm=TRUE)
    
    
    idx_results[k,1] = modelName; idx_results[k,2] = "wet_fraction"; 
    idx_results[k,3] = bYear; idx_results[k,4] = eYear ;
    idx_results[k,5] = wf_model[i]; k = k + 1
    
    idx_results[k,1] = modelName; idx_results[k,2] = "wet_spell_length"
    idx_results[k,3] = bYear; idx_results[k,4] = eYear  
    idx_results[k,5] = wet_spell_model[i]; k = k + 1
    
    idx_results[k,1] = modelName; idx_results[k,2] = "dry_spell_length"
    idx_results[k,3] = bYear; idx_results[k,4] = eYear  
    idx_results[k,5] = dry_spell_model[i]; k = k + 1
    
    idx_results[k,1] = modelName; idx_results[k,2] = "Percentile_5"
    idx_results[k,3] = bYear; idx_results[k,4] = eYear  
    idx_results[k,5] = p_05_model[i]; k = k + 1
    
    idx_results[k,1] = modelName; idx_results[k,2] = "bias"
    idx_results[k,3] = bYear; idx_results[k,4] = eYear  
    idx_results[k,5] = b_model[i]; k = k + 1
    
    idx_results[k,1] = modelName; idx_results[k,2] = "mean"
    idx_results[k,3] = bYear; idx_results[k,4] = eYear  
    idx_results[k,5] = mean_model[i]; k = k + 1
  }
  
  #FUTURE GCM
  for (i in 1:ceiling((maxGCMYear - min(ensemble_year_f$time_year))/30)){
    
    bYear = min(ensemble_year_f$time_year) + (i - 1) * 30
    eYear = min(ensemble_year_f$time_year)  + (i) * 30
    
    aym_b = which(ensemble_day_f$time_year >= bYear)
    aym_b = aym_b[1]
    aym_e = which(ensemble_day_f$time_year < eYear)
    aym_e = tail(aym_e, n=1)
    
    aym = aym_b:aym_e
    
    wf_model[i] = wet_fraction(ensemble_day_f$pr_E1[aym],1)
    
    spellLengths_model = wet_dry_spell_length(ensemble_day_f$pr_E1[aym],1)
    wet_spell_model[i] = spellLengths_model[6]
    dry_spell_model[i] = spellLengths_model[5]
    mean_model[i] = mean(ensemble_day_f$pr_E1[aym],na.rm=TRUE)
    
    p_05_model[i] = percentile(ensemble_day_f$pr_E1[aym],0.05)
    b_model[i] = bias(data_d$avg_value[ayo],ensemble_day_f$pr_E1[aym])
    
    idx_results[k,1] = modelName; idx_results[k,2] = "wet_fraction"; 
    idx_results[k,3] = bYear; idx_results[k,4] = eYear ;
    idx_results[k,5] = wf_model[i]; k = k + 1
    
    idx_results[k,1] = modelName; idx_results[k,2] = "wet_spell_length"
    idx_results[k,3] = bYear; idx_results[k,4] = eYear  
    idx_results[k,5] = wet_spell_model[i]; k = k + 1
    
    idx_results[k,1] = modelName; idx_results[k,2] = "dry_spell_length"
    idx_results[k,3] = bYear; idx_results[k,4] = eYear  
    idx_results[k,5] = dry_spell_model[i]; k = k + 1
    
    idx_results[k,1] = modelName; idx_results[k,2] = "Percentile_5"
    idx_results[k,3] = bYear; idx_results[k,4] = eYear  
    idx_results[k,5] = p_05_model[i]; k = k + 1
    
    idx_results[k,1] = modelName; idx_results[k,2] = "bias"
    idx_results[k,3] = bYear; idx_results[k,4] = eYear  
    idx_results[k,5] = b_model[i]; k = k + 1
    
    idx_results[k,1] = modelName; idx_results[k,2] = "mean"
    idx_results[k,3] = bYear; idx_results[k,4] = eYear  
    idx_results[k,5] = mean_model[i]; k = k + 1
  }
  
  write.table(idx_results,paste(fPath,paste(modelName,fConc,varName,"stats.csv",sep="_"),sep=""), sep=",",row.names=FALSE)
}

# ---------------------------------------
# IV. KNN BOOTSRAP FUNCTION 
# ---------------------------------------

knn_bootstrap = function(refDate = c(month = 1, day = 1, year = 1850),
                         main_folder = "F:/GCM_ClimateTool",
                         results_folder = "/Results/test",
                         modelName = "CCSM4",
                         futures = c("rcp85"),
                         varName = "pr",
                         varLabel = "Precipitation (mm)",
                         minObsYear = 1980,
                         maxObsYear = 2016,
                         nearWindow = 15,
                         minGCMYear = 2015, 
                         maxGCMYear = 2040,
                         expNumber = 1,
                         JPmode = "Window", # Options: "Yearly" | "Monthly" | "Window"
                         alignHistYears = TRUE,
                         HistRepro = FALSE) {
## override params
     #   refDate = c(month = df$BaseMonth[rcp_pr_ind], day = df$BaseDay[rcp_pr_ind], year = df$BaseYear[rcp_pr_ind])
     # main_folder = mFolder
     # results_folder = rFolder
     # modelName = model
     # futures = experiment2
     # varName = "pr"
     # varLabel = "Precipitation [mm]"
     # minObsYear = 1980
     # maxObsYear = 2016
     # nearWindow = 15
     # minGCMYear = 2040
     # maxGCMYear = 2060
     # expNumber = 1
     # JPmode = "Window"
     # alignHistYears = TRUE
     # HistRepro = FALSE

  bootstrapLength <- maxGCMYear - minGCMYear + 1
  YearlyThreshReCalc <- TRUE # experiment trigger for only the Yearly JPmode scheme for activating re-calc of wet-state thresholds for each daily brick
  MonthlyThreshPar <- TRUE
  
  # LOADS HISTORICAL GCM DATA
  # fPath <- paste(main_folder,results_folder,'/',modelName,'/',sep="")
  fPath <- paste(main_folder,results_folder,'/',sep="")
  fConc = futures
  
  fName = paste(fPath,paste(modelName,'historical',varName,"ensemble_day.Rda",sep="_"),sep="")
  load(fName,verbose=TRUE)    
  
  ensemble_day$month_idx = (ensemble_day$time_year - refDate[3]) * 12 + ensemble_day$time_month
  ensemble_day_h <- ensemble_day[ensemble_day$time_year >= minObsYear & ensemble_day$time_year <= maxObsYear,]
  ensemble_day_h$julian_day <- julian(ensemble_day_h$time_month,ensemble_day_h$time_day,1850,origin = c(month = 1, day = 1, 1850))
  
  fName = paste(fPath,paste(modelName,'historical',varName,"ensemble_year.Rda",sep="_"),sep="")
  load(fName,verbose=TRUE)   
  ensemble_year_h <- ensemble_year[ensemble_year$time_year >= minObsYear & ensemble_year$time_year <= maxObsYear,]
  
  fName = paste(fPath,paste(modelName,'historical',varName,"ensemble_month.Rda",sep="_"),sep="")
  load(fName)   
  ensemble_month_h <- ensemble_month[ensemble_month$time_year >= minObsYear & ensemble_month$time_year <= maxObsYear,]
  
  # LOADS FUTURE GCM DATA
  fConc = futures
  
  fName = paste(fPath,paste(modelName,fConc,varName,"ensemble_day.Rda",sep="_"),sep="")
  load(fName,verbose=TRUE)    
  ensemble_day$month_idx = (ensemble_day$time_year - refDate[3]) * 12 + ensemble_day$time_month
  ensemble_day_f <- ensemble_day[ensemble_day$time_year >= minGCMYear & ensemble_day$time_year <= maxGCMYear,]
  ensemble_day_f$julian_day <- julian(ensemble_day_f$time_month,ensemble_day_f$time_day,1850,origin = c(month = 1, day = 1, 1850))

  fName = paste(fPath,paste(modelName,fConc,varName,"ensemble_year.Rda",sep="_"),sep="")
  load(fName,verbose=TRUE)   
  ensemble_year_f <- ensemble_year[ensemble_year$time_year >= minGCMYear & ensemble_year$time_year <= maxGCMYear,]

  fName = paste(fPath,paste(modelName,fConc,varName,"ensemble_month.Rda",sep="_"),sep="")
  load(fName)   
  ensemble_month_f <- ensemble_month[ensemble_month$time_year >= minGCMYear & ensemble_month$time_year <= maxGCMYear,]
  
  # observed historic

  fName = paste(fPath,paste(modelName, "observed", varName,"day.Rda",sep="_"),sep="")
  load(fName,verbose=TRUE)
  fName = paste(fPath,paste(modelName, "observed", varName,"year.Rda",sep="_"),sep="")
  load(fName,verbose=TRUE)
  fName = paste(fPath,paste(modelName, "observed", varName,"month.Rda",sep="_"),sep="")
  load(fName,verbose=TRUE)
  
  #------------------------------------------------------------------------------------------------------------
  # Filters Historic Years of GCM and Observed Data to Coincide, and Future GCM by minGCMYear, Bootstrap Length
  #------------------------------------------------------------------------------------------------------------
  
  if(alignHistYears == TRUE)
  {
    maxObsYear <- min(max(data_d$Year), max(ensemble_day_h$time_year))
  } else {
    maxObsYear <- max(data_d$Year)
  }
  
  maxGCMfYear <- minGCMYear + bootstrapLength - 1
  
  data_d <- data_d[data_d$Year >= minObsYear & data_d$Year <= maxObsYear,]
  data_m <- data_m[data_m$Year >= minObsYear & data_m$Year <= maxObsYear,]
  data_y <- data_y[data_y$Year >= minObsYear & data_y$Year <= maxObsYear,]
  
  ensemble_day_h <- ensemble_day_h[ensemble_day_h$time_year >= minObsYear & ensemble_day_h$time_year <= maxObsYear,]
  ensemble_month_h <- ensemble_month_h[ensemble_month_h$time_year >= minObsYear & ensemble_month_h$time_year <= maxObsYear,]
  ensemble_year_h <- ensemble_year_h[ensemble_year_h$time_year >= minObsYear & ensemble_year_h$time_year <= maxObsYear,]
  
  ensemble_day_f <- ensemble_day_f[ensemble_day_f$time_year >= minGCMYear & ensemble_day_f$time_year <= maxGCMfYear,]
  ensemble_month_f <- ensemble_month_f[ensemble_month_f$time_year >= minGCMYear & ensemble_month_f$time_year <= maxGCMfYear,]
  ensemble_year_f <- ensemble_year_f[ensemble_year_f$time_year >= minGCMYear & ensemble_year_f$time_year <= maxGCMfYear,]
  
  #-----------------------------------------------------------------
  # Bootstraping daily
  #-----------------------------------------------------------------
  
  # a. Puts the data in a more efficient data_frame
  a = c("timeStamp","Year","Month","Day","avg_value")  
  data_int <- subset(data_d[a], data_d$Year >= minObsYear)  
  
  # adds a new column to data frame with the value of the next_day
  nData = dim(data_int)[1]                               
  nextDay = data_int$avg_value[2:nData]
  data_int$avg_value_f = c(nextDay,0)                 
  data_int$julian_day <- julian(data_int$Month,data_int$Day,1970,origin = c(month = 1, day = 1, 1970))
  
  # b. Produces the CURRENT CLIMATE transition matrix: Probability(wet_state_t|wet_state_t-1) using the markov chain transition model informed by the GCM
  
  Threshold = 1  # dry/wet threshold [mm] used for the yearly scheme and for monthly scheme with constant thresholds
  extreme_wet_percentile = 0.75
  
  if (JPmode == "Yearly")
  {
    extreme_p = percentile(data_int$avg_value,1 - extreme_wet_percentile)  # wet/extremely wet thresh for observed historical data   
    extreme_l = max(Threshold,percentile(data_int$avg_value,0.99)) # dry/wet thresh for obs hist data
    
    extreme_p_gcm = percentile(ensemble_day$pr_E1,1 - extreme_wet_percentile)  # wet/extremely wet thresh for gcm historical data   
    extreme_l_gcm = max(Threshold,percentile(ensemble_day$pr_E1,0.99)) # dry/wet thresh for gcm hist data
    
    gcm_ref_idx = which(ensemble_day$time_year > minObsYear)
    
    jp = tri_state_joint_probability(data_int$avg_value,extreme_l,extreme_p) 
    
    jp_gcm_avg = tri_state_joint_probability(ensemble_day$pr_E1[gcm_ref_idx], extreme_l_gcm, extreme_p_gcm)
    
  } else if (JPmode == "Monthly"){
    
    gcm_ref_idx = which(ensemble_day$time_year > minObsYear)
    
    extreme_p = sapply(1:12, function(x) percentile(data_int$avg_value[data_int$Month == x],1 - extreme_wet_percentile))  # extremely wet    
    extreme_l = sapply(1:12, function(x) max(Threshold,percentile(data_int$avg_value[data_int$Month == x],0.99))) # wet
    
    jp = tri_state_joint_prob_monthly(data_d[,c("Year", "Month", "avg_value")],extreme_l,extreme_p, monthly_thresholds = MonthlyThreshPar)   
    
    extreme_p_gcm = sapply(1:12, function(x) percentile(ensemble_day$pr_E1[ensemble_day$time_month == x],1 - extreme_wet_percentile))  # extremely wet    
    extreme_l_gcm = sapply(1:12, function(x) max(Threshold,percentile(ensemble_day$pr_E1[ensemble_day$time_month == x],0.99))) # wet
    
    jp_gcm_avg = tri_state_joint_prob_monthly(ensemble_day[gcm_ref_idx,c("time_year", "time_month", "pr_E1")],extreme_l_gcm,extreme_p_gcm, monthly_thresholds = MonthlyThreshPar)
    
  } else if (JPmode == "Window") {
    
    extreme_p = sapply(0:364, function(x) {
      jdays <- (x - nearWindow):(x + nearWindow)
      jdays[jdays < 0] <- jdays[jdays < 0] + 365
      percentile(data_int$avg_value[data_int$julian_day %in% jdays],1 - extreme_wet_percentile)  # wet/extremely wet thresh for observed historical data  
    }) # extremely wet    
    extreme_l = sapply(0:364, function(x) {
      jdays <- (x - nearWindow):(x + nearWindow)
      jdays[jdays < 0] <- jdays[jdays < 0] + 365
      val <- max(Threshold,percentile(data_int$avg_value[data_int$julian_day %in% jdays],0.99)) # dry/wet thresh for obs hist data
      
      if(val >= extreme_p[x+1] | val == 0)
      {
        val <- extreme_p[x+1]/2
      }
      val
    }) # wet
    
    jp <- matrix(NA, ncol = 3*365, nrow = 3)
    
    for(i in 0:364)
    {
      jdays <- (i - nearWindow):(i + nearWindow)
      jdays[jdays < 0] <- jdays[jdays < 0] + 365
      
      jp[,(i*3 + 1):(i*3 + 3)] <- tri_state_joint_probability(data_int$avg_value[data_int$julian_day %in% jdays],extreme_l[i+1],extreme_p[i+1])
    }
    
    gcm_ref_idx = which(ensemble_day_h$time_year > minObsYear)
    
    ensemble_day_filt <- ensemble_day_h
    
    extreme_p_gcm = sapply(0:364, function(x) {
      jdays <- (x- nearWindow):(x + nearWindow)
      jdays[jdays < 0] <- jdays[jdays < 0] + 365
      percentile(ensemble_day_filt$pr_E1[ensemble_day_filt$julian_day %in% jdays],1 - extreme_wet_percentile) # wet/extremely wet thresh for gcm historical data   
    }) 
    # extremely wet    
    extreme_l_gcm = sapply(0:364, function(x) {
      jdays <- (x - nearWindow):(x + nearWindow)
      jdays[jdays < 0] <- jdays[jdays < 0] + 365
      val <- max(Threshold,percentile(ensemble_day_filt$pr_E1[ensemble_day_filt$julian_day %in% jdays],0.99)) # dry/wet thresh for gcm hist data
      
      if(val >= extreme_p_gcm[x+1] | val == 0)
      {
        val <- extreme_p_gcm[x+1]/2
      }
      val
    }) # wet
    
    jp_gcm_avg <- matrix(NA, ncol = 3*365, nrow = 3)
    for(i in 0:364)
    {
      jdays <- (i - nearWindow):(i + nearWindow)
      jdays[jdays < 0] <- jdays[jdays < 0] + 365
      
      jp_gcm_avg[,(i*3 + 1):(i*3 + 3)] <- tri_state_joint_probability(ensemble_day_filt$pr_E1[ensemble_day_filt$julian_day %in% jdays],extreme_l_gcm[i+1],extreme_p_gcm[i+1])
    }
  }
  
  simYear <-NULL
  simMonth <-NULL
  simDay <-NULL
  
  Year <- NULL
  Month <- NULL
  Day <- NULL
  bs_value <- NULL
  bs_wet_state <- NULL
  bs_brick_size <- NULL
  
  bs_result = data.frame(simYear,simMonth,simDay,Year,Month,Day,bs_value, bs_wet_state, bs_brick_size) # data.frame with the bootstrapping results
  nYears = bootstrapLength
  
  for (y in 1:nYears) {
    # Generates the JOINT PROBABILITY SHIFT from a data_set informed by the GCM. Sum(jp_shift) must be equal to 0)
    print(c(y,"of",nYears))
    if(JPmode == "Yearly")
    {
      gcm_ref_idx_f = which(ensemble_day_f$time_year == minGCMYear + y - 1) 
      jp_gcm_f = tri_state_joint_probability(ensemble_day_f$pr_E1[gcm_ref_idx_f],extreme_l_gcm,extreme_p_gcm)
      
      # Simulation SHIFTED JOINT PROBABILITY from a data_set informed by the GCM
      jp_shift = jp_gcm_f - jp_gcm_avg
      jp_sim = jp + jp_shift
      jp_sim = replace(jp_sim,jp_sim < 0,0)
      
      if(HistRepro == TRUE)
      {
        jp_sim <- jp # for assessing historical re-production performance
      }
      
    } else if (JPmode == "Monthly") 
    {
      gcm_ref_idx_f = which(ensemble_day_f$time_year == minGCMYear + y - 1) 
      jp_gcm_f = tri_state_joint_prob_monthly(ensemble_day_f[gcm_ref_idx_f, c("time_year", "time_month", "pr_E1")],extreme_l_gcm,extreme_p_gcm, monthly_thresholds = MonthlyThreshPar)
      
      # Simulation SHIFTED JOINT PROBABILITY from a data_set informed by the GCM
      jp_shift = jp_gcm_f - jp_gcm_avg
      jp_sim = jp + jp_shift
      jp_sim = replace(jp_sim,jp_sim < 0,0)
      
      if(HistRepro == TRUE)
      {
        jp_sim <- jp # for assessing historical re-production performance
      }
    }
    
    if(JPmode == "Yearly")
    {
      jp_sim = jp_sim / sum(jp_sim)
      
      # Calcultates the Markov TRANSITION MATRIX from a data_set: P(future|present).
      tr_sim = mat.or.vec(3, 3)
      
      tr_sim[1,] = jp_sim[1,]/sum(jp_sim[1,])
      tr_sim[2,] = jp_sim[2,]/sum(jp_sim[2,])
      tr_sim[3,] = jp_sim[3,]/sum(jp_sim[3,])
      
      tr_sim[is.nan(tr_sim)] <- 0
      
      tr_sim_mar <- NULL
      tr_sim_mar[1] = sum(jp_sim[,1])  # marginal probability of wet state 1: dry
      tr_sim_mar[2] = sum(jp_sim[,2])  #                         wet state 2: wet
      tr_sim_mar[3] = sum(jp_sim[,3])  #                         wet state 3: extremely wet
      
    } else if (JPmode == "Monthly"){
      
      tr_sim = as.data.frame(mat.or.vec(3, 3*12))
      row.names(tr_sim) <- row.names(jp_sim)
      colnames(tr_sim) <- colnames(jp_sim)
      
      tr_sim_mar <- NULL
      
      for(m in 1:12)
      {
        jp_sim[,(m*3 - 2):(m*3)] <- round(jp_sim[,(m*3 - 2):(m*3)]/sum(jp_sim[,(m*3 - 2):(m*3)]),5)
        
        # Calcultates the Markov TRANSITION MATRIX from a data_set: P(future|present).
        
        tr_sim[1,(m*3 - 2):(m*3)] <- jp_sim[1,(m*3 - 2):(m*3)]/sum(jp_sim[1,(m*3 - 2):(m*3)])
        tr_sim[2,(m*3 - 2):(m*3)] <- jp_sim[2,(m*3 - 2):(m*3)]/sum(jp_sim[2,(m*3 - 2):(m*3)])
        tr_sim[3,(m*3 - 2):(m*3)] <- jp_sim[3,(m*3 - 2):(m*3)]/sum(jp_sim[3,(m*3 - 2):(m*3)])
        
        tr_sim[tr_sim == "NaN"] <- 0 # removes NaNs resulting from a 0,0,0 jp row (i.e. Feb has no days above p95)
        
        tr_sim_mar[m*3 - 2] <- sum(jp_sim[,m*3 - 2])  # marginal probability of wet state 1: dry
        tr_sim_mar[m*3 - 1] <- sum(jp_sim[,m*3 - 1])  #                         wet state 2: wet
        tr_sim_mar[m*3] <- sum(jp_sim[,m*3])  #                         wet state 3: extremely wet
      }
    } 
    
    simYear <-1
    simMonth <-1
    simDay <-1
    Year <- 1
    Month <- 1
    Day <- 1
    bs_value <- 1
    bs_wet_state <- 1
    bs_brick_size <- 1
    
    bs_sequence = data.frame(simYear,simMonth,simDay,Year,Month,Day,bs_value, bs_wet_state, bs_brick_size)  # initializes a data.frame to store the results for each year
    
    for (i in 0:364) {   # 364 days
      # a. makes the "brick": the subset of "local" data 
      if (i < nearWindow) {
        brick <- subset(data_int, (data_int$julian_day <= i + nearWindow) | (data_int$julian_day > 364 + i - nearWindow))
        ens_brick <- subset(ensemble_day_filt, (ensemble_day_filt$julian_day <= i + nearWindow) | (ensemble_day_filt$julian_day > 364 + i - nearWindow))
        ens_brick_f <- subset(ensemble_day_f, (ensemble_day_f$julian_day <= i + nearWindow) | (ensemble_day_f$julian_day > 364 + i - nearWindow))
        
      } else if(i > 364 - nearWindow)
      {
        brick <- subset(data_int, (data_int$julian_day >= i - nearWindow) | (data_int$julian_day < i - 364 + nearWindow))
        ens_brick <- subset(ensemble_day_filt, (ensemble_day_filt$julian_day >= i - nearWindow) | (ensemble_day_filt$julian_day < i - 364 + nearWindow))
        ens_brick_f <- subset(ensemble_day_f, (ensemble_day_f$julian_day >= i - nearWindow) | (ensemble_day_f$julian_day < i - 364 + nearWindow))
      } else {
        brick <- subset(data_int, (data_int$julian_day <= i + nearWindow) & (data_int$julian_day >= i - nearWindow))
        ens_brick <- subset(ensemble_day_filt, (ensemble_day_filt$julian_day <= i + nearWindow) & (ensemble_day_filt$julian_day >= i - nearWindow))
        ens_brick_f <- subset(ensemble_day_f, (ensemble_day_f$julian_day <= i + nearWindow) & (ensemble_day_f$julian_day >= i - nearWindow))
      }
      
      if(JPmode == "Monthly")
      {
        tr_thresholds_var <- c(extreme_l ,extreme_p)
        mons <- unique(brick$Month)
        mon_weights <- sapply(mons, function(x) sum(brick$Month == x)/dim(brick)[1])
        l_thresh_mweight <- sum(sapply(1:length(mons), function(x) tr_thresholds_var[mons[x]]*mon_weights[x]))
        h_thresh_mweight <- sum(sapply(1:length(mons), function(x) tr_thresholds_var[mons[x] + 12]*mon_weights[x]))
      }
      
      if(JPmode == "Window")
      {
        brickdays <- (i - nearWindow):(i + nearWindow)
        brickdays[brickdays < 0] <- brickdays[brickdays < 0] + 365
        brickdays[brickdays > 364] <- brickdays[brickdays > 364] - 365
        
        if(i < nearWindow & y > 1)
        {
          gcm_ref_idx_f = c(which(ensemble_day_f$time_year == minGCMYear + y - 2 & ensemble_day_f$julian_day %in% brickdays[brickdays > 364 - nearWindow]), which(ensemble_day_f$time_year == minGCMYear + y - 1 & ensemble_day_f$julian_day %in% brickdays[brickdays <= nearWindow]))
          
        } else if(i > 364 - nearWindow & y < nYears)
        {
          gcm_ref_idx_f = c(which(ensemble_day_f$time_year == minGCMYear + y - 1 & ensemble_day_f$julian_day %in% brickdays[brickdays > 364 - nearWindow]), which(ensemble_day_f$time_year == minGCMYear + y & ensemble_day_f$julian_day %in% brickdays[brickdays <= nearWindow]))
          
        } else {
          
          gcm_ref_idx_f = which(ensemble_day_f$time_year == minGCMYear + y - 1 & ensemble_day_f$julian_day %in% brickdays)
        }
        
        jp_gcm_f = tri_state_joint_probability(ensemble_day_f$pr_E1[gcm_ref_idx_f],extreme_l_gcm[i+1],extreme_p_gcm[i+1])
        
        jp_shift = jp_gcm_f - jp_gcm_avg[,(i*3 + 1):(i*3 + 3)]
        
        gcm_f_avg <- mean(ensemble_day_f$pr_E1[gcm_ref_idx_f])
        
        # Simulation SHIFTED JOINT PROBABILITY from a data_set informed by the GCM
        jp_sim = jp[,(i*3 + 1):(i*3 + 3)] + jp_shift
        jp_sim = replace(jp_sim,jp_sim < 0,0)
        
        if(HistRepro == TRUE)
        {
          jp_sim <- jp[,(i*3 + 1):(i*3 + 3)] # for assessing historical re-production performance
        }
        
        jp_sim = jp_sim / sum(jp_sim)
        
        # Calcultates the Markov TRANSITION MATRIX from a data_set: P(future|present).
        tr_sim = mat.or.vec(3, 3)
        
        tr_sim[1,] = jp_sim[1,]/sum(jp_sim[1,])
        tr_sim[2,] = jp_sim[2,]/sum(jp_sim[2,])
        tr_sim[3,] = jp_sim[3,]/sum(jp_sim[3,])
        
        tr_sim[is.nan(tr_sim)] <- 0
        
        tr_sim_mar <- NULL
        tr_sim_mar[1] = sum(jp_sim[,1])  # marginal probability of wet state 1: dry
        tr_sim_mar[2] = sum(jp_sim[,2])  #                         wet state 2: wet
        tr_sim_mar[3] = sum(jp_sim[,3])  #                         wet state 3: extremely wet
        
      } # end JPmode == "Window" logic
      
      if(YearlyThreshReCalc == TRUE & JPmode == "Yearly")
      {
        extreme_p = percentile(brick$avg_value,1 - extreme_wet_percentile)  # wet/extremely wet thresh for observed historical data   
        extreme_l = max(Threshold,percentile(brick$avg_value,0.99))
        tr_thresholds_var <- c(extreme_l,extreme_p)
      } else if (JPmode == "Window")
      {
        tr_thresholds_var <- c(extreme_l[i+1],extreme_p[i+1])
      } else {
        tr_thresholds_var <- c(extreme_l,extreme_p)# the low and high thresholds for the year or brick, and the 12 monthly low thresholds and 12 monthly high thresholds in the monthly scheme (1:12 and 13:24, resp.)
      }
      
      
      # corrects an error, in the case of extremely dry historical spells in the near window, that results in extremeWet PR therhsold to be 0.
      if (tr_thresholds_var[2] <= tr_thresholds_var[1]){
        tr_thresholds_var[2] = tr_thresholds_var[1] + 1 
      }
      
      if (i==0) {   # required to initialize the synthetic generator. is based on the marginal distribution
        
        a = runif(1) # a random UNIFORM, i = 0 will always be Jan 1, so no need to alter based on month, tr_sim_mar[1] is for Jan
        if(JPmode == "Yearly" | JPmode == "Window")
        {
          if (a <= tr_sim_mar[1]) {
            c_wet_state = 1                # dry
            new_brick <- subset(brick, brick$avg_value <= tr_thresholds_var[1]) 
          } else if (a > tr_sim_mar[2] + tr_sim_mar[1]){  # extremely wet
            c_wet_state = 3     
            new_brick <- subset(brick, brick$avg_value > tr_thresholds_var[2])        
          } else {
            c_wet_state = 2                # wet
            new_brick <- subset(brick, brick$avg_value > tr_thresholds_var[1] & brick$avg_value <= tr_thresholds_var[2]) 
          }
        } else if(JPmode == "Monthly")
        {
          
          if (a <= tr_sim_mar[1]) {
            c_wet_state = 1                # dry
            new_brick <- subset(brick, brick$avg_value <= l_thresh_mweight) 
          } else if (a > tr_sim_mar[2] + tr_sim_mar[1]){  # extremely wet
            c_wet_state = 3     
            new_brick <- subset(brick, brick$avg_value > h_thresh_mweight)        
          } else {
            c_wet_state = 2                # wet
            new_brick <- subset(brick, brick$avg_value > l_thresh_mweight & brick$avg_value <= h_thresh_mweight) 
          }
        }
        
        
        #first step: no Nearest neighbors
        nDat <- dim(new_brick)[1] 
        idx <- ceiling(runif(1) * nDat)
        
##      dffilt <- df[df$Model == model & df$Experiment == " hist-rcp85" & df$Variable == varName,]
        dffilt <- df      #[df$Model == model & df$Experiment == " hist-rcp85" & df$Variable == varName,]
        
        #HistBaseDate <- c(month = dffilt$BaseMonth, day = dffilt$BaseDay, year = dffilt$BaseYear)
        HistBaseDate <- c(month = 1, day = 1, year = 2006) # HA temp test
        
        dum_date = julian(1,1,minGCMYear + y - 1, origin = HistBaseDate)
        simDate = month.day.year(dum_date + i, origin = HistBaseDate)
        
        bs_sequence$simYear[i+1] <- simDate$year
        bs_sequence$simMonth[i+1] <- simDate$month
        bs_sequence$simDay[i+1] <- simDate$day      
        
        bs_sequence$Year[i+1] <- new_brick$Year[idx]      # i+1 because the julian_day indexing starts at 0 
        bs_sequence$Month[i+1] <- new_brick$Month[idx]
        bs_sequence$Day[i+1] <- new_brick$Day[idx]
        bs_sequence$bs_value[i+1] <- new_brick$avg_value[idx]  
        bs_sequence$bs_brick_size[i+1] <- dim(new_brick)[1]
        bs_sequence$bs_wet_state[i+1] <- c_wet_state 
        
      } else { # The rest of the days
        
        # rolls the UNIFORM dice to determine future wet_state given current wet_state AND
        # filters the brick to the dates with BOTH current and future wet states 
        
        if(JPmode == "Yearly" | JPmode == "Window")
        {
          if (c_wet_state == 1) {                                               # dry
            new_brick <- subset(brick, brick$avg_value <= tr_thresholds_var[1])        
          } else if (c_wet_state == 2) {                                         # wet
            new_brick <- subset(brick, brick$avg_value > tr_thresholds_var[1] & brick$avg_value <= tr_thresholds_var[2]) 
          } else {                                           # extremely wet
            new_brick <- subset(brick, brick$avg_value > tr_thresholds_var[2])          
          }
        } else if(JPmode == "Monthly")
        {
          if (c_wet_state == 1) {                                               # dry
            new_brick <- subset(brick, brick$avg_value <= l_thresh_mweight)        
          } else if (c_wet_state == 2) {                                         # wet
            new_brick <- subset(brick, brick$avg_value > l_thresh_mweight & brick$avg_value <= h_thresh_mweight) 
          } else {                                           # extremely wet
            new_brick <- subset(brick, brick$avg_value > h_thresh_mweight)          
          }
        }
        
        if(YearlyThreshReCalc == TRUE & JPmode == "Yearly")
        {
          extreme_p = percentile(new_brick$avg_value_f, 1 - extreme_wet_percentile)    # recalculates the wet, extremely wet boudaries for local data 

          if(length(extreme_p) == 0){      #sometimes, the number of days that fulfill the condition is very small < 20, to calculate correctly the percentile
            tr_thresholds_var = c(Threshold, max(new_brick$avg_value_f)) 
            #print(paste("i: ",i,",", tr_thresholds_var,"Fixed!"))
          } else {
            if(extreme_p > Threshold) {
              tr_thresholds_var = c(Threshold,extreme_p) 
            }else{
              tr_thresholds_var = c(Threshold, max(new_brick$avg_value_f)) 
            }
            #print(paste("i: ",i,",", tr_thresholds_var))
          }
        }
        
        a = runif(1)
        
        if(JPmode == "Yearly" | JPmode == "Window")
        {
          if (a <= tr_sim[c_wet_state,1]) {    
            f_wet_state = 1 
            new_brick_f <- subset(new_brick, new_brick$avg_value_f <= tr_thresholds_var[1])   
            
            # during very wet spells, dryest days are even higher that 0.3 (Threshold). This codes fixes that
            if (dim(new_brick_f)[1] == 0) {
              new_brick_f <- new_brick[which.min(new_brick$avg_value_f),] 
            }
            
          } else if (a > tr_sim[c_wet_state,1] + tr_sim[c_wet_state,2]){  
            f_wet_state = 3     
            new_brick_f <- subset(new_brick, new_brick$avg_value_f > tr_thresholds_var[2])
            # during very dry spells, wettest days don't show up in the brick. This codes fixes that
            if (dim(new_brick_f)[1] == 0) {
              new_brick_f <- new_brick[which.max(new_brick$avg_value_f),]  
            }
            
          } else {
            f_wet_state = 2                                               
            new_brick_f <- subset(new_brick, new_brick$avg_value_f > tr_thresholds_var[1] & new_brick$avg_value_f < tr_thresholds_var[2])
            if (dim(new_brick_f)[1] == 0) {
              new_brick_f <- new_brick[which.min(abs(new_brick$avg_value_f - mean(new_brick$avg_value_f))),]  
            }
          }
        } else if(JPmode == "Monthly") {
          
          mons <- unique(new_brick$Month)
          mon_weights <- sapply(mons, function(x) sum(new_brick$Month == x)/dim(new_brick)[1])
          
          if (a <= sum(sapply(1:length(mons), function(x) tr_sim[c_wet_state, mons[x]*3 - 2]*mon_weights[x]))) {    
            f_wet_state = 1 
            new_brick_f <- subset(new_brick, new_brick$avg_value_f <= l_thresh_mweight)   
            # during very wet spells, dryest days are even higher that 0.3 (Threshold). This codes fixes that
            if (dim(new_brick_f)[1] == 0) {
              minval = min(new_brick$avg_value_f) 
              new_brick_f <- subset(new_brick, new_brick$avg_value_f <= minval)  
            }
            
          } else if (a > sum(sapply(1:length(mons), function(x) tr_sim[c_wet_state, mons[x]*3 - 2]*mon_weights[x])) + sum(sapply(1:length(mons), function(x) tr_sim[c_wet_state, mons[x]*3 - 1]*mon_weights[x]))){  
            f_wet_state = 3     
            new_brick_f <- subset(new_brick, new_brick$avg_value_f >= h_thresh_mweight)
            # during very dry spells, wetest days don't show up in the brick. This codes fixes that
            if (dim(new_brick_f)[1] == 0) {
              minval = max(new_brick$avg_value_f) 
              new_brick_f <- subset(new_brick, new_brick$avg_value_f >= minval)  
            }
            
          } else {
            f_wet_state = 2                                               
            new_brick_f <- subset(new_brick, new_brick$avg_value_f > l_thresh_mweight & new_brick$avg_value_f <= h_thresh_mweight) 
            if (dim(new_brick_f)[1] == 0) {
              minval = max(new_brick$avg_value_f) 
              new_brick_f <- subset(new_brick, new_brick$avg_value_f >= minval)  
            }
          }
        } # end monthly new_brick_f if routine
        
        #--------------------------------------------------------------------------      
        # nearest neighbors
        #--------------------------------------------------------------------------
        
        new_brick_f$distance = abs(new_brick_f$avg_value - (bs_sequence$bs_value[i]))
        new_brick_f$weight = 1 - (new_brick_f$distance / (max(new_brick_f$distance) + 1)) 
        new_brick_f$prob_acum =  cumsum(new_brick_f$weight) / sum(new_brick_f$weight)
        
        # rolls the dice to choose a given neighbor
        a = runif(1)
        idx <- which(new_brick_f$prob_acum > a)[1] 
        
        if (nrow(new_brick_f) > 0) {
        
          f_date = month.day.year(if(leap.year(new_brick_f$Year[idx]) == TRUE){new_brick_f$julian_day[idx] + 1} else {new_brick_f$julian_day[idx]} , origin = c(month = 1, day = 1, new_brick_f$Year[idx]))
          
          bs_sequence[i+1,] = bs_sequence[i,]
          
          dum_date = julian(1,1,minGCMYear + y - 1, origin = HistBaseDate)
          if(leap.year(minGCMYear + y - 1) == TRUE & i >= 59) # keeps Feb 29 out of output series since day loop (i) only has 365 days
          {
            simDate = month.day.year(dum_date + i + 1, origin = HistBaseDate)
          } else {
            simDate = month.day.year(dum_date + i, origin = HistBaseDate)
          }
        } else {
          simDate = month.day.year(dum_date, origin = HistBaseDate) # ********* provisional. para ver si el problema el la ocurrencia de ciertas trancisiones en el registro hist?rico **********
        }

        bs_sequence$simYear[i+1] <- simDate$year
        bs_sequence$simMonth[i+1] <- simDate$month
        bs_sequence$simDay[i+1] <- simDate$day 
        
        bs_sequence$Year[i+1] <- f_date$year
        bs_sequence$Month[i+1] <- f_date$month
        bs_sequence$Day[i+1] <- f_date$day
        bs_sequence$bs_value[i+1] <- new_brick_f$avg_value_f[idx]  
        bs_sequence$bs_brick_size[i+1] <- dim(new_brick_f)[1]
        bs_sequence$bs_wet_state[i+1] <- f_wet_state 
        
        # updates the current wet state variable
        c_wet_state = f_wet_state  
      }
    } # next day, i
    bs_result = rbind(bs_result,bs_sequence) 
  } # next year, y
  
  
#------------------------------------------------------------------------------------------#
  #### EXTREMES CORRECTION OF BS OUTPUT #### 
  if(HistRepro != TRUE)
  {
    ######################################################################
    ### Future Boot-strapped Data ###
    bs_result_orig <- bs_result # preserve a copy of the observations
    
    exThreshold <- 2
    repeat {
      GPDbs <- try(gpd(sort(bs_result$bs_value),threshold=exThreshold),silent = FALSE)
      if (substr(GPDbs[1],1,5) != "Error" && GPDbs$n.exceed < 200) { break }
      exThreshold <- exThreshold + 2
    } 
    
    bs_pctl <- val2pctl(exThreshold, bs_result$bs_value)
    bs_extr <- rgpd(10000, GPDbs$par.ests[1], GPDbs$threshold, GPDbs$par.ests[2])
    
    idx.ext.bs <- which(bs_result$bs_value >= GPDbs$threshold) #the index of the extreme precip of the observations
    
    GPDbs <- gpd(sort(bs_result$bs_value),threshold=exThreshold)
    bs_extr <- rgpd(10000,GPDbs$par.ests[1], GPDbs$threshold, GPDbs$par.ests[2])
    
    ### For historic GCM
    pr_ensemble_day <- ensemble_day[is.finite(ensemble_day$pr_E1) & ensemble_day$time_year >= minObsYear & ensemble_day$time_year < minGCMYear,]
    exThreshold <- percentile(pr_ensemble_day$pr_E1, 1 - bs_pctl)
    
    GPDhist <- gpd(sort(pr_ensemble_day$pr_E1),threshold=exThreshold)
    
    gcmh_pctl <- bs_pctl
    gcmh_extr <- rgpd(10000,GPDhist$par.ests[1], GPDhist$threshold, GPDhist$par.ests[2])
    
    #########################################################################################
    
    ens  <- subset(ensemble_day_f[,c(2,3,5)],time_year >= minGCMYear & time_year <= minGCMYear + bootstrapLength)
    p.ens <- ens[ens$pr_E1 >1,]
    ### For Future GCM
    exThreshold <- percentile(p.ens$pr_E1, 1 - bs_pctl)
    GPDens <- gpd(sort(p.ens$pr_E1),threshold=exThreshold)
    
    gcmf_pctl <- bs_pctl
    gcmf_extr <- rgpd(10000, GPDens$par.ests[1], GPDens$threshold, GPDens$par.ests[2])
    
    #### update the parameters of the Synth Bootstrapped data. #############################################
    GPDbeta  <- GPDbs$par.ests[2]* (GPDens$par.ests[2]/GPDhist$par.ests[2])
    GPDalpha <-  GPDbs$par.ests[1] + (GPDens$par.ests[1]-GPDhist$par.ests[1]) 
    GPDThresh <-  GPDbs$threshold * (GPDens$threshold / GPDhist$threshold)
    
    xbs_new <- rgpd(10000, GPDalpha, GPDThresh, GPDbeta) #find the new values using the new extreme values
    
    bsorig <- sort(bs_extr)
    bsnew <- sort(xbs_new)
    
    bs_extr_new <- bs_result$bs_value[bs_result$bs_value >= GPDbs$threshold]
    bs_extr_new <- sapply(1:length(bs_extr_new), function(a) bsnew[which.min(abs(bsorig - bs_extr_new[a]))])
    
    xbs_new <- sort(xbs_new)
    xbs_new <- sapply(1:100, function(a) { 
      val <- a*length(bs_extr)/100 
      linapprox(val, 1:length(xbs_new), xbs_new) })
    
    xbs <- sort(bs_extr)
    xbs <- sapply(1:100, function(a) { 
      val <- a*length(xbs)/100 
      linapprox(val, 1:length(xbs), xbs) })
    
    gcmh <- sort(gcmh_extr)
    gcmh <- sapply(1:100, function(a) linapprox(a*(length(gcmh)/100), 1:length(gcmh), gcmh))
    
    gcmf <- sort(gcmf_extr)
    gcmf <- sapply(1:100, function(a) linapprox(a*(length(gcmf)/100), 1:length(gcmf), gcmf))
    
    bsold <- bs_result
    #subsitute those values into the observations.. now have a new set of obs.
    bs_result$bs_value[idx.ext.bs] <- bs_extr_new
    bs_Ratios <- bs_result$bs_value[idx.ext.bs]/bsold$bs_value[idx.ext.bs]
    
    bs_excorr <- c(bs_pctl, (GPDens$par.ests[1]-GPDhist$par.ests[1]), (GPDens$threshold / GPDhist$threshold), (GPDens$par.ests[2]/GPDhist$par.ests[2]))
    bs_excorr <- as.data.frame(t(bs_excorr))
    colnames(bs_excorr) <- c("bs_pctl", "GPDalphaDiff", "GPDThreshRatio", "GPDbetaRatio")
    
    save(bs_excorr,file=paste(fPath,paste(modelName,"bs", "pr",futures,"day_ExtremesCorrected_key.Rda",sep="_"),sep=""))
  }
  
  # saves a CSV
  write.table(bs_result,paste(fPath,paste(modelName,fConc,varName,expNumber,"bootstrap.csv",sep="_"),sep=""), sep=",",row.names=FALSE)
  
  # Saves the data frame with the bootstraped signal
  save(bs_result,file=paste(fPath,paste(modelName,fConc,varName,expNumber,"bootstrap.Rda",sep="_"),sep=""))
  
  #----------------------------------------------------------------------
  # Plots to show the seasonal patterns in the observation and in the GCM
  #----------------------------------------------------------------------
  
  pdf(file=paste0(fPath,modelName,'_',fConc,"_",varName,'_Bootstrap_Boxplots.pdf',sep=""), height=8, width=14, onefile=TRUE, family='Helvetica', paper="a4r", pointsize=12) 
  
  par(mfrow = c(2,1), mar = c(2,4,4,1),xpd = TRUE)
  
  # shows some results:
  maxy <- 1.03*max(bs_result$bs_value[is.finite(bs_result$bs_value)], data_int$avg_value[is.finite(data_int$avg_value)])
  
  boxplot(bs_value~simMonth,data=bs_result, main = modelName, col="coral1", boxwex  = 0.2, at = 1:12 + 0.2,yaxs = "i", ylim = c(0,maxy), ylab = "Daily Avg Precip by Month (mm)")
  boxplot(avg_value~Month,data=data_int,col="cyan", boxwex  = 0.2, at = 1:12 - 0.2,add= T)
  legend(x="top",inset = c(0, -0.11), legend = c("Observed", "Synthetic"),fill= c("cyan","coral1"), horiz = T, cex = 0.75)
  
  bs_ym_sum <- bs_result[,c("simMonth", "simYear", "bs_value")] %>% group_by(simMonth, simYear) %>% summarise_all(sum)
  bs_ym_sum <- bs_ym_sum[order(bs_ym_sum$simYear),]
  
  # a. box plot of stationality
  maxy <- 1.03*max(ensemble_month$pr_E1[is.finite(ensemble_month$pr_E1)], data_m$avg_value[is.finite(data_m$avg_value)], ensemble_month_f$pr_E1[is.finite(ensemble_month_f$pr_E1)], bs_ym_sum$bs_value[is.finite(bs_ym_sum$bs_value)])
  
  ay = which(data_m$Year >= minObsYear) #1970)
  nam = c("","","","","","","","","","","","")
  
  boxplot(pr_E1~time_month,data=ensemble_month, main='', 
          xlab="Month", ylab="Monthly Avg Precip (mm)",col="green", boxwex  = 0.15, at = 1:12 - 0.0,yaxs = "i", ylim = c(0, maxy))   # GCM historico
  
  boxplot(avg_value[ay]~Month[ay],data=data_m, 
          xlab="Month",names = nam, ylab="Monthly Avg Precip (mm)",col="cyan",  boxwex  = 0.15, at = 1:12 - 0.2,add=TRUE) # Observado  
  
  boxplot(pr_E1~time_month,data=ensemble_month_f, 
          xlab="Month",names = nam, ylab="Monthly Avg Precip (mm)",col="orange", boxwex  = 0.15, at = 1:12 + 0.2,add=TRUE) #GCM Future
  
  boxplot(bs_value~simMonth,data=bs_ym_sum, xlab="Month",names = nam, ylab="Monthly Avg Precip (mm)",col="coral1", boxwex  = 0.15, at = 1:12 + 0.4,add=TRUE) #GCM Future
  
  legend(x="top", inset = c(0,-0.11), legend = c(paste0("Observed ", minObsYear, "-", maxObsYear), paste0("GCM ", minObsYear, "-", maxObsYear), paste0("GCM ",minGCMYear,"-",maxGCMfYear), paste0("Synth ",minGCMYear,"-",maxGCMfYear)),fill= c("cyan","green", "orange", "coral1"), horiz = T, cex = 0.75)
  
  dev.off()
  
  #---------------------------------------------------------
  # Plot Upper/Lower Thresh Comp Plots ####
  #--------------------------------------------------------
  
  pdf(file=paste0(fPath,modelName,'_',fConc,"_",varName,"_Threshold_Comp_Plots.pdf"), paper = 'a4r', height=8.5, width=11)
  par(mfrow = c(1, 1), xpd = TRUE, mar = c(2, 4, 5, 2))
  
  if(JPmode == "Window")
  {
    extreme_p_gcm_f = sapply(0:364, function(x) {
      jdays <- (x- nearWindow):(x + nearWindow)
      jdays[jdays < 0] <- jdays[jdays < 0] + 365
      percentile(ensemble_day_f$pr_E1[ensemble_day_f$julian_day %in% jdays],1 - extreme_wet_percentile) # wet/extremely wet thresh for gcm historical data   
    }) 
    # extremely wet    
    extreme_l_gcm_f = sapply(0:364, function(x) {
      jdays <- (x - nearWindow):(x + nearWindow)
      jdays[jdays < 0] <- jdays[jdays < 0] + 365
      val <- max(Threshold,percentile(ensemble_day_f$pr_E1[ensemble_day_f$julian_day %in% jdays],0.99)) # dry/wet thresh for gcm hist data
      if(val >= extreme_p_gcm_f[x+1] | val == 0)
      {
        val <- extreme_p_gcm_f[x+1]/2
      }
      val
    })
  } else if(JPmode == "Yearly")
  {
    extreme_p_gcm_f <- percentile(ensemble_day_f$pr_E1,1 - extreme_wet_percentile)
    extreme_l_gcm_f <- max(Threshold,percentile(ensemble_day_f$pr_E1,0.99))
  } else if(JPmode == "Monthly")
  {
    extreme_p_gcm_f = sapply(1:12, function(x) percentile(ensemble_day_f$pr_E1[ensemble_day_f$time_month == x],1 - extreme_wet_percentile))  # extremely wet    
    extreme_l_gcm_f = sapply(1:12, function(x) max(Threshold,percentile(ensemble_day_f$pr_E1[ensemble_day_f$time_month == x],0.99)))
  }
  
  plot(extreme_p, type = 'l', lwd = 2, ylim = c(min(extreme_p, extreme_p_gcm, extreme_p_gcm_f),max(extreme_p, extreme_p_gcm, extreme_p_gcm_f)), ylab = "mm", xlab = '', main = paste0("Upper Threshold Values (p",gsub("0.","",extreme_wet_percentile),")"))
  lines(extreme_p_gcm, lwd = 1.8, lty = "dashed", col = 4)
  lines(extreme_p_gcm_f, lwd = 1.8, lty = "dashed", col = 2)
  
  legend(x = "top", inset = c(0,-0.05), legend = c("Observed", "Model Historical  ", "Model Future"), horiz = TRUE, col = c("black", "blue", "red"), lty = c(1,2,2), pch = c(NA,NA,NA), cex = 0.75)
  
  plot(extreme_l, type = 'l', lwd = 2, ylim = c(min(extreme_l, extreme_l_gcm, extreme_l_gcm_f),max(extreme_l, extreme_l_gcm, extreme_l_gcm_f)), ylab = "mm", xlab = '', main = paste0("Lower Threshold Values (",Threshold," mm or p01)"))
  lines(extreme_l_gcm, lwd = 1.8, lty = "dashed", col = 4)
  lines(extreme_l_gcm_f, lwd = 1.8, lty = "dashed", col = 2)
  
  legend(x = "top", inset = c(0,-0.05), legend = c("Observed", "Model Historical  ", "Model Future"), horiz = TRUE, col = c("black", "blue", "red"), lty = c(1,2,2), pch = c(NA,NA,NA), cex = 0.75)
  
  dev.off()
  #-------------------------------------------------------------------------------------------------
  
  ## Plot Seasonal Joint Probability Charts for Obersved Historical and Historical GCM Daily Data
  
  pdf(file=paste0(fPath,modelName,'_',fConc,"_",varName,"_JP_Monthly_plots.pdf"))
  par(mfrow = c(3, 3), mar = c(2,3,5.5,2))
  
  high_thr = percentile(data_int$avg_value,1 - extreme_wet_percentile)  # wet/extremely wet thresh for observed historical data   
  lower_thr = max(Threshold,percentile(data_int$avg_value,0.99)) # dry/wet thresh for obs hist data
  
  jp <- tri_state_joint_prob_monthly(data_int[,c("Year", "Month", "avg_value")], lower_thr, high_thr, monthly_thresholds = FALSE)
  
  jp_bs <- tri_state_joint_prob_monthly(bs_result[,c("simYear", "simMonth", "bs_value")], lower_thr, high_thr, monthly_thresholds = FALSE)
  
  high_thr_gcm = percentile(ensemble_day$pr_E1,1 - extreme_wet_percentile)  # wet/extremely wet thresh for observed historical data   
  lower_thr_gcm = max(Threshold,percentile(ensemble_day$pr_E1,0.99)) # dry/wet thresh for obs hist data
  
  jp_gcm_hist <- tri_state_joint_prob_monthly(ensemble_day[,c("time_year", "time_month", "pr_E1")], lower_thr_gcm, high_thr_gcm, monthly_thresholds = FALSE)
  
  high_thr_gcm_f = percentile(ensemble_day_f$pr_E1,1 - extreme_wet_percentile)  # wet/extremely wet thresh for observed historical data   
  lower_thr_gcm_f = max(Threshold,percentile(ensemble_day_f$pr_E1,0.99)) # dry/wet thresh for obs hist data
  
  jp_gcm_f <- tri_state_joint_prob_monthly(ensemble_day_f[,c("time_year", "time_month", "pr_E1")], lower_thr_gcm_f, high_thr_gcm_f, monthly_thresholds = FALSE)
  
  for(i in 1:3)
  {
    plot(1:12,jp[i,seq(1,34,3)], type = 'l', xlab = "", ylab = "Joint Probability", ylim = c(min(jp[i,seq(1,34,3)], jp_gcm_hist[i,seq(1,34,3)],jp_gcm_f[i,seq(1,34,3)], jp_bs[i, seq(1,34,3)]), max(jp[i,seq(1,34,3)], jp_gcm_hist[i,seq(1,34,3)],jp_gcm_f[i,seq(1,34,3)], jp_bs[i, seq(1,34,3)])))
    title(paste(row.names(jp)[i], "toLow"), line = 1)
    lines(1:12,jp_gcm_hist[i,seq(1,34,3)], lty = "dashed", col = "blue")
    lines(1:12,jp_gcm_f[i,seq(1,34,3)], lty = "dashed", col = "red")
    lines(1:12, jp_bs[i, seq(1,34,3)], lty = "dashed", col = "green4")
    
    plot(1:12,jp[i,seq(2,35,3)], type = 'l', xlab = "", ylab = "Joint Probability", ylim = c(min(jp[i,seq(2,35,3)], jp_gcm_hist[i,seq(2,35,3)],jp_gcm_f[i,seq(2,35,3)], jp_bs[i,seq(2,35,3)]), max(jp[i,seq(2,35,3)], jp_gcm_hist[i,seq(2,35,3)],jp_gcm_f[i,seq(2,35,3)], jp_bs[i,seq(2,35,3)]))) 
    title(paste(row.names(jp)[i], "toMid"), line = 1)
    lines(1:12,jp_gcm_hist[i,seq(2,35,3)], lty = "dashed", col = "blue")
    lines(1:12,jp_gcm_f[i,seq(2,35,3)], lty = "dashed", col = "red")
    lines(1:12,jp_bs[i,seq(2,35,3)], lty = "dashed", col = "green4")
    
    plot(1:12,jp[i,seq(3,36,3)], type = 'l', xlab = "", ylab = "Joint Probability", ylim = c(min(jp[i,seq(3,36,3)], jp_gcm_hist[i,seq(3,36,3)],jp_gcm_f[i,seq(3,36,3)], jp_bs[i,seq(3,36,3)]), max(jp[i,seq(3,36,3)], jp_gcm_hist[i,seq(3,36,3)],jp_gcm_f[i,seq(3,36,3)], jp_bs[i,seq(3,36,3)])))
    title(paste(row.names(jp)[i], "toHigh"), line = 1)
    lines(1:12,jp_gcm_hist[i,seq(3,36,3)], lty = "dashed", col = "blue")
    lines(1:12,jp_gcm_f[i,seq(3,36,3)], lty = "dashed", col = "red")
    lines(1:12,jp_bs[i,seq(3,36,3)], lty = "dashed", col = "green4")
  }
  
  mtext(paste("Obs Low Threshold =", round(lower_thr,2),"mm",    "   Obs Upper Threshold (p", extreme_wet_percentile*100,") =", round(high_thr,2),"mm (Black)"), outer = TRUE, cex = 0.5, line = -1)
  mtext(paste("BS Synth Low Threshold =", round(lower_thr,2),"mm",    "   BS Synth Upper Threshold (p", extreme_wet_percentile*100,") =", round(high_thr,2),"mm (Green)"), outer = TRUE, cex = 0.5, line = -1.7)
  mtext(paste("GCM Hist Low Threshold =", round(lower_thr_gcm,2),"mm",    "   GCM Hist Upper Threshold (p", extreme_wet_percentile*100,") =", round(high_thr_gcm,2),"mm (Blue)"), outer = TRUE, cex = 0.5, line = -2.4)
  mtext(paste("GCM Future Low Threshold =", round(lower_thr_gcm,2),"mm",    "   GCM Future Upper Threshold (p", extreme_wet_percentile*100,") =", round(high_thr_gcm,2),"mm (Red)"), outer = TRUE, cex = 0.5, line = -3.1)
  
  dev.off()
  
#----------------------------------------------------------------------------------
  ##### PERCENTILE COMPARISON CHARTS ####
  #pdf(file=paste0(main_folder,results_folder,modelName,'_',fConc,"_",varName,"_Percentiles_Comp.pdf"), paper = 'a4r', height=8.5, width=11)
  #par(mfrow = c(2, 1), mar = c(2, 4, 5, 2))
  
  # p0 to p99
  #pctls_obs <- sapply(1:99, function(x) percentile(data_int$avg_value, (1 - x/100)))
  #pctls_gcmh <- sapply(1:99, function(x) percentile(ensemble_day$pr_E1, (1 - x/100)))
  #pctls_gcmf <- sapply(1:99, function(x) percentile(ensemble_day_f$pr_E1, (1 - x/100)))
  #pctls_bs <- sapply(1:99, function(x) percentile(bs_result$bs_value, (1 - x/100)))
  
  #plot(seq(0.01,0.99,0.01),pctls_obs, xlab = 'Percentile', ylab = 'mm', ylim = c(min(pctls_gcmh,pctls_gcmf, pctls_bs, pctls_obs), max(pctls_gcmh,pctls_gcmf, pctls_bs, pctls_obs)), type = 'o')
  #lines(seq(0.01,0.99,0.01),pctls_gcmh, col = 'blue', type = 'o')
  #lines(seq(0.01,0.99,0.01),pctls_gcmf, col = 'red', type = 'o')
  #lines(seq(0.01,0.99,0.01),pctls_bs, col = 'green', type = 'o')
  
  #par(xpd = TRUE)
  #legend('top', legend = c("Obs Hist", "GCM Hist", "GCM Future", "Synth Bootstrap"), col = c(1, "blue", "red", "green"), pch = c(1,1,1,1), horiz = T, inset = c(0, -0.3))
  #par(xpd = FALSE)
  
  # p90 to p95 #
  
  #pctls_obs <- sapply(1:99, function(x) percentile(data_int$avg_value, (0.1 - x/1000)))
  #pctls_gcmh <- sapply(1:99, function(x) percentile(ensemble_day$pr_E1, (0.1 - x/1000)))
  #pctls_gcmf <- sapply(1:99, function(x) percentile(ensemble_day_f$pr_E1, (0.1 - x/1000)))
  #pctls_bs <- sapply(1:99, function(x) percentile(bs_result$bs_value, (0.1 - x/1000)))
  
  #plot(seq(0.901,0.999,0.001),pctls_obs, xlab = 'Percentile', ylab = 'mm', ylim = c(min(pctls_gcmh,pctls_gcmf, pctls_bs, pctls_obs), max(pctls_gcmh,pctls_gcmf, pctls_bs, pctls_obs)), type = 'o')
  #lines(seq(0.901,0.999,0.001),pctls_gcmh, col = 'blue', type = 'o')
  #lines(seq(0.901,0.999,0.001),pctls_gcmf, col = 'red', type = 'o')
  #lines(seq(0.901,0.999,0.001),pctls_bs, col = 'green', type = 'o')
  
  #par(xpd = TRUE)
  #legend('top', legend = c("Obs Hist", "GCM Hist", "GCM Future", "Synth Bootstrap"), col = c(1, "blue", "red", "green"), pch = c(1,1,1,1), horiz = T, inset = c(0, -0.3))
  #par(xpd = FALSE)
  
  #grid()
  
  #dev.off()
}


# ---------------------------------------
# IV. HEAT SIGNAL FUNCTION 
# ---------------------------------------

heatSignal = function(refDate = c(month = 1, day = 1, year = 1850),
                      main_folder = "G:/GCM_ClimateTool",
                      results_folder = "/Results/test",
                      modelName = "CCSM4",
                      futures = c("rcp85"),
                      varName = "tas",
                      varLabel = "Temperature [c]",
                      minObsYear = 1980,
                      maxObsYear = 2016,
                      maxGCMYear = 2100,
                      minGCMYear = 2015) {
  
 # refDate = c(month = df$BaseMonth[rcp_tas_ind], day = df$BaseDay[rcp_tas_ind], year = df$BaseYear[rcp_tas_ind])
 # main_folder = mFolder
 # results_folder = rFolder
 # modelName = model
 # futures = experiment2
 # varName = "tas"
 # varLabel = "Temperature [c]"
 # minObsYear = 1980
 # maxObsYear = 2015
 # minGCMYear = 2040
 # maxGCMYear = 2060
  
  heatSignalLength <- maxGCMYear - minGCMYear + 1 
  # LOADS DATA
  #fPath <- paste(main_folder,results_folder,'/',modelName,'/',sep="")
  fPath <- paste(main_folder,results_folder,'/',sep="")
  fConc = futures
  
  fName = paste(fPath,paste(modelName,fConc,varName,"ensemble_day.Rda",sep="_"),sep="")
  load(fName,verbose=TRUE)
  ensemble_day$month_idx = (ensemble_day$time_year - refDate[3]) * 12 + ensemble_day$time_month
  ensemble_day_f <- ensemble_day[ensemble_day$time_year >= minGCMYear & ensemble_day$time_year <= maxGCMYear,]
  ensemble_day_h <- ensemble_day[ensemble_day$time_year >= minObsYear & ensemble_day$time_year <= maxObsYear,]
  
  fName = paste(fPath,paste(modelName,fConc,varName,"ensemble_year.Rda",sep="_"),sep="")
  load(fName)   
  ensemble_year_f <- ensemble_year[ensemble_year$time_year >= minGCMYear & ensemble_year$time_year <= maxGCMYear,]
  ensemble_year_h <- ensemble_year[ensemble_year$time_year >= minObsYear & ensemble_year$time_year <= maxObsYear,]
  
  fName = paste(fPath,paste(modelName,fConc,varName,"ensemble_month.Rda",sep="_"),sep="")
  load(fName)   
  ensemble_month_f <- ensemble_month[ensemble_month$time_year >= minGCMYear & ensemble_month$time_year <= maxGCMYear,]
  ensemble_month_h <- ensemble_month[ensemble_month$time_year >= minObsYear & ensemble_month$time_year <= maxObsYear,]
  
  # observed historic
  if(varName == "pr")
  {
    #fName = paste(fPath,paste(modelName, "observed", varName,futures,"day_ExtremesCorrected.Rda",sep="_"),sep="")
    fName = paste(fPath,paste(modelName, "observed", varName,"day.Rda",sep="_"),sep="")
  } else {
    fName = paste(fPath,paste(modelName, "observed", varName,"day.Rda",sep="_"),sep="")
  }

  load(fName,verbose=TRUE)
  fName = paste(fPath,paste(modelName, "observed", varName,"year.Rda",sep="_"),sep="")
  load(fName,verbose=TRUE)
  fName = paste(fPath,paste(modelName, "observed", varName,"month.Rda",sep="_"),sep="")
  load(fName,verbose=TRUE)
  
  data_d <- data_d[data_d$Year >= minObsYear & data_d$Year <= max(ensemble_year$time_year),]
  data_m <- data_m[data_m$Year >= minObsYear & data_m$Year <= max(ensemble_year$time_year),]
  data_y <- data_y[data_y$Year >= minObsYear & data_y$Year <= max(ensemble_year$time_year),]
  
#  ensemble_day_h <- ensemble_day_h[ensemble_day_h$time_year >= minObsYear,]
#  ensemble_month_h <- ensemble_month_h[ensemble_month_h$time_year >= minObsYear,]
#  ensemble_year_h <- ensemble_year_h[ensemble_year_h$time_year >= minObsYear,]
  
  minObsYear <- min(data_d$Year)
  maxObsYear <- max(data_d$Year)
  
  # Creates a symmetrically smoothed (weekly moving average of signal):
  # average of current sample, 3 future samples, and 3 past samples 
  f7 <- rep(1/7,7)
  ensemble_day_f$wk_ma_signal <- filter(ensemble_day_f$tas_E1,f7,sides=2)
  ensemble_day$wk_ma_signal <- filter(ensemble_day$tas_E1,f7,sides=2)
  data_d$wk_ma_signal <- filter(data_d$avg_value,f7,sides=2)
  
  #### Assigns the first and last three days or the series to the nearest value (4th and 4th to last) so as to not have 3 NAs leading and ending the series ###
  ensemble_day_f$wk_ma_signal[1:3] <- ensemble_day_f$wk_ma_signal[4] 
  ensemble_day_f$wk_ma_signal[(nrow(ensemble_day_f) - 3):nrow(ensemble_day_f)] <- ensemble_day_f$wk_ma_signal[nrow(ensemble_day_f) - 4]
  ensemble_day$wk_ma_signal[1:3] <- ensemble_day$wk_ma_signal[4]
  ensemble_day$wk_ma_signal[(nrow(ensemble_day) - 3):nrow(ensemble_day)] <- ensemble_day$wk_ma_signal[nrow(ensemble_day) - 4]
  data_d$wk_ma_signal[1:3] <- data_d$wk_ma_signal[4]
  data_d$wk_ma_signal[(nrow(data_d) - 3):nrow(data_d)] <- data_d$wk_ma_signal[nrow(data_d) - 4]
  
  # adds a julian day column for indexing purposes  
  ensemble_day_f$julian_day = julian(ensemble_day_f$time_month,ensemble_day_f$time_day,2005,
                                     origin = c(month = 1, day = 1, year = 2005))
  
  ensemble_day$julian_day = julian(ensemble_day$time_month,ensemble_day$time_day,2005,
                                   origin = c(month = 1, day = 1, year = 2005))
  
  data_d$julian_day = julian(data_d$Month,data_d$Day,2005,
                             origin = c(month = 1, day = 1, year = 2005))
  
  
  # creates a new data frame for the reference period and then a dataframe of average ta signal for the reference period 
  
  ensemble_day_ref <- subset(ensemble_day, ensemble_day$time_year >= minObsYear & ensemble_day$time_year <= maxObsYear)
  ma_weekly_ref <- aggregate(ensemble_day_ref, by=list(ensemble_day_ref$julian_day), FUN=mean,na.rm=TRUE)
  
  obs_day_ref <- subset(data_d, data_d$Year >= minObsYear & data_d$Year <= maxObsYear)
  obs_ma_weekly_ref <- aggregate(obs_day_ref, by=list(obs_day_ref$julian_day), FUN=mean,na.rm=TRUE)
  
  # creates a new data frame for the future period 
  
  ensemble_day_f$wk_ma_ref = ensemble_day_f$julian_day
  
  for (d in 0:364) {
    n = ensemble_day_f$julian_day %in% d
    ensemble_day_f$wk_ma_ref[t(n)] = ma_weekly_ref$wk_ma_signal[d+1]
    ensemble_day_f$obs_wk_ma_ref[t(n)] <- obs_ma_weekly_ref$wk_ma_signal[d+1]
  }
  
  ensemble_day_f$heatSignal = ensemble_day_f$wk_ma_signal - ensemble_day_f$wk_ma_ref
  
  df_heat_signal <- subset(ensemble_day_f, ensemble_day_f$time_year >= minGCMYear & ensemble_day_f$time_year < minGCMYear + heatSignalLength)
  
  write.table(df_heat_signal,paste(fPath,paste(modelName,fConc,varName,"heatSignal.csv",sep="_"),sep=""), sep=",",row.names=FALSE)
  
  # Saves the data frame with the heat signal
  save(df_heat_signal,file=paste(main_folder,results_folder,"/",paste(modelName,"heatSignal",".Rda",sep="_"),sep=""))
  
  dummyDate = julian(df_heat_signal$time_month,df_heat_signal$time_day,df_heat_signal$time_year,origin = c(month =1, day = 1,year = 1900))
  df_heat_signal$timeIdxs = format(as.Date(dummyDate, origin = "1900-01-01"), "%d-%m-%Y") #chron(dummyDate, origin = c(month =1, day = 1,year = 1900),format = "d-m-Y")
  drops <- c("time_year","time_month","time_day","tas_E1","month_idx","wk_ma_signal","julian_day","wk_ma_ref")
  df_heat_signal_WP = df_heat_signal[,!(names(df_heat_signal) %in% drops)]
  write.table(df_heat_signal_WP,paste(fPath,paste(modelName,fConc,"TS_HeatSignal_WP.csv",sep="_"),sep=""), sep=",",row.names=FALSE, na = "-9999")
  
  #-------------------------------------------------------------------------
  # PLOTS OF TEMPERATURE COMPARISON OBS, GCM HIST, GCM FUTURE
  #-------------------------------------------------------------------------
  
  pdf(file=paste(fPath,modelName,'_',fConc,"_",varName,'_Temp_Med_Comp_BoxPlots.pdf',sep=""), height=8.5, width=11, onefile=TRUE, family='Helvetica', paper="a4r", pointsize=12) 
  
  par(mfrow = c(2,1), mar = c(2,4,4.5,1), xpd = TRUE)
  
  maxy <- max(ensemble_day$tas_E1[is.finite(ensemble_day$tas_E1)], data_d$avg_value[is.finite(data_d$avg_value)], ensemble_day_f$tas_E1[is.finite(ensemble_day_f$tas_E1)])
  miny <- min(ensemble_day$tas_E1[is.finite(ensemble_day$tas_E1)], data_d$avg_value[is.finite(data_d$avg_value)], ensemble_day_f$tas_E1[is.finite(ensemble_day_f$tas_E1)])
  
  nam = c("","","","","","","","","","","","")
  
  boxplot(avg_value~Month,data=data_d, xlab="Month", ylab="Monthly Avg Temp Media [C]",col="cyan", boxwex  = 0.15, at = 1:12, ylim = c(miny, maxy), main = paste0("Avg. Temperature (",model,")")) #GCM Future
  
  boxplot(tas_E1~time_month,data=ensemble_day, xlab="Month",col="green", boxwex  = 0.15, at = 1:12 + 0.2, add = T, names = nam) #GCM Future
  
  boxplot(tas_E1~time_month,data=ensemble_day_f, xlab="Month",col="orange", boxwex  = 0.15, at = 1:12 + 0.4, add = T, names = nam) #GCM Future
  
  legend(x="top", inset = c(0,-0.12), legend = c(paste0("Observed ", minObsYear, "-", maxObsYear), paste0("GCM ", minObsYear, "-", maxObsYear), paste0("GCM ",minGCMYear,"-",maxGCMYear)),fill= c("cyan","green", "orange"), horiz = T, cex = 0.75)
  
  #------------------------------------------------------
  
  #Raw Histograms of Obs, GCM Hist, GCm Future
  
  tempstats <- data.frame(ID = rep(NA, 3), Mean = rep(NA,3), SD = rep(NA,3))
  
  tempstats$ID <- c("Observed", "GCM Hist", "GCM Future")
  tempstats$Mean <- round(c(mean(data_d$avg_value), mean(ensemble_day$tas_E1), mean(ensemble_day_f$tas_E1)),2)
  tempstats$SD <- round(c(sd(data_d$avg_value), sd(ensemble_day$tas_E1), sd(ensemble_day_f$tas_E1)),2)
  
  par(mar = c(4,4,3,2), xpd = FALSE)
  
  multhist(list(data_d$avg_value, ensemble_day$tas_E1, ensemble_day_f$tas_E1), beside = TRUE, breaks = 50, col = c("blue","green3","orange"), border = NA, xlab = "Avg. Daily Temp [C]", main = "", ylab = "Density", freq = FALSE)
  
  legend(x="topright", legend = c(paste0("Observados ", minObsYear, "-", maxObsYear,"\n mean = ",tempstats$Mean[1],"   sd = ", tempstats$SD[1],"\n"), paste0("GCM ", minObsYear, "-", maxObsYear,"\n mean = ",tempstats$Mean[2],"   sd = ", tempstats$SD[2]), paste0("GCM ",minGCMYear,"-",maxGCMYear,"\n mean = ",tempstats$Mean[3],"   sd = ", tempstats$SD[3])), fill= c("blue","green3", "orange"), cex = 0.75, bty = 'n')
  
  dev.off()
  
  return(df_heat_signal)
}

# --------------------------------------
# Apply a bootstrapped series
# --------------------------------------
sinteticSeries = function(refDate = c(month = 1, day = 1, year = 1850),
                          modelName = "MPI-ESM-MR",
                          futures = c("rcp85"),
                          main_folder = "G:/GCM_ClimateTool",
                          results_folder = "/Results/test",
                          expNumber = 1,
                          fullObsCatalog = TRUE, # T/F to apply the bs date series to full obs station catalog or not
                          HistRepro = TRUE,
                          minObsYear = 1980,
                          maxObsYear = 2015) {

  ## experimentation  
#  refDate = c(month = df$BaseMonth[rcp_pr_ind], day = df$BaseDay[rcp_pr_ind], year = df$BaseYear[rcp_pr_ind])
  #  modelName = model
  #  futures = experiment
  #main_folder = mFolder
  #results_folder = rFolder
  #expNumber = n
  #fullObsCatalog = FALSE
  #HistRepro = FALSE
  #minObsYear = 1980
  #maxObsYear = 2015
  
  extrcorr_regiones <- c("La Paz")
  fullCatalogFilePr = paste0("/OBSERVATIONS/Bolivia/","SERIES_pr_DIARIO_v4.Rda")
  fullCatalogFileTemp = "/OBSERVATIONS/Bolivia/SERIES_tas_DIARIO_v4.Rda"
  
  # LOADS DATA
  #fPath <- paste(main_folder,results_folder,'/',modelName,'/',sep="")
  fPath <- paste(main_folder,results_folder,'/',sep="")
  fConc = futures
  
  # GCM precip futures
  fName = paste(fPath,paste(modelName,fConc,"pr","ensemble_day.Rda",sep="_"),sep="")
  load(fName) # ensemble_day   
  pr_ensemble_day_f <- ensemble_day
  
  # GCM precip historic
  fName = paste(fPath,paste(modelName,"hist-rcp85","pr","ensemble_day.Rda",sep="_"),sep="")
  
  load(fName) # ensemble_day
  pr_ensemble_day <- ensemble_day
  
  # GCM temperature futures
  fName = paste(fPath,paste(modelName,fConc,"tas","ensemble_day.Rda",sep="_"),sep="")
  load(fName) # ensemble_day    
  tas_ensemble_day_f <- ensemble_day
  
  # GCM temperature historic
  fName = paste(fPath,paste(modelName,"hist-rcp85","tas","ensemble_day.Rda",sep="_"),sep="")
  load(fName) # ensemble_day
  tas_ensemble_day <- ensemble_day
  
  # observed precip historic
  fName = paste(fPath,paste(modelName, "observed", "pr_day.Rda",sep="_"),sep="")
  load(fName) # data_d
  pr_data_d <- data_d
  
  if(HistRepro != TRUE)
  {
    # Correct Obs Hist Data for Extremes #
    fName = paste(fPath,paste(modelName, "bs_pr",experiment, "day_ExtremesCorrected_key.Rda",sep="_"),sep="")
    load(fName)
    
    sta_ind <- grep("pr_", colnames(pr_data_d))
    for(s in sta_ind)
    {
      GPDobs <- gpd(sort(pr_data_d[,s]),threshold=percentile(pr_data_d[,s], 1 - bs_excorr$bs_pctl))
      GPDbeta  <- GPDobs$par.ests[2]*bs_excorr$GPDbetaRatio
      GPDalpha <-  GPDobs$par.ests[1] + bs_excorr$GPDalphaDiff #^(GPDhist$par.ests[2]/GPDobs$par.ests[2])
      GPDThresh <-  GPDobs$threshold * bs_excorr$GPDThreshRatio
      obs_extr <- rgpd(10000, GPDobs$par.ests[1], GPDobs$threshold, GPDobs$par.ests[2])
      xObs_new <- rgpd(10000, GPDalpha, GPDThresh, GPDbeta) #find the new values using the new extreme values
      obsorig <- sort(obs_extr)
      obsnew <- sort(xObs_new)
      
      obs_extr_new <- pr_data_d[,s][which(pr_data_d[,s] >= GPDobs$threshold)]
      obs_extr_new <- sapply(1:length(obs_extr_new), function(a) obsnew[which.min(abs(obsorig - obs_extr_new[a]))])
      
      idx.ext.obs <- which(pr_data_d[,s] >= percentile(pr_data_d[,s], 1 - bs_excorr$bs_pctl))
      pr_data_d[,s][idx.ext.obs] <- obs_extr_new
    }
    pr_data_d$avg_value <- rowMeans(pr_data_d[,sta_ind], na.rm = T) 
  }
 
  # observed tas historic
  fName = paste(fPath,paste(modelName, "observed", "tas","day.Rda",sep="_"),sep="")
  load(fName,verbose=TRUE) # data_d
  tas_data_d <- data_d
  
  # bootstraped signal: bs_result
  fName = paste(fPath,paste(modelName,fConc,"pr",expNumber,"bootstrap.Rda",sep="_"),sep="")
  load(fName,verbose=TRUE) # bs_result
  
  # heat signal: df_heat_signal
  fName = paste(fPath,paste(modelName,"heatSignal_.Rda",sep="_"),sep="")
  load(fName,verbose=TRUE) # df_heat_signal
  
  # creates a new data-frame wih the boostraped signal for local observations 
  pr_sintetic <- bs_result
  
  # adds an index to the bs_result
  lnt = dim(pr_sintetic)[1]
  pr_sintetic$index = seq(1,lnt)
  
  pr_sintetic$bs_timeStamp = julian(pr_sintetic$Month,pr_sintetic$Day,pr_sintetic$Year,origin = c(month =1, day = 1,year = 1900))
  pr_data_d$bs_timeStamp = julian(pr_data_d$Month,pr_data_d$Day,pr_data_d$Year,origin = c(month =1, day = 1,year = 1900))
  
  if(fullObsCatalog == TRUE)
  {
    # observed precip historic whole region ALL stations
    fName = paste0(mFolder, fullCatalogFilePr)
    load(fName) # pr_data_d_ALL
    
    fName = paste0(mFolder, fullCatalogFileTemp)
    load(fName) # tas_data_d_ALL
    
    pr_data_d_ALL$bs_timeStamp = julian(pr_data_d_ALL$Month,pr_data_d_ALL$Day,pr_data_d_ALL$Year,origin = c(month =1, day = 1,year = 1900))
    bs_pr_local_data_ALL <- (merge(pr_sintetic,pr_data_d_ALL, by = 'bs_timeStamp'))
    bs_pr_local_data_ALL <- bs_pr_local_data_ALL[order(bs_pr_local_data_ALL$index),]
    dummyDate = julian(bs_pr_local_data_ALL$simMonth,bs_pr_local_data_ALL$simDay,bs_pr_local_data_ALL$simYear,origin = c(month =1, day = 1,year = 1900))
    bs_pr_local_data_ALL$bs_timeStamp = format(as.Date(dummyDate, origin = "1900-01-01"), "%m/%d/%Y")
    
    # if(HistRepro != TRUE)
    # {
    # for(r in extrcorr_regiones)
    # {
    #   fPath <- paste0(main_folder,"/Results/Magdalena-Cauca-Test/",r,"/")
    #   fName <- paste(fPath,paste(modelName, "bs_pr",futures, "day_ExtremesCorrected_key.Rda",sep="_"),sep="")
    #   load(fName)
    #   fName <- paste(fPath,paste(modelName, "observed_pr_year.Rda", sep = "_"),sep="")
    #   load(fName)
    #   
    #   sta_reg <- colnames(data_y)[grep("pr_", colnames(data_y))]
    #   sta_reg <- sta_reg[which(sta_reg %in% colnames(bs_pr_local_data_ALL))]
    #   sta_ind <- sapply(sta_reg, function(x) which(colnames(bs_pr_local_data_ALL) == x))
    #   
    #   for(s in sta_ind)
    #   {
    #     GPDobs <- try(gpd(sort(bs_pr_local_data_ALL[,s]),threshold=percentile(bs_pr_local_data_ALL[,s], 1 - bs_excorr$bs_pctl)), silent = TRUE)
    #     if(class(GPDobs) == "try-error")
    #     {
    #       Threshold <- 2
    #       repeat {
    #         GPDobs <- try(gpd(sort(bs_pr_local_data_ALL[,s]),threshold=Threshold),silent = TRUE)
    #         if (substr(GPDobs[1],1,5) != "Error" && GPDobs$n.exceed < 200) { break }
    #         Threshold <- Threshold + 2
    #       }
    #       rm(Threshold)
    #     }
    #     GPDbeta  <- GPDobs$par.ests[2]*bs_excorr$GPDbetaRatio
    #     GPDalpha <-  GPDobs$par.ests[1] + bs_excorr$GPDalphaDiff #^(GPDhist$par.ests[2]/GPDobs$par.ests[2])
    #     GPDThresh <-  GPDobs$threshold * bs_excorr$GPDThreshRatio
    #     obs_extr <- rgpd(10000, GPDobs$par.ests[1], GPDobs$threshold, GPDobs$par.ests[2])
    #     xObs_new <- rgpd(10000, GPDalpha, GPDThresh, GPDbeta) #find the new values using the new extreme values
    #     obsorig <- sort(obs_extr)
    #     obsnew <- sort(xObs_new)
    #     
    #     obs_extr_new <- bs_pr_local_data_ALL[,s][which(bs_pr_local_data_ALL[,s] >= GPDobs$threshold)]
    #     obs_extr_new <- sapply(1:length(obs_extr_new), function(a) obsnew[which.min(abs(obsorig - obs_extr_new[a]))])
    #     
    #     idx.ext.obs <- which(bs_pr_local_data_ALL[,s] >= GPDobs$threshold)
    #     bs_pr_local_data_ALL[,s][idx.ext.obs] <- obs_extr_new
    #   } # next sta
    # } # next region
    # } # end HistRepro logic
    #PRECIP
   
    #fPath = paste(main_folder,results_folder,"/",sep="")
    write.table(bs_pr_local_data_ALL,paste(fPath,paste(modelName,futures,"pr_Sintetic_Series_ALL_Stations",expNumber,".csv",sep="_"),sep=""), sep=",",row.names=FALSE)
    
    save(bs_pr_local_data_ALL,file = paste0(fPath,paste(modelName,futures,"pr_Sintetic_Series_ALL_Stations",expNumber,".Rda",sep="_")))
    
    drops <- c("simYear","simMonth","simDay","Year.x","Month.x","Day.x","bs_value",  "bs_wet_state",	"bs_brick_size","index","TimeStamp","Year.y","Month.y","Day.y")
    bs_local_WPformat_ALL <- bs_pr_local_data_ALL[,!(names(bs_pr_local_data_ALL) %in% drops)]
    
    #adds the date strings and rounds the numbers for a neater file
    write.table(cbind(bs_local_WPformat_ALL$bs_timeStamp,round(bs_local_WPformat_ALL[,!(names(bs_local_WPformat_ALL) %in% "bs_timeStamp")],2)),
                paste(fPath,paste(modelName,futures,"pr_Sintetic_Series_ALL_Stations",expNumber,"WP.csv",sep="_"),sep=""), 
                sep=",",row.names=FALSE,col.names=FALSE,na = "-9999",quote=FALSE)
    
    #TEMP
    heatsig <- read.csv(paste0(fPath,paste(modelName,futures,"tas","heatSignal.csv",sep="_")))
    
    tas_data_d_ALL$bs_timeStamp = julian(tas_data_d_ALL$Month,tas_data_d_ALL$Day,tas_data_d_ALL$Year,origin = c(month =1, day = 1,year = 1900))
    tas_local_data_ALL <- (merge(pr_sintetic,tas_data_d_ALL,  by = 'bs_timeStamp'))
    tas_local_data_ALL <- tas_local_data_ALL[order(tas_local_data_ALL$index),]
    
    nd = dim(tas_local_data_ALL)
    tas_local_data_ALL[,grep("tas", colnames(tas_local_data_ALL))[1]:nd[2]] = tas_local_data_ALL[,grep("tas", colnames(tas_local_data_ALL))[1]:nd[2]] + df_heat_signal$heatSignal[1:nd[1]]
    
    dummyDate = julian(tas_local_data_ALL$simMonth,tas_local_data_ALL$simDay,tas_local_data_ALL$simYear,origin = c(month =1, day = 1,year = 1900))
    tas_local_data_ALL$bs_timeStamp = format(as.Date(dummyDate, origin = "1900-01-01"), "%m/%d/%Y")
    tas_local_data_ALL[,c("timeStamp", "Year.y", "Month.y", "Day.y", "index")] <- NULL
    
    ALL_heatsig <- heatsig
    if(length(which(ALL_heatsig$time_month == 2 & ALL_heatsig$time_day == 29)) > 0)
    {
      ALL_heatsig <- ALL_heatsig[-which(ALL_heatsig$time_month == 2 & ALL_heatsig$time_day == 29),]
    }

    tas_local_data_ALL$heatSignal <- ALL_heatsig$heatSignal
    
    write.table(tas_local_data_ALL,paste(fPath,paste(modelName,futures,"tas_Sintetic_Series_ALL_stations",expNumber,".csv",sep="_"),sep=""), sep=",",row.names=FALSE)
    
    save(tas_local_data_ALL,file = paste0(fPath,paste(modelName,futures,"tas_Sintetic_Series_ALL_stations",expNumber,".Rda",sep="_"),sep=""))
    
    tas_local_WPformat_ALL = tas_local_data_ALL[,!(names(tas_local_data_ALL) %in% drops)]
    write.table(cbind(tas_local_WPformat_ALL$bs_timeStamp,round(tas_local_WPformat_ALL[,!(names(tas_local_WPformat_ALL) %in% "bs_timeStamp")],1)),
                paste(fPath,paste(modelName,fConc,"tas_Sintetic_Series_ALL_stations",expNumber,"WP.csv",sep="_"),sep=""), 
                sep=",",row.names=FALSE,col.names=FALSE, na = "-9999",quote=FALSE)
  } # end fullObsCatalog logic
  ####
  
  tas_data_d$bs_timeStamp = julian(tas_data_d$Month,tas_data_d$Day,tas_data_d$Year,origin = c(month =1, day = 1,year = 1900))
  
  bs_pr_local_data <- (merge(pr_sintetic,pr_data_d, by = 'bs_timeStamp'))
  tas_local_data <- (merge(pr_sintetic,tas_data_d,  by = 'bs_timeStamp'))
  
  bs_pr_local_data <- bs_pr_local_data[order(bs_pr_local_data$index),]
  tas_local_data <- tas_local_data[order(tas_local_data$index),]
  
  # adds heat
  nd = dim(tas_local_data)
  tas_local_data[,grep("tas", colnames(tas_local_data))[1]:nd[2]] = tas_local_data[,grep("tas", colnames(tas_local_data))[1]:nd[2]] + df_heat_signal$heatSignal[1:nd[1]]

  
  dummyDate = julian(bs_pr_local_data$simMonth,bs_pr_local_data$simDay,bs_pr_local_data$simYear,origin = c(month =1, day = 1,year = 1900))
  bs_pr_local_data$bs_timeStamp = format(as.Date(dummyDate, origin = "1900-01-01"), "%m/%d/%Y") #chron(dummyDate, origin = c(month =1, day = 1,year = 1900),format = "m/d/y")
  
  dummyDate = julian(tas_local_data$simMonth,tas_local_data$simDay,tas_local_data$simYear,origin = c(month =1, day = 1,year = 1900))
  tas_local_data$bs_timeStamp = format(as.Date(dummyDate, origin = "1900-01-01"), "%m/%d/%Y") #chron(dummyDate, origin = c(month =1, day = 1,year = 1900),format = "m/d/y")
  
  write.table(bs_pr_local_data,paste(fPath,paste(modelName,futures,"pr_Sintetic_Series",expNumber,".csv",sep="_"),sep=""), sep=",",row.names=FALSE)
  
  save(bs_pr_local_data,file = paste0(fPath,paste(modelName,futures,"pr_Sintetic_Series",expNumber,".Rda",sep="_"),sep=""))
  
  tas_local_data[,c("timeStamp", "Year.y", "Month.y", "Day.y", "index")] <- NULL

  heatsig <- read.csv(paste0(mFolder, rFolder,paste(modelName,futures,"tas","heatSignal.csv",sep="_")))
  
  if(nrow(heatsig) != nrow(tas_local_data) & length(which(heatsig$time_month == 2 & heatsig$time_day == 29)) > 0)
  {
    heatsig <- heatsig[-which(heatsig$time_month == 2 & heatsig$time_day == 29),]
  }

  tas_local_data$heatSignal <- heatsig$heatSignal
  
  write.table(tas_local_data,paste(fPath,paste(modelName,futures,"tas_Sintetic_Series",expNumber,".csv",sep="_"),sep=""), sep=",",row.names=FALSE)
  
  save(tas_local_data, file = paste0(fPath,paste(modelName,futures,"tas_Sintetic_Series",expNumber,".Rda",sep="_"),sep=""))
  
  drops <- c("simYear","simMonth","simDay","Year.x","Month.x","Day.x","bs_value",  "bs_wet_state",	"bs_brick_size","index","TimeStamp","Year.y","Month.y","Day.y")
  bs_local_WPformat = bs_pr_local_data[,!(names(bs_pr_local_data) %in% drops)]
  
  #adds the date strings and rounds the numbers for a neater file
  write.table(cbind(bs_local_WPformat$bs_timeStamp,round(bs_local_WPformat[,!(names(bs_local_WPformat) %in% "bs_timeStamp")],2)),
              paste(fPath,paste(modelName,futures,"pr_Sintetic_Series",expNumber,"WP.csv",sep="_"),sep=""), 
              sep=",",row.names=FALSE,col.names=TRUE,na = "-9999",quote=FALSE)
  
  tas_local_WPformat = tas_local_data[,!(names(tas_local_data) %in% drops)]
  write.table(cbind(tas_local_WPformat$bs_timeStamp,round(tas_local_WPformat[,!(names(tas_local_WPformat) %in% "bs_timeStamp")],1)),
              paste(fPath,paste(modelName,fConc,"tas_Sintetic_Series",expNumber,"WP.csv",sep="_"),sep=""), 
              sep=",",row.names=FALSE,col.names=TRUE, na = "-9999",quote=FALSE)
  
} 

# ---------------------------------------
# IV. EXTREME VALUES SHIFT
# --------------------------------------
extremePR_GPD = function(modelName = "MPI-ESM-MR",
                         futures = c("rcp85"),
                         main_folder = "G:/GCM_ClimateTool/",
                         results_folder = "/Results/test",
                         minObsYear = 1980,
                         minGCMYear = 2015,
                         maxGCMYear = 2050) {
  
 # modelName = model
 # futures = experiment2
 # main_folder = mFolder
 # results_folder = rFolder
 # minObsYear = 1980
 # minGCMYear = 2040
 # maxGCMYear = 2060
  
  # Observed Historic Record
  # fPath <- paste(main_folder,results_folder,'/',modelName,'/',sep="")
  fPath <- paste(main_folder,results_folder,'/',sep="")
  fName <- paste(fPath,paste(modelName, "observed", "pr","day.Rda",sep="_"),sep="") # Historic Observed 
  load(fName,verbose=TRUE)
  
  # GCM precip futures
  fName = paste(fPath,paste(modelName,futures,"pr","ensemble_day.Rda",sep="_"),sep="")
  load(fName,verbose=TRUE)    
  pr_ensemble_day_f <- ensemble_day
  #pr_ensemble_day_f$month_idx = (ensemble_day_f$time_year - refDate[3]) * 12 + ensemble_day_f$time_month
  
  # GCM precip historic
  fName = paste(fPath,paste(modelName,"historical","pr","ensemble_day.Rda",sep="_"),sep="")
  load(fName,verbose=TRUE) 
  pr_ensemble_day <- ensemble_day
  
      pdf(file=paste(fPath,"/",modelName,'_',futures,'_Extremes.pdf',sep=""), 
          height=12, width=8, onefile=TRUE, family='Helvetica', paper="a4r", pointsize=12) 
      par(mfrow = c(2,2))
      
      ######################################################################
  ### For observations
  data_d <- data_d[which(is.finite(data_d$avg_value) & data_d$Year >= minObsYear & data_d$Year < minGCMYear),]
  data_d_obs <- data_d # preserve a copy of the observations
  
  Threshold <- 2
  repeat {
    GPDobs <- try(gpd(sort(data_d$avg_value),threshold=Threshold),silent = TRUE)
    if (substr(GPDobs[1],1,5) != "Error" && GPDobs$n.exceed < 200) { break }
    Threshold <- Threshold + 2
  } 
  
  obs_pctl <- val2pctl(Threshold, data_d$avg_value)
  obs_extr <- rgpd(10000, GPDobs$par.ests[1], GPDobs$threshold, GPDobs$par.ests[2])
  
  ii  <- order(obs_extr)
  gval <- rbind(obs_extr)[,ii]
  mnx <- min(obs_extr)
  mxx <- max(obs_extr)
  
 # hist(obs_extr,freq=FALSE,xlim=c(mnx,mxx), xlab= "Pcp (mm/day)", main = paste("OBS Pcp Above ", round(GPDobs$threshold,2),"mm (p",round(100*obs_pctl,2),")", " - ",round(GPDobs$par.ests[1],2),round(GPDobs$par.ests[2],2)))
  hist(obs_extr,freq=FALSE,xlim=c(mnx,mxx), xlab= "Pcp (mm/day)", main = paste("OBS Pcp Above ", round(GPDobs$threshold,2),"mm (p",round(100*obs_pctl,2),")", "(par.ests ",round(GPDobs$par.ests[1],2),round(GPDobs$par.ests[2],2),")"),
       cex.lab=0.7, cex.axis=0.7, cex.main=0.7, cex.sub=0.7)
  
  
  x <- obs_extr
  
  curve(dgpd(x, xi = GPDobs$par.ests[1], GPDobs$threshold, beta = GPDobs$par.ests[2]), col = "red", lty = 2, lwd = 3, add = TRUE)
  
  box()
  
  idx.ext.obs <- which(data_d$avg_value >= GPDobs$threshold) #the index of the extreme precip of the observations
  
  ### For historic GCM
  pr_ensemble_day <- pr_ensemble_day[is.finite(pr_ensemble_day$pr_E1) & pr_ensemble_day$time_year >= minObsYear & pr_ensemble_day$time_year < minGCMYear,]
  Threshold <- percentile(pr_ensemble_day$pr_E1, 1 - obs_pctl)
  
  GPDhist <- gpd(sort(pr_ensemble_day$pr_E1),threshold=Threshold)

  gcmh_pctl <- obs_pctl
  gcmh_extr <- rgpd(10000,GPDhist$par.ests[1], GPDhist$threshold, GPDhist$par.ests[2])
  
  ii  <- order(gcmh_extr)
  gval <- rbind(gcmh_extr)[,ii]
  mnx <- min(gcmh_extr)
  mxx <- max(gcmh_extr)
  
  #hist(gcmh_extr,freq=FALSE, xlim=c(mnx,mxx), xlab= "Precip (mm/day)",main = paste("ENS Hist ", modelName, "Pcp Above ", round(Threshold,2),"mm (p",round(100*gcmh_pctl,2),")", " - ",round(GPDhist$par.ests[1],2),round(GPDhist$par.ests[2],2)))
  hist(gcmh_extr,freq=FALSE, xlim=c(mnx,mxx), xlab= "Precip (mm/day)",main = paste("ENS Hist ", modelName, "Pcp Above ", round(Threshold,2),"mm (p",round(100*gcmh_pctl,2),")","(par.ests ",round(GPDhist$par.ests[1],2),round(GPDhist$par.ests[2],2),")"),
               cex.lab=0.7, cex.axis=0.7, cex.main=0.7, cex.sub=0.7)
  
  
  x <- gcmh_extr
  curve(dgpd(x, xi = GPDhist$par.ests[1], Threshold, beta = GPDhist$par.ests[2]), col = 2, lty = 2, lwd = 3, add = TRUE)
  box()
  
  #########################################################################################

  ens  <- subset(pr_ensemble_day_f[,c(2,3,5)],time_year >= minGCMYear & time_year <= maxGCMYear)
  p.ens <- ens[ens$pr_E1 >1,]
  ### For Future GCM
  Threshold <- percentile(p.ens$pr_E1, 1 - obs_pctl)
  GPDens <- gpd(sort(p.ens$pr_E1),threshold=Threshold)
  
  gcmf_pctl <- obs_pctl
  gcmf_extr <- rgpd(10000, GPDens$par.ests[1], GPDens$threshold, GPDens$par.ests[2])
  
  ii  <- order(gcmf_extr)
  gval <- rbind(gcmf_extr)[,ii]
  mnx <- min(gcmf_extr)
  mxx <- max(gcmf_extr)
  
  #hist(gcmf_extr,freq=FALSE,xlim=c(mnx,mxx), xlab= "Pcp (mm/day)", main = paste("ENS Future Pcp Above ",round(GPDens$threshold,2),"mm (p",round(100*gcmf_pctl,2),")", " - ",round(GPDens$par.ests[1],2),round(GPDens$par.ests[2],2)))
  hist(gcmf_extr,freq=FALSE,xlim=c(mnx,mxx), xlab= "Pcp (mm/day)", main = paste("ENS Future Pcp Above ",round(GPDens$threshold,2),"mm (p",round(100*gcmf_pctl,2),")","(par.ests ",round(GPDens$par.ests[1],2),round(GPDens$par.ests[2],2),")"),
       cex.lab=0.7, cex.axis=0.7, cex.main=0.7, cex.sub=0.7)
  
  
  x <- gcmf_extr
  curve(dgpd(x, xi = GPDens$par.ests[1], GPDens$threshold, 
             beta = GPDens$par.ests[2]), col = 2, lty = 2, lwd = 3, add = TRUE)
  box()
  
  #### update the parameters of the Synth Bootstrapped data. #############################################
  GPDbeta  <- GPDobs$par.ests[2]* (GPDens$par.ests[2]/GPDhist$par.ests[2])
  GPDalpha <-  GPDobs$par.ests[1] + (GPDens$par.ests[1]-GPDhist$par.ests[1]) #^(GPDhist$par.ests[2]/GPDobs$par.ests[2])
  GPDThresh <-  GPDobs$threshold * (GPDens$threshold / GPDhist$threshold)
  
  xObs_new <- rgpd(10000, GPDalpha, GPDThresh, GPDbeta) #find the new values using the new extreme values
  obsorig <- sort(obs_extr)
  obsnew <- sort(xObs_new)
  
  obs_extr_new <- data_d_obs$avg_value[data_d_obs$avg_value >= GPDobs$threshold]
  obs_extr_new <- sapply(1:length(obs_extr_new), function(a) obsnew[which.min(abs(obsorig - obs_extr_new[a]))])
  
  xObs_new <- sort(xObs_new)
  xObs_new <- sapply(1:100, function(a) { 
    val <- a*length(obs_extr)/100 
    linapprox(val, 1:length(xObs_new), xObs_new) })
  
  xObs <- sort(obs_extr)
  xObs <- sapply(1:100, function(a) { 
    val <- a*length(xObs)/100 
    linapprox(val, 1:length(xObs), xObs) })
  
  gcmh <- sort(gcmh_extr)
  gcmh <- sapply(1:100, function(a) linapprox(a*(length(gcmh)/100), 1:length(gcmh), gcmh))

  gcmf <- sort(gcmf_extr)
  gcmf <- sapply(1:100, function(a) linapprox(a*(length(gcmf)/100), 1:length(gcmf), gcmf))
  
  plot(xObs,type="l",lty=1,main="Extremes Corrected gpd Update",col="gray40",lwd=3,ylim=c(min(xObs_new,xObs, gcmf, gcmh),max(xObs_new,xObs, gcmf, gcmh)), ylab = "Precip (mm)")
  lines(xObs_new,type="l",lty=1,col="red",lwd=3)
  lines(gcmh,col="gray10")
  lines(gcmf,col="red",lty=2)

  grid()
  legend("topleft",legend=c("obs","New Obs","Hist","Ens"),
         lty=c(1,1,1,2),lwd=c(3,3,1,1),col=c("gray40","red","gray10","red"))
  
  obsOld <- data_d
  data_d$avg_value[idx.ext.obs] <- obs_extr_new #update the avg_value fields
  
  Ratios <- data_d$avg_value[idx.ext.obs]/data_d_obs$avg_value[idx.ext.obs] #find the ratios of the historic and new
  sta_ind <- grep("pr_",colnames(data_d))
  OtherStations <- data_d[idx.ext.obs,sta_ind]  #use the ratio to update the other station values
  data_d[idx.ext.obs,sta_ind] <- Ratios * OtherStations  #replace those new values,mutliply Ratio_Vector * data_table 
  
  obsY <- aggregate(data_d$avg_value, by = list(data_d$Year), sum,na.rm=TRUE)
  obsYold <- aggregate(obsOld$avg_value, by = list(obsOld$Year), sum,na.rm=TRUE)
  
  #overwrite the observed values 
  save(data_d,file=paste(main_folder,results_folder,"/",paste(modelName,"observed", "pr",futures,"day_ExtremesCorrected_OldMethod.Rda",sep="_"),sep=""))
  dev.off()
  
}  #end extreme gpd
