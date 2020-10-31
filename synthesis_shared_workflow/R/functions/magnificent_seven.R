#===============================================================================
#Functions for calculating the 'magnificent 7' from Archfield et al. (2015) An objective and...
#needs the lmomco library for calculating L-moments
#Created 12/4/2017
#===============================================================================
#------------------------------------------------- 
#Calculating amplitude and phase
#------------------------------------------------- 
  get_seasonality <- function(timeseries, var, standardize){
    #Calculating the decimal year
      decimal_year <- as.numeric(timeseries[, "Year"]) + (timeseries[, "DOY"] / 365.25)
      
    #Standardize 
      ifelse(standardize == "yes", standard <- scale(timeseries[, var], center = TRUE, scale = TRUE), 
        standard <- timeseries[, var])
        
    #3) Use linear model to fit 
      x_mat <- cbind(1, sin(2 * pi * decimal_year), cos(2 * pi * decimal_year))
      seasonfit <- .lm.fit(x_mat, standard)
      b1 <- as.vector(coef(seasonfit)[2])
      b2 <- as.vector(coef(seasonfit)[3]) 
      
    #Now compute the amplitude and phase of the seasonal signal
      amplitude <- round(sqrt((b2^2) + (b1^2)), digits = 2)
      
    #phase<-round(atan((-seasonB)/seasonA),digits=2)
      MaxDay <- function(b1, b2){
        version1 <- 365.25 * ((pi / 2) - atan(b2 / b1)) / (2 * pi)
        version2 <- 365.25 * ((pi/2)-pi-atan(b2/b1))/(2*pi)
        MaxDay <- if(b1 > 0) version1 else 365.25 + version2
        MaxDay <- if(b1 == 0 & b2 > 0) 365.25 else MaxDay
        MaxDay <- if(b1 == 0 & b2 < 0) 365.25/2 else MaxDay
        MaxDay <- if(b1 == 0 & b2 == 0) NA else MaxDay
        return(MaxDay)
      } #End MaxDay function
        
      phase <- MaxDay(b1,b2)
      
    #Combining them all together
      get_seasonalityv <- cbind(amplitude, phase)
      
    return(get_seasonalityv)
  } #End get seasonality function

#-------------------------------------------------
#Calculating the AR(1) correlation coefficient 
#-------------------------------------------------
  ar_fun <- function(timeseries, var){
    #Adding a column for the month
      timeseries$Month <- as.numeric(format( as.Date(timeseries[, "DOY"] - 1, origin = 
        paste(timeseries[, "Year"], "-01-01", sep = "")) , format = "%m"))   
  
    #"deseasonalize" the data
      #Calculating the long term monthly mean	
        monmeans <- aggregate(timeseries[, var], list(timeseries[, "Month"]), 
          FUN = mean, na.rm = TRUE)
        
      #Merging together the timeseries and long term monthly means
        mon_merge <- merge(timeseries, monmeans, by.x = "Month", by.y = "Group.1")
      
      #Subtracting long term monthly mean from each day
        mon_merge$deseason <-  mon_merge[, var] - mon_merge[, "x"]
        
      #Ordering the merged data
        ordered <- mon_merge[order(mon_merge$Year, mon_merge$DOY), ]
        
    #Standardizing to have 0 mean and unit variance
      stand <- scale(ordered[, "deseason"], center = TRUE, scale = TRUE)
  
    #Fitting the autoregressive model  
      ar_mod <- ar(stand, aic = FALSE, order.max = 1, method = "yule-walker")
      
    #Getting the AR(1) correlation coefficient
      ar_cor <- round(ar_mod$ar, digits = 3)
      
    return(ar_cor)
  } #End function for calculating AR(1) 

#-------------------------------------------------
#Calculating the "Magnificent 7"
#-------------------------------------------------
  mag7_fun <- function(timeseries, var, standardize){
    #Calculating L-moments 
      l_mom <- lmomco::lmom.ub(timeseries[, var])
      
      lam1 <- round(lmomco::l_mom$L1, digits = 2) #Mean
      tau2 <- round(lmomco::l_mom$LCV, digits = 2) #Coefficient of L-variation
      tau3 <- round(lmomco::l_mom$TAU3, digits = 2) #L-skew
      tau4 <- round(lmomco::l_mom$TAU4, digits = 2) #L-kurtosis
      
    #Calculating phase and amplitude
      seasonality <- get_seasonality(timeseries, var, standardize)
      
    #Calculating the autoregressive lag 1 correlation coefficient
      ar1 <- ar_fun(timeseries, var)
      
    #Combining them together
      mag7 <-data.frame(lam1, tau2, tau3, tau4, seasonality, ar1)
      
    return(mag7)
      
  }#End mag7_fun    

