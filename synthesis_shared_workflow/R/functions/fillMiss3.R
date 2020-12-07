#' Function for filling missing values
#' @description This function fills missing values and is only a minor adaptation
#' of the fillMiss function from the waterData package. The description of the
#' parameters from the original function are taken directly from the waterData package.
#'
#' @param dataset The data frame containing the desired variable
#' @param var The column name of the variable to be gap-filled
#' @param block The size of the largest block of missing data that the function will fill-in
#' @param pmiss The maximum amount of the missing data that can be missing in the dataset for fill-in procedure to be performed.
#' @param model The type of structural time series model, see StructTS. The default value is trend. If level is used, the results of fillMiss, which by default applies a fixed-interval smoothing to the time series, tsSmooth, will be very close to linear interpolation.
#' @param smooth A logical that indicates whether or not to apply tsSmooth to the structured time series.
#' @param plot A logical that indicates whether or not to plot the data
#'
#' @return Returns gap-filled data for the selected variable
#'
#' @export
#===============================================================================
#Fill missing discharge data
#Based on my fillMiss2 function (4/30/2019) that I used for my SFS 2019 analysis,
#which is a minor adaptation of the fillMiss function from the waterData package.
#NOTE: COME BACK AND MAKE ERROR MESSAGES WORK
#Created 8/23/2019
#===============================================================================
fillMiss3 <- function(dataset, var, block = 30, pmiss = 40, model = "trend",
  smooth = TRUE, plot = TRUE, ... ){
  #Find missing values in the dataset
    pck <- is.na(dataset[, var])

  #Catch if the dataset contains only NA values
  if(all(pck == TRUE) == FALSE){
    #-------------------------------------------------
    #Calculate information about the number and length of gaps
    #-------------------------------------------------
      if(sum(pck) > 0){
        #Calculate the percentage of missing values in the dataset
          percent <- round((sum(pck) / length(dataset[, var])) * 100,
            digits = 2)

        # #Output message to the user for the number of missing values
        #   message(paste("There are ", sum(pck), "missing values for",
        #     dataset, "or", percent, "percent of the data.", sep=" "))

        #Calculate the largest consecutive set of missing data
          rles <- rle(is.na(dataset[, var]))
          max.mis <- max(rles$lengths[rles$values])

        # #Output message for largest consecutive set of missing data
        #   message(paste("The maximum block of missing", dataset$staid[1],
        #     "data is", max.mis, "days long.", sep=" "))

      } else {
          # message(paste("No missing values for", dataset, sep = " "))
          percent <- 0
          max.mis <- 0
        }

    #-------------------------------------------------
    #Filling gaps in the dataset
    #-------------------------------------------------
      #Check to see if there are too many missing values for gap filling
        if(percent >= pmiss | max.mis >= block) {
          # message(paste("Too much missing data for", dataset,
          #   "Cannot fill in missing values.", sep=" "))
        } #End if statement

      #If there are not too many missing values then proceed
        if(percent < pmiss & max.mis < block){
          my.series <- window(dataset[, var])

          #First value can't be NA for StructTS, replace NA with nearest value P.S. 2019
            if(is.na(my.series)[1]){my.series[1] <- zoo::na.locf(my.series,
              option = "nocb", na.remaining = "rev")[1]}

          #Fit a structural time series model (includes error check to see
          #if the model can be fit. If not, use type = level which I think
          #approximates linear interpolation)
            struct_try <- try(StructTS(my.series, type = model), silent = TRUE)

            #If there is no error fit the default model = "trend" option
              if(class(struct_try)[1] != 'try-error'){
                my.struct <- StructTS(my.series, type = model)
              } #End if statement

              if(class(struct_try)[1] == 'try-error'){
                message(paste0("type = level was used for ", var))
                my.struct <- StructTS(my.series, type = "level")
              } #End if statement

          #Perform smoothing (if chosen)
            if(smooth) fit <- tsSmooth(my.struct) else fit <- fitted(my.struct)

          #Replace fit values that exceed the minimum and maximum values for the dataset
            data_min <- min(dataset[, var], na.rm = TRUE)
            data_max <- max(dataset[, var], na.rm = TRUE)            
            
            fit[fit[, 1] < data_min, 1] <- data_min
            fit[fit[, 1] > data_max, 1] <- data_max
            
          #-------------------------------------------------
          #Plotting
          #-------------------------------------------------
            if(plot == TRUE){
              plot(my.series, typ = "l", lwd = 4, xlab = "Observation",
                ylab = "Observed and estimated times series", ...)

              #Add fitted data
                lines(fit[, 1], col = "green")

              #Add legend
                leg.txt <- c("Observed values", "New time series")
                legend("topleft", leg.txt, col = c("black", "green"), lwd = c(4, 1),
                  bty = "n", ncol = 2, cex = 0.8)
            } #End if statement

          #Gap-filling
            dataset$filled <- dataset[, var]
            dataset$filled[pck] <- fit[pck, 1]

          # #Output message of how many gaps were filled
          #   message(paste("Filled in", sum(pck), "values for", dataset, sep = " "))

        } #End if statement

    return(dataset[, "filled"])
  } #End if statement

  #If the data contains only NA values
    if(all(pck == TRUE) == TRUE){
      message("Input variable contained only NA values so only NA values were returned")
      return(rep(NA, nrow(dataset)))

    } #End if statement

} #End fillMiss3 function
