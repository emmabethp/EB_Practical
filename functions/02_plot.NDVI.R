plot.NDVI <- function(dat, fit = NA, ...){
  
  if(!is.null(dat)){ # Begin if/else statement
    
    # Base plot of the data points
    plot(dat$age, dat$NDVI, ylab = "NDVI", xlab = "Postfire age (Years)")
    
    if(!is.na(fit[1])){ #Begin inner if statement
      
      # Overlay the fitted model on the plot
      lines(dat$age, predict(fit, list(x = dat$age)), col = 'skyblue', lwd = 3)
      
    } # End inner if statement
    
  } else {
    print("plot.NDVI: input data not provided or invalid")
  } # End if/else statement
}