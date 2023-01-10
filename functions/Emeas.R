
#### Emeas function ####


Emeasfun <- function(numtrials, numsam, length, dMeas, r1, r2, r3, la){

  for (j in 1:numtrials) {
    #vector on ones
    Isoslope <- rep(1,numsam)
    
    for (n in 1:numsam) {
      if(n == 1) Isoslope[n] <- (((abs(dMeas[n+1]-dMeas[n]))/(0.5*(length[n+1]+length[n]))))
      else
        if(n == numsam) Isoslope[n] <- (((abs(dMeas[n]-dMeas[n-1]))/(0.5*(length[n]+length[n-1]))))
        else
          Isoslope[n]=((((( (abs(dMeas[n+1] - dMeas[n])) / ( 0.5*(length[n+1]+length[n]) ))   
                          +   ( (abs(dMeas[n] - dMeas[n-1])) /(0.5* (length[n]+length[n-1]))))/2)))
    }
    
    RElength = r2*rnorm(numsam,mean = 0, sd = 1);		#generates vector of random length errors from normal distribution
    LengthError = Isoslope*RElength 	#calculates isotope error caused by length errors (eq **.**)  
    LEdMeas = LengthError + dMeas   	#adds LengthError to the measured data vector    
    
    #The following loops check to see if the 'error' isotope value is above or below 
    #the highest or lowest adjacent values in dMeas. 
    #as this would be a physically impossible result of sample length error
    
    for (n in 1:numsam) {
      
      if(n == 1){
        if(LEdMeas[n] < min(dMeas[n:n+1])){
          LEdMeas[n] <- (min(dMeas[n:n+1]))
        } else if(LEdMeas[n] > max(dMeas[n:n+1])){
          LEdMeas[n] <- (max(dMeas[n:n+1]))
        } else LEdMeas[n] <- LEdMeas[n]
      } else if(n == numsam){
        if(LEdMeas[n] > max(dMeas[n:n-1]))
          LEdMeas[n] <- min(dMeas[n:n-1])
        else if(LEdMeas[n] > max(dMeas[n:n-1])){
          LEdMeas[n] <- max(dMeas[n:n-1])
        } else LEdMeas[n] <- LEdMeas[n]
      } else
        if(LEdMeas[n] < min(dMeas[n-1:n+1])){
          LEdMeas[n] <- min(dMeas[n-1:n+1])
        } else if(LEdMeas[n] > max(dMeas[n-1:n+1])){
          LEdMeas[n] <- max(dMeas[n-1:n+1])
        } else LEdMeas[n] <- LEdMeas[n]
        
    }
    
    # The following section calculates the depth-dependent isotope error.  It fits a 
    # cubic spline to the the measured data dMeas, 
    # and creates a vector of interpolated delta values for each unit length on the x-axis. 
    # A new vector is created that is shifted la 
    # units, and the two vectors are subtracted to give DELTA-delta values reflecting the 
    # difference between the isotope ratio at the outside 
    # enamel surface and the enamel-dentine junction.  These values are then multiplied by 
    # the sampling depth uncertainty to give a per mil uncertainty.  
    
    totallength = rep(1,numsam)
    totallength[1] = length[1]
    
    for (n in 2:numsam) {
      totallength[n] <- totallength[n-1] + length[n]
    }
    
    addbefore = dMeas[1]*rep(1,la)
    
    addafter = dMeas[numsam]*rep(1,la)
    
    xx <- 0:totallength[numsam]
    
    dInterp = spline(totallength,dMeas,xout = xx)
    dInterp <- dInterp$y
    dInterp
    
    dInterpShift = c(addbefore,dInterp)  
    
    dInterp = c(dInterp, addafter)  
    
    Dd = dInterpShift-dInterp  
    
    Deltadelta = rep(1, numsam) 
    
    for (n in 1:numsam) {
      Deltadelta[n] = Dd[totallength[n]]
    }
    
    
    REdepth = (r3*rnorm(numsam, mean = 0, sd = 1))/la          #generates vector of random depth errors
    
    DepthError = Deltadelta*REdepth
    
    
    #end of depth error section   
    
    AnalysisError = r1*rnorm(numsam, mean = 0, sd = 1) 
    length(AnalysisError)
    
    SqSumError = (LengthError^2 + DepthError^2 + AnalysisError^2)^0.5  
    length(SqSumError)
    E = SqSumError %*% SqSumError 
    
    Edist[j] <<- E  
    
    allTrials[j] = SqSumError  
    
    dMeasError = dMeas + SqSumError  
    
  }   
  
  
  totallength = rep(1,numsam) 
  totallength[1] = length[1]  
  
  for (n in 2:numsam) {
    totallength[n] = totallength[n-1]+length[n]
  }
  
  output <- as.data.frame(cbind(totallength, allTrials, dMeas, dMeasError))
  return(output)
}

