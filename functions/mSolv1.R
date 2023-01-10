
PasseyInverse<-function(Length,dMeas,depth,la,lm,maxlength,minlength,mindepth,df,nsolxns, finit, openindx, avelength, r1, r2, r3, maxratio, minratio, stdev, numsam){
  MEST<-list()
  DPE<-c()
  S<-list()
  for(n in 1:nsolxns){
    print(n)

    #### SECTION 2 ####
    # addition of random error and definitions
    # add random error to length vector, by matrix multiplying the length precision
    # with a vector of the length numsam of random numbers between 0 and 2
    
    rlength <- r2*rnorm(numsam)	
    Length <- round(rlength) + Length #round to integers and add error
    
    #check that lengths are between min and max threshold specified previously
    for(a in 1:numsam){
      if(Length[a]>maxlength){
        Length[a]<-maxlength
      }else{
        if(Length[a]<minlength){
          Length[a]<-minlength
        }else{ 
          Length[a]<-Length[a]
        }
      }
    }
    rdepth <- r3*rnorm(numsam)
    rdepth<-round(rdepth)
    depth2 <- round(rdepth) + depth
    
    for(z in 1:numsam){
      if(depth2[z]>la){
        depth2[z]<-la
      }
      if(depth2[z]<mindepth){
        depth2[z]<-mindepth
      }
    }
    depth<-depth2
    
    #Determines the number of m's distal and proximal to those that directly correspond with d's
    #numbefore reflects the m's that are sampled into at the beginning of the profile because of sampling depth
    #numafter reflects the m's that contribute to the isotope values of the
    #final samples in open-ended cases.
    
    #round to nearest higher integer
    numbefore<-ceiling(la/avelength)
    numafter<-ceiling((lm-openindx)/avelength) +1
    lengthbefore<-avelength*rep(1,numbefore)
    lengthafter<-avelength*rep(1,numafter)
    
    #more definitions
    fmat<- 1 - finit
    numcol <<- numbefore + numsam + numafter
    Length2<-c(lengthbefore,Length,lengthafter,lm)
    depth<-c(depth,lengthbefore)
    
    #### SECTION 3.1 BINARY MATRIX ####
    # constructing the averaging matrix

    B<-matrix(nrow=Length2[1],ncol=numcol,0) # make matrix (fill,nrows,ncols) full of zeroes
    B[,1]<-1 # fill first column with ones
    
    for(m in 2:numcol){
      F<-matrix(nrow=Length2[m],ncol=numcol,0)
      F[,m]<-1
      B1<-rbind(B,F)
      B<-B1
    }
    
    # expand with additional non contributing rows
    F<-matrix(nrow=(((sum(Length2)-openindx)+1):sum(Length2))-nrow(B),ncol=numcol,0)
    B<-rbind(B,F)
    B[((sum(Length2)-openindx)+1):sum(Length2),]<-c(0)
    B<-B[,1:(numcol-1)] # remove last column
    
    #### SECTION 3.2 MATURATION AVERAGE OF BINARY MATRIX ####
    o<-1
    AB <- c(finit*(colMeans(B[o:(o+lm-1),])) + fmat*(colMeans(B[o:(o+lm-1),])))
    AB <- AB/(sum(B[o,])) # make seed row for AB matrix, will be deleted later
    for(o in 1:(sum(Length2)-lm)){
      p<-c(finit*(B[o,]) + fmat*(colMeans(B[o:(o+lm-1),])))
      if(sum(p)==0){
        p<-p
      }else{
        p<-(p/sum(p)) # makes all rows sum to one
      }    
      AB<-rbind(AB,p)
    }
    AB<-AB[2:nrow(AB),] # remove seed row
    
    #### CUMULATIVE LENGTH VECTOR ####
    
    clength<-c(Length2[1])
    for(q in 2:numcol){
      cl<-Length2[q] + clength[q-1]
      clength<-c(clength, cl)
    }
    #### FINAL CALCULATION OF A ####
    
    A<-AB[1,]
    for(k in (numbefore+1):(numsam+numbefore)){
      E<-AB[1,]
      for (j in 0:(depth[k-numbefore] - 1)){
        e<-colMeans(AB[(clength[k-1]-j+1):(clength[k] - j),])
        E<-rbind(E,e)
      }
      E<-E[2:nrow(E),]
      meanE<-colMeans(E)
      A <- rbind(A,meanE)
    }
    A<-A[2:nrow(A),]
    
    #### INVERSION ####
    
    I<-diag(numsam)
    dMeasr<-dMeas + r1*rnorm(numsam) #check other rnorm above
    NB<-numbefore
    Na<-numafter - 1
    mm <- matrix(nrow=numsam+numbefore+numafter-1,ncol=1,1)
    mm[1,]<-((maxratio-minratio)*runif(1))+minratio
    
    for(x in 2:(numsam+numbefore+numafter-1)){
      mm[x,]<-mm[x-1,]+stdev*rnorm(1)
      while (mm[x,]>maxratio){
        mm[x,]<-mm[x-1,]+stdev*rnorm(1)
      }
      while(mm[x,]<minratio){
        mm[x,]<-mm[x-1,]+stdev*rnorm(1)
      }
    }
    
    # end of reference vector determination
    
    AA<-A%*%Conj(t(A))
    epsilon<-df*I
    AAep<-AA+epsilon
    test<-matrix.power(AAep,-1)
    mEst<-mm + (Conj(t(A))%*%(test))%*%(dMeasr - A%*%mm)
    MEST[[n]]<-mEst 
    dPred<-A%*%mEst	
    dpe<-dPred - dMeas
    DPE[n] <<- Conj(t(dpe))%*%dpe
    S[[n]]<- dpe           
  }
  
  dim <-nrow(MEST[[1]])
  
  # write objects for function output
  
  totallength<-matrix(1,nrow=numsam,ncol=1) 
  totallength[1,]<-Length2[1]  
  
  xx <- avelength*rep(1,numbefore)
  zz <- avelength*rep(1,numafter-1)
  totallength <- c(xx,totallength,zz)
  
  vec1 <- dMeasr[1]*rep(1,numbefore)
  vec2 <- dMeasr[numsam]*rep(1,numafter-1)
  dMeasd <- c(vec1,dMeasr,vec2)
  
  for(n in 2:(numsam+numbefore+numafter-1)){
    totallength[n] <- c(totallength[n-1]+Length2[n])
  }

  MEST_df <- do.call(cbind, MEST)
  
  solvout <- cbind(totallength, dMeasd, dMeas, MEST_df)
  
  ntrials=1:nsolxns
  tcols=paste("trial",ntrials,sep="")
  
  colnames(solvout) = c("totallength", "dMeasd","dMeas", tcols)
  solvout <<- as.data.frame(solvout)
  
}

