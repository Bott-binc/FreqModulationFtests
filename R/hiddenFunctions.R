library("multitaper")

# Sine Related -------------------------

#' Single Sinusoidal Taper
#'
#' This returns only the one sin taper that you are interested in at the single point
#'
#' @param N Length of time series
#' @param k Number of tapers
#' @param t Specific time point that you are interested in t starts at zero
#'
#' @return Only the one sin taper that you are interested in at the single point
sinTaperSingle <- function(N, k, t){
  return(sqrt(2/(N+1)) * sin(((k+1)*pi*(t + 1))/(N+1)))
}



#' Sinusoidal Taper Matrix
#'
#' This is the same as the one Ben and Yatharth added to multitaper but is returned in the same orientation as the one used in Percival and Walden
#'
#' @param N Length of Time series
#' @param k Number of tapers starting at 0
#'
#' @return Matrix of tapers for the given N and k values, N columns starting at n = 0

sineTaperMatrix <- function(N, k){
  stopifnot(N >= 8, k >= 1)

  taper <- sin(outer(1:N, 1:k)*pi/(N+1))
  taper <-  taper*sqrt(2/(N+1))

  out <- as.matrix(taper)

  return(out)
}

#' Non Centered eigenCoef for Sine Tapers
#'
#' @param n Total number of observations
#' @param k Number of tapers
#' @param Xt Time series
#' @param f Specific frequency you are looking at
#' @param deltat default is 1 unless otherwise inputted
#' @param sineRet = FALSE , allows for the tapers to be returned for computational efficiency
#'
#' @return Vector of the eigenCoef for each k 0, ..., K-1
eigenCoefSine <- function(n, k, Xt, f, deltat = 1, sineRet = FALSE){

  EigenCoef <- vector(length = k)
  tapers <- sineTaperMatrix(N = n, k = k)
  for(i in 1:k){
    vec <- t(tapers[,i]) %*% as.matrix(Xt * exp(-1i*2*pi*f*0:(n-1)*deltat))
    EigenCoef[i] <- sum(vec)
  }
  if(sineRet){
    return(list(EigenCoef = EigenCoef, Sine = tapers))
  }
  else{
    return(EigenCoef)
  }
}



#' Non Centered eigenCoef Using DPSS and fft for speed and zero padding ability
#'
#' @param N Total number of observations
#' @param k Number of tapers
#' @param Xt Time series
#' @param passInTaper leave as null unless you are passing in a taper matrix for sine
#' @param deltat default is one, changes the nyquist if altered
#' @param returnSineMat if you want to use Sine taper matrix outside function can return if TRUE
#' @param pad if true, will automatically zero padd the fft
#' @param penalty 1 is no penalty , 0.2  would give seq(from =  1, to =  1/(0.2*k, length.out = k)
#' penalty to each respective taper
#'
#' @return Matrix of $EigenCoef for each frequency, as well as the frequencies under $Freq in the list
#' if returnSineMat = TRUE, will also return $SineMat
eigenCoefSineFFT <- function(N, k, Xt, deltat = 1, passInTaper = NULL,
                             returnSineMat = FALSE, pad = TRUE,
                             penalty = 1, penaltyType = "ScaledExp"){
  if(penaltyType %in% c("ScaledExp", "Cos", "Clip", "CustomE", "CustomO", "Custom")){
    #do nothing and continue
  }else{
    stop("No penalty type that is available was selected, please check documentation")
  }
  if(class(penalty) == "character"){

    stop("penalty is reserved for the weighting, please put the type of penalty in penaltyType arg")
  }
  EigenCoef <- matrix(nrow = N,ncol = k)
  if(is.null(passInTaper)){
    sine <- FreqModulationFtests:::sineTaperMatrix(N, k)
  }
  else{
    sine = passInTaper
  }

  if(!pad){

    for(i in 1:k){
      vec <- sine[,i] * as.matrix(Xt)
      EigenCoef[,i] <- fft(vec)
    }
    if(penalty != 1){
      EigenCoef <- t(apply(EigenCoef, MARGIN = 1, FUN = function(x){x/seq(from = 1, to = penalty*k, length.out = k)}))
    }
    ret <- data.frame(EigenCoef = EigenCoef, Freq = c(seq(from = 0, to = 1/(2*deltat),
                                                          by = 1/(N*deltat)),
                                                      seq(from = (-1/(2*deltat) + 1/(N*deltat)),
                                                          to = (-1/(N*deltat)),
                                                          by = 1/(N*deltat))))
  }else{


    nFFT <- 2^ceiling(log2(2*N))
    nextPowerOfTwo <- 2^ceiling(log2(2*N)-1)
    taper <- sine * Xt
    pad <- rbind(taper, matrix(0, nrow = nFFT - N, ncol = k))
    EigenCoef <- mvfft(pad)[1:(nextPowerOfTwo + 1), ,drop = FALSE]
    if(penaltyType == "Cos"){
      weightMat <- matrix((1/2*cos((pi*1:k)/k) + 1/2), nrow = nrow(EigenCoef), ncol = k, byrow = TRUE)
      EigenCoef <- EigenCoef * weightMat
    }
    else if(penaltyType == "Clip"){
      clipPoint <- penalty#2/3
      weight <- c(rep(1, length.out = floor(k*clipPoint)), rep(0, length.out = k - floor(k*clipPoint)))#c(rep(1, length.out = penalty), rep(0, length.out = k - penalty)) #
      weightMat <- matrix(weight, nrow = nrow(EigenCoef), ncol = k, byrow = TRUE)
      EigenCoef <- EigenCoef * weightMat
    }
    else if(penaltyType == "CustomO"){
      # this is coded to remove the even right now
      odd <- rep(c(1,0), length.out = k) # this will keep all the odd ones and remove even ones
      clipPoint <- penalty#2/3
      weight <- c(rep(1, length.out = floor(k*clipPoint)), rep(0, length.out = k - floor(k*clipPoint)))
      weightMat <- matrix(odd * weight, nrow = nrow(EigenCoef), ncol = k, byrow = TRUE)
      EigenCoef <- EigenCoef * weightMat
    }
    else if(penaltyType == "CustomE"){
      # this is coded to remove the even right now
      odd <- rep(c(0,1), length.out = k) # this will keep all the odd ones and remove even ones
      clipPoint <- penalty#2/3
      weight <- c(rep(1, length.out = floor(k*clipPoint)), rep(0, length.out = k - floor(k*clipPoint)))
      weightMat <- matrix(odd * weight, nrow = nrow(EigenCoef), ncol = k, byrow = TRUE)
      EigenCoef <- EigenCoef * weightMat
    }
    else if(penaltyType == "Custom"){
      # this is coded to remove the even right now
      #weight <- rep(c(1,0), length.out = k) # this will keep all the odd ones and remove even ones
      #weightMat <- matrix(weight, nrow = nrow(EigenCoef), ncol = k, byrow = TRUE)
      weight <- c(rep(0, penalty), rep(1, k-penalty))
      weightMat <- matrix(weight, nrow = nrow(EigenCoef), ncol = k, byrow = TRUE)
      EigenCoef <- EigenCoef*weight #* weightMat

    }
    else if(penaltyType == "ScaledExp"){
      if(penalty == 1){
        #do nothing as no penalty is selected
      }else{
        #EigenCoef <- t(apply(EigenCoef, MARGIN = 1, FUN = function(x){x/seq(from = 1, to = penalty*k, length.out = k)}))
        weight <- 1/seq(from = 1, to = penalty*k, length.out = k)
        EigenCoef <- EigenCoef * matrix(weight, nrow = nrow(EigenCoef), ncol = k, byrow = TRUE)
      }
    }
    ret <- data.frame(EigenCoef = EigenCoef, Freq = seq(from = 0, to = 1/(2*deltat),
                                                        by = 1/(2*nextPowerOfTwo*deltat)))
  }

  ret <- ret[order(ret$Freq, decreasing = FALSE),]
  eSpecReorder <- as.matrix(ret[,1:(ncol(ret)-1)])
  rownames(eSpecReorder) <- ret$Freq
  colnames(eSpecReorder) <- 1:ncol(eSpecReorder)

  if(!returnSineMat){

    return(list(EigenCoef = eSpecReorder, Freq = ret$Freq))
  }
  else{

    return(list(EigenCoef = eSpecReorder, Freq = ret$Freq,SineMat = sine))
  }
}

#' Standard Inversion for Polynomial Projection Filter Sine for a freq f
#'
#' @param xt time series
#' @param N Total number of observations
#' @param k Number of tapers
#' @param f specified frequency you want to look at
#' @param FFT if you want to use fft to speed up or not have delta t = 1 with zero padding
#' @param deltat  = 1 if not specified, but is only used if not using FFT
#'
#' @return Standard Inverse vector for frequency f
standardInverseSineSingleFreq <- function(xt, N, k, f, deltat = 1, FFT = FALSE){

  v <- sineTaperMatrix(N, k)
  if(!FFT){
    # Note this could be maybe sped up with eigenCoefSine if needed in future
    Xf <- xt * exp(-1i*2*pi*f*0:(N-1)*deltat)
    stdInverse <- v %*% t(v) %*% as.matrix(Xf)
  }
  else{

    Y <- eigenCoefSineFFT(N, k, xt, deltat = deltat)
    stdInverse <- v %*% as.matrix(Y$EigenCoef[which(Y$Freq == f),])
  }
  return(stdInverse)
}



#' Standard Inverse for Sine tapers
#'
#' @param xt time series
#' @param N Total number of observations
#' @param k Number of tapers
#' @param retSineTapers Will return calculated sine taper matrix for future use
#' @param deltat default is 1 unless using different nyquist
#' @param passInSineMat If sine matrix was calculated already
#' @param passInSineUnder If sine undersampling matrix was calculated already
#' @param penalty 1 is no penalty , 0.2  would give seq(from =  1, to =  1/(0.2*k, length.out = k)
#' penalty to each respective taper
#' @param penaltyOnTapersStdInv applys the penalty to all tapers in the calculation instead of just weighting the
#' eigenCoef's
#'
#' @return returns matrix of standard inverse along with frequencies corresponding to the columns
standardInverseSine <- function(xt, N, k, deltat = 1, passInSineMat = NULL,
                                retSineTapers = FALSE,
                                passInSineUnder = NULL,
                                penalty = 1, penaltyType = "ScaledExp",
                                penaltyOnTapersStdInv = FALSE){
  if(is.null(passInSineMat)){
    v <- sineTaperMatrix(N, k)
  }
  else{
    v <- passInSineMat
  }

  if(is.null(passInSineUnder)){ #using the penalty is the same here as passing in the tapers with the penalty on them instead
    Y <- eigenCoefSineFFT(N, k, xt, deltat = deltat, passInTaper = v,
                          pad = FALSE, penalty = penalty, penaltyType = penaltyType)
  }else{
    Y <- eigenCoefSineFFT(N, k, xt, deltat = deltat, passInTaper = v,
                          pad = TRUE, penalty = penalty, penaltyType = penaltyType)
  }

  if(penaltyOnTapersStdInv){ # we add the penalty and use this new set of tapers for the standard inverse
    if(is.null(passInSineUnder)){
        penaltyV <- v %*% diag(1/seq(from = 1, to = penalty*k, length.out = k))
        stdInverse <- tcrossprod(penaltyV, Y$EigenCoef)
      }else{#reduced matrix of tapers used here
        passInSineUnder <- passInSineUnder %*% diag(1/seq(from = 1, to = penalty*k, length.out = k))
        stdInverse <- tcrossprod(passInSineUnder, Y$EigenCoef)
      }
  }else{
    if(is.null(passInSineUnder)){
      stdInverse <- tcrossprod(v, Y$EigenCoef)
    }else{#reduced matrix of tapers used here
      stdInverse <- tcrossprod(passInSineUnder, Y$EigenCoef)
    }
  }


  colnames(stdInverse) <- Y$Freq
  if(!retSineTapers){
    return(list(StdInverse = stdInverse, Freq = Y$Freq))
  }
  else{
    return(list(StdInverse = stdInverse,Freq = Y$Freq, sineTapers = v))
  }
}


#' Standard Inversion first Derivative for Polynomial Projection Filter Sine
#'
#'NOTE: this could be sped up by only calculating the fft once!! for all f
#'
#' @param xt time series
#' @param N Total number of observations
#' @param k Number of tapers
#' @param deltat default is 1 unless using different nyquist
#' @param returnSineMat will return the sine taper matrix back to the parent function
#' @param passInSineMat If sine matrix was calculated already
#' @param passInSineUnder If sine undersampling matrix was calculated already
#' @param penalty 1 is no penalty , 0.2  would give seq(from =  1, to =  1/(0.2*k, length.out = k)
#' penalty to each respective taper
#'
#' @return list of $StdInverse and $Freq which are the columns of StdInverse,
#' if returnSineMat = TRUE then will also return $SineMat
standardInverseSineDer <- function(xt, N, k, deltat = 1, passInSineMat = NULL,
                                   returnSineMat = FALSE, passInSineUnder = NULL,
                                   penalty = 1, penaltyType = "ScaledExp"){
  if(is.null(passInSineUnder)){
    FirstDir <- FirstDerSineTaper(N, k)
    if(!returnSineMat){
      if(is.null(passInSineMat)){ # no tapers will be passed to EigenCoeffft
        Y <- eigenCoefSineFFT(N, k, xt, deltat = deltat, pad = FALSE,
                              penalty = penalty, penaltyType = penaltyType)
      }
      else{ # passing in sine tapers from outside function
        Y <- eigenCoefSineFFT(N, k, xt, deltat = deltat, passInTaper = passInSineMat,
                              pad = FALSE, penalty = penalty, penaltyType = penaltyType)
      }
    }
    else{
      if(is.null(passInSineMat)){ # no tapers will be passed to EigenCoeffft
        Y <- eigenCoefSineFFT(N, k, xt, deltat = deltat, returnSineMat = TRUE,
                              penalty = penalty, penaltyType = penaltyType)
      }
      else{ # passing in sine tapers from outside function
        Y <- eigenCoefSineFFT(N, k, xt, deltat = deltat, passInTaper = passInSineMat,
                                  returnSineMat = TRUE, penalty = penalty,
                              penaltyType = penaltyType)
      }
    }
    stdInverse <- tcrossprod(FirstDir, Y$EigenCoef)
  }
  else{ # reduction
    FirstDir <- FirstDerSineTaper(N = nrow(passInSineUnder), k)

    if(!returnSineMat){
      if(is.null(passInSineMat)){ # no tapers will be passed to EigenCoeffft
        Y <- eigenCoefSineFFT(N = nrow(passInSineMat), k, xt, deltat = deltat,
                              pad = TRUE, penalty = penalty, penaltyType = penaltyType)
      }
      else{ # passing in sine tapers from outside function
        #this is the main one

        Y <- eigenCoefSineFFT(N = nrow(passInSineMat), k, xt, deltat = deltat,
                                  passInTaper = passInSineMat, pad = TRUE,
                              penalty = 0.8)
      }
    }
    else{
      if(is.null(passInSineMat)){ # no tapers will be passed to EigenCoeffft
        Y <- eigenCoefSineFFT(N = nrow(passInSineMat), k, xt, deltat = deltat,
                                  returnSineMat = TRUE, pad = TRUE, penalty = penalty
                              ,penaltyType = penaltyType)
      }
      else{ # passing in sine tapers from outside function
        Y <- eigenCoefSineFFT(N = nrow(passInSineMat), k, xt, deltat = deltat,
                                  passInTaper = passInSineMat,
                                  returnSineMat = TRUE, pad = TRUE, penalty = penalty,
                              penaltyType = penaltyType)
      }
    }

    stdInverse <- tcrossprod(FirstDir, Y$EigenCoef)
  }




  colnames(stdInverse) <- Y$Freq
  if(!returnSineMat){
    return(list(StdInverse = stdInverse, Freq = Y$Freq))
  }
  else{
    return(list(StdInverse = stdInverse, Freq = Y$Freq, SineMat = Y$SineMat))
  }
}


#' Composite inverse for Sine Tapers using Gram schmidt hat matrix
#'
#' @param N Total number of observations
#' @param k Number of tapers
#' @param xt Time series
#' @param f Specified frequency
#' @param p Highest degree polynomial
#' @param polynomialPart If you only want the polynomial part GH^tY == G_pC^hat_p used in ftest
#' @param returnHG Will return H and G along with the Z/includes in list if poly is true
#' @param deltat  = 1 if not specified
#'
#' @return Vector of compInverse for freq f, if polynomialOnly, it will return
#' list of cHat and polynomialPart of composit inverse which is G_p*cHat
compInverseSine <- function(N, k, xt, f, p, deltat = 1, polynomialPart = FALSE, returnHG = FALSE){

  Y <- eigenCoefSine(N, k, Xt = xt, f = f, sineRet = TRUE, deltat = deltat)
  v <- Y$Sine

  SineGram <- HatMatGMatSine(N = N, k = k, p = p)
  H <- SineGram$H
  G <- SineGram$G

  if(!returnHG){
    if(!polynomialPart){
      zCompSine <-  ((v %*% (diag(k) - H %*% t(H))) + G %*% t(H)) %*% Y$EigenCoef
      return(zCompSine)
    }
    else{
      cHat <- t(H) %*% Y$EigenCoef
      polyPart <- G %*% cHat

      return(list(chat = cHat, polyPart = polyPart))
    }
  }
  else{
    if(!polynomialPart){
      zCompSine <-  ((v %*% (diag(k) - H %*% t(H))) + G %*% t(H)) %*% Y$EigenCoef
      return(list(z = zCompSine, H = H, G = G))
    }
    else{
      cHat <- t(H) %*% Y$EigenCoef
      polyPart <- G %*% cHat

      return(list(cHat = cHat, polyPart = polyPart, H = H, G = G))
    }
  }

}



#' U_kp used in the regression for sine tapers
#'
#' @param N Total number of observations
#' @param k the number of tapers 0 to k-1
#' @param p Highest degree der used/ highest degree polynomial 0 to p
#' @param round default 15 decimals to obtain exactly 0 in the non mod 2 points,
#' I believe it has to do with the rounding on pi for why its not exactly zero.
#' @param passInSineTapers used for speeding up the calculation.
#'
#' @return matrix
USine <- function(N, k, p, round = 15, passInSineTapers = NULL){
  n <- 1:N
  if(is.null(passInSineTapers)){
    v <- sineTaperMatrix(N, k)
  }else{
    v <- passInSineTapers
  }

  uMat <- matrix(nrow = k, ncol = p + 1)
  for(i in 0:p){
    uMat[,i+1] <- as.matrix(apply(v *
                                    ((2/(N-1)) * (n-1) - 1)^i, 2, FUN = sum))
  }
  uMat <- round(uMat, round)

  return(uMat)
}

#' GramShmidt Orthonormalized Hat and Polynomial Matrices for Sine
#'
#' @param N Total Number of Observations
#' @param k Number of tapers
#' @param p Highest degree polynomial
#' @param passInSineTapers used for speeding up the calculation.
#'
#' @return List of H hat matrix and G polynomial matrix
HatMatGMatSine <- function(N, k, p, passInSineTapers = NULL){

  Umat <- USine(N, k, p, passInSineTapers = passInSineTapers)
  Rmat <- RnpMat(N = N,P = p)
  sineGram <- GramSchmidtMod(uMat = Umat, rMat = Rmat)

  return(list(H = sineGram$H, G = sineGram$G))
}

#' First Derivative of Sine Tapers
#'
#' @param N Total number of Observations
#' @param k Total number of tapers
#'
#' @return Matrix of evaluated derivatives for sine taper
FirstDerSineTaper <- function(N,k){
  stopifnot(N >= 8, k >= 1)

  taper <- cos(outer(1:N, 1:k)*pi/(N+1))
  taper <-  taper*sqrt(2/(N+1)) * (pi/(N+1))
  taper <- as.matrix(taper)*matrix(1:k, nrow = N, ncol = k, byrow = TRUE)

  return(taper)
}

#' Instantaneous Frequency for Sine
#'
#' @param xt time series
#' @param N Total number of observations
#' @param k Number of tapers
#' @param deltat default is 1, used for nyquist calculation
#' @param passInSineMat if sine matrix was calculated somewhere else
#' @param returnSineMat if you want to return SineTapers for use elsewhere
#' @param passInSineUnder If sine undersampling matrix was calculated already
#' @param penalty 1 is no penalty , 0.2  would give seq(from =  1, to =  1/(0.2*k, length.out = k)
#' penalty to each respective taper
#' @param penaltyOnTapersStdInv applys the penalty to all tapers in the calculation instead of just weighting the
#' eigenCoef's
#'
#' @return list of $InstFreq and $Freq which is the column of instFreq
#' if returnSineMat = TRUE then list of $InstFreq, $Freq, and $SineMat
instFreqSine <- function(xt, N, k, deltat, passInSineMat = NULL, returnSineMat = FALSE,
                         passInSineUnder = NULL, penalty = 1, penaltyType = "ScaledExp",
                         penaltyOnTapersStdInv = FALSE){
  if(is.null(passInSineUnder)){
    if(is.null(passInSineMat)){ # need to create sine taper matrix in stdInverse then pass into stdInvDer
      stdInv <- standardInverseSine(xt = xt, N = N, k = k, deltat = deltat,
                                    retSineTapers = TRUE, penalty = penalty,
                                    penaltyType = penaltyType,
                                    penaltyOnTapersStdInv = penaltyOnTapersStdInv)
      stdInvDer <- standardInverseSineDer(xt = xt, N = N, k = k,
                                          deltat = deltat, passInSineMat = stdInv$sineTapers
                                          , penalty = penalty, penaltyType = penaltyType)
      sineMat <- stdInv$sineTapers
    }
    else{ # we already have sine taper matrix so we will pass into stdinv, then into stdInvDer
      stdInv <- standardInverseSine(xt = xt, N = N, k = k, deltat = deltat,
                                    retSineTapers = FALSE, passInSineMat = passInSineMat
                                    , penalty = penalty, penaltyType = penaltyType,
                                    penaltyOnTapersStdInv = penaltyOnTapersStdInv)
      stdInvDer <- standardInverseSineDer(xt = xt, N = N, k = k,
                                          deltat = deltat, passInSineMat = passInSineMat
                                          , penalty = penalty, penaltyType = penaltyType)
      sineMat <- passInSineMat
    }
  }else{ #need to do some undersampling!
    if(is.null(passInSineMat)){ # need to create sine taper matrix in stdInverse then pass into stdInvDer
      stdInv <- standardInverseSine(xt = xt, N = N, k = k, deltat = deltat,
                                    retSineTapers = TRUE, passInSineMat = FALSE,
                                    passInSineUnder = passInSineUnder, penalty = penalty,
                                    penaltyType = penaltyType,
                                    penaltyOnTapersStdInv = penaltyOnTapersStdInv)
      stdInvDer <- standardInverseSineDer(xt = xt, N = N, k = k,
                                          deltat = deltat, passInSineMat = stdInv$sineTapers,
                                          passInSineUnder = passInSineUnder, penalty = penalty,
                                          penaltyType = penaltyType)
      sineMat <- stdInv$sineTapers
    }
    else{ # we already have sine taper matrix so we will pass into stdinv, then into stdInvDer
      stdInv <- standardInverseSine(xt = xt, N = N, k = k, deltat = deltat,
                                    retSineTapers = FALSE, passInSineMat = passInSineMat,
                                    passInSineUnder = passInSineUnder, penalty = penalty,
                                    penaltyType = penaltyType,
                                    penaltyOnTapersStdInv = penaltyOnTapersStdInv)
      stdInvDer <- standardInverseSineDer(xt = xt, N = N, k = k,
                                          deltat = deltat, passInSineMat = passInSineMat,
                                          passInSineUnder = passInSineUnder, penalty = penalty,
                                          penaltyType = penaltyType)
      sineMat <- passInSineMat
    }
  }

  A <- Re(stdInv$StdInverse)
  B <- Im(stdInv$StdInverse)
  ADot <- Re(stdInvDer$StdInverse)
  BDot <- Im(stdInvDer$StdInverse)
  rSqrd <- A^2 + B^2

  instFreq <- (A*BDot - B*ADot)/(2*pi*rSqrd)

  if(!returnSineMat){
    return(list(InstFreq = instFreq, Freq = stdInv$Freq))
  }
  else{
    return(list(InstFreq = instFreq, Freq = stdInv$Freq, SineMat = sineMat))
  }
}

#' Eigen coefficients for Instantainous Frequency Sine Tapers
#'
#' @param xt time series
#' @param N Total number of observations
#' @param k Number of tapers
#' @param deltat default is 1, used for nyquist calculation
#' @param returnSineMat if you want to return SineTapers for use elsewhere
#' @param passInSineTapers  if sine matrix was calculated somewhere else
#' @param passInSineUnder If sine undersampling matrix was calculated already
#' @param penalty 1 is no penalty , 0.2  would give seq(from =  1, to =  1/(0.2*k, length.out = k)
#' penalty to each respective taper
#' @param penaltyOnTapersStdInv applys the penalty to all tapers in the calculation instead of just weighting the
#' eigenCoef's
#'
#' @return returns $PSI kxn matrix, $Freq which are the columns of PSI.  If
#' returnSineMat = TRUE, will also return tapers for future use
eigenCoefSineInstFrequency <- function(xt, N, k, deltat,
                                           returnSineMat = FALSE,
                                           passInSineTapers = NULL,
                                           passInSineUnder = NULL,
                                           penalty = 1,
                                           penaltyOnTapersStdInv = FALSE,
                                           penaltyType = "ScaledExp"){
  if(is.null(passInSineUnder)){
    if(is.null(passInSineTapers)){ # need instFreq to create the sineTapers
      instFreq <- instFreqSine(xt = xt, N = N, k = k, deltat = deltat,
                               returnSineMat = TRUE, penalty = penalty,
                               penaltyType = penaltyType,
                               penaltyOnTapersStdInv = penaltyOnTapersStdInv)
      sineTaper <- instFreq$SineMat
    }
    else{ # will pass in sinetapers to instFreq
      instFreq <- instFreqSine(xt = xt, N = N, k = k, deltat = deltat,
                               returnSineMat = FALSE, passInSineMat = passInSineTapers,
                               penalty = penalty, penaltyType = penaltyType,
                               penaltyOnTapersStdInv = penaltyOnTapersStdInv)
      sineTaper <- passInSineTapers
    }
    v <- sineTaper



    PSI <- crossprod(v, instFreq$InstFreq)#t(v) %*% instFreq$InstFreq
  }
  else{ # we want under sampling so we need to chnge a few things
    if(is.null(passInSineTapers)){ # need instFreq to create the sineTapers
      instFreq <- instFreqSine(xt = xt, N = N, k = k, deltat = deltat,
                               returnSineMat = TRUE,
                               passInSineMat = passInSineTapers,
                               passInSineUnder = passInSineUnder,
                               penalty = penalty, penaltyType = penaltyType,
                               penaltyOnTapersStdInv = penaltyOnTapersStdInv)
      sineTaper <- instFreq$SineMat
    }
    else{ # will pass in sinetapers to instFreq
      instFreq <- instFreqSine(xt = xt, N = N, k = k, deltat = deltat,
                               returnSineMat = FALSE,
                               passInSineMat = passInSineTapers,
                               passInSineUnder = passInSineUnder,
                               penalty = penalty, penaltyType = penaltyType,
                               penaltyOnTapersStdInv = penaltyOnTapersStdInv)
      sineTaper <- passInSineTapers
    }
    v <- passInSineUnder



    PSI <- crossprod(v, instFreq$InstFreq)#t(v) %*% instFreq$InstFreq

  }



  if(!returnSineMat){
    return(list(PSI = PSI, Freq = instFreq$Freq))
  }
  else{ # will return the DPSS
    return(list(PSI = PSI, Freq = instFreq$Freq,SineTaper = sineTaper))
  }

}

#' regression for Sine Tapers using Gram schmidt hat matrix specifically for Inst Freq
#'
#' @param N Total number of observations
#' @param k Number of tapers
#' @param p Highest degree polynomial
#' @param passInSineMat if you have already computed the matrix of tapers
#' @param returnSineTapers if you need the tapers for other functions
#' @param instFreqEigen inst frequency eigen spectrum
#' @param returnRp if you want residuals to be returned and calculated as well
#' @param withoutzeroPoly if you dont want zeroth order polynomials to be in test statistics
#' @param returnGCHatp will return GCHatp to the user if needed
#'
#' @return returnlist of $cHat and polynomialPart of composit inverse $GCHatp, $H and $G, \
#'if returnSineTapers = TRUE then will return tapers as well $sineTapers
regressionSineInstFreq <- function(N, k, instFreqEigen, p,
                                   passInSineMat = NULL,
                                   returnSineTapers = FALSE,
                                   returnRp = FALSE,
                                   withoutzeroPoly = FALSE,
                                   returnGCHatp = FALSE){

  if(is.null(passInSineMat)){
    v <- sineTaperMatrix(N = N, k = k)
  }
  else{
    v <- passInSineMat
  }

  if(returnGCHatp){
    SineGram <- HatMatGMatSine(N = N, k = k, p = p, passInSineTapers = v)
    if(withoutzeroPoly){
      H <- SineGram$H[,-1] # removes the 0th order column
      G <- SineGram$G[,-1]
      cHat <- t(H) %*% instFreqEigen
      polyPart <- G %*% cHat
    }else{ # still has the zeroth order polynomial in it
      H <- SineGram$H
      G <- SineGram$G
      cHat <- t(H) %*% instFreqEigen
      polyPart <- G %*% cHat
    }

    if(returnRp){
      rp <- (diag(k) - H %*% t(H)) %*% instFreqEigen

      if(!returnSineTapers){
        return(list(cHat = cHat, GCHatp =  polyPart, rp = rp, H = H, G = G))
      }
      else{
        return(list(cHat = cHat, GCHatp =  polyPart, rp = rp, H = H, G = G,
                    sineTapers = v))
      }
    }
    else{
      if(!returnSineTapers){
        return(list(cHat = cHat, GCHatp =  polyPart, H = H, G = G))
      }
      else{
        return(list(cHat = cHat, GCHatp =  polyPart, H = H, G = G,
                    sineTapers = v))
      }
    }
  }else{
    SineGram <- HatMatGMatSine(N = N, k = k, p = p, passInSineTapers = v)
    if(withoutzeroPoly){
      H <- SineGram$H[,-1] # removes the 0th order column
      G <- SineGram$G[,-1]
      cHat <- t(H) %*% instFreqEigen

    }else{ # still has the zeroth order polynomial in it
      H <- SineGram$H
      G <- SineGram$G
      cHat <- t(H) %*% instFreqEigen

    }
    if(returnRp){
      rp <- (diag(k) - H %*% t(H)) %*% instFreqEigen

      if(!returnSineTapers){
        return(list(cHat = cHat, rp = rp, H = H, G = G))
      }
      else{
        return(list(cHat = cHat, rp = rp, H = H, G = G,
                    sineTapers = v))
      }
    }
    else{
      if(!returnSineTapers){
        return(list(cHat = cHat, H = H, G = G))
      }
      else{
        return(list(cHat = cHat, H = H, G = G,
                    sineTapers = v))
      }
    }
  }


}


# Polynomial Related --------------------------------------


#' Matrix whos pth column is the unscaled monomial of degree p on the interval -1 to 1
#'
#' @param n the nth point in the series 0 to N-1
#' @param p the degree of the polynomial 0 to P
#' @param N The total number of observations
#'
#' @return signle n,pth element from the matrix
Rnp <- function(n,p,N){

  return(((2/(N-1) * n) - 1)^p)
}




#' Polynomial matrix generation
#'
#' @param N total number of observations
#' @param P highest degree polynomial
#'
#' @return matrix of all polynomial numbers
RnpMat <- function(N,P){
  mat <- matrix(nrow = N, ncol = P+1)
  for(i in 0:(N-1)){
    for(j in 0:P){
      mat[i+1,j+1] <- Rnp(i,j,N)
    }
  }
  return(mat)
}



# Slepian Related -----------------------------------


#' Matrix of Ukp's defined by 2.36 in Kians thesis for the DPSS
#'
#' @param N Total number of points in your series including zero pad
#' @param k number of tapers used 0 to k-1
#' @param p highest derivative used/highest degree polynomial starts at 0
#' @param w bandwidth parameter, usually found with Shannons number
#' @param round if 16 you will obtain not zero in k != p mod 2
#'
#' @return matrix
UDpss <- function(N, k, p, w, round = 14, passInDPSS = NULL){
  n <- 1:N
  if(is.null(passInDPSS)){
    v <- multitaper::dpss(N, k, nw = N*w)$v
  }else{
    v <- passInDPSS
  }
  uMat <- matrix(nrow = k, ncol = p+1)
  for(i in 0:p){
    uMat[,i+1] <- as.matrix(apply(v *
                                    ((2/(N-1)) * (n-1) - 1)^i, 2, FUN = sum))
  }
  uMat <- round(uMat, round)

  return(uMat)
}


#' Non Centered eigenCoef Using DPSS
#'
#' @param n Total number of observations
#' @param k Number of tapers
#' @param w Width that gets passed into dpss
#' @param Xt Time series
#' @param f Specific frequency you are looking at
#' @param deltat = 1 unless otherwise inputted
#' @param dpssRet  = FALSE, will not return dpss use true to help not
#' re-computing dpss multiple times in program
#'
#' @return Vector of the eigenCoef for each k 0, ..., K-1
eigenCoefdpss <- function(n, k, w, Xt, f, deltat = 1, dpssRet = FALSE){

  EigenCoef <- vector(length = k)
  dpFull <- multitaper::dpss(n, k, n*w)
  dp <- dpFull$v
  for(i in 1:k){
    vec <- t(dp[,i]) %*% as.matrix(Xt * exp(-1i*2*pi*f*0:(n-1)*deltat))
    EigenCoef[i] <- sum(vec)
  }
  if(dpssRet){
    return(list(EigenCoef = EigenCoef,
                Dpss = dpFull))
  }
  else{
    return(EigenCoef)
  }
}

#' Non Centered eigenCoef Using DPSS and fft for speed
#'
#' @param n Total number of observations
#' @param k Number of tapers
#' @param w Width that gets passed into dpss
#' @param Xt Time series
#' @param passInDPSSMat For limiting the number of calculations of DPSS this is dpss$v
#' @param deltat = 1 by default
#' @param pad if true will zero padd automatically the Xt that is inputted into the function
#' @param penalty 1 is no penalty , 0.2  would give seq(from =  1, to =  1/(0.2*k, length.out = k)
#' penalty to each respective taper
#'
#' @return Matrix of EigenCoef for each frequency, as well as the frequencies under Freq in the list
eigenCoefdpssFFT <- function(n, k, w, Xt, deltat = 1, passInDPSSMat = NULL,
                             pad = FALSE, penalty = 1, penaltyType = "ScaledExp"){
  EigenCoef <- matrix(nrow = n,ncol = k)
  if(is.null(passInDPSSMat)){
    dp <- multitaper::dpss(n, k, n*w)$v
  }
  else{
    dp <- passInDPSSMat
  }
  if(!pad){

    for(i in 1:k){
      vec <- dp[,i] * as.matrix(Xt)
      EigenCoef[,i] <- fft(vec)
    }

    if(penalty != 1){
      EigenCoef <- t(apply(EigenCoef, MARGIN = 1, FUN = function(x){x/seq(from = 1, to = penalty*k, length.out = k)}))
    }
    ret <- data.frame(EigenCoef = EigenCoef, Freq = c(seq(from = 0, to = 1/(2*deltat),
                                                          by = 1/(n*deltat)),
                                                      seq(from = (-1/(2*deltat) + 1/(n*deltat)),
                                                          to = (-1/(n*deltat)),
                                                          by = 1/(n*deltat))))
  }else{
    nFFT <- 2^ceiling(log2(2*n))
    nextPowerOfTwo <- 2^ceiling(log2(2*n)-1)
    taper <- dp * Xt
    pad <- rbind(taper, matrix(0, nrow = nFFT - n, ncol = k))
    EigenCoef <- mvfft(pad)[1:(nextPowerOfTwo + 1), ,drop = FALSE]


    if(penaltyType == "mtm"){
      spec <- multitaper::spec.mtm(ts(Xt, deltat = deltat), nw = n*w,k = k, adaptiveWeighting = TRUE,
                                   returnInternals = TRUE, plot = FALSE)
      EigenCoef <- EigenCoef * spec$mtm$eigenCoefWt
    }
    else if(penaltyType == "ScaledExp"){
      if(penalty == 1){
        # do nothing as that is no penalty
      }else{
        EigenCoef <- t(apply(EigenCoef, MARGIN = 1, FUN = function(x){x/seq(from = 1, to = penalty*k, length.out = k)}))
      }
    }
    else if(penaltyType == "Clip"){
      clipPoint <- penalty#2/3
      weight <- c(rep(1, length.out = floor(k*clipPoint)), rep(0, length.out = k - floor(k*clipPoint)))
      weightMat <- matrix(weight, nrow = nrow(EigenCoef), ncol = k, byrow = TRUE)
      EigenCoef <- EigenCoef * weightMat
    }
    ret <- data.frame(EigenCoef = EigenCoef, Freq = seq(from = 0, to = 1/(2*deltat),
                                                        by = 1/(2*nextPowerOfTwo*deltat)))
  }
  ret <- ret[order(ret$Freq, decreasing = FALSE),]
  eSpecReorder <- as.matrix(ret[,1:(ncol(ret)-1)])
  rownames(eSpecReorder) <- ret$Freq
  colnames(eSpecReorder) <- 1:ncol(eSpecReorder)
  return(list(EigenCoef = eSpecReorder, Freq = ret$Freq))
}



#' Standard Inversion for Polynomial Projection Filter DPSS at freq f
#'
#' @param xt time series
#' @param N Total number of observations
#' @param k Number of tapers
#' @param w Width that gets passed into dpss
#' @param f specified frequency you want to look at
#' @param FFT if you want to use fft to speed up or not have delta t = 1 with zero padding
#' @param deltat  = 1 if not specified, but is only used if not using FFT
#' @param retDPSS if TRUE, returns calculated DPSS for pass in to other functions
#'
#' @return Standard Inverse vector for frequency f
standardInverseDPSSSingleFreq <- function(xt, N, w, k, f, deltat = 1, FFT = FALSE, retDPSS = FALSE){

  dp <- multitaper::dpss(N, k, N*w)
  v <- dp$v
  if(!FFT){
    Xf <- xt * exp(-1i*2*pi*f*0:(N-1)*deltat)
    stdInverse <- v %*% t(v) %*% as.matrix(Xf)
  }
  else{
    Y <- eigenCoefdpssFFT(N, k, w, xt, deltat = deltat, passInDPSSMat = v)
    stdInverse <- v %*% as.matrix(Y$EigenCoef[which(Y$Freq == f),])
  }

  if(!retDPSS){
    return(stdInverse)
  }
  else{
    return(list(StdInverse = stdInverse, DPSSv = dp$v, DPSSEigen = dp$eigen))
  }
}

#' Standard Inverse with fft for all frequencies using DPSS
#'
#' @param xt time series
#' @param N Total number of observations
#' @param k Number of tapers
#' @param w Width that gets passed into dpss
#' @param retDPSS If true, will return the dpss for pass into other functions
#' @param deltat default is 1, used for nyquist calculation
#' @param passInDPSS if dpss has already been calculated somewhere else this is the full
#' object of dpss both $v and $eigen
#' @param passInDPSSReduced if dpss undersampling was calculated already
#' @param penalty 1 is no penalty , 0.2  would give seq(from =  1, to =  1/(0.2*k, length.out = k)
#' penalty to each respective taper
#'
#' @return returns standard inverse along with frequencies that are the columns
standardInverseDPSS <- function(xt, N, k, w, deltat = 1, retDPSS = FALSE,
                                passInDPSS = NULL, passInDPSSReduced = NULL,
                                penalty = 1, penaltyType = "ScaledExp"){

  if(is.null(passInDPSS)){
    dp <- multitaper::dpss(N, k, N*w)
    v <- dp$v
  }
  else{
    dp <- passInDPSS
    v <- dp$v
  }

  if(is.null(passInDPSSReduced)){ # we don't zeropadd
    Y <- eigenCoefdpssFFT(N, k, w, xt, deltat = deltat, passInDPSSMat = v,
                          pad = FALSE, penalty = penalty, penaltyType = penaltyType)
  }else{
    Y <- eigenCoefdpssFFT(N, k, w, xt, deltat = deltat, passInDPSSMat = v, pad = TRUE,
                          penalty = penalty, penaltyType = penaltyType)
  }
  #stdInverse <- matrix(nrow = N, ncol = length(Y$Freq))
  #index <- 0
  # for(i in 1:length(Y$Freq)){
  #   #index <- index + 1
  #   stdInverse[,i] <- v %*% as.matrix(Y$EigenCoef[i,])
  # }

  if(is.null(passInDPSSReduced)){
    stdInverse <- tcrossprod(v,Y$EigenCoef)
  }else{
    stdInverse <- tcrossprod(passInDPSSReduced$v,Y$EigenCoef)
  }
  colnames(stdInverse) <- Y$Freq
  if(!retDPSS){
    return(list(StdInverse = stdInverse, Freq = Y$Freq))
  }
  else{
    return(list(StdInverse = stdInverse,Freq = Y$Freq, DPSS = dp))
  }
}

#' Standard Inversion for Polynomial Projection Filter First Derivative DPSS
#'
#' @param xt time series
#' @param N Total number of observations
#' @param k Number of tapers
#' @param w Width that gets passed into dpss
#' @param deltat default is 1 used for nyquist
#' @param passInDPSS used for speed if you have already calculated the dpss object somewhere else
#' this includes both $v and $eigen
#' @param returnDPSS returns dpss if needed outside function
#' @param passInDPSSReduced if dpss undersampling was calculated already
#' @param penalty 1 is no penalty , 0.2  would give seq(from =  1, to =  1/(0.2*k, length.out = k)
#' penalty to each respective taper
#'
#' @return Standard Inverse Derivative Z' vector $StdInverse as well as $Freq which are the columns of stdInverse
#' if returnDPSS = TRUE will also return $DPSS which is the dpss object
standardInverseDPSSFirstDir <- function(xt, N, w, k, deltat = 1, passInDPSS = NULL,
                                        returnDPSS = FALSE, passInDPSSReduced = NULL,
                                        penalty = 1, penaltyType = "ScaledExp"){#, FirstDir = NULL){
  if(is.null(passInDPSSReduced)){
    if(is.null(passInDPSS)){
      FirstDir <- FirstDerTimeDomSlepians(N = N, w = w, k = k, returnDPSS = TRUE)
      dpfull <- FirstDir$DPSS
      dp <- dpfull$v
    }
    else{

      FirstDir <- FirstDerTimeDomSlepians(N = N, w = w, k = k,
                                          passInFullDPSS = passInDPSS,
                                          returnDPSS = FALSE)
      dpfull <- passInDPSS
      dp <- dpfull$v
    }
    # if(is.null(FirstDir)){
    #FirstDir <- FirstDerTimeDomSlepians(N = N, w = w, k = k, returnDPSS = TRUE)
    #   dp <- FirstDir$DPSSV
    # }
    # else{
    #   dp <- dpss(N, k, N*W)$v
    # }


    Y <- eigenCoefdpssFFT(N, k, w, xt, deltat = deltat, passInDPSSMat = dp,
                          penalty = penalty, penaltyType = penaltyType)

    stdInverse <- matrix(nrow = N, ncol = length(Y$Freq))
    #index <- 0

    stdInverse <- tcrossprod(FirstDir$Vdot,Y$EigenCoef)
  }else{ #undersamp so
    if(is.null(passInDPSS)){
      stop("you must pass in DPSS and reduced DPSS to use this function")
    }
    else{

      FirstDir <- FirstDerTimeDomSlepians(N = nrow(passInDPSSReduced$v),
                                          w = nrow(passInDPSS$v)*w/nrow(passInDPSSReduced$v),  #scaled w for under sampling
                                          k = k,
                                          passInFullDPSS = passInDPSSReduced,
                                          returnDPSS = FALSE)
      dpfull <- passInDPSS
      dp <- dpfull$v
    }


    Y <- eigenCoefdpssFFT(nrow(dp), k, w, xt, deltat = deltat, passInDPSSMat = dp, pad = TRUE, penalty = penalty)

    stdInverse <- matrix(nrow = nrow(dp), ncol = length(Y$Freq))
    #index <- 0

    stdInverse <- tcrossprod(FirstDir$Vdot,Y$EigenCoef)
  }

  # for(i in 1:length(Y$Freq)){
  #   #index <- index + 1
  #   stdInverse[,i] <- FirstDir$Vdot %*% as.matrix(Y$EigenCoef[i,])
  # }
  colnames(stdInverse) <- Y$Freq
  if(!returnDPSS){
    return(list(StdInverse = stdInverse, Freq = Y$Freq))
  }
  else{
    return(list(StdInverse = stdInverse, Freq = Y$Freq,
                DPSS = dpfull))
  }
}


#' Composite inverse for Dpss using from schmidt hat matrix
#'
#' @param N Total number of observations
#' @param k Number of tapers
#' @param w Width that gets passed into dpss
#' @param xt Time series
#' @param f Specified frequency
#' @param p Highest degree polynomial
#' @param deltat  = 1 if not specified
#' @param polynomialOnly If you only want the polynomial part GH^tY == G_pC^hat_p used in ftest
#' @param returnHG Will return H and G along with the Z/includes in list if poly is true
#'
#' @return Vector of compInverse for freq f, if polynomialOnly, it will return
#' list of cHat and polynomialPart of composit inverse which is G_p*cHat
compInverseDpss <- function(N, k, w, xt, f, p, deltat = 1, polynomialOnly = FALSE, returnHG = FALSE){

  Y <- eigenCoefdpss(N, k, w, xt, f = f, dpssRet = TRUE, deltat = deltat)
  v <- Y$Dpss$v

  dpssGram <- HatMatGMatDpss(N = N, k = k, p = p, w = w)
  H <- dpssGram$H
  G <- dpssGram$G

  if(!returnHG){
    if(!polynomialOnly){ # just z
      zCompDp <-  ((v %*% (diag(k) - H %*% t(H))) + G %*% t(H)) %*% Y$EigenCoef
      return(zCompDp)
    }
    else{ # poly part only
      cHat <- t(H) %*% Y$EigenCoef
      polyPart <- G %*% cHat

      return(list(cHat = cHat, polyPart = polyPart))
    }
  }
  else{ # returning G and H as well
    if(!polynomialOnly){
      zCompDp <-  ((v %*% (diag(k) - H %*% t(H))) + G %*% t(H)) %*% Y$EigenCoef
      return(list(z = zCompDp, H = H, G = G))
    }
    else{ # poly part and G and H
      cHat <- t(H) %*% Y$EigenCoef
      polyPart <- G %*% cHat

      return(list(cHat = cHat, polyPart = polyPart, H = H, G = G))
    }
  }


}


#' Derivative function that creates matrix of pth derivative columns and k th taper rows #'for the pth derivative with respect to f.
#'
#' May want to use the optimization of this in the gegenbauerPolynomialDerivation.Rmd
#'
#' @param N Total number of points in your series including zero pad
#' @param k number of tapers used
#' @param p highest derivative used
#' @param w bandwidth parameter, usually found with Shannons number
#' @param f Frequency you are looking at
#'
#' @return A matrix

eigenFunctionDerDpss <- function(N, k, p, w, f){
  Der <- matrix(nrow = k, ncol = p)
  Dpss <- multitaper::dpss(N, k = k, nw = N*w)$v
  for(i in 1:k){
    for(j in 1:p){
      sumVec <- vector(length = N)
      for(n in 1:N){
        sumVec[n] <- Dpss[n,i] * ((n - 1)-(N-1)/2)^j * exp(1i*2*pi*f*((n-1) - (N-1)/2))
      }
      if(i %% 2 == 0){

        Der[i,j] <- (1i*2*pi)^j * sum(sumVec)
      }
      else{
        Der[i,j] <- 1i * (1i*2*pi)^j * sum(sumVec)
      }

    }
  }
  return(Der)
}

#' GramShmidt Orthonormalized Hat and Polynomial Matrices for Dpss
#'
#' @param N Total Number of Observations
#' @param k Number of tapers
#' @param p Highest degree polynomial
#' @param w Band Width parameter passed into dpss
#'
#' @return List of H hat matrix and G polynomial matrix
HatMatGMatDpss <- function(N, k, p, w, passInDPSS = NULL){

  Umat <- UDpss(N = N, k = k, p = p, w = w, passInDPSS = passInDPSS)
  Rmat <- RnpMat(N = N,P = p)
  dpssGram <- GramSchmidtMod(uMat = Umat, rMat = Rmat)

  return(list(H = dpssGram$H, G = dpssGram$G))
}

#'  Time domain slepians **Untested**
#'
#'  Calculates the time domain version of the slepians, it is not used in any of
#'  the functions but could be usefull to know in the future.
#'
#' @param N Total number of observations
#' @param w bandwidth parameter
#' @param k number of tapers
#' @param returnDPSS returnDPSS for limiting calculations of dpss
#' @param passInFullDPSS this is the object returned by dpss.  both $v and $eigen
#'
#' @return matrix of der for each t = 0, ..., N-1 by k = 0, ..., k-1

TimeDomSlepians <- function(N, w, k, passInFullDPSS = NULL, returnDPSS = FALSE){

  if(is.null(passInFullDPSS)){
    dpMat <- matrix(nrow = N, ncol = k)
    dp <- multitaper::dpss(N, k, N*w)
    dpMat <- dp$v
  }
  else{
    dp <- passInFullDPSS
    dpMat <- dp$v
  }

  diff <- 1:(N - 1)
  #
  # #Generating rho matrix
  firstCol <- c(2*w, sin(2*pi*w*diff)/(pi*diff))
  #Generating rho matrix
  # rho <- matrix(nrow = N, ncol = N)
  # for(i in 1:N){
  #   for(j in 1:N){
  #     diff <- i-j
  #     if(diff == 0){
  #       rho[i,j] <- 2*w
  #     }
  #     else{
  #       rho[i,j] <-  sin(2*pi*w*diff)/(pi*diff)
  #     }
  #   }
  # }

  # due to the symmetry of the matrix generated with the double for loop above,
  #we can do it with only a vector and a single for loop making it much faster!
  rho <- matrix(nrow = N, ncol = N)
  rho[,1] <- firstCol
  for(i in 2:N){
    rho[seq(from = i, to = N, by = 1), i] <- firstCol[1:(N-(i-1))]
    rho[seq(from = 1, to = (i-1), by = 1), i] <- firstCol[i:2]#((N-(i-1))+1):N]
  }

  #(2.32) in kians thesis
  V <- t(apply((rho %*% dpMat), MARGIN = 1, FUN = function(x) {x/dp$eigen}))#matrix(dp$eigen, nrow = N, ncol = k, byrow = TRUE)
  if(!returnDPSS){
    return(list(V = V))
  }
  else{
    return(list(V = V, DPSS = dp))
  }
}


#' First derivative of time domain slepians
#'
#' @param N Total number of observations
#' @param w bandwidth parameter
#' @param k number of tapers
#' @param returnDPSS returnDPSS for limiting calculations of dpss
#' @param passInFullDPSS this is the object returned by dpss.  both $v and $eigen
#'
#' @return matrix of der for each t = 0, ..., N-1 by k = 0, ..., k-1

FirstDerTimeDomSlepians <- function(N, w, k, passInFullDPSS = NULL, returnDPSS = FALSE){

  if(is.null(passInFullDPSS)){
    dpMat <- matrix(nrow = N, ncol = k)
    dp <- multitaper::dpss(N, k, N*w)
    dpMat <- dp$v
  }
  else{
    dp <- passInFullDPSS
    dpMat <- dp$v
  }

  diff <- 1:(N - 1)
  #
  # #Generating rho matrix
  a <- c(0, 2*w*cos(2*pi*w*diff)/diff)
  b <- c(0, sin(2*pi*w*diff)/(pi*diff^2))
  firstCol <- a - b
  #Generating rho matrix
  # rho <- matrix(nrow = N, ncol = N)
  # for(i in 1:N){
  #   for(j in 1:N){
  #     diff <- i-j
  #     if(diff == 0){
  #       rho[i,j] <- 0
  #     }
  #     else{
  #       rho[i,j] <- 2*w * cos(2*pi*w*diff)/diff - sin(2*pi*w*diff)/(pi*diff^2)
  #     }
  #   }
  # }

  # due to the symmetry of the matrix generated with the double for loop above,
  #we can do it with only a vector and a single for loop making it much faster!
  rho <- matrix(nrow = N, ncol = N)
  rho[,1] <- firstCol
  for(i in 2:N){
    rho[seq(from = i, to = N, by = 1), i] <- firstCol[1:(N-(i-1))]
    rho[seq(from = 1, to = (i-1), by = 1), i] <- -1*firstCol[i:2]#((N-(i-1))+1):N]
  }

  #(2.32) in kians thesis
  Vdot <- t(apply((rho %*% dpMat), MARGIN = 1, FUN = function(x) {x/dp$eigen}))#matrix(dp$eigen, nrow = N, ncol = k, byrow = TRUE)
  if(!returnDPSS){
    return(list(Vdot = Vdot))
  }
  else{
    return(list(Vdot = Vdot, DPSS = dp))
  }
}

#' Instantaneous Frequency for DPSS
#'
#' @param xt time series
#' @param N Total number of observations
#' @param k Number of tapers
#' @param w Width that gets passed into dpss
#' @param deltat default is 1, used for nyquist calculation
#' @param passInDPSS if dpss has already been calculated somewhere else this is the full
#' object of dpss both $v and $eigen
#' @param passInDPSSUnder if dpss undersampling was calculated already
#' @param returnDPSS if you want to return dpss for use elsewhere
#' @param penalty 1 is no penalty , 0.2  would give seq(from =  1, to =  1/(0.2*k, length.out = k)
#' penalty to each respective taper
#'
#' @return list of $InstFreq and $Freq which is the column of instFreq
#' if returnDPSS = TRUE then list of $InstFreq, $Freq, and $DPSS
instFreqDPSS <- function(xt, N, k, w, deltat, passInDPSS = NULL, returnDPSS = FALSE,
                         passInDPSSUnder = NULL, penalty = 1, penaltyType = "ScaledExp"){
  if(is.null(passInDPSSUnder)){ # no undersampling will take place
    if(is.null(passInDPSS)){ # need to create dpss in stdInverse then pass into stdInvDer


      stdInv <- standardInverseDPSS(xt = xt, N = N, k = k, w = w, deltat = deltat,
                                    retDPSS = TRUE, penalty = penalty, penaltyType = penaltyType)
      stdInvDer <- standardInverseDPSSFirstDir(xt = xt, N = N, w = w, k = k,
                                               deltat = deltat, passInDPSS = stdInv$DPSS,
                                               penalty = penalty, penaltyType = penaltyType)

      dp <- stdInv$DPSS
    }
    else{ # we already have dpss so we will pass into std, then into stdInvDer
      stdInv <- standardInverseDPSS(xt = xt, N = N, k = k, w = w, deltat = deltat,
                                    retDPSS = FALSE, passInDPSS = passInDPSS,
                                    penalty = penalty, penaltyType = penaltyType)
      stdInvDer <- standardInverseDPSSFirstDir(xt = xt, N = N, w = w, k = k,
                                               deltat = deltat, passInDPSS = passInDPSS,
                                               penalty = penalty, penaltyType = penaltyType)
      dp <- passInDPSS
    }
  }else{ #need to do some undersampling!
    if(is.null(passInDPSS)){ # need to create dpss in stdInverse then pass into stdInvDer


      stdInv <- standardInverseDPSS(xt = xt, N = N, k = k, w = w, deltat = deltat,
                                    retDPSS = TRUE, passInDPSSReduced = passInDPSSUnder,
                                    penalty = penalty, penaltyType = penaltyType)
      stdInvDer <- standardInverseDPSSFirstDir(xt = xt, N = N, w = w, k = k,
                                               deltat = deltat, passInDPSS = stdInv$DPSS,
                                               passInDPSSReduced = passInDPSSUnder,
                                               penalty = penalty, penaltyType = penaltyType)

      dp <- stdInv$DPSS
    }
    else{ # we already have dpss so we will pass into std, then into stdInvDer
      stdInv <- standardInverseDPSS(xt = xt, N = N, k = k, w = w, deltat = deltat,
                                    retDPSS = TRUE, passInDPSS = passInDPSS,
                                    passInDPSSReduced = passInDPSSUnder,
                                    penalty = penalty, penaltyType = penaltyType)
      stdInvDer <- standardInverseDPSSFirstDir(xt = xt, N = N, w = w, k = k,
                                               deltat = deltat, passInDPSS = stdInv$DPSS,
                                               passInDPSSReduced = passInDPSSUnder,
                                               penalty = penalty, penaltyType = penaltyType)
      dp <- passInDPSS
    }
  }


  A <- Re(stdInv$StdInverse)
  B <- Im(stdInv$StdInverse)
  ADot <- Re(stdInvDer$StdInverse)
  BDot <- Im(stdInvDer$StdInverse)
  rSqrd <- A^2 + B^2

  instFreq <- (A*BDot - B*ADot)/(2*pi*rSqrd)

  if(!returnDPSS){
    return(list(InstFreq = instFreq, Freq = stdInv$Freq))
  }
  else{
    return(list(InstFreq = instFreq, Freq = stdInv$Freq, DPSS = dp))
  }
}

#' Eigen Coefficients for DPSS Instantanious Frequency
#'
#' @param xt time series
#' @param N Total number of observations
#' @param k Number of tapers
#' @param w Width that gets passed into dpss
#' @param deltat default is 1, used for nyquist calculation
#' @param nu note nu = 0 is only working right now, implimentation is there for nu != 0 but not running parallel
#' @param returnDPSS if you want to return dpssfor use elsewhere
#' @param passInDPSS  if dpss has already been calculated somewhere else this is the full
#' object of dpss both $v and $eigen
#' @param passInDPSSUnder if dpss undersampling was calculated already
#' @param penalty 1 is no penalty , 0.2  would give seq(from =  1, to =  1/(0.2*k, length.out = k)
#' penalty to each respective taper
#'
#' @return returns $PSI kxn matrix, $Freq which are the columns of PSI.  If
#' returnDPSS = TRUE, will also return tapers $v and $eigen for future use
eigenCoefDPSSInstFrequency <- function(xt, N, k, w, deltat, nu = 0, returnDPSS = FALSE,
                                           passInDPSS = NULL, passInDPSSUnder = NULL,
                                       penalty = 1, penaltyType = "ScaledExp"){
  if(is.null(passInDPSSUnder)){#runs with no under sampling
    if(is.null(passInDPSS)){ # need instFreq to create the dpss
      instFreq <- instFreqDPSS(xt = xt, N = N, k = k, w = w, deltat = deltat,
                               returnDPSS = TRUE, penalty = penalty, penaltyType = penaltyType)
      dp <- instFreq$DPSS
    }
    else{ # will pass DPSS inot instFreq
      instFreq <- instFreqDPSS(xt = xt, N = N, k = k, w = w, deltat = deltat,
                               returnDPSS = FALSE, passInDPSS = passInDPSS,
                               penalty = penalty, penaltyType = penaltyType)
      dp <- passInDPSS
    }
    #need to use the undersampling dpss here
    v <- dp$v

    #if(nu == 0){
    PSI <- crossprod(v, instFreq$InstFreq)#t(v) %*% instFreq$InstFreq
    #}
  }else{ #we want under sampling so we need to change a few things
    if(is.null(passInDPSS)){ # need instFreq to create the dpss
      instFreq <- instFreqDPSS(xt = xt, N = N, k = k, w = w, deltat = deltat,
                               returnDPSS = TRUE, passInDPSSUnder = passInDPSSUnder,
                               penalty = penalty, penaltyType = penaltyType)
      dp <- instFreq$DPSS
    }
    else{ # will pass DPSS inot instFreq
      instFreq <- instFreqDPSS(xt = xt, N = N, k = k, w = w, deltat = deltat,
                               returnDPSS = FALSE, passInDPSS = passInDPSS,
                               passInDPSSUnder = passInDPSSUnder, penalty = penalty,
                               penaltyType = penaltyType)
      dp <- passInDPSS
    }

    v <- passInDPSSUnder$v

    #if(nu == 0){
    PSI <- crossprod(v, instFreq$InstFreq)#t(v) %*% instFreq$InstFreq
    #}
  }


  if(!returnDPSS){
    return(list(PSI = PSI, Freq = instFreq$Freq))
  }
  else{ # will return the DPSS
    return(list(PSI = PSI, Freq = instFreq$Freq,DPSS = dp))
  }

}

#' regression for DPSS Tapers using Gram schmidt hat matrix specifically for Inst Freq
#'
#' @param N Total number of observations
#' @param k Number of tapers
#' @param p Highest degree polynomial
#' @param passInDPSS if you have already computed dpss() with $v and $eigen
#' @param returnDPSS if you need the tapers for other functions
#' @param instFreqEigen inst frequency eigen spectrum
#' @param withoutzeroPoly used if using modified f tests without 0th order polynomial
#' @param w bandwidth for the multitaper
#' @param returnRp will return the residuals if needed
#' @param returnGCHatp will return GCHatP if needed
#'
#' @return returnlist of $cHat and polynomialPart of composit inverse $GCHatp, $H and $G, \
#'if returnDPSS = TRUE then will return tapers as well $DPSS
regressionDPSSInstFreq <- function(N, k, w, instFreqEigen, p,
                                   passInDPSS = NULL,
                                   returnDPSS = FALSE,
                                   returnRp = FALSE,
                                   withoutzeroPoly = FALSE,
                                   returnGCHatp = FALSE){

  if(is.null(passInDPSS)){
    dp <- multitaper::dpss(n = N, k = k, nw = N*w)
    v <- dp$v
  }
  else{
    dp <- passInDPSS
    v <- dp$v
  }

  if(returnGCHatp){
    dpGram <- HatMatGMatDpss(N = N, k = k, p = p, w = w, passInDPSS = v)
    if(withoutzeroPoly){
      H <- dpGram$H[,-1] # removes the 0th order column
      G <- dpGram$G[,-1]
      cHat <- crossprod(H, instFreqEigen)#t(H) %*% instFreqEigen
      polyPart <- G %*% cHat
    }else{ # still has the zeroth order polynomial in it
      H <- dpGram$H
      G <- dpGram$G
      cHat <- crossprod(H, instFreqEigen)#t(H) %*% instFreqEigen
      polyPart <- G %*% cHat
    }


    if(returnRp){
      rp <- (diag(k) - H %*% t(H))%*% instFreqEigen

      if(!returnDPSS){
        return(list(cHat = cHat, GCHatp =  polyPart, rp = rp, H = H, G = G))
      }
      else{
        return(list(cHat = cHat, GCHatp =  polyPart, rp = rp, H = H, G = G,
                    DPSS = dp))
      }
    }
    else{
      if(!returnDPSS){
        return(list(cHat = cHat, GCHatp =  polyPart, H = H, G = G))
      }
      else{
        return(list(cHat = cHat, GCHatp =  polyPart, H = H, G = G,
                    DPSS = dp))
      }
    }
  }else{
    dpGram <- HatMatGMatDpss(N = N, k = k, p = p, w = w, passInDPSS = v)
    if(withoutzeroPoly){
      H <- dpGram$H[,-1] # removes the 0th order column
      G <- dpGram$G[,-1]
      cHat <- crossprod(H, instFreqEigen)#t(H) %*% instFreqEigen
      #polyPart <- G %*% cHat
    }else{ # still has the zeroth order polynomial in it
      H <- dpGram$H
      G <- dpGram$G
      cHat <- crossprod(H, instFreqEigen)#t(H) %*% instFreqEigen
      #polyPart <- G %*% cHat
    }


    if(returnRp){
      rp <-  (diag(k) - H %*% t(H))%*% instFreqEigen

      if(!returnDPSS){
        return(list(cHat = cHat, rp = rp, H = H, G = G))
      }
      else{
        return(list(cHat = cHat, rp = rp, H = H, G = G,
                    DPSS = dp))
      }
    }
    else{
      if(!returnDPSS){
        return(list(cHat = cHat, H = H, G = G))
      }
      else{
        return(list(cHat = cHat, H = H, G = G,
                    DPSS = dp))
      }
    }
  }
}

# Gram Schmidt ---------------------------------


#' modified gram schmidt from 2.47 and .48 and Kian Thesis
#'
#' @param uMat U matrix used for the regression
#' @param rMat Matrix of polynomials used in regression
#'
#' @return orthonormalized U and orthogonal r
GramSchmidtMod <- function(uMat, rMat){
  for(p in 1:ncol(uMat)){
    if(p == 1 | p == 2){
      innerPrdt <- as.numeric(1/sqrt(t(uMat[,p]) %*% uMat[,p]))
      uMat[,p] <- innerPrdt * uMat[,p]
      rMat[,p] <- innerPrdt * rMat[,p]
    }
    else{
      kMin <- 2 #p = 4,6,8,...
      if(p %% 2 == 1){
        kMin <- 1 #p = 3,5,7,9,...
      }
      # either start on the first row or the second row depending on the pth column
      for(l in seq(from = kMin, to = (p-1), by = 2 )){# skipping every other one
        innerPrdt <- as.numeric(t(uMat[,p]) %*% uMat[,l])
        uMat[,p] <- uMat[,p] - innerPrdt * uMat[,l]
        rMat[,p] <- rMat[,p] - innerPrdt * rMat[,l]
      }
      normConst <- as.numeric(1/(sqrt(t(uMat[,p])%*%uMat[,p])))
      uMat[,p] <- normConst * uMat[,p] #normalize
      rMat[,p] <- normConst * rMat[,p]
    }
  }
  return(list(H = uMat, G = rMat))
}

# Inner Ftest Functions -------------------------


#' Single iteration of F4Mod for Parallel applications
#'
#' can be used inside mclapply.  note that there is no checks in here as they are assumed
#' to be ran in the parent user function
#'
#' @param xt vector of time series observations
#' @param k single k that will be used for the f test
#' @param p largest degree of the polynomial modulation suspected
#' @param deltat  = 1 by default
#' @param w multitaper bandwidth parameter
#' @param dpss = FALSE uses sine instead
#' @param undersampleNumber = 100 by default  this is the number used in undersampling the tapers
#' @param N length of Xt
#'
#' @return $cHat chat used in caluclation of f , $PSI capital PSI , $Freq vector of frequencies , $k number of tapers used
singleIterationForParallel4 <- function(xt, N, k, p, deltat = 1, w = NULL, dpss = FALSE, undersampleNumber = 100){


  if(is.null(undersampleNumber)){
    stop("need to set undersample amount")
  }
  if(dpss){ #Use DPSS taper
    if(is.null(w)){
      stop("need to set w for dpss")
    }
    dp <- multitaper::dpss(n = N, k = k, nw = N*w)
    dpUnder <- multitaper::dpss(n = undersampleNumber, k = k, nw = N*w)
    instFreqEigen <- eigenCoefDPSSInstFrequency(xt = xt, N = N, k = k, w = w,
                                                    deltat = deltat,
                                                    returnDPSS = FALSE, passInDPSS = dp,
                                                    passInDPSSUnder = dpUnder)
    fStuff <- regressionDPSSInstFreq(N = N, k = k, w = w, instFreqEigen = instFreqEigen$PSI,
                                     p = p, passInDPSS = dpUnder ,returnDPSS = FALSE,
                                     returnRp = FALSE, withoutzeroPoly = TRUE)
  }
  else{ #Sine Tapers are used
    sine <- sineTaperMatrix(N = N, k = k)
    sineUnder <- sineTaperMatrix(N = undersampleNumber, k = k)
    instFreqEigen <- eigenCoefSineInstFrequency(xt = xt, N = N, k = k,deltat = deltat,
                                                    returnSineMat = FALSE, passInSineTapers = sine,
                                                    passInSineUnder = sineUnder)

    fStuff <- regressionSineInstFreq(N = N, k = k, instFreqEigen = instFreqEigen$PSI,
                                     p = p, returnSineTapers = FALSE,
                                     passInSineMat = sineUnder,
                                     returnRp = FALSE, withoutzeroPoly = TRUE)
  }


  #removingzero and nyquist frequencies
  zeroNyquist <- c(length(instFreqEigen$Freq),which(instFreqEigen$Freq == 0))
  instFreqEigen$Freq <- instFreqEigen$Freq[-zeroNyquist]
  Freq <- instFreqEigen$Freq
  instFreqEigen$PSI <- instFreqEigen$PSI[,-zeroNyquist]
  fStuff$cHat <- fStuff$cHat[,-zeroNyquist]


  return(list(cHat = fStuff$cHat, PSI = instFreqEigen$PSI, Freq = Freq, k = k))
}




#' Single iteration of F3Mod for Parallel applications
#'
#' can be used inside mclapply.  note that there is no checks in here as they are assumed
#' to be ran in the parent user function
#'
#' @param xt vector of time series observations
#' @param k single k that will be used for the f test
#' @param p largest degree of the polynomial modulation suspected
#' @param deltat  = 1 by default
#' @param w multitaper bandwidth parameter
#' @param dpss = FALSE uses sine instead
#' @param undersampleNumber = 100 by default  this is the number used in undersampling the tapers
#' @param confLevel level of confidence used in F3Mod, default is 1-1/N
#' @param returnFTestVars = FALSE.  Used for when more information about inner Ftests is needed. NOTE: this will be very memory intensive
#' @param penalty 1 is no penalty , 0.2  would give seq(from =  1, to =  1/(0.2*k, length.out = k)
#' penalty to each respective taper
#' @param penaltyOnTapersStdInv applys the penalty to all tapers in the calculation instead of just weighting the
#' eigenCoef's
#'
#' @return $F3Mod, $Freq, $significantFreq, $k

singleIterationForParallel <- function(xt, k, p, deltat = 1, w = NULL, dpss = FALSE,
                                       undersampleNumber = 100,
                                       confLevel = (1-(1/length(xt))),
                                       # altSig = FALSE,
                                       returnFTestVars = FALSE,
                                       penalty = 1, penaltyType = "ScaledExp",
                                       penaltyOnTapersStdInv = FALSE){

  N = length(xt)
  if(penaltyType == "mtm" & dpss == FALSE){
    stop("Adaptive weighting only works for dpss, set dpss = TRUE")
  }
  if(is.null(undersampleNumber)){
    stop("need to set undersample amount")
  }
  if(dpss){ #Use DPSS taper
    if(is.null(w)){
      stop("need to set w for dpss")
    }
    dp <- multitaper::dpss(n = N, k = k, nw = N*w)
    dpUnder <- multitaper::dpss(n = undersampleNumber, k = k, nw = N*w)
    instFreqEigen <- eigenCoefDPSSInstFrequency(xt = xt, N = N, k = k, w = w,
                                                    deltat = deltat,
                                                    returnDPSS = FALSE, passInDPSS = dp,
                                                    passInDPSSUnder = dpUnder, penalty = penalty,
                                                penaltyType = penaltyType)
    fStuff <- regressionDPSSInstFreq(N = N, k = k, w = w, instFreqEigen = instFreqEigen$PSI,
                                     p = p, passInDPSS = dp ,returnDPSS = FALSE,
                                     returnRp = FALSE, withoutzeroPoly = TRUE)
  }
  else{ #Sine Tapers are used

    sine <- sineTaperMatrix(N = N, k = k)
    sineUnder <- sineTaperMatrix(N = undersampleNumber, k = k)


    instFreqEigen <- eigenCoefSineInstFrequency(xt = xt, N = N, k = k,deltat = deltat,
                                                    returnSineMat = FALSE, passInSineTapers = sine,
                                                    passInSineUnder = sineUnder, penalty = penalty,
                                                penaltyType = penaltyType,
                                                penaltyOnTapersStdInv = penaltyOnTapersStdInv)

    fStuff <- regressionSineInstFreq(N = N, k = k, instFreqEigen = instFreqEigen$PSI,
                                     p = p, returnSineTapers = FALSE,
                                     passInSineMat = sine,
                                     returnRp = FALSE, withoutzeroPoly = TRUE)
  }


  #removingzero and nyquist frequencies
  zeroNyquist <- c(length(instFreqEigen$Freq),which(instFreqEigen$Freq == 0))
  instFreqEigen$Freq <- instFreqEigen$Freq[-zeroNyquist]
  Freq <- instFreqEigen$Freq
  instFreqEigen$PSI <- instFreqEigen$PSI[,-zeroNyquist]
  fStuff$cHat <- fStuff$cHat[,-zeroNyquist]

  normPhiSq <- colSums(instFreqEigen$PSI^2)
  normcHatWOutZeroSq <- 0


 if(!returnFTestVars){

     if(p ==1){
       F3 <-  matrix(nrow = 1, ncol = length(instFreqEigen$Freq))
       colnames(F3) <- Freq
       significantFreq <- list()

       normcHatWOutZeroSq <- normcHatWOutZeroSq + fStuff$cHat^2
       F3[p,] <- (fStuff$cHat)^2/
         ((normPhiSq - normcHatWOutZeroSq)/(k - p))
       FcutOff <- qf(confLevel, df1 = 1, df2 = (k-p), lower.tail = TRUE)
       significantFreq[[p]] <- Freq[which(F3[p,] >= FcutOff)]

     }else{
       F3 <-  matrix(nrow = nrow(fStuff$cHat), ncol = length(instFreqEigen$Freq))
       colnames(F3) <- Freq
       significantFreq <- list()
       for(P in 1:nrow(fStuff$cHat)){ # this is 1:p as we are removing zero so P-1 is actually P

         normcHatWOutZeroSq <- normcHatWOutZeroSq + fStuff$cHat[P,]^2
         F3[P,] <- (fStuff$cHat[P,])^2/
           ((normPhiSq - normcHatWOutZeroSq)/(k - P))
         FcutOff <- qf(confLevel, df1 = 1, df2 = (k-P), lower.tail = TRUE)
         significantFreq[[P]] <- Freq[which(F3[P,] >= FcutOff)]


       }
     }
     return(list(F3Mod = F3, Freq = Freq, significantFreq = significantFreq,
                 k = k))

 }else{# if the Ftest variables needed to be returned
   if(dpss){
     tapers <- dpUnder
   }else{
     tapers <- sineUnder
   }

     if(p ==1){
       F3 <-  matrix(nrow = 1, ncol = length(instFreqEigen$Freq))
       colnames(F3) <- Freq
       significantFreq <- list()

       normcHatWOutZeroSq <- normcHatWOutZeroSq + fStuff$cHat^2
       F3[p,] <- (fStuff$cHat)^2/
         ((normPhiSq - normcHatWOutZeroSq)/(k - p))
       FcutOff <- qf(confLevel, df1 = 1, df2 = (k-p-1), lower.tail = TRUE)
       significantFreq[[p]] <- Freq[which(F3[p,] >= FcutOff)]

     }else{
       F3 <-  matrix(nrow = nrow(fStuff$cHat), ncol = length(instFreqEigen$Freq))
       colnames(F3) <- Freq
       significantFreq <- list()
       for(P in 1:nrow(fStuff$cHat)){ # this is 1:p as we are removing zero so P-1 is actually P

         normcHatWOutZeroSq <- normcHatWOutZeroSq + fStuff$cHat[P,]^2
         F3[P,] <- (fStuff$cHat[P,])^2/
           ((normPhiSq - normcHatWOutZeroSq)/(k - P))
         FcutOff <- qf(confLevel, df1 = 1, df2 = (k-P-1), lower.tail = TRUE)
         significantFreq[[P]] <- Freq[which(F3[P,] >= FcutOff)]


       }
     }
     return(list(F3Mod = F3, Freq = Freq, significantFreq = significantFreq,
                 k = k, ftestvars = list(instFreqEigen = instFreqEigen,
                                         fStuff = fStuff, tapers = tapers)))
 }



}


