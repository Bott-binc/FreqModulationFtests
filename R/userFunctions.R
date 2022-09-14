# Data Generation ----------------------------------------

#' Modulation Generation for grumble distribution
#'
#' @param N Length of series you want to create
#' @param P Highest degree polynomial will use all degrees less than P other than 0
#' @param BLinear desired bandwidth for linear modulation
#' @param linConstCoef linear constant term
#' @param linCoef linear term coefficient
#' @param BQuadratic desired bandwidth for quadratic modulation
#' @param quadConstCoef quadratic constant term
#' @param quadLinCoef quadratic linear term coefficient
#' @param quadCoef quadratic term coefficient
#' @param BCubic desired bandwidth for cubic modulation
#' @param cubeConstCoef cubic constant term
#' @param cubeLinCoef cubic linear term coefficient
#' @param cubeQuadCoef cubic quadratic term coefficient
#' @param cubeCoef cubic term coefficient
#' @param quartConstCoef quartic constant term
#' @param quartLinCoef quartic linear term coefficient
#' @param quartQuadCoef quartic quadratic term coefficient
#' @param quartCubeCoef quartic cubic term coefficient
#' @param quartCoef quartic term coefficient
#' @param BQuartic desired bandwidth for quartic modulation
#' @param AmpLinear amplitude of linear modulation
#' @param AmpQuad amplitude of quadratic modulation
#' @param AmpCube amplitide of cubic modulation
#' @param AmpQuart amplitude of quartic modulation
#' @param fLin desired linear modulation carrier frequency
#' @param fQuad desired quadratic modulation carrier frequency
#' @param fCube desired cubic modulation carrier frequency
#' @param fQuart desired quartic modulation carrier frequency
#' @param checkBandWidth def = FALSE  if TRUE will print out the bandwidths to check
#' @param ar  =  0.5, 0.3, -0.1 AR coef for noise generation
#' @param ma  = 0.6 MA coef for noise generation
#' @param noiseScale = 1*6/pi^2 ratio of noise to pure signal
#' @param plotXt def = FALSE, will plot xt with noise if needed
#'
#' @return $xt for noisy data and $xtNoNoise for signal without noise and $noise for just the pure noise that was used
#' @export
grumbelModulationGeneration <- function(N,P,
                                        BLinear, linConstCoef, linCoef,
                                        AmpLinear,
                                        BQuadratic = 0, quadConstCoef = 0, quadLinCoef = 0, quadCoef = 0,
                                        AmpQuad = 0,
                                        BCubic = 0, cubeConstCoef = 0, cubeLinCoef = 0, cubeQuadCoef = 0, cubeCoef = 0,
                                        AmpCube = 0,
                                        BQuartic = 0, quartConstCoef = 0, quartLinCoef = 0, quartQuadCoef = 0, quartCubeCoef = 0, quartCoef = 0,
                                        AmpQuart = 0,
                                        fLin = 0.1, fQuad = 0.3, fCube = 0.31, fQuart = 0.4,
                                        checkBandWidth = FALSE, ar = c(0.5, 0.3, -0.1),
                                        ma = c(0.6), noiseScale = 1*6/pi^2, plotXt = FALSE, seed = NULL

                                        ){

  if(P %in% 1:4){

    n <- 0:(N-1)
    nFFT <- 2^ceiling(log2(2*N))
    freq <- seq(1/nFFT,0.5, by = 1/nFFT)

    # To ensure we are choosing a frequency that will be contained in the series
    f1 <- freq[which.min(abs(freq-fLin))]
    f2 <- freq[which.min(abs(freq-fQuad))]
    f3 <- freq[which.min(abs(freq-fCube))]
    f4 <- freq[which.min(abs(freq-fQuart))]

    # creating time indexes
    tstep <- 2.0/(N-1) # creates the step for the specific N so we end up with t in -1 to 1
    tt <- n * tstep - 1.0  # this runs from -1 to 1

    #modulating functions
    Linear <-  (linConstCoef + linCoef * tt)   #linear
    Quadratic <-  (quadConstCoef + quadLinCoef * tt + quadCoef * tt^2)     #quadratic
    Cubic <-  (cubeConstCoef + cubeLinCoef*tt + cubeQuadCoef*tt^2 + cubeCoef*tt^3) #cubic
    Quartic <-  (quartConstCoef + quartLinCoef*tt + quartQuadCoef*tt^2 +
                           quartCubeCoef*tt^3 + quartCoef*tt^4) # quartic

    bwLin <- max(abs(Linear))
    bwQuad <- max(abs(Quadratic))
    bwCube <- max(abs(Cubic))
    bwQuart <- max(abs(Quartic))

    if(P == 1){
      correctionLinearbw <- BLinear/(ceiling((bwLin)*100)/100) # finds closest correction factor to three digits below the desired bandwidth

      FMLinear <- Linear*correctionLinearbw

      if(checkBandWidth){
        print(paste0("linear = " ,max(abs(FMLinear)))) # just smaller than the w0
      }

      #then computing the 'integrals'

      modulationLinear <- cumsum(FMLinear)*2*pi
      InnerCosLin <- 2*pi*f1*n + modulationLinear

      modulation <- AmpLinear*cos(InnerCosLin)
    }
    else if(P == 2){

      correctionLinearbw <- BLinear/(ceiling((bwLin)*100)/100) # finds closest correction factor to three digits below the desired bandwidth
      correctionQuadbw <- BQuadratic/(ceiling((bwQuad)*100)/100)

      FMLinear <- Linear*correctionLinearbw
      FMQuadratic <- Quadratic*correctionQuadbw


      if(checkBandWidth){
        print(paste0("linear = " ,max(abs(FMLinear)))) # just smaller than the w0
        print(paste0("Quadratic = ",max(abs(FMQuadratic))))

      }

      #then computing the 'integrals'

      modulationLinear <- cumsum(FMLinear)*2*pi
      modulationQuadratic <- cumsum(FMQuadratic)*2*pi

      InnerCosLin <- 2*pi*f1*n + modulationLinear
      InnerCosQuad <- 2*pi*f2*n + modulationQuadratic

      modulation <- AmpLinear*cos(InnerCosLin) + AmpQuad*cos(InnerCosQuad)

    }
    else if(P == 3){

      correctionLinearbw <- BLinear/(ceiling((bwLin)*100)/100) # finds closest correction factor to three digits below the desired bandwidth
      correctionQuadbw <- BQuadratic/(ceiling((bwQuad)*100)/100)
      correctionCubebw <- BCubic/(ceiling((bwCube)*100)/100)

      FMLinear <- Linear*correctionLinearbw
      FMQuadratic <- Quadratic*correctionQuadbw
      FMCubic <- Cubic*correctionCubebw


      if(checkBandWidth){
        print(paste0("linear = " ,max(abs(FMLinear)))) # just smaller than the w0
        print(paste0("Quadratic = ",max(abs(FMQuadratic))))
        print(paste0("Cubic = ",max(abs(FMCubic))))

      }

      #then computing the 'integrals'

      modulationLinear <- cumsum(FMLinear)*2*pi
      modulationQuadratic <- cumsum(FMQuadratic)*2*pi
      modulationCubic <- cumsum(FMCubic)*2*pi


      InnerCosLin <- 2*pi*f1*n + modulationLinear
      InnerCosQuad <- 2*pi*f2*n + modulationQuadratic
      InnerCosCube <- 2*pi*f3*n + modulationCubic

      modulation <- AmpLinear*cos(InnerCosLin) + AmpQuad*cos(InnerCosQuad) +
                    AmpCube*cos(InnerCosCube)
    }
    else{
      correctionLinearbw <- BLinear/(ceiling((bwLin)*100)/100) # finds closest correction factor to three digits below the desired bandwidth
      correctionQuadbw <- BQuadratic/(ceiling((bwQuad)*100)/100)
      correctionCubebw <- BCubic/(ceiling((bwCube)*100)/100)
      correctionQuartbw <- BQuartic/(ceiling((bwQuart)*100)/100)

      FMLinear <- Linear*correctionLinearbw
      FMQuadratic <- Quadratic*correctionQuadbw
      FMCubic <- Cubic*correctionCubebw
      FMQuartic <- Quartic*correctionQuartbw

      if(checkBandWidth){
        print(paste0("linear = " ,max(abs(FMLinear)))) # just smaller than the w0
        print(paste0("Quadratic = ",max(abs(FMQuadratic))))
        print(paste0("Cubic = ",max(abs(FMCubic))))
        print(paste0("Quartic = ",max(abs(FMQuartic))))
      }

      #then computing the 'integrals'

      modulationLinear <- cumsum(FMLinear)*2*pi
      modulationQuadratic <- cumsum(FMQuadratic)*2*pi
      modulationCubic <- cumsum(FMCubic)*2*pi
      modulationQuartic <- cumsum(FMQuartic)*2*pi

      InnerCosLin <- 2*pi*f1*n + modulationLinear
      InnerCosQuad <- 2*pi*f2*n + modulationQuadratic
      InnerCosCube <- 2*pi*f3*n + modulationCubic
      InnerCosQuart <- 2*pi*f4*n + modulationQuartic

      modulation <- AmpLinear*cos(InnerCosLin) + AmpQuad*cos(InnerCosQuad) +
        AmpCube*cos(InnerCosCube) + AmpQuart*cos(InnerCosQuart)
    }

    if(!is.null(seed)){ # allows the user to keep the noise generation the same
      set.seed(seed)
    }
    ARMA <- list(ar = ar, ma = ma)
    noiseInnov <- VGAM::rgumbel(N, scale = noiseScale)
    noise <- stats::arima.sim(model = ARMA, n = N, innov = noiseInnov)
    xt <- modulation + as.numeric(noise)

  }else{
    stop("P can only be up to degree 4")
  }

  if(plotXt){
    plot(xt, x = 1:N, type = "l")
  }

  return(list(xt = xt, xtNoNoise = modulation, noise = noise))
}



# F Tests -------------------------------------------------




#' F1 test statistic
#'
#' @param xt time series
#' @param N Total number of observations
#' @param k Number of tapers
#' @param p Highest degree polynomial you want to test for
#' @param deltat Time interval between each observation
#' @param w Only needed if dpss = TRUE
#' @param dpss  = FALSE unless you want to use dpss, it will do sine tapers by delfault
#' @param returnInstFreqAndRegression  = FALSE, this speeds up the other f tests so you can pass in information
#' @param withoutZeroDegree TRUE if wanting modified test statistic that does not use zero degree polynomials
#'
#' @return $F1testStat and $Freq corresponding to the f1 test statistic, if reutrnInstFreqAndRegression = TRUE
#' it will also return $necessaryTestStuff$instFreqEigen and $necessaryTestStuff$regressionInstFreq
#'
#' @export
F1Test <- function(xt, N, k, p, deltat = 1, w = NULL, dpss = FALSE, returnInstFreqAndRegression = FALSE,
                   withoutZeroDegree = FALSE){

  if(!withoutZeroDegree){
    if(dpss){ # dpss is used
      if(is.null(w)){
        stop("need to set w for dpss")
      }
      instFreqEigen <- eigenSpectrumDPSSInstFrequency(xt = xt, N = N, k = k, w = w, deltat = deltat,
                                                      returnDPSS = TRUE)
      fStuff <- regressionDPSSInstFreq(N = N, k = k, w = w, instFreqEigen = instFreqEigen$PHI,
                                       p = p, passInDPSS = instFreqEigen$DPSS,returnDPSS = FALSE,
                                       returnRp = FALSE)
    }
    else{ # the sine tapers are used
      instFreqEigen <- eigenSpectrumSineInstFrequency(xt = xt, N = N, k = k,deltat = deltat,
                                                      returnSineMat = TRUE)

      fStuff <- regressionSineInstFreq(N = N, k = k, instFreqEigen = instFreqEigen$PHI,
                                       p = p, returnSineTapers = FALSE,
                                       passInSineMat = instFreqEigen$SineTaper,
                                       returnRp = FALSE)
    }
    # removing nyquist and the zero as they aren't complex (no -nyquist from FFt)
    Freq <- instFreqEigen$Freq[-c(length(instFreqEigen$Freq),which(instFreqEigen$Freq == 0))]
    F1 <- matrix(nrow = nrow(fStuff$cHat), ncol = length(instFreqEigen$Freq) - 2)
    colnames(F1) <- Freq
    for(P in 1:nrow(fStuff$cHat)){ # this P is 0 through P hense the P-1
      for(f in Freq){
        indexF1 <- which(colnames(F1) == f)
        index <- which(colnames(fStuff$cHat) == f)
        F1[P,indexF1] <- ((norm(fStuff$cHat[1:P,index], type = "2"))^2/((P - 1) + 1))/
          (((norm(instFreqEigen$PHI[,index], type = "2"))^2 - (norm(fStuff$cHat[1:P,index], type = "2"))^2)/(k - (P - 1) - 1))
      }
    }
  }else{ # without zero degree
    if(dpss){ # dpss is used
      if(is.null(w)){
        stop("need to set w for dpss")
      }
      instFreqEigen <- eigenSpectrumDPSSInstFrequency(xt = xt, N = N, k = k, w = w, deltat = deltat,
                                                      returnDPSS = TRUE)
      fStuff <- regressionDPSSInstFreq(N = N, k = k, w = w, instFreqEigen = instFreqEigen$PHI,
                                       p = p, passInDPSS = instFreqEigen$DPSS,returnDPSS = FALSE,
                                       returnRp = FALSE, withoutzeroPoly =
                                         TRUE)
    }
    else{ # the sine tapers are used
      instFreqEigen <- eigenSpectrumSineInstFrequency(xt = xt, N = N, k = k,deltat = deltat,
                                                      returnSineMat = TRUE)

      fStuff <- regressionSineInstFreq(N = N, k = k, instFreqEigen = instFreqEigen$PHI,
                                       p = p, returnSineTapers = FALSE,
                                       passInSineMat = instFreqEigen$SineTaper,
                                       returnRp = FALSE, withoutzeroPoly = TRUE)
    }

    # removing nyquist and the zero as they aren't complex (no -nyquist from FFt)
    Freq <- instFreqEigen$Freq[-c(length(instFreqEigen$Freq),which(instFreqEigen$Freq == 0))]
    F1 <- matrix(nrow = nrow(fStuff$cHat), ncol = length(instFreqEigen$Freq) - 2)
    colnames(F1) <- Freq
    for(P in 1:nrow(fStuff$cHat)){ # this P is 1 through p
      for(f in Freq){
        indexF1 <- which(colnames(F1) == f)
        index <- which(colnames(fStuff$cHat) == f)
        F1[P,indexF1] <- ((norm(fStuff$cHat[1:P,index], type = "2"))^2/(P))/
          (((norm(instFreqEigen$PHI[,index], type = "2"))^2 - (norm(fStuff$cHat[1:P,index], type = "2"))^2)/(k - (P)))
      }
    }
  }




  #making the return
  if(!returnInstFreqAndRegression){
    return(list(F1testStat = F1, Freq = Freq))
  }
  else{
    return(list(F1testStat = F1, Freq = Freq, necessaryTestStuff = list(instFreqEigen = instFreqEigen,
                                                                        regressionInstFreq = fStuff)))
  }
}


#' F2 test statistic
#'
#' @param xt time series
#' @param N Total number of observations
#' @param k Number of tapers
#' @param p Highest degree polynomial you want to test for
#' @param deltat Time interval between each observation
#' @param w Only needed if dpss = TRUE
#' @param dpss  = FALSE unless you want to use dpss, it will do sine tapers by delfault
#' @param returnInstFreqAndRegression  = FALSE, this speeds up the other f tests so you can pass in information
#' @param passInInstFreqAndRegression leave null unless trying to simultainously do all test statistics
#'
#' @return $F2testStat $Freq and $necessaryTestStuff used to pass into another ftest function
#'
#' @export
F2Test <- function(xt, N, k, p, deltat = 1, w = NULL, dpss = FALSE, passInInstFreqAndRegression = NULL,
                   returnInstFreqAndRegression = FALSE){

  if(dpss){ # dpss is used
    if(is.null(w)){
      stop("need to set w for dpss")
    }
    if(is.null(passInInstFreqAndRegression)){
      instFreqEigen <- eigenSpectrumDPSSInstFrequency(xt = xt, N = N, k = k, w = w, deltat = deltat,
                                                      returnDPSS = TRUE)
      fStuff <- regressionDPSSInstFreq(N = N, k = k, w = w, instFreqEigen = instFreqEigen$PHI,
                                       p = p, passInDPSS = instFreqEigen$DPSS,returnDPSS = FALSE,
                                       returnRp = FALSE, returnGCHatp = TRUE)
    }else{
      instFreqEigen <- passInInstFreqAndRegression$instFreqEigen
      fStuff <- passInInstFreqAndRegression$regressionInstFreq
    }
  }
  else{ # the sine tapers are used
    if(is.null(passInInstFreqAndRegression)){
      instFreqEigen <- eigenSpectrumSineInstFrequency(xt = xt, N = N, k = k,deltat = deltat,
                                                      returnSineMat = TRUE)

      fStuff <- regressionSineInstFreq(N = N, k = k, instFreqEigen = instFreqEigen$PHI,
                                       p = p, returnSineTapers = FALSE,
                                       passInSineMat = instFreqEigen$SineTaper,
                                       returnRp = FALSE, returnGCHatp = TRUE)
    }else{
      instFreqEigen <- passInInstFreqAndRegression$instFreqEigen
      fStuff <- passInInstFreqAndRegression$regressionInstFreq
    }
  }


  # removing nyquist and the zero as they aren't complex (no -nyquist from FFt)
  Freq <- instFreqEigen$Freq[-c(length(instFreqEigen$Freq),which(instFreqEigen$Freq == 0))]
  F2 <- matrix(nrow = nrow(fStuff$cHat), ncol = length(instFreqEigen$Freq) - 2)
  colnames(F2) <- Freq
  for(P in 1:nrow(fStuff$cHat)){
    for(f in Freq){
      indexF2 <- which(colnames(F2) == f)
      index <- which(colnames(fStuff$cHat) == f)
      F2[P,indexF2] <- ((norm(fStuff$G[,1:P] %*% as.matrix(fStuff$cHat[1:P,index]), type = "2"))^2/((P-1)+1))/
        (((norm(instFreqEigen$PHI[,index], type = "2"))^2 - (norm(fStuff$cHat[1:P,index], type = "2"))^2)/(k - (P-1) - 1))
    }
  }

  #making the return
  if(!returnInstFreqAndRegression){
    return(list(F2testStat = F2, Freq = Freq))
  }
  else{
    return(list(F2testStat = F2, Freq = Freq, necessaryTestStuff = list(instFreqEigen = instFreqEigen,
                                                                        regressionInstFreq = fStuff)))
  }
}


#' F3 test statistic
#'
#' @param xt time series
#' @param N Total number of observations
#' @param k Number of tapers
#' @param p Highest degree polynomial you want to test for
#' @param deltat Time interval between each observation
#' @param w Only needed if dpss = TRUE
#' @param dpss  = FALSE unless you want to use dpss, it will do sine tapers by default
#' @param returnInstFreqAndRegression  = FALSE, this speeds up the other f tests so you can pass in information
#' @param withoutZeroDegree TRUE if wanting modified test statistic that does not use zero degree polynomials (no undersampling if FALSE)
#' @param undersample True or FALSE, allows for faster run time while maintaining most accuracy, note that this will
#' also cause zero padding to take place.  the user DOES NOT need to manually zero pad
#' @param undersampleNumber A numeric of the number the user wants to undersample, usually 100 is a good start
#'
#' @return $F3testStat, $Freq, $necessaryTestStuff used to pass into another ftest function
#'
#' @export
F3Test <- function(xt, N, k, p, deltat = 1, w = NULL, dpss = FALSE,
                   returnInstFreqAndRegression = FALSE, withoutZeroDegree = TRUE,
                   undersample = FALSE, undersampleNumber = NULL){

  if(!withoutZeroDegree){
    if(dpss){ # dpss is used
      if(is.null(w)){
        stop("need to set w for dpss")
      }
        instFreqEigen <- eigenSpectrumDPSSInstFrequency(xt = xt, N = N, k = k, w = w, deltat = deltat,
                                                        returnDPSS = TRUE)
        fStuff <- regressionDPSSInstFreq(N = N, k = k, w = w, instFreqEigen = instFreqEigen$PHI,
                                         p = p, passInDPSS = instFreqEigen$DPSS,returnDPSS = FALSE,
                                         returnRp = FALSE)
    }
    else{ # the sine tapers are used

        instFreqEigen <- eigenSpectrumSineInstFrequency(xt = xt, N = N, k = k,deltat = deltat,
                                                        returnSineMat = TRUE)

        fStuff <- regressionSineInstFreq(N = N, k = k, instFreqEigen = instFreqEigen$PHI,
                                         p = p, returnSineTapers = FALSE,
                                         passInSineMat = instFreqEigen$SineTaper,
                                         returnRp = FALSE)

    }

    #removingzero and nyquist frequencies
    zeroNyquist <- c(length(instFreqEigen$Freq),which(instFreqEigen$Freq == 0))
    instFreqEigen$Freq <- instFreqEigen$Freq[-zeroNyquist]
    Freq <- instFreqEigen$Freq
    instFreqEigen$PHI <- instFreqEigen$PHI[,-zeroNyquist]
    fStuff$cHat <- fStuff$cHat[,-zeroNyquist]


    normPhiSq <- colSums(instFreqEigen$PHI^2)
    normcHatwithZeroSq <- 0
    # removing nyquist and the zero as they aren't complex (no -nyquist from FFt)

    F3 <- matrix(nrow = nrow(fStuff$cHat), ncol = length(instFreqEigen$Freq))
    colnames(F3) <- Freq
    for(P in 1:nrow(fStuff$cHat)){ # this is 0 to p

        normcHatwithZeroSq <- normcHatwithZeroSq + fStuff$cHat[P,]^2
        F3[P,] <- (fStuff$cHat[P,])^2/
          (((normPhiSq - normcHatwithZeroSq)/(k - (P-1) - 1)))

    }
  }else{#not using zeroth degree in the test. (best test statistic at this point)
    if(!undersample){
      if(dpss){ #Use DPSS taper
        if(is.null(w)){
          stop("need to set w for dpss")
        }

        instFreqEigen <- eigenSpectrumDPSSInstFrequency(xt = xt, N = N, k = k, w = w, deltat = deltat,
                                                        returnDPSS = TRUE,)
        fStuff <- regressionDPSSInstFreq(N = N, k = k, w = w, instFreqEigen = instFreqEigen$PHI,
                                         p = p, passInDPSS = instFreqEigen$DPSS,returnDPSS = FALSE,
                                         returnRp = FALSE,  withoutzeroPoly = TRUE)
      }
      else{ #Sine Tapers are used
        instFreqEigen <- eigenSpectrumSineInstFrequency(xt = xt, N = N, k = k,deltat = deltat,
                                                        returnSineMat = TRUE)

        fStuff <- regressionSineInstFreq(N = N, k = k, instFreqEigen = instFreqEigen$PHI,
                                         p = p, returnSineTapers = FALSE,
                                         passInSineMat = instFreqEigen$SineTaper,
                                         returnRp = FALSE, withoutzeroPoly = TRUE)
      }
    }else{
      if(is.null(undersampleNumber)){
        stop("need to set undersample amount")
      }
      if(dpss){ #Use DPSS taper
        if(is.null(w)){
          stop("need to set w for dpss")
        }
        dp <- multitaper::dpss(n = N, k = k, nw = N*w)
        dpUnder <- multitaper::dpss(n = undersampleNumber, k = k, nw = N*w)
        instFreqEigen <- eigenSpectrumDPSSInstFrequency(xt = xt, N = N, k = k, w = w, deltat = deltat,
                                                        returnDPSS = FALSE, passInDPSS = dp,
                                                        passInDPSSUnder = dpUnder)
        fStuff <- regressionDPSSInstFreq(N = N, k = k, w = w, instFreqEigen = instFreqEigen$PHI,
                                         p = p, passInDPSS = dpUnder ,returnDPSS = FALSE,
                                         returnRp = FALSE, withoutzeroPoly = TRUE)
      }
      else{ #Sine Tapers are used
        sine <- sineTaperMatrix(N = N, k = k)
        sineUnder <- sineTaperMatrix(N = undersampleNumber, k = k)
        instFreqEigen <- eigenSpectrumSineInstFrequency(xt = xt, N = N, k = k,deltat = deltat,
                                                        returnSineMat = FALSE, passInSineTapers = sine,
                                                        passInSineUnder = sineUnder)

        fStuff <- regressionSineInstFreq(N = N, k = k, instFreqEigen = instFreqEigen$PHI,
                                         p = p, returnSineTapers = FALSE,
                                         passInSineMat = sineUnder,
                                         returnRp = FALSE, withoutzeroPoly = TRUE)
      }
    }

    #removingzero and nyquist frequencies
    zeroNyquist <- c(length(instFreqEigen$Freq),which(instFreqEigen$Freq == 0))
    instFreqEigen$Freq <- instFreqEigen$Freq[-zeroNyquist]
    Freq <- instFreqEigen$Freq
    instFreqEigen$PHI <- instFreqEigen$PHI[,-zeroNyquist]
    fStuff$cHat <- fStuff$cHat[,-zeroNyquist]

    normPhiSq <- colSums(instFreqEigen$PHI^2)
    normcHatWOutZeroSq <- 0

    F3 <-  matrix(nrow = nrow(fStuff$cHat), ncol = length(instFreqEigen$Freq))
    colnames(F3) <- Freq
    for(P in 1:nrow(fStuff$cHat)){ # this is 1:p as we are removing zero so P-1 is actually P

        normcHatWOutZeroSq <- normcHatWOutZeroSq + fStuff$cHat[P,]^2
        F3[P,] <- (fStuff$cHat[P,])^2/
          ((normPhiSq - normcHatWOutZeroSq)/(k - P))

    }
  }

  #making the return
  if(!returnInstFreqAndRegression){
    return(list(F3testStat = F3, Freq = Freq))
  }
  else{
    return(list(F3testStat = F3, Freq = Freq, necessaryTestStuff = list(instFreqEigen = instFreqEigen,
                                                                        regressionInstFreq = fStuff)))
  }
}


#' A combined F test that returns F1 modified, F3 modified and non modified modified referring to removing the
#' zero polynomial test  It is the fastest version if you want all three tests at the same time
#'
#' @param xt time series
#' @param N Total number of observations
#' @param k Number of tapers
#' @param p Highest degree polynomial you want to test for
#' @param deltat Time interval between each observation
#' @param w Only needed if dpss = TRUE
#' @param dpss  = FALSE unless you want to use dpss, it will do sine tapers by default
#' @param returnInstFreqAndRegression  = FALSE, this speeds up the other f tests so you can pass in information
#' @param undersample True or FALSE, allows for faster run time while maintaining most accuracy, note that this will
#' also cause zero padding to take place.  the user DOES NOT need to manually zero pad
#' @param undersampleNumber A numeric of the number the user wants to undersample, usually 100 is a good start
#'
#' @return $F1Mod $F3Mod without 0Poly and $F3 with 0 poly, as well as the frequencies.  if returnInstFreq = TRUE
#' it will also reuturn all things needed to run the other test stats.
#'
#' @export
FtestCombined <- function(xt, N, k, p, deltat = 1, w = NULL, dpss = FALSE,
                          returnInstFreqAndRegression = FALSE,
                          undersample = FALSE, undersampleNumber = NULL){

  if(!returnInstFreqAndRegression){
    if(!undersample){
      if(dpss){ #Use DPSS taper
        if(is.null(w)){
          stop("need to set w for dpss")
        }

        instFreqEigen <- eigenSpectrumDPSSInstFrequency(xt = xt, N = N, k = k, w = w, deltat = deltat,
                                                        returnDPSS = TRUE)
        fStuff <- regressionDPSSInstFreq(N = N, k = k, w = w, instFreqEigen = instFreqEigen$PHI,
                                         p = p, passInDPSS = instFreqEigen$DPSS,returnDPSS = FALSE,
                                         returnRp = FALSE)
      }
      else{ #Sine Tapers are used
        instFreqEigen <- eigenSpectrumSineInstFrequency(xt = xt, N = N, k = k,deltat = deltat,
                                                        returnSineMat = TRUE)

        fStuff <- regressionSineInstFreq(N = N, k = k, instFreqEigen = instFreqEigen$PHI,
                                         p = p, returnSineTapers = FALSE,
                                         passInSineMat = instFreqEigen$SineTaper,
                                         returnRp = FALSE, withoutzeroPoly = FALSE)
      }
    }else{
      if(is.null(undersampleNumber)){
        stop("need to set undersample amount")
      }
      if(dpss){ #Use DPSS taper
        if(is.null(w)){
          stop("need to set w for dpss")
        }
        dp <- multitaper::dpss(n = N, k = k, nw = N*w)
        dpUnder <- multitaper::dpss(n = undersampleNumber, k = k, nw = N*w)
        instFreqEigen <- eigenSpectrumDPSSInstFrequency(xt = xt, N = N, k = k, w = w, deltat = deltat,
                                                        returnDPSS = FALSE, passInDPSS = dp,
                                                        passInDPSSUnder = dpUnder)
        fStuff <- regressionDPSSInstFreq(N = N, k = k, w = w, instFreqEigen = instFreqEigen$PHI,
                                         p = p, passInDPSS = dpUnder ,returnDPSS = FALSE,
                                         returnRp = FALSE)
      }
      else{ #Sine Tapers are used
        sine <- sineTaperMatrix(N = N, k = k)
        sineUnder <- sineTaperMatrix(N = undersampleNumber, k = k)
        instFreqEigen <- eigenSpectrumSineInstFrequency(xt = xt, N = N, k = k,deltat = deltat,
                                                        returnSineMat = FALSE, passInSineTapers = sine,
                                                        passInSineUnder = sineUnder)

        fStuff <- regressionSineInstFreq(N = N, k = k, instFreqEigen = instFreqEigen$PHI,
                                         p = p, returnSineTapers = FALSE,
                                         passInSineMat = sineUnder,
                                         returnRp = FALSE, withoutzeroPoly = FALSE)
      }
    }
  }else{ #we want to return and calculate the Rp too
    if(!undersample){
      if(dpss){ #Use DPSS taper
        if(is.null(w)){
          stop("need to set w for dpss")
        }

        instFreqEigen <- eigenSpectrumDPSSInstFrequency(xt = xt, N = N, k = k, w = w, deltat = deltat,
                                                        returnDPSS = TRUE)
        fStuff <- regressionDPSSInstFreq(N = N, k = k, w = w, instFreqEigen = instFreqEigen$PHI,
                                         p = p, passInDPSS = instFreqEigen$DPSS,returnDPSS = FALSE,
                                         returnRp = TRUE)
        taperMat <- instFreqEigen$DPSS$v
      }
      else{ #Sine Tapers are used
        instFreqEigen <- eigenSpectrumSineInstFrequency(xt = xt, N = N, k = k,deltat = deltat,
                                                        returnSineMat = TRUE)

        fStuff <- regressionSineInstFreq(N = N, k = k, instFreqEigen = instFreqEigen$PHI,
                                         p = p, returnSineTapers = FALSE,
                                         passInSineMat = instFreqEigen$SineTaper,
                                         returnRp = TRUE, withoutzeroPoly = FALSE)
        taperMat <- instFreqEigen$SineTaper
      }
    }else{
      if(is.null(undersampleNumber)){
        stop("need to set undersample amount")
      }
      if(dpss){ #Use DPSS taper
        if(is.null(w)){
          stop("need to set w for dpss")
        }
        dp <- multitaper::dpss(n = N, k = k, nw = N*w)
        dpUnder <- multitaper::dpss(n = undersampleNumber, k = k, nw = N*w)
        instFreqEigen <- eigenSpectrumDPSSInstFrequency(xt = xt, N = N, k = k, w = w, deltat = deltat,
                                                        returnDPSS = FALSE, passInDPSS = dp,
                                                        passInDPSSUnder = dpUnder)
        fStuff <- regressionDPSSInstFreq(N = N, k = k, w = w, instFreqEigen = instFreqEigen$PHI,
                                         p = p, passInDPSS = dpUnder ,returnDPSS = FALSE,
                                         returnRp = TRUE)
        taperMat <- dpUnder$v
      }
      else{ #Sine Tapers are used
        sine <- sineTaperMatrix(N = N, k = k)
        sineUnder <- sineTaperMatrix(N = undersampleNumber, k = k)
        instFreqEigen <- eigenSpectrumSineInstFrequency(xt = xt, N = N, k = k,deltat = deltat,
                                                        returnSineMat = FALSE, passInSineTapers = sine,
                                                        passInSineUnder = sineUnder)

        fStuff <- regressionSineInstFreq(N = N, k = k, instFreqEigen = instFreqEigen$PHI,
                                         p = p, returnSineTapers = FALSE,
                                         passInSineMat = sineUnder,
                                         returnRp = TRUE, withoutzeroPoly = FALSE)
        taperMat <- sineUnder
      }
    }
  }




  #removingzero and nyquist frequencies
  zeroNyquist <- c(length(instFreqEigen$Freq),which(instFreqEigen$Freq == 0))
  instFreqEigen$Freq <- instFreqEigen$Freq[-zeroNyquist]
  Freq <- instFreqEigen$Freq
  instFreqEigen$PHI <- instFreqEigen$PHI[,-zeroNyquist]
  fStuff$cHat <- fStuff$cHat[,-zeroNyquist]

  #need to calculate the H and Chat for the removed zero HERE
  HWoutZero <- fStuff$H[,-1] # removes the 0th order column
  #GWOutZero <- fStuff$G[,-1] dont think we need this for f1 and f3
  cHatWOutZero <- t(HWoutZero) %*% instFreqEigen$PHI
  #cHatWOutZero <- as.matrix(cHatWOutZero[, -zeroNyquist])

  F1Reduced <- matrix(nrow = nrow(cHatWOutZero), ncol = length(instFreqEigen$Freq))
  F3Reduced <- matrix(nrow = nrow(cHatWOutZero), ncol = length(instFreqEigen$Freq))
  F3 <- matrix(nrow = nrow(fStuff$cHat), ncol = length(instFreqEigen$Freq))
  colnames(F1Reduced) <- Freq
  colnames(F3Reduced) <- Freq
  colnames(F3) <- Freq

  normPhiSq <- colSums(instFreqEigen$PHI^2)
  normcHatWOutZeroSq <- 0
  normcHatwithZeroSq <- 0
  rpmodSqrdWithoutZero <- matrix(nrow = (nrow(fStuff$cHat)-1), ncol = length(Freq))
  colnames(rpmodSqrdWithoutZero) <- Freq
  for(P in 1:nrow(fStuff$cHat)){# this P is 1 through p
    normcHatwithZeroSq <- normcHatwithZeroSq + fStuff$cHat[P,]^2
    if(P != nrow(fStuff$cHat)){ # will do P = 1, ...P for the reduced and P = 0 to P for F3
      normcHatWOutZeroSq <- normcHatWOutZeroSq + cHatWOutZero[P,]^2
      F1Reduced[P,] <- (normcHatWOutZeroSq/(P))/
        ((normPhiSq - normcHatWOutZeroSq)/(k - (P)))
      F3Reduced[P,] <- (cHatWOutZero[P,])^2/
        ((normPhiSq - normcHatWOutZeroSq)/(k - P))
      F3[P,] <- (fStuff$cHat[P,])^2/
        (((normPhiSq - normcHatwithZeroSq)/(k - (P-1) - 1)))
      rpmodSqrdWithoutZero[P,] <- normPhiSq - normcHatwithZeroSq
    }else{
      F3[P,] <- (fStuff$cHat[P,])^2/
        (((normPhiSq - normcHatwithZeroSq)/(k - (P-1) - 1)))
    }
  }

  if(!returnInstFreqAndRegression){
    return(list(F1Mod = F1Reduced,
                F3Mod = F3Reduced,
                F3 = F3,
                Freq = Freq))
  }else{ #if they want the regression stuff to be returned as well
    return(list(F1Mod = F1Reduced,
                F3Mod = F3Reduced,
                F3 = F3,
                Freq = Freq,
                necessaryTestStuff = list(instFreqEigen = instFreqEigen,
                                          regressionInstFreq = fStuff,
                                          cHatWOutZero = cHatWOutZero,
                                          rpmodSqrdWithoutZero = rpmodSqrdWithoutZero,
                                          rp = fStuff$rp)))#(taperMat %*% (diag(k) - HWoutZero %*% t(HWoutZero)))%*% instFreqEigen$PHI )))
  }
}

#' F4Test
#'
#' @param xt time series
#' @param N Total number of observations
#' @param k Number of tapers
#' @param p Highest degree polynomial you want to test for
#' @param deltat Time interval between each observation
#' @param w Only needed if dpss = TRUE
#' @param dpss  = FALSE unless you want to use dpss, it will do sine tapers by default
#' @param undersampleNumber A numeric of the number the user wants to undersample, usually 100 is a good start
#'
#' @return $F4testStat, $Freq, $F14TestStat this is the f1 not the f3 mod
#'
#' @export
F4Test <- function(xt, N, k, p, deltat = 1, w = NULL, dpss = FALSE, undersampleNumber = 100){
  initial <- TRUE
  loopNum <- 0
  for(K in k){
    loopNum <- loopNum + 1

    if(is.null(undersampleNumber)){
      stop("need to set undersample amount")
    }
    if(dpss){ #Use DPSS taper
      if(is.null(w)){
        stop("need to set w for dpss")
      }
      dp <- multitaper::dpss(n = N, k = K, nw = N*w[loopNum])
      dpUnder <- multitaper::dpss(n = undersampleNumber, k = K, nw = N*w[loopNum])
      instFreqEigen <- eigenSpectrumDPSSInstFrequency(xt = xt, N = N, k = K, w = w[loopNum],
                                                      deltat = deltat,
                                                      returnDPSS = FALSE, passInDPSS = dp,
                                                      passInDPSSUnder = dpUnder)
      fStuff <- regressionDPSSInstFreq(N = N, k = K, w = w[loopNum], instFreqEigen = instFreqEigen$PHI,
                                       p = p, passInDPSS = dpUnder ,returnDPSS = FALSE,
                                       returnRp = FALSE, withoutzeroPoly = TRUE)
    }
    else{ #Sine Tapers are used
      sine <- sineTaperMatrix(N = N, k = K)
      sineUnder <- sineTaperMatrix(N = undersampleNumber, k = K)
      instFreqEigen <- eigenSpectrumSineInstFrequency(xt = xt, N = N, k = K,deltat = deltat,
                                                      returnSineMat = FALSE, passInSineTapers = sine,
                                                      passInSineUnder = sineUnder)

      fStuff <- regressionSineInstFreq(N = N, k = K, instFreqEigen = instFreqEigen$PHI,
                                       p = p, returnSineTapers = FALSE,
                                       passInSineMat = sineUnder,
                                       returnRp = FALSE, withoutzeroPoly = TRUE)
    }


    #removingzero and nyquist frequencies
    zeroNyquist <- c(length(instFreqEigen$Freq),which(instFreqEigen$Freq == 0))
    instFreqEigen$Freq <- instFreqEigen$Freq[-zeroNyquist]
    Freq <- instFreqEigen$Freq
    instFreqEigen$PHI <- instFreqEigen$PHI[,-zeroNyquist]
    fStuff$cHat <- fStuff$cHat[,-zeroNyquist]


    if(initial){                            # K iterations X P polynomials X zeroPadd size
      normcHatWOutZeroSq <- array(0, dim = c(length(k), nrow(fStuff$cHat), ncol(fStuff$cHat)))
      normPhiSq <- matrix(nrow = length(k), ncol = ncol(fStuff$cHat))
      cHat <- array(0, dim=c(length(k), nrow(fStuff$cHat), ncol(fStuff$cHat)))
      initial <- FALSE
    }
    normPhiSq[loopNum,] <- colSums(instFreqEigen$PHI^2)

    #p = 1 th for loop iteration
    normcHatWOutZeroSq[loopNum,1,] <-  fStuff$cHat[1,]^2
    cHat[loopNum,1,] <- fStuff$cHat[1,]^2
    # F3[1,] <- (fStuff$cHat[1,])^2/
    #   ((normPhiSq - normcHatWOutZeroSq[k,1,])/(k - 1))
    for(P in 2:nrow(fStuff$cHat)){ # this is 1:p as we are removing zero so P-1 is actually P
      cHat[loopNum,P,] <- fStuff$cHat[P,]^2
      normcHatWOutZeroSq[loopNum,P,] <- normcHatWOutZeroSq[loopNum,(P-1),] + fStuff$cHat[P,]^2
      # F3[P,] <- (fStuff$cHat[P,])^2/
      #   ((normPhiSq - normcHatWOutZeroSq[k,P,])/(k - P))

    }
  }
  F3 <- F1 <- matrix(nrow = nrow(fStuff$cHat), ncol = length(instFreqEigen$Freq))
  colnames(F3) <- Freq

  for(P in 1:nrow(fStuff$cHat)){
    F1[P,] <- (colSums(normcHatWOutZeroSq[,P,])/sum(rep(P, times = length(k))))/
      (((colSums(normPhiSq - normcHatWOutZeroSq[,P,])))/(sum((k - rep(P, times=length(k))))))
    F3[P,] <- (((colSums(cHat[,P,]^2)))/(length(k)))/
      (((colSums(normPhiSq - normcHatWOutZeroSq[,P,])))/(sum((k - rep(P, times=length(k))))))
  }
  #making the return
  return(list(F4testStat = F3, Freq = Freq,
              F14testStat = F1 ))
}



#' F4Test
#'
#' w is chosen by shannons number based on k
#'
#' @param xt time series
#' @param N Total number of observations
#' @param p Highest degree polynomial you want to test for
#' @param deltat Time interval between each observation
#' @param dpss  = FALSE unless you want to use dpss, it will do sine tapers by default
#' @param undersampleNumber A numeric of the number the user wants to undersample, usually 100 is a good start
#' @param k vector of tapers used in the f test
#' @param cores must be 1 if on windows, number of cores used for parallelization
#'
#' @return $F4testStat, $Freq, $F14TestStat this is the f1 not the f3 mod
#'
#' @export
F4Testpar <- function(xt, N, k, p, deltat = 1, dpss = FALSE, undersampleNumber = 100, cores = 1){

  if(is.null(undersampleNumber)){
    stop("need to set undersample amount")
  }

  if(dpss){
    fullDat <- parallel::mclapply(X = k,FUN = function(x){
      return(singleIterationForParallel4(xt = xt, N = N, k = x, w = ((x+1)/(2*length(xt))), p = p, deltat = deltat,
                               undersampleNumber = undersampleNumber, dpss = TRUE))
    }, mc.cores = cores, mc.cleanup = TRUE)
  }else{
    fullDat <- parallel::mclapply(X = k,FUN = function(x){
      return(singleIterationForParallel4(xt = xt, N = N, k = x, p = p, deltat = deltat,
                               undersampleNumber = undersampleNumber, dpss = FALSE))
    }, mc.cores = cores, mc.cleanup = TRUE)
  }

  if(p != 1){
                               # K iterations X P polynomials X zeroPadd size
  normcHatWOutZeroSq <- array(0, dim = c(length(k), nrow(fullDat[[1]]$cHat), ncol(fullDat[[1]]$cHat)))
  normPhiSq <- matrix(nrow = length(k), ncol = ncol(fullDat[[1]]$cHat))
  cHat <- array(0, dim=c(length(k), nrow(fullDat[[1]]$cHat), ncol(fullDat[[1]]$cHat)))
  }else{
    normcHatWOutZeroSq <- array(0, dim = c(length(k), length(fullDat[[1]]$cHat),1))
    normPhiSq <- matrix(nrow = length(k), ncol = 1)
    cHat <- array(0, dim=c(length(k), length(fullDat[[1]]$cHat), 1))
  }

  for(loopNum in 1:length(k)){

    normPhiSq[loopNum,] <- colSums(fullDat[[loopNum]]$PHI^2)

    #p = 1 th for loop iteration
    normcHatWOutZeroSq[loopNum,1,] <-  fullDat[[loopNum]]$cHat[1,]^2
    cHat[loopNum,1,] <- fullDat[[loopNum]]$cHat[1,]^2
    # F3[1,] <- (fStuff$cHat[1,])^2/
    #   ((normPhiSq - normcHatWOutZeroSq[k,1,])/(k - 1))
    for(P in 2:nrow(fullDat[[loopNum]]$cHat)){ # this is 1:p as we are removing zero so P-1 is actually P
      cHat[loopNum,P,] <- fullDat[[loopNum]]$cHat[P,]^2
      normcHatWOutZeroSq[loopNum,P,] <- normcHatWOutZeroSq[loopNum,(P-1),] + fullDat[[loopNum]]$cHat[P,]^2
      # F3[P,] <- (fStuff$cHat[P,])^2/
      #   ((normPhiSq - normcHatWOutZeroSq[k,P,])/(k - P))

    }
  }
  Freq <- fullDat[[1]]$Freq
  F3 <- F1 <- matrix(nrow = nrow(fullDat[[1]]$cHat), ncol = length(Freq))
  colnames(F3) <- Freq

  for(P in 1:nrow(fullDat[[1]]$cHat)){
    F1[P,] <- (colSums(normcHatWOutZeroSq[,P,])/sum(rep(P, times = length(k))))/
      (((colSums(normPhiSq - normcHatWOutZeroSq[,P,])))/(sum((k - rep(P, times=length(k))))))
    F3[P,] <- (((colSums(cHat[,P,]^2)))/(length(k)))/
      (((colSums(normPhiSq - normcHatWOutZeroSq[,P,])))/(sum((k - rep(P, times=length(k))))))
  }
  #making the return
  return(list(F4testStat = F3, Freq = Freq,
              F14testStat = F1 ))
}

#' F3TestParallel
#'
#' w is chosen by shannons number based on k
#'
#' @param xt time series
#' @param N = length(xt)Total number of observations
#' @param p Highest degree polynomial you want to test for
#' @param deltat Time interval between each observation
#' @param dpss  = FALSE unless you want to use dpss, it will do sine tapers by default
#' @param undersampleNumber A numeric of the number the user wants to undersample, usually 100 is a good start
#' @param k vector of tapers used in the f test
#' @param cores must be 1 if on windows, number of cores used for parallelization
#' @param confLevel default is 1-1/N, level of confidence used in the Ftest
#'
#' @return $F3testStat, $Freq, $significantFrequencies
#'
#' @export
F3Testpar <- function(xt, k, p, N = length(xt), deltat = 1, dpss = FALSE, undersampleNumber = 100, cores = 1,
                      confLevel = (1-(1/length(xt)))){

  if(is.null(undersampleNumber)){
    stop("need to set undersample amount")
  }

  if(dpss){
    fullDat <- parallel::mclapply(X = k,FUN = function(x){
      return(singleIterationForParallel(xt = xt, k = x, w = ((x+1)/(2*length(xt))), p = p, deltat = deltat,
                                        undersampleNumber = undersampleNumber, dpss = TRUE,
                                        confLevel = (1-(1/length(xt)))))
    }, mc.cores = cores, mc.cleanup = TRUE)
  }else{
    fullDat <- parallel::mclapply(X = k,FUN = function(x){
      return(singleIterationForParallel(xt = xt, k = x, p = p, deltat = deltat,
                                        undersampleNumber = undersampleNumber, dpss = FALSE,
                                        confLevel = (1-(1/length(xt)))))
    }, mc.cores = cores, mc.cleanup = TRUE)
  }
  Freq = fullDat[[1]]$Freq

  significantFrequencies<- matrix(0,nrow = p, ncol = length(Freq))
  for(i in 1:length(k)){
    for(j in 1:p){
      if(length(as.vector(fullDat[[i]]$significantFreq[[j]])) == 0){

      }else{
        #apparently you have to check each sig freq separatly
        indexesOfSigFreq <- match(fullDat[[i]]$significantFreq[[j]], Freq)
        significantFrequencies[j,indexesOfSigFreq] =
          significantFrequencies[j,indexesOfSigFreq] + 1
      }

    }

  }
  prop <- significantFrequencies/length(k)



  #making the return
  return(list(F3TestStat = fullDat, Freq = Freq, sigFreq = significantFrequencies,
              proportionSig = prop))
}
