###########################################################################
###########################################################################
###########################################################################

dwt.2d <- function(x, wf, J=4, boundary="periodic")
{

  m <- dim(x)[1]
  storage.mode(m) <- "integer"
  n <- dim(x)[2]
  storage.mode(n) <- "integer"

  dict <- wave.filter(wf)
  L <- dict$length
  storage.mode(L) <- "integer"
  h <- dict$hpf
  storage.mode(h) <- "double"
  g <- dict$lpf
  storage.mode(g) <- "double"

  z <- matrix(0, m/2, n/2)
  storage.mode(z) <- "double"

  x.wt <- vector("list", 3*J+1)
  x.names <- NULL
  for(j in 1:J) {
    out <- .C("two_D_dwt", "Image"=as.double(x), "Rows"=m, "Cols"=n, 
                "filter.length"=L, "hpf"=h, "lpf"=g, "LL"=z, "LH"=z,
                "HL"=z, "HH"=z)[7:10]
    if(j < J) {
      index <- (3*j-2):(3*j)
      x.wt[index] <- out[-1]
      x.names <- c(x.names, sapply(names(out)[-1], paste, j, sep=""))
      x <- out[[1]]
      m <- dim(x)[1]
      storage.mode(m) <- "integer"
      n <- dim(x)[2]
      storage.mode(n) <- "integer"
      z <- matrix(0, m/2, n/2)
      storage.mode(z) <- "double"
    }
    else {
      index <- (3*j):(3*(j+1)) - 2
      x.wt[index] <- out[c(2:4,1)]
      x.names <- c(x.names, sapply(names(out)[c(2:4,1)], paste, j, sep=""))
    }
  }

  names(x.wt) <- x.names
  attr(x.wt, "J") <- J
  attr(x.wt, "wavelet") <- wf
  attr(x.wt, "boundary") <- boundary
  x.wt
}

###########################################################################
###########################################################################
###########################################################################

idwt.2d <- function(y)
{
  J <- attributes(y)$J

  dict <- wave.filter(attributes(y)$wavelet)
  L <- dict$length
  storage.mode(L) <- "integer"
  h <- dict$hpf
  storage.mode(h) <- "double"
  g <- dict$lpf
  storage.mode(g) <- "double"

  LL <- paste("LL", J, sep="")
  y.in <- y[[LL]]

  for(j in J:1) {
    LH <- paste("LH", j, sep="")
    HL <- paste("HL", j, sep="")
    HH <- paste("HH", j, sep="")
    
    m <- dim(y.in)[1]
    storage.mode(m) <- "integer"
    n <- dim(y.in)[2]
    storage.mode(n) <- "integer"
    x <- matrix(0, 2*m, 2*n)
    storage.mode(x) <- "double"

    out <- .C("two_D_idwt", as.double(y.in), as.double(y[[LH]]),
              as.double(y[[HL]]), as.double(y[[HH]]), m, n, L, h, g,
              "Y"=x)
    y.in <- out$Y
  }
  zapsmall(y.in)
}

###########################################################################
###########################################################################
###########################################################################

modwt.2d <- function(x, wf, J=4, boundary="periodic")
{
  m <- dim(x)[1]
  storage.mode(m) <- "integer"
  n <- dim(x)[2]
  storage.mode(n) <- "integer"

  dict <- wave.filter(wf)
  L <- dict$length
  storage.mode(L) <- "integer"
  h <- dict$hpf / sqrt(2)
  storage.mode(h) <- "double"
  g <- dict$lpf / sqrt(2)
  storage.mode(g) <- "double"

  z <- matrix(0, m, n)
  storage.mode(z) <- "double"

  x.wt <- vector("list", 3*J+1)
  x.names <- NULL
  for(j in 1:J) {
    out <- .C("two_D_modwt", "Image"=as.double(x), "Rows"=m, "Cols"=n,
              "Level"=j, "filter.length"=L, "hpf"=h, "lpf"=g, "LL"=z,
              "LH"=z, "HL"=z, "HH"=z)[8:11]
    if(j < J) {
      index <- (3*j-2):(3*j)
      x.wt[index] <- out[-1]
      x.names <- c(x.names, sapply(names(out)[-1], paste, j, sep=""))
      x <- out$LL
    }
    else {
      index <- (3*j):(3*(j+1)) - 2
      x.wt[index] <- out[c(2:4,1)]
      x.names <- c(x.names, sapply(names(out)[c(2:4,1)], paste, j, sep=""))
    }
  }

  names(x.wt) <- x.names
  attr(x.wt, "J") <- J
  attr(x.wt, "wavelet") <- wf
  attr(x.wt, "boundary") <- boundary
  x.wt
}

###########################################################################
###########################################################################
###########################################################################

imodwt.2d <- function(y)
{
  J <- attributes(y)$J

  dict <- wave.filter(attributes(y)$wavelet)
  L <- dict$length
  storage.mode(L) <- "integer"
  h <- dict$hpf / sqrt(2)
  storage.mode(h) <- "double"
  g <- dict$lpf / sqrt(2)
  storage.mode(g) <- "double"

  LL <- paste("LL", J, sep="")
  y.in <- y[[LL]]

  for(j in J:1) {
    LH <- paste("LH", j, sep="")
    HL <- paste("HL", j, sep="")
    HH <- paste("HH", j, sep="")
    
    m <- dim(y.in)[1]
    storage.mode(m) <- "integer"
    n <- dim(y.in)[2]
    storage.mode(n) <- "integer"
    x <- matrix(0, m, n)
    storage.mode(x) <- "double"

    out <- .C("two_D_imodwt", as.double(y.in), as.double(y[[LH]]),
              as.double(y[[HL]]), as.double(y[[HH]]), m, n, j, L,
              h, g, "Y"=x)
    y.in <- out$Y
  }
  zapsmall(y.in)
}

###########################################################################
###########################################################################
###########################################################################

plot.dwt.2d <- function(x, cex.axis=1, plot=TRUE, ...)
{
  J <- attributes(x)$J
  X <- x[[paste("LL", J, sep="")]]
  for(j in J:1) {
    x.names <- sapply(c("LH","HL","HH"), paste, j, sep="")
    X <- rbind(cbind(X, x[[x.names[2]]]),
               cbind(x[[x.names[1]]], x[[x.names[3]]]))
  }
  M <- dim(X)[1]; N <- dim(X)[2]
  if(plot) {
    image(1:M, 1:N, X, col=rainbow(128), axes=FALSE, xlab="", ylab="", ...)
    x.label <- NULL
    lines(c(0,N,N,0,0) + 0.5, c(0,0,M,M,0) + 0.5)
    for(j in J:1) {
      lines(c(M/2^j,M/2^j) + 0.5, 2*c(0,N/2^j) + 0.5)
      lines(2*c(0,M/2^j) + 0.5, c(N/2^j,N/2^j) + 0.5)
    }
    at <- c((3*N+2)/2^(1:J+1),(N+2)/2^(J+1))
    labs <- c(paste("H",1:J,sep=""), paste("L",J,sep=""))
    axis(side=1, at=at, label=labs, tick=FALSE, cex.axis=cex.axis)
    axis(side=2, at=at, label=labs, tick=FALSE, cex.axis=cex.axis)
  }
  else
    return(X)
  invisible()
}

###########################################################################
###########################################################################
###########################################################################

denoise.dwt.2d <- function(x, wf = "la8", J = 4, method = "universal", 
                           H = 0.5, noise.dir = 3, rule = "hard")
{
  soft <- function(x, delta) sign(x) * pmax(abs(x) - delta, 0)
  hard <- function(x, delta) ifelse(abs(x) > delta, x, 0)

  n <- length(x)
  x.dwt <- dwt.2d(x, wf, J)
  if(noise.dir == 3)
    sigma.mad <- list(HH = mad(x.dwt$HH1), HL = mad(x.dwt$HL1), 
  		      LH = mad(x.dwt$LH1))
  else {
    noise <- x.dwt$jj
    sigma.mad <- list(HH = mad(noise), HL = mad(noise), LH = mad(noise))
  }
    thresh <- list(HH = rep(sqrt(2 * sigma.mad$HH^2 * log(n)), J), 
		   HL = rep(sqrt(2 * sigma.mad$HL^2 * log(n)), J),
		   LH = rep(sqrt(2 * sigma.mad$LH^2 * log(n)), J))

  if(method == "long-memory")
    thresh <- lapply(thresh, function(x,J,H) 2^(0:(J-1)*(H-1/2))*x, J=J, H=H)
  for(j in 1:J) {
    jj <- paste("HL", j, sep = "")
    if(rule == "hard")
      x.dwt[[jj]] <- hard(x.dwt[[jj]], thresh$HL[j])
    else 
      x.dwt[[jj]] <- soft(x.dwt[[jj]], thresh$HL[j])
    jj <- paste("LH", j, sep = "")
    if(rule == "hard")
      x.dwt[[jj]] <- hard(x.dwt[[jj]], thresh$LH[j])
    else 
      x.dwt[[jj]] <- soft(x.dwt[[jj]], thresh$LH[j])
    jj <- paste("HH", j, sep = "")
    if(rule == "hard")
      x.dwt[[jj]] <- hard(x.dwt[[jj]], thresh$HH[j])
    else 
      x.dwt[[jj]] <- soft(x.dwt[[jj]], thresh$HH[j])
  }
  idwt.2d(x.dwt)
}

###########################################################################
###########################################################################
###########################################################################

denoise.modwt.2d <- function(x, wf = "la8", J = 4, method = "universal", 
  H = 0.5, rule = "hard")
{
  soft <- function(x, delta) sign(x) * pmax(abs(x) - delta, 0)
  hard <- function(x, delta) ifelse(abs(x) > delta, x, 0)
  n <- length(x)
  x.modwt <- modwt.2d(x, wf, J)
  sigma.mad <- list(HH = sqrt(2) * mad(x.modwt$HH1),
		    HL = sqrt(2) * mad(x.modwt$HL1),
		    LH = sqrt(2) * mad(x.modwt$LH1))
    thresh <- list(HH = rep(sqrt(2 * sigma.mad$HH^2 * log(n))/2^(1:J), J), 
		   HL = rep(sqrt(2 * sigma.mad$HL^2 * log(n))/2^(1:J), J), 
		   LH = rep(sqrt(2 * sigma.mad$LH^2 * log(n))/2^(1:J), J))
  if(method == "long-memory")
    thresh <- lapply(thresh, function(x,J,H) 2^(0:(J-1)*(H-1/2))*x, J=J, H=H)
  for(j in 1:J) {
    jj <- paste("HL", j, sep = "")
    if(rule == "hard")
      x.modwt[[jj]] <- hard(x.modwt[[jj]], thresh$HL[j])
    else 
      x.modwt[[jj]] <- soft(x.modwt[[jj]], thresh$HL[j])
    jj <- paste("LH", j, sep = "")
    if(rule == "hard")
      x.modwt[[jj]] <- hard(x.modwt[[jj]], thresh$LH[j])
    else 
      x.modwt[[jj]] <- soft(x.modwt[[jj]], thresh$LH[j])
    jj <- paste("HH", j, sep = "")
    if(rule == "hard")
      x.modwt[[jj]] <- hard(x.modwt[[jj]], thresh$HH[j])
    else 
     x.modwt[[jj]] <- soft(x.modwt[[jj]], thresh$HH[j])
  }
  imodwt.2d(x.modwt)
}
