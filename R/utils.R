a2g <- function(a1,a2)
{
  i <- ifelse(a1 < a2,a2,a1)
  j <- ifelse(a1 < a2,a1,a2)
  genocode <- ifelse (j==0, 0, i*(i-1)/2 + j)
  invisible (genocode)

}

g2a.c <- function (s)
{
    d <- 1 + 8 * (s - 1)
    u <- 1 + ((1 + (sqrt(d) - 1) - 1) / 2)
    u <- ceiling(u)
    l = s - u * (u - 1) / 2
    list (l=l,u=u)
}

g2a <- function (x)
{
    i <- 1 + floor((sqrt(8 * x + 1) - 1)/2)
    j <- x - i * (i - 1)/2
    i <- ifelse(j == 0, i - 1, i)
    j <- ifelse(j == 0, i, j)
    return(cbind(j,i))
}

is.miss <- function(data,data.int,miss.val=0)
{
   if (data.int==2) # genotype
   {
      id <- array(FALSE, length(data))
      for (i in 1:length(miss.val))
          id <- id | data==miss.val[i]
   }
   else # allele
   {
      id <- array(FALSE, length(data[,1]))  
      for (i in 1:length(miss.val))
          id <- id | apply(data==miss.val[i],1,any)
   }
   return (id)  
}

# adapted from gcp.c on 24/9/2004
# 17-6-2004 JH Zhao

revhap <- function(loci,hapid)
{
   nloci <- length(loci)
   nalleles <- vector("numeric",nloci)
   nalleles[nloci] <- 1
   m <- nloci
   for(k in nloci:1)
   {
      m <- m - 1
      nalleles[m] <- nalleles[m+1] * loci[k]
   }
  n <- length(hapid)
  hap <- matrix(0,n,nloci)
  for (i in 1:n)
  {
    l <- hapid[i] - 1
    for (j in 1:nloci)
    {
      hap[i,j] <- floor(l/nalleles[j])
      if (j==nloci) hap[i,j] <- l
      else l <- l%%nalleles[j]
      hap[i,j] <- hap[i,j] + 1
    }
  }
  invisible(hap)

}

revhap.i <- function(loci,hapid)
{
   nloci <- length(loci)
   nalleles <- vector("numeric",nloci)
   nalleles[nloci] <- 1
   m <- nloci
   for(k in nloci:1)
   {
      m <- m - 1
      nalleles[m] <- nalleles[m+1] * loci[k]
   }

  hap <- vector("numeric",nloci)
  l <- hapid - 1
  for (j in 1:nloci)
  {
    hap[j] <- floor(l/nalleles[j])
    if (j==nloci) hap[j] <- l
    else l <- l%%nalleles[j]
    hap[j] <- hap[j] + 1
  }
  invisible(hap)

}

gcode <- function(a1,a2) {

  i <- ifelse(a1 < a2,a2,a1)
  j <- ifelse(a1 < a2,a1,a2)
  genocode <- i*(i-1)/2 + j
  return(genocode)

}

ungcode <- function(x) {

  i <- 1 + floor((sqrt(8*x+1)-1)/2)
  j <- x - i*(i-1)/2
  
  # the following 2 lines were added as a patch to make this ungcode work:
  i <- ifelse(j==0,i-1,i)
  j <- ifelse(j==0,i,j)

  return(cbind(j,i))
}

grec2g <- function (h, n, t)
{
  hh <- h
  for (i in 1:n)
  {
    hh[,i] <- t$alist[[i]]$allele[h[,i]]
  }
  invisible(hh)
}

# a simple scheme to represent SNPs
# similar to a2g() and not unlike Mike Weale's Matlab function
# 10-2-2005 JHZ
m2plem <- function(a1,a2)
{
  a <- vector("numeric")
# 1=A, 2=a
  for (i in 1:length(a1))
  {
      if(a1[i]==1&a2[i]==2) {
         a[i] <- 0
      }
      else if (a1[i]==1&a2[i]==1) {
         a[i] <- 1
      }
      else if (a1[i]==2&a2[i]==2) {
         a[i] <- 2
      }
  }
  invisible(a)
}

plem2m <- function(a)
{
  a1 <- vector("numeric")
  a2 <- vector("numeric")
# 1=A, 2=a
  for (i in 1:length(a))
  {
      if(a[i]==0) {
         a1[i] <- 1
         a2[i] <- 2
      }
      else if (a[i]==1) {
         a1[i] <- 1
         a2[i] <- 1
      }
      else if (a[i]==2) {
         a1[i] <- 2
         a2[i] <- 2
      }
  }
  invisible(list(a1,a2))
}

micombine <- function (est, std.err, confidence = 0.95)
{
    qstar <- est[[1]]
    for (i in 2:length(est)) {
        qstar <- cbind(qstar, est[[i]])
    }
    qbar <- apply(qstar, 1, mean)
    u <- std.err[[1]]
    for (i in 2:length(std.err)) {
        u <- cbind(u, std.err[[i]])
    }
    if (!is.null(dimnames(qstar)[[1]]))
        dimnames(u)[[1]] <- dimnames(qstar)[[1]]
    u <- u^2
    ubar <- apply(u, 1, mean)
    bm <- apply(qstar, 1, var)
    m <- dim(qstar)[2]
    tm <- ubar + ((1 + (1/m)) * bm)
    rem <- (1 + (1/m)) * bm/ubar
    nu <- (m - 1) * (1 + (1/rem))^2
    alpha <- 1 - (1 - confidence)/2
    low <- qbar - qt(alpha, nu) * sqrt(tm)
    up <- qbar + qt(alpha, nu) * sqrt(tm)
    pval <- 2 * (1 - pt(abs(qbar/sqrt(tm)), nu))
    fminf <- (rem + 2/(nu + 3))/(rem + 1)
    result <- list(est = qbar, std.err = sqrt(tm), df = nu, signif = pval,
        lower = low, upper = up, r = rem, fminf = fminf)
    result
}
