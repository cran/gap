# 08-2-2004 Working but miss.value needs to be added later
# 09-2-2004 Add using genotype counts
# 11-2-2004 Works ok in all three modes
# 17-2-2004 Add is.miss function

hwe <- function(data, is.count=FALSE, is.genotype=FALSE, yates.correct=FALSE, miss.val=0)
{
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
      list(l=j,u=i)
  }
  is.miss <- function(data,is.genotype,miss.val=0)
  {
     if (is.genotype)
     {
        id <- array(FALSE, length(data))
        for (i in 1:length(miss.val))
            id <- id | data==miss.val[i]
     }
     else
     {
        id <- array(FALSE, length(data[,1]))  
        for (i in 1:length(miss.val))
            id <- id | apply(data==miss.val[i],1,any)
     }
     return (id)  
  }
  if (!is.count)
  {
     if (!is.genotype)
     {
        data <- data[!is.miss(data,is.genotype),]
        n.obs <- length(data[,1])
        genotype <- array(0,n.obs)
        n.allele <- max(data)
        a1 <- data[,1]
        a2 <- data[,2]
        for (i in 1:n.obs)
        {
            l <- min(a1[i],a2[i])
            u <- max(a1[i],a2[i])
            genotype[i] <- l + u * (u - 1) / 2
            geno.table <- table(genotype)
        }
     }
     else
     {
        data <- data[!is.miss(data,is.genotype)]
        geno.table <- table(data)
        n.obs <- length(data)
        a1 <- a2 <- array(0,n.obs)
        genotype <- data
        s <- max(data)
        z <- g2a(s)
        n.allele <- z$u
        for (i in 1:n.obs)
        {
            z <- g2a(data[i])
            a1[i] <- z$l
            a2[i] <- z$u
        }
     }
     n.obs2 <- 2 * n.obs
     allele.freq <- table(c(a1,a2)) / n.obs2
     n.geno <- length(geno.table)
     geno.name <- as.numeric(names(geno.table))

     x2 <- lrt <- 0
     n.genotype <- n.allele * (n.allele + 1) / 2
     all.genotype <- array(0,n.genotype)

     for (i in 1:n.geno)
     {
         j <- geno.name[i]
         all.genotype[j] <- geno.table[i]
     }
     to.warn <- FALSE
     for (i in 1:n.genotype)
     {
         o <- all.genotype[i]
         z <- g2a(i)
         a1 <- z$l
         a2 <- z$u
         e <- ifelse(a1==a2,1,2) * allele.freq[a1] * allele.freq[a2] * n.obs
         if (e < 0.5) to.warn <- TRUE
         if (yates.correct)
            x2 <- x2 + (abs (o - e) - 0.5)^2 / e
         else
            x2 <- x2 + (o - e)^2 / e
         if (o>0) lrt <- lrt + o * log(o / e)
     }
  }
  else
  {
      n.genotype <- length(data)
      n.obs <- sum(data)
      e <- array(0,n.genotype)
      n.allele <- (sqrt(1 + 8 * n.genotype) - 1) / 2
      allele.freq <- array(0,n.allele)
      to.warn <- FALSE
      k <- 0
      for (i in 1:n.allele)
      {
          for (j in 1:i)
          {
              k <- k + 1
              allele.freq[i] <- allele.freq[i] + data[k]
              allele.freq[j] <- allele.freq[j] + data[k]
              if(e[k] < 0.5) to.warn <- TRUE
          }
      }
      for (i in 1:n.allele) allele.freq[i] <- allele.freq[i] / n.obs / 2
      k <- 0
      for (i in 1:n.allele)
          for (j in 1:i)
          {
              k <- k + 1
              e[k] <- ifelse(i==j,1,2) * allele.freq[i] * allele.freq[j] * n.obs
          }
      x2 <- lrt <- 0
      for (i in 1:n.genotype)
      {
          o <- data[i]
          x2 <- x2 + (o - e[i])^2 / e[i]
          if (yates.correct) x2 <- x2 + (abs(o - e[i]) - 0.5)^2 / e[i]
          if (o > 0) lrt <- lrt + o * log(o / e[i])
      }
  }
  if (to.warn) cat("there is at least one cell with expected value < 0.5\n")
  df <- n.genotype - n.allele
  rho <- x2 / n.obs
  cat("Pearson x2=",x2,"df=",df,"p=",1-pchisq(x2,df),sep="\t")
  cat("\n")
  list (allele.freq=allele.freq,x2=x2, p.x2=1-pchisq(x2,df),
        lrt=lrt, p.lrt=1-pchisq(x2,df), df=df, rho=rho)
}
