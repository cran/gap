# This function fits a Bradley-Terry Model to a squared contingency table
# JH Zhao 5/1/2004

bt <- function(x)
{
  n <- dim(x)[1]
  a <- - diag(n)
  a[,1] <- 1
  allele <- a
  for (i in 2:(n-1))
  {
    a <- - diag(n)
    a[,i] <- 1
    allele <- rbind(allele,a)
  }
  count <- rep(0,n*(n-1))
  y <- count + 1
  k <- 1
  for (i in 1:n)
  {
    for (j in 1:n) if(i!=j)
    {
      count[k] <- x[i,j]
      k <- k + 1
    }
  }
  # This generates format as required by ETDT
  # adapted from bt.sas (16/05/1999)
  toETDT <- function(a)
  {
     n <- dim(a)[1]
     nn1 <- n*(n-1)/2;
     C<-tr<-i1<-i2<-rep(0,nn1)
     ij <- matrix(rep(0,n*nn1),nrow=nn1)
     vi <- rep(0,n)
     l <- 1
     for(i in 1:(n-1))
     {
        for(j in (i+1):n)
        {
          C[l] <- a[i,j]+a[j,i]
          tr[l] <- a[i,j]
          i1[l] <- i
          i2[l] <- j
          ij[l,i] <- 1
          ij[l,j] <- -1
          l <- l + 1
        }
        vi[i] <- i
     }
     cbind(i1,i2,C,tr,ij)
  }

  bt.glm<-glm(y~-1+allele,weights=count,family="binomial")
  list(y=y,count=count,allele=allele,bt.glm=bt.glm,etdt.dat=toETDT(x))
}
