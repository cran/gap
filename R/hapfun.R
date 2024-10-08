#' `HapDesign` and `HapFreqSE` both accept a [haplo.stats::haplo.em] object to derieve a design/dosage
#' matrix and standard error of haplotype frequency estimates. The former is appropriate for haplotype trend
#' regression (HTR), e.g., within the generalized linear model (GLM) framework to be equivllant to a formal
#' approach as implemented in the package haplo.stats and hap.score. However, they are expected to be compatible
#' with objects from [`gc.em`] and [`hap.em`]. The two functions are
#' included as courtesy of Prof Andrea Foulkes from the useR!2008 tutorial.
#' @noRd

##########################################################################
# Description: This function creates a design matrix with i,j
# element equal to the conditional expectation
# of the number of copies of haplotype j for
# individual i based on the output from haplo.em()
# Input: HaploEM (object resulting from haplo.em()), Output: XmatHap
##########################################################################
HapDesign <- function(HaploEM)
{
   Nobs <- length(unique(HaploEM$indx.subj))
   Nhap <- length(HaploEM$hap.prob)
   XmatHap <- matrix(data=0,nrow=Nobs,ncol=Nhap)
   for (i in 1:Nobs)
   {
       IDSeq <- seq(1:sum(HaploEM$nreps))[HaploEM$indx.subj==i]
       for (j in 1:length(IDSeq))
       {
          XmatHap[i,HaploEM$hap1code[IDSeq][j]] <- XmatHap[i,HaploEM$hap1code[IDSeq][j]] + HaploEM$post[IDSeq][j]
          XmatHap[i,HaploEM$hap2code[IDSeq][j]] <- XmatHap[i,HaploEM$hap2code[IDSeq][j]] + HaploEM$post[IDSeq][j]
       }
    }
    return(XmatHap)
}

##########################################################################
# Description: This function creates a vector with jth element
# equal to the standard error of haplotype j
# based on the output from haplo.em()
# Input: HaploEM (object resulting from haplo.em()), Output: HapSE
##########################################################################
HapFreqSE <- function(HaploEM)
{
   HapMat <- HapDesign(HaploEM)
   Nobs <- length(unique(HaploEM$indx.subj))
   Nhap <- length(HaploEM$hap.prob)
   S.Full<-matrix(data=0, nrow=Nobs, ncol=Nhap-1)
   for(i in 1:Nobs)
   {
      for(k in 1:(Nhap-1))
      {
         S.Full[i,k] <- HapMat[i,k]/HaploEM$hap.prob[k]-HapMat[i,Nhap]/HaploEM$hap.prob[Nhap]
      }
   }
   Score <- t(S.Full)%*%S.Full
   invScore <- solve(Score)
   HapSE <- c(sqrt(diag(invScore)),sqrt(t(rep(1,Nhap-1))%*%invScore%*%rep(1,Nhap-1)))
   return(HapSE)
}
