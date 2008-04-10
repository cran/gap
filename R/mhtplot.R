# 25-3-2008 MRC-Epid JHZ

mhtplot <- function(data,logscale=TRUE, base=10, cutoffs=c(3,5,7,9), color=NULL,labels=paste(1:22,sep=""),...)
{
  chr <- data[,1]
  pos <- data[,2]
  p <- data[,3]
  affy <- as.vector(table(chr))
  CM <- cumsum(affy)
  n.markers <- sum(affy)
  n.chr <- length(affy)
  id <- 1:n.chr
  eps <- .Machine$double.eps
  dp <- seq(eps,1,length=n.markers)
  if (logscale) y <- -log(dp,base)
  else y <- dp
  par(xaxt="n",yaxt="n")
  plot(pos,y,type="n",xlab="",ylab="",axes=FALSE,...)
  axis(1,tick=FALSE)
  axis(2,tick=FALSE)
  colorlist <- colors()
  if(is.null(color)) color <- sample(colorlist,22)
  for(i in 1:n.chr)
  {
     u <- CM[i]
     l <- CM[i]-affy[i]+1
     cat("Plotting points ",l,"-",u,"\n");
     chr <- l:u
     if (logscale) y <- -log(p[chr])
     else y <- p[chr]
     points(pos[chr],y,col=color[i],...)
  }
  par(xaxt="s",yaxt="s")
  axis(1,at=c(0,CM[-n.chr]),labels=labels)
  if (logscale)
  {
     abline(h=cutoffs)
     mtext(eval(expression(cutoffs)),2,at=cutoffs)
     # axis(2,at=cutoffs,labels=eval(expression(10^(-cutoffs))))
  }
}
