mhtplot <- function(data, usepos=FALSE, logscale=TRUE, base=10, cutoffs=c(4,6,8), colors=NULL, labels=NULL, gap=NULL, ...)
{
  chr <- data[,1]
  pos <- newpos <- data[,2]
  p <- data[,3]
  tablechr <- table(chr)
  allchr <- as.vector(tablechr)
  n.chr <- length(allchr)
  colorlist <- colors()
  if(is.null(colors)) colors <- sample(colorlist,n.chr)
  if(is.null(labels)) labels <- names(tablechr)
  if(is.null(gap)) gap <- 1
  CMindex <- cumsum(allchr)
  for(i in 1:n.chr)
  {
     u <- CMindex[i]
     l <- CMindex[i]-allchr[i]+1
     chr <- l:u
     if (usepos) d <- diff(pos[chr])
     else d <- rep(1,allchr[i])
     newpos[chr] <- c(gap,d)
  }
  CM <- cumsum(as.numeric(newpos))
  dp <- seq(base^min(-cutoffs),1,length=sum(allchr))
  if (logscale) y <- -log(dp,base)
  else y <- dp
  par(xaxt="n",yaxt="n")
  plot(CM,y,type="n",xlab="",ylab="",axes=FALSE,...)
  axis(1,tick=FALSE)
  axis(2,tick=FALSE)
  par(xaxt="s",yaxt="s")
  for(i in 1:n.chr)
  {
     u <- CMindex[i]
     l <- CMindex[i]-allchr[i]+1
     chr <- l:u
     cat("Plotting points ",l,"-",u,"\n");
     if (logscale) y <- -log(p[chr],base)
     else y <- p[chr]
     points(CM[chr],y,col=colors[i],...)
     axis(1,at=ifelse(i==1,0,CM[l]),labels=labels[i],...)
  }
  abline(h=cutoffs)
  mtext(paste(cutoffs," "),2,at=cutoffs)
  mtext(paste("-log",base,"(Observed value)",sep=""),2,line=2.5,las=0)
  mtext("Chromosome",1,line=2.5,las=0)
}

#20-2-2009 MRC-Epid JHZ
