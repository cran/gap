mhtplot <- 
function(data, usepos=FALSE, logscale=TRUE, base=10, cutoffs=c(4,6,8), colors=NULL, labels=NULL, xlab=NULL, gap=NULL, ...)
{
  data2 <- data[!apply(is.na(data),1,any),]
  chr <- data2[,1]
  pos <- newpos <- data2[,2]
  p <- data2[,3]
  tablechr <- table(chr)
  allchr <- as.vector(tablechr)
  n.chr <- length(allchr)
  colorlist <- colors()
  if(is.null(colors)) colors <- sample(colorlist,n.chr)
  if(is.null(labels)) labels <- names(tablechr)
  if(is.null(gap)) gap <- 0
  CMindex <- cumsum(allchr)
  for(i in 1:n.chr)
  {
     u <- CMindex[i]
     l <- CMindex[i]-allchr[i]+1
     chr <- l:u
     if (usepos) d <- diff(pos[chr]) else d <- rep(1,allchr[i]-1)
     newpos[chr] <- c(gap,d)
  }
  CM <- cumsum(as.numeric(newpos))
  dp <- seq(min(p),max(p),length=sum(allchr))
  if (logscale) y <- -log(dp,base) else y <- dp
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
     if (logscale) y <- -log(p[chr],base) else y <- p[chr]
     points(CM[chr],y,col=colors[i],...)
     axis(1,at=ifelse(i==1,CM[1],CM[l]),labels=labels[i],...)
  }
  abline(h=cutoffs)
  axis(2,at=cutoffs,lwd=0)
  mtext(ifelse(logscale,paste("-log",base,"(Observed value)",sep=""),"Observed value"),2,line=2.5,las=0)
  if (!is.null(xlab)) xlabel <- xlab else xlabel <- ifelse(is.null(names(chr)),"Chromosome",names(chr))
  mtext(xlabel,1,line=2.5,las=0)
}

#20-11-2009 MRC-Epid JHZ
