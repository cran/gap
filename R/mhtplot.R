mht.control <- function(usepos=FALSE, logscale=TRUE, base=10, cutoffs=NULL, colors=NULL, labels=NULL, gap=NULL, cex=0.4, yline=3, xline=3)
{
   list(usepos=usepos, logscale=logscale, base=base, cutoffs=cutoffs, colors=colors, labels=labels, gap=gap, cex=cex, yline=yline, xline=3)
}

hmht.control <- function(data=NULL, colors=NULL, yoffset=0.25, cex=1.5, boxed=FALSE)
{
   list(data=data,colors=colors,yoffset=yoffset,cex=cex,boxed=boxed)
}

mhtplot <- function(data, control=mht.control(), hcontrol=hmht.control(), ...)
{
  require(grid)
  data2 <- data[!apply(is.na(data),1,any),]
  chr <- data2[,1]
  pos <- newpos <- data2[,2]
  p <- data2[,3]
  tablechr <- table(chr)
  allchr <- as.vector(tablechr)
  n.chr <- length(allchr)
  usepos <- control$usepos
  logscale <- control$logscale
  base <- control$base
  cutoffs <- control$cutoffs
  colors <- control$colors
  labels <- control$labels
  gap <- control$gap
  pcex <- control$cex
  yline <- control$yline
  xline <- control$xline
  colorlist <- colors()
  if(is.null(colors)) colors <- sample(colorlist,n.chr)
  if(is.null(labels)) labels <- names(tablechr)
  if(is.null(gap)) gap <- 0
  if (!is.null(hcontrol$data))
  {
     hdata <- hcontrol$data
     hdata2 <- hdata[!apply(is.na(hdata),1,any),]
     hchr <- hdata2[,1]
     hpos <- hnewpos <- hdata2[,2]
     hp <- hdata2[,3]
     hname <- hdata2[,4]
     hyoffs <- hcontrol$yoffs
     hboxed <- hcontrol$boxed
     hcex <- hcontrol$cex
     htablechr <- table(hchr)
     hallchr <- as.vector(htablechr)
     hn.chr <- length(hallchr)
     htablename <- table(hname)
     htablechrname <- table(hchr,hname)
     if(is.null(colors)) hcolors <- rep("red",length(htablename))
     else hcolors <- hcontrol$colors
  }
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
  args <- list(...)
  if ("ylim"%in%names(args)) dp <- seq(args$ylim[1],args$ylim[2],length=sum(allchr))
  else dp <- seq(min(p),max(p),length=sum(allchr))
  if (logscale) y <- -log(dp,base) else y <- dp
  y1 <- min(y)
  par(xaxt="n",yaxt="n")
  plot(CM,y,type="n",ann=FALSE,axes=FALSE,...)
  axis(1)
  axis(2)
  par(xaxt="s",yaxt="s")
  for(i in 1:n.chr)
  {
     u <- CMindex[i]
     l <- CMindex[i]-allchr[i]+1
     chr <- l:u
     cat("Plotting points ",l,"-",u,"\n");
     if (logscale) y <- -log(p[chr],base) else y <- p[chr]
     col.chr <- colors[i]
     points(CM[chr],y,col=col.chr,cex=pcex,...)
     text(ifelse(i==1,CM[1],CM[l]),y1,pos=1,offset=1.5,labels[i],srt=45,...)
  }
  j <- 1
  for(i in 1:n.chr)
  {
     u <- CMindex[i]
     l <- CMindex[i]-allchr[i]+1
     chr <- l:u
     if (logscale) y <- -log(p[chr],base) else y <- p[chr]
     col.chr <- colors[i]
     if (!is.null(hcontrol$data))
     {
        chrs <- names(htablechr)[as.character(i)==names(htablechr)]
        if (length(chrs)>=1)
        {
           hchrname <- colnames(htablechrname)[htablechrname[j,]>0]
           for(k in 1:length(hchrname))
           {
              hregion <- hpos[hchr==i&hname==hchrname[k]]
              hl <- chr[pos[chr]==hregion[1]]
              hu <- chr[pos[chr]==hregion[length(hregion)]]
              cat("   ... highlighting",hl,"-",hu,hchrname[k],"\n")
              l1 <- hl-l+1
              l2 <- hu-l+1
              col.chr[l1:l2] <- hcolors[j]
              if (hboxed)
              {
                 tg <- textGrob(hchrname[k])
                 rg <- rectGrob(x=CM[chr][l1],y=max(y[l1:l2])+hyoffs,
                                width=1.1*grobWidth(tg),height=1.3*grobHeight(tg),
                                gp=gpar(col="black",lwd=2.5))
                 boxedText <- gTree(children=gList(tg,rg))
                 grid.draw(boxedText)
              } else text(CM[chr][l1],max(y[l1:l2]+hyoffs),hchrname[k],cex=hcex)
              points(CM[l+(l1:l2)],y[l1:l2],col=col.chr[l1:l2],cex=pcex,...)
           }
           j <- j+1
        }
     }
  }
  if(!is.null(cutoffs)) abline(h=cutoffs)
  mtext(ifelse(logscale,paste("-log",base,"(Observed value)",sep=""),"Observed value"),2,line=yline,las=0)
  if ("xlab"%in%names(args)) xlabel <- xlab else xlabel <- ifelse(is.null(names(chr)),"Chromosome",names(chr))
  mtext(xlabel,1,line=xline,las=0)
}

#17-11-2010 MRC-Epid JHZ
