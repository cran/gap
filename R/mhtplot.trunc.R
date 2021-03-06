mhtplot.trunc <- function (x, chr = "CHR", bp = "BP", p = NULL, log10p = NULL, z = NULL, snp = "SNP",
                           col = c("gray10", "gray60"),
                           chrlabs = NULL, suggestiveline = -log10(1e-05),
                           genomewideline = -log10(5e-08), highlight = NULL,
                           annotatelog10P = NULL, annotateTop = FALSE, cex.mtext=1.5, cex.text=0.7,
                           mtext.line = 2, cex.y = 1, y.ax.space = 5, y.brk1, y.brk2, delta=0.05, ...)
{
  for (q in c("calibrate","plotrix")) {
     if (length(grep(paste("^package:", q, "$", sep=""), search())) == 0) {
        if (!requireNamespace(q, quietly = TRUE))
        warning(paste("mhtplot.trunc needs package `", q, "' to be fully functional; please install", sep=""))
     }
  }
  CHR <- BP <- BP.x <- BP.y <- SNP <- log10P <- index <- NULL
  if (!(chr %in% names(x))) stop(paste("Column", chr, "not found!"))
  if (!is.numeric(x[[chr]])) 
     stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
  if (!(bp %in% names(x)))  stop(paste("Column", bp, "not found!"))
  if (!is.null(p)) log10P <- -log10pvalue(x[[p]])
  if (is.null(p) & !is.null(log10p)) log10P <- -as.numeric(x[[log10p]])
  if (is.null(p) & is.null(log10p) & !is.null(z)) log10P <- -log10p(as.numeric(x[[z]]))
  if(is.null(p) & is.null(log10p) & is.null(z)) stop("At least one of p, log10p, or z (priority given in that order) is needed")
  if (y.brk2 <= y.brk1) stop("y.brk2 must be larger than y.brk1")
  if (!(snp %in% names(x))) warning(paste("No SNP column found. OK unless you're trying to highlight."))
  d <- data.frame(CHR = x[[chr]], BP = as.integer(x[[bp]]), log10P = log10P)
  d <- subset(d, !is.na(CHR) & !is.na(BP) & !is.na(log10P))
  if (!is.null(x[[snp]])) d <- transform(d, SNP = x[[snp]])
  d.order <- with(d,order(CHR, BP))
  d <- within(d[d.order,], {pos <- NA; index <- NA})
  ind <- 0
  for (i in unique(with(d,CHR))) {
    ind <- ind + 1
    d[with(d,CHR) == i, "index"] <- ind
  }
  nchr <- length(unique(with(d,CHR)))
  if (nchr == 1) {
    d["pos"] <- d["BP"]
    ticks <- floor(length(d["pos"]))/2 + 1
    xlabel <- paste("Chromosome", unique(with(d,CHR)), "position")
    labs <- ticks
  } else {
    lastbase <- 0
    ticks <- NULL
    for (i in unique(with(d,index))) {
      if (i == 1) d[with(d,index) == i, "pos"] <- d[with(d,index) == i, "BP"]
      else {
        lastbase = lastbase + tail(with(subset(d, index == i - 1),BP), 1)
        d[with(d,index) == i, "pos"] = d[with(d,index) == i, "BP"] + lastbase
      }
      ticks <- c(ticks, (min(d[with(d,index) == i, "pos"]) + max(d[with(d,index) == i, "pos"]))/2 + 1)
    }
    xlabel <- "Chromosome"
    labs <- unique(with(d,CHR))
  }
  xmax <- ceiling(max(with(d,pos)) * 1.03)
  xmin <- floor(max(with(d,pos)) * -0.03)
  max.y <- ceiling(max(with(d,log10P), na.rm=TRUE))
  if (y.brk2 > max.y ){
      message(paste("max.y is", max.y))
      stop("User error: Upper breakpoint must be lower than maximum -log10(P-value)")
  }
  offset <- y.brk2-y.brk1
  d <- within(d, {
    gapped <- log10P > y.brk1 & log10P < y.brk2
    above <- log10P > y.brk2
    log10P[gapped] <- NA
    log10P[above] <- log10P[above] - offset
  })
  def_args <- list(xaxt = "n", yaxt="n", bty = "n", xaxs = "i",
                   las = 1, pch = 20, xlim = c(xmin, xmax),
                   ylim = c(0, ceiling(max(with(d,log10P), na.rm=TRUE))),
                   xlab = "", ylab = "")
  dotargs <- list(...)
  do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))
  mtext(text = xlabel, side = 1, line = mtext.line, cex = cex.mtext, font=2)
  mtext(text = expression(-log[10](italic(p))), side=2, line = mtext.line, cex = cex.mtext, font=2)
  y.lab.tick.pos <- seq(from = 0, by = y.ax.space, to = ceiling(max.y) - offset + y.ax.space / 5)
  pre.brk.labs <- seq(from = 0, by = y.ax.space, to = y.brk1-y.ax.space)
  post.brk.labs <- seq(from = y.brk2, by=y.ax.space, to = max(y.lab.tick.pos))
  y.labels <- c(pre.brk.labs, seq(from=y.brk2, by=y.ax.space, length.out=length(y.lab.tick.pos)-length(pre.brk.labs)))
  axis(side=2, at=y.lab.tick.pos, labels=y.labels, cex.axis=cex.y, las=1)
  plotrix::axis.break(axis = 2, breakpos = y.brk1, style = "slash")
  if (!is.null(chrlabs)) {
    if (is.character(chrlabs)) {
       if (length(chrlabs) == length(labs)) labs <- chrlabs
       else warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
    }
    else warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
  }
  if (nchr == 1) axis(1, ...) else axis(1, at = ticks, labels = labs)
  col <- rep(col, max(with(d,CHR)))
  if (nchr == 1) with(d, points(pos, log10P, pch = 20, col = col[1], ...))
  else {
    icol = 1
    for (i in unique(with(d,index))) {
      with(d[with(d,index) == unique(with(d,index))[i], ], points(pos, log10P, col = col[icol], pch = 20, ...))
      icol = icol + 1
    }
  }
  if (suggestiveline) abline(h = suggestiveline, col = "blue")
  if (genomewideline) abline(h = genomewideline, col = "red")
  if (!is.null(highlight)) {
    if (any(!(highlight %in% with(d,SNP)))) warning("You're trying to highlight SNPs that don't exist in your results.")
    d.highlight = d[which(with(d,SNP) %in% highlight), ]
    with(d.highlight, points(pos, log10P, col = "red", pch = 20, ...))
    d.column <- subset(merge(d,d.highlight[c("CHR","BP")],by=c("CHR")),BP.x>(1-delta)*BP.y & BP.x<(1+delta)*BP.y)
    print(nrow(d.column))
    with(d.column,points(pos, log10P, col = "red", pch = 20, ...))
  }
  if (!is.null(annotatelog10P)) {
    topHits = subset(d, log10P >= annotatelog10P)
    if (!annotateTop) {
      with(subset(topHits,SNP %in% highlight),
           calibrate::textxy(pos, log10P, offset = 0.625, pos = 3, labs = SNP, cex = cex.text, font = 4), ...)
    }
    else {
      topHits <- topHits[order(with(topHits,log10P)), ]
      topSNPs <- NULL
      for (i in unique(with(topHits,CHR))) {
        chrSNPs <- topHits[with(topHits,CHR) == i, ]
        topSNPs <- rbind(topSNPs, chrSNPs[1, ])
      }
      with(topSNPs,calibrate::textxy(pos, log10P, offset = 0.625, pos = 3, labs = SNP, cex = cex.text, font = 4),...)
    }
  }
}
