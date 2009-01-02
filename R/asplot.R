asplot <- function (snp, locusname, chr, locus, gmap, glist, best.pval=NULL, sf=c(3,8), logpmax=9, pch=23)
{
    hit <- locus[snp, ]
    if (is.null(best.pval)) best.pval <- hit$PVAL
    lb <- min(locus$POS) - 1e3
    ub <- max(locus$POS) + 1e3
    lu <- ub - lb
    center <- lb + lu / 2
    center.100kb <- round(center/1e5) * 1e5
    offset.100kb <- round((lu/3)/1e5) * 1e5
    center.three <- c(center.100kb - offset.100kb, center.100kb, center.100kb + offset.100kb)
    offset <- logpmax / sf[1]
    ylim <- logpmax + offset
    yadj <- -offset + ylim / sf[2]
    keep <- subset(gmap, gmap[, 1] > lb & gmap[, 1] < ub)
    genes <- subset(glist, (glist$START > lb & glist$START < ub) | (glist$STOP > lb & glist$STOP < ub))
    print(genes)
    par(mar = c(4, 4, 3, 4))
    plot(keep[, 1], yadj + ((keep[, 2]/60) * (3 * ylim / 4)), type = "l", col = "lightblue", lwd = 1, xlim = c(lb, ub), ylim = c(-offset, logpmax), xlab = "", ylab = "", main = locusname, axes = F)
    mtext(paste("Chromosome", chr, "position (kb)", sep = " "), side = 1, line = 2.5)
    axis(1, at = center.three, labels = center.three / 1e3, las = 1)
    axis(2, at = seq(0, logpmax, 2), labels = seq(0, logpmax, 2), las = 1)
    mtext("-log10(Observed p)", side = 2, at = logpmax / 2, line = 2)
    axis(4, at = c(yadj, yadj + ylim / 4, yadj + ylim / 2, yadj + 3 * ylim / 4), labels = c("0", "20", "40", "60"), las = 1)
    mtext("Recombination rate (cM/Mb)", side = 4, at = (-offset + ylim / 2), line = 2)
    box()
    lines(c(lb, ub), c(0, 0), lty = "dotted", lwd = 1, col = "black")
    points(hit$POS, -(log10(hit$PVAL)), pch = pch, cex = 2.5, bg = "red")
    text(hit$POS, -(log10(hit$PVAL)), labels = c(row.names(hit)), pos = 3, offset = 1)
    if (-(log10(best.pval)) < logpmax) {
        points(hit$POS, -(log10(best.pval)), pch = pch, cex = 2.5, bg = "blue")
        text(hit$POS, -(log10(best.pval)), labels = c(paste("P=", best.pval, sep = "")), pos = 4, offset = 2)
    } else {
        points(hit$POS, logpmax, pch = pch, cex = 2.5, bg = "blue")
        text(hit$POS, logpmax, labels = c(paste("P=", best.pval, sep = "")), pos = 4, offset = 1)
    }
    colors <- c("red","orange","yellow","white","grey")
    markers.in.strong.ld <- subset(locus, (row.names(locus) != snp & locus$RSQR >= 0.8 & locus$TYPE == "typed"))
    markers.in.moderate.ld <- subset(locus, (row.names(locus) != snp & locus$RSQR >= 0.5 & locus$RSQR < 0.8 & locus$TYPE == "typed"))
    markers.in.weak.ld <- subset(locus, (row.names(locus) != snp & locus$RSQR >= 0.2 & locus$RSQR < 0.5 & locus$TYPE == "typed"))
    markers.not.in.ld <- subset(locus, (row.names(locus) != snp & locus$RSQR < 0.2 & locus$TYPE == "typed"))
    points(markers.in.strong.ld$POS, -log10(markers.in.strong.ld$PVAL), pch = pch, cex = 1.25, bg = colors[1])
    points(markers.in.moderate.ld$POS, -log10(markers.in.moderate.ld$PVAL), pch = pch, cex = 1.25, bg = colors[2])
    points(markers.in.weak.ld$POS, -log10(markers.in.weak.ld$PVAL), pch = pch, cex = 1.25, bg = colors[3])
    points(markers.not.in.ld$POS, -log10(markers.not.in.ld$PVAL), pch = pch, cex = 1, bg = colors[4])
    imputed <- subset(locus, row.names(locus) != snp & locus$TYPE == "imputed")
    points(imputed$POS, -log10(imputed$PVAL), pch = pch, cex = 1, bg = colors[5])
    n.genes <- nrow(genes)
    for (i in 1:n.genes) {
        adj <- -offset+2*(i%%4)/3
        if (genes[i, ]$STRAND == "+") arrows(max(genes[i, ]$START, lb), adj, min(genes[i, ]$STOP, ub), adj, length = 0.05, lwd = 2, code = 2, lty = "solid", col = "darkgreen")
        else arrows(max(genes[i, ]$START, lb), adj, min(genes[i, ]$STOP, ub), adj, length = 0.05, lwd = 2, code = 1, lty = "solid", col = "darkgreen")
        if (!is.na(genes[i, ]$GENE)) text(genes[i, ]$START + (genes[i, ]$SIZE/2), adj +  ylim / 40, labels = genes[i, ]$GENE, cex = 0.7)
    }
    ltext <- rbind("0.8-","0.5-","0.2-","0.0-","imp.")
    legend("topleft",legend=ltext,title="LD (r^2)",fill=colors)
}
