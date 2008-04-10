asplot <- function (snp, locusname, chr, locus, range, best.pval, pch, recomb, genelist)
{
    hit <- locus[snp, ]
    min.pos <- min(locus$POS) - 1e+03
    max.pos <- max(locus$POS) + 1e+03
    size.pos <- max.pos - min.pos
    center.pos <- min.pos + (size.pos/2)
    center.100kb.pos <- round(center.pos/1e+05) * 1e+05
    offset.100kb.pos <- round((size.pos/3)/1e+05) * 1e+05
    offset <- (range * 4/3) - range
    big.range <- range + offset
    ystart.gene <- -offset
    ystart.recomb <- -offset + (big.range/8)
    keep.recomb <- subset(recomb, recomb[, 1] > min.pos & recomb[, 
        1] < max.pos)
    genes.in.locus <- subset(genelist, (genelist$START > min.pos & 
        genelist$START < max.pos) | (genelist$STOP > min.pos & 
        genelist$STOP < max.pos))
    print(genes.in.locus)
    markers.in.strong.ld <- subset(locus, (row.names(locus) != 
        snp & locus$RSQR >= 0.8 & locus$TYPE == "typed"))
    markers.in.moderate.ld <- subset(locus, (row.names(locus) != 
        snp & locus$RSQR >= 0.5 & locus$RSQR < 0.8 & locus$TYPE == 
        "typed"))
    markers.in.weak.ld <- subset(locus, (row.names(locus) != 
        snp & locus$RSQR >= 0.2 & locus$RSQR < 0.5 & locus$TYPE == 
        "typed"))
    markers.not.in.ld <- subset(locus, (row.names(locus) != snp & 
        locus$RSQR < 0.2 & locus$TYPE == "typed"))
    imputed.in.strong.ld <- subset(locus, (row.names(locus) != 
        snp & locus$RSQR >= 0.8 & locus$TYPE == "imputed"))
    imputed.in.moderate.ld <- subset(locus, (row.names(locus) != 
        snp & locus$RSQR >= 0.5 & locus$RSQR < 0.8 & locus$TYPE == 
        "imputed"))
    imputed.in.weak.ld <- subset(locus, (row.names(locus) != 
        snp & locus$RSQR >= 0.2 & locus$RSQR < 0.5 & locus$TYPE == 
        "imputed"))
    imputed.not.in.ld <- subset(locus, (row.names(locus) != snp & 
        locus$RSQR < 0.2 & locus$TYPE == "imputed"))
    par(mar = c(4, 4, 3, 4))
    plot(keep.recomb[, 1], ystart.recomb + ((keep.recomb[, 2]/60) * 
        (6 * big.range/8)), type = "l", col = "lightblue", lwd = 1, 
        xlim = c(min.pos, max.pos), ylim = c(-offset, range), 
        xlab = "", ylab = "", main = locusname, axes = F)
    mtext(paste("Chromosome", chr, "position (kb)", sep = " "), 
        side = 1, line = 2.5)
    axis(1, at = c(center.100kb.pos - offset.100kb.pos, center.100kb.pos, 
        center.100kb.pos + offset.100kb.pos), labels = c((center.100kb.pos - 
        offset.100kb.pos)/1000, center.100kb.pos/1000, (center.100kb.pos + 
        offset.100kb.pos)/1000), las = 1)
    axis(2, at = seq(0, range, 2), labels = seq(0, range, 2), 
        las = 1)
    mtext("-log10(Observed p)", side = 2, at = (range/2), line = 2)
    axis(4, at = c(ystart.recomb, ystart.recomb + (big.range/4), 
        ystart.recomb + (2 * big.range/4), ystart.recomb + (3 * 
            big.range/4)), labels = c("0", "20", "40", "60"), 
        las = 1)
    mtext("Recombination rate (cM/Mb)", side = 4, at = (-offset + 
        big.range/2), line = 2)
    box()
    lines(c(min.pos, max.pos), c(0, 0), lty = "dotted", lwd = 1, 
        col = "black")
    points(hit$POS, -(log10(hit$PVAL)), pch = pch, cex = 2.5, 
        bg = "red")
    text(hit$POS, -(log10(hit$PVAL)), labels = c(row.names(hit)), 
        pos = 3, offset = 1)
    if (-(log10(best.pval)) < range) {
        points(hit$POS, -(log10(best.pval)), pch = pch, cex = 2.5, 
            bg = "blue")
        text(hit$POS, -(log10(best.pval)), labels = c(paste("P=", 
            best.pval, sep = "")), pos = 4, offset = 2)
    }
    else {
        points(hit$POS, range, pch = pch, cex = 2.5, bg = "blue")
        text(hit$POS, range, labels = c(paste("P=", best.pval, 
            sep = "")), pos = 4, offset = 1)
    }
    points(markers.not.in.ld$POS, -(log10(markers.not.in.ld$PVAL)), 
        pch = pch, cex = 1, bg = "white")
    points(markers.in.weak.ld$POS, -(log10(markers.in.weak.ld$PVAL)), 
        pch = pch, cex = 1.25, bg = "yellow")
    points(markers.in.moderate.ld$POS, -(log10(markers.in.moderate.ld$PVAL)), 
        pch = pch, cex = 1.25, bg = "orange")
    points(markers.in.strong.ld$POS, -(log10(markers.in.strong.ld$PVAL)), 
        pch = pch, cex = 1.25, bg = "red")
    points(imputed.not.in.ld$POS, -(log10(imputed.not.in.ld$PVAL)), 
        pch = pch, cex = 1, bg = "grey")
    points(imputed.in.weak.ld$POS, -(log10(imputed.in.weak.ld$PVAL)), 
        pch = pch, cex = 1, bg = "grey")
    points(imputed.in.moderate.ld$POS, -(log10(imputed.in.moderate.ld$PVAL)), 
        pch = pch, cex = 1, bg = "grey")
    points(imputed.in.strong.ld$POS, -(log10(imputed.in.strong.ld$PVAL)), 
        pch = pch, cex = 1, bg = "grey")
    for (i in 1:nrow(genes.in.locus)) {
        if (genes.in.locus[i, ]$STRAND == "+") {
            arrows(max(genes.in.locus[i, ]$START, min.pos), -offset, 
                min(genes.in.locus[i, ]$STOP, max.pos), -offset, 
                length = 0.05, lwd = 2, code = 2, lty = "solid", 
                col = "darkgreen")
        }
        else {
            arrows(max(genes.in.locus[i, ]$START, min.pos), -offset, 
                min(genes.in.locus[i, ]$STOP, max.pos), -offset, 
                length = 0.05, lwd = 2, code = 1, lty = "solid", 
                col = "darkgreen")
        }
        if (!is.na(genes.in.locus[i, ]$GENE)) {
            text(genes.in.locus[i, ]$START + (genes.in.locus[i, 
                ]$SIZE/2), -offset + (big.range/20), labels = genes.in.locus[i, 
                ]$GENE, cex = 0.8, srt = 45)
        }
    }
    ltext <- rbind("0.8-","0.5-","0.2-","0.0-","imp.")
    lcol <- c("red","orange","yellow","white","grey")
    legend("topleft",legend=ltext,title="LD (r^2)",fill=lcol)
}
