a2g <- function(a1,a2)
{
  i <- ifelse(a1 < a2,a2,a1)
  j <- ifelse(a1 < a2,a1,a2)
  genocode <- ifelse (j==0, 0, i*(i-1)/2 + j)
  invisible (genocode)

}

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
    return(cbind(j,i))
}

is.miss <- function(data,data.int,miss.val=0)
{
   if (data.int==2) # genotype
   {
      id <- array(FALSE, length(data))
      for (i in 1:length(miss.val))
          id <- id | data==miss.val[i]
   }
   else # allele
   {
      id <- array(FALSE, length(data[,1]))  
      for (i in 1:length(miss.val))
          id <- id | apply(data==miss.val[i],1,any)
   }
   return (id)  
}

# adapted from gcp.c on 24/9/2004
# 17-6-2004 JH Zhao

revhap <- function(loci,hapid)
{
   nloci <- length(loci)
   nalleles <- vector("numeric",nloci)
   nalleles[nloci] <- 1
   m <- nloci
   for(k in nloci:1)
   {
      m <- m - 1
      nalleles[m] <- nalleles[m+1] * loci[k]
   }
  n <- length(hapid)
  hap <- matrix(0,n,nloci)
  for (i in 1:n)
  {
    l <- hapid[i] - 1
    for (j in 1:nloci)
    {
      hap[i,j] <- floor(l/nalleles[j])
      if (j==nloci) hap[i,j] <- l
      else l <- l%%nalleles[j]
      hap[i,j] <- hap[i,j] + 1
    }
  }
  invisible(hap)

}

revhap.i <- function(loci,hapid)
{
   nloci <- length(loci)
   nalleles <- vector("numeric",nloci)
   nalleles[nloci] <- 1
   m <- nloci
   for(k in nloci:1)
   {
      m <- m - 1
      nalleles[m] <- nalleles[m+1] * loci[k]
   }

  hap <- vector("numeric",nloci)
  l <- hapid - 1
  for (j in 1:nloci)
  {
    hap[j] <- floor(l/nalleles[j])
    if (j==nloci) hap[j] <- l
    else l <- l%%nalleles[j]
    hap[j] <- hap[j] + 1
  }
  invisible(hap)

}

# adapted from haplo.score 8/11/2004

Ginv<-function(x)
{
savesvd<-svd(x)                                  
U.svd<-savesvd$u                                   
V.svd<-savesvd$v
d.svd<-savesvd$d
eps<-.000001                                       
maxd<-max(d.svd)                                  
w<-ifelse((d.svd/maxd) < eps, rep(0,length(d.svd)), 1/d.svd)
df<-sum(d.svd/maxd >= eps)
Ginv<-V.svd %*% diag(w) %*% t(U.svd)

list(Ginv=Ginv,rank=df)
}

geno.recode <- function(geno, miss.val=0){

n.loci <- ncol(geno)/2
alist <- vector("list",n.loci)
grec <- NULL

for(i in 1:n.loci){
   t <- (i-1)*2 + 1
   tmp <- allele.recode(geno[,t],geno[,t+1], miss.val=miss.val)
   grec <- cbind(grec,tmp$a1,tmp$a2)
   alist[[i]] <- list(allele=tmp$allele.label)
}

return(list(grec=grec, alist=alist))

}

allele.recode <- function (a1, a2, miss.val = NA)
{
    n <- length(a1)
    if (is.factor(a1))
        a1 <- as.character(a1)
    if (is.factor(a2))
        a2 <- as.character(a2)
    is.ch <- is.character(a1) | is.character(a2)
    if (is.ch) {
        t <- factor(c(a1, a2), exclude = miss.val)
    }
    if (!is.ch) {
        lev <- sort(unique(c(a1, a2)))
        t <- factor(c(a1, a2), levels = lev, exclude = miss.val)
    }
    allele.label <- levels(t)
    t <- as.numeric(t)
    a1 <- t[1:n]
    a2 <- t[(n + 1):(2 * n)]
    return(list(a1 = a1, a2 = a2, allele.label = allele.label))
}

gcode <- function(a1,a2) {

  i <- ifelse(a1 < a2,a2,a1)
  j <- ifelse(a1 < a2,a1,a2)
  genocode <- i*(i-1)/2 + j
  return(genocode)

}

ungcode <- function(x) {

  i <- 1 + floor((sqrt(8*x+1)-1)/2)
  j <- x - i*(i-1)/2
  
  # the following 2 lines were added as a patch to make this ungcode work:
  i <- ifelse(j==0,i-1,i)
  j <- ifelse(j==0,i,j)

  return(cbind(j,i))
}

grec2g <- function (h, n, t)
{
  hh <- h
  for (i in 1:n)
  {
    hh[,i] <- t$alist[[i]]$allele[h[,i]]
  }
  invisible(hh)
}

# compositeLD adapted from Jason Sinnwell

"compositeLD" <-
function (a1, a2, b1, b2) 
{
    len1 <- length(a1)
    if (length(a2) != len1 | length(b1) != len1 | length(b2) != 
        len1) {
        stop("Input allele vectors must all be the same length")
    }
    zed <- apply(!is.na(cbind(a1, a2, b1, b2)), 1, all)
    if (sum(!zed) == length(zed)) {
        stop("all obs have at least one missing allele")
    }
    a1 <- a1[zed]
    a2 <- a2[zed]
    b1 <- b1[zed]
    b2 <- b2[zed]
    tmp.geno <- geno.recode(cbind(a1, a2, b1, b2))
    a1 <- tmp.geno$grec[, 1]
    a2 <- tmp.geno$grec[, 2]
    b1 <- tmp.geno$grec[, 3]
    b2 <- tmp.geno$grec[, 4]
    a.allele <- tmp.geno$alist[[1]]$allele
    b.allele <- tmp.geno$alist[[2]]$allele
    locus.info <- function(a1, a2) {
        t1 <- ifelse(a1 < a2, a1, a2)
        t2 <- ifelse(a2 > a1, a2, a1)
        a1 <- t1
        a2 <- t2
        afreq <- table(c(a1, a2))
        afreq <- afreq/sum(afreq)
        a <- sort(unique(c(a1, a2)))
        tmp <- expand.grid(a, a)
        tmp <- cbind(tmp[, 2], tmp[, 1])
        tmp <- tmp[tmp[, 1] <= tmp[, 2], ]
        u1 <- tmp[, 1]
        u2 <- tmp[, 2]
        ugeno <- gcode(u1, u2)
        geno <- factor(gcode(a1, a2), levels = ugeno)
        n.geno <- table(geno)
        exp.hwe <- afreq[u1] * afreq[u2]
        exp.hwe <- ifelse(u1 != u2, 2 * exp.hwe, exp.hwe)
        obs <- n.geno/sum(n.geno)
        D.hw <- (exp.hwe - obs)/2
        hom <- u1 == u2
        het <- !hom
        indx.hom <- (1:length(hom))[hom]
        for (i in indx.hom) {
            v1 <- u1[i]
            D.hw[i] <- sum(D.hw[(u1 == v1 | u2 == v1) & het])
        }
        ugeno.df = data.frame(a1 = tmp[, 1], a2 = tmp[, 2], geno = ugeno, 
            n = n.geno, probs.obs = obs, prob.hwe = exp.hwe, 
            D.hw = D.hw)
        geno.df = data.frame(a1 = a1, a2 = a2, geno = geno)
        return(list(afreq = afreq, geno.df = geno.df, ugeno.df = ugeno.df))
    }
    setup.countmat <- function(n.a) {
        single.mat <- matrix(rep(diag(0.5, n.a), n.a), ncol = n.a, 
            byrow = T)
        block.mat <- matrix(apply(diag(0.5, n.a), 1, rep, times = rep(n.a, 
            n.a)), ncol = n.a)
        mat <- single.mat + block.mat
        rm.indx <- unlist(sapply(1:(n.a - 1), indx <- function(x, 
            n) {
            return(x * n + 1:x)
        }, n = n.a))
        mat <- mat[-rm.indx, ]
        return(mat)
    }
    D.hw <- function(a.count, b.count, c.count, geno.vec, n) {
        pa <- sum(a.count * geno.vec)/n
        pb <- sum(b.count * geno.vec)/n
        D.ab <- sum(c.count * geno.vec)/n - 2 * pa * pb
        return(D.ab)
    }
    tmp.a <- locus.info(a1, a2)
    tmp.b <- locus.info(b1, b2)
    tbl.count <- table(tmp.a$geno.df$geno, tmp.b$geno.df$geno)
    geno.vec <- as.vector(tbl.count)
    prgeno.obs <- geno.vec/sum(geno.vec)
    n <- sum(geno.vec)
    d.tmp <- ifelse(tmp.a$ugeno.df$a1 == tmp.a$ugeno.df$a2, tmp.a$ugeno$D.hw.Freq, 
        -2 * tmp.a$ugeno$D.hw.Freq)
    marg.a <- tmp.a$ugeno.df$prob.hwe + d.tmp
    d.tmp <- ifelse(tmp.b$ugeno.df$a1 == tmp.b$ugeno.df$a2, tmp.b$ugeno$D.hw.Freq, 
        -2 * tmp.b$ugeno$D.hw.Freq)
    marg.b <- tmp.b$ugeno.df$prob.hwe + d.tmp
    prgeno.null <- as.vector(marg.a %o% marg.b)
    prgenoHWE <- as.vector(tmp.a$ugeno.df$prob.hwe %o% tmp.b$ugeno.df$prob.hwe)
    n.a <- length(tmp.a$afreq)
    n.b <- length(tmp.b$afreq)
    a.sub <- setup.countmat(n.a)
    b.sub <- setup.countmat(n.b)
    a.counts <- matrix(apply(a.sub, 2, rep, times = nrow(b.sub)), 
        ncol = n.a)
    b.counts <- matrix(apply(b.sub, 1, rep, times = nrow(a.sub)), 
        ncol = n.b, byrow = T)
    D.df <- NULL
    c.counts <- NULL
    for (a in 1:(n.a - 1)) {
        for (b in 1:(n.b - 1)) {
            c.vec <- ifelse(a.counts[, a] == 0 | b.counts[, b] == 
                0, 0, trunc(a.counts[, a] + b.counts[, b]))
            c.vec <- ifelse(a.counts[, a] == 0.5 & b.counts[, 
                b] == 0.5, 0.5, c.vec)
            D.df <- rbind(D.df, c(a, b, D.hw(a.counts[, a], b.counts[, 
                b], c.vec, geno.vec, n)))
            c.counts <- cbind(c.counts, c.vec)
        }
    }
    pa.vec <- t(a.counts) %*% prgeno.null
    pb.vec <- t(b.counts) %*% prgeno.null
    D.df <- data.frame(D.df)
    names(D.df) <- c("a", "b", "delta")
    nr <- nrow(D.df)
    delta.cov <- matrix(rep(0, nr^2), ncol = nr)
    delta.covHWE <- matrix(rep(0, nr^2), ncol = nr)
    mtx <- cbind(D.df$a, D.df$b)
    indx <- expand.grid(1:nr, 1:nr)
    indx <- indx[!(indx[, 1] < indx[, 2]), ]
    indx <- matrix(c(indx[, 1], indx[, 2]), byrow = F, ncol = 2)
    indx2 <- cbind(indx[, 2], indx[, 1])
    nr.counts <- nrow(a.counts)
    part1 <- 2 * (matrix(rep(pa.vec[mtx[indx[, 1], 1]], nr.counts), 
        nrow = nr.counts, byrow = T) * b.counts[, mtx[indx[, 
        1], 2]] + matrix(rep(pb.vec[mtx[indx[, 1], 2]], nr.counts), 
        nrow = nr.counts, byrow = T) * a.counts[, mtx[indx[, 
        1], 1]])
    part2 <- 2 * (matrix(rep(pa.vec[mtx[indx[, 2], 1]], nr.counts), 
        nrow = nr.counts, byrow = T) * b.counts[, mtx[indx[, 
        2], 2]] + matrix(rep(pb.vec[mtx[indx[, 2], 2]], nr.counts), 
        nrow = nr.counts, byrow = T) * a.counts[, mtx[indx[, 
        2], 1]])
    tmat1 <- (c.counts[, indx[, 1]] - part1) * (c.counts[, indx[, 
        2]] - part2)
    tmat <- tmat1 * prgeno.null
    tmatHWE <- tmat1 * prgenoHWE
    t.tot <- (-t(prgeno.null) %*% c.counts[, indx[, 1]] + 4 * 
        pa.vec[mtx[indx[, 1], 1]] * pb.vec[mtx[indx[, 1], 2]]) * 
        (-t(prgeno.null) %*% c.counts[, indx[, 2]] + 4 * pa.vec[mtx[indx[, 
            2], 1]] * pb.vec[mtx[indx[, 2], 2]])
    t.totHWE <- (-t(prgenoHWE) %*% c.counts[, indx[, 1]] + 4 * 
        pa.vec[mtx[indx[, 1], 1]] * pb.vec[mtx[indx[, 1], 2]]) * 
        (-t(prgenoHWE) %*% c.counts[, indx[, 2]] + 4 * pa.vec[mtx[indx[, 
            2], 1]] * pb.vec[mtx[indx[, 1], 2]])
    sum.null <- unlist(apply(tmat, 2, sum))
    sum.HWE <- unlist(apply(tmatHWE, 2, sum))
    delta.cov[indx2] <- delta.cov[indx] <- (sum.null - as.vector(t.tot))/n
    delta.covHWE[indx2] <- delta.covHWE[indx] <- (sum.HWE - as.vector(t.totHWE))/n
    stat.delta <- D.df$delta^2/diag(delta.cov)
    stat.deltaHWE <- D.df$delta^2/diag(delta.covHWE)
    if (length(D.df$delta) == 1) {
        stat.global <- stat.delta
        stat.globalHWE <- stat.deltaHWE
        df <- 1
        dfHWE <- 1
    }
    else {
        tmp <- Ginv(delta.cov)
        inv <- tmp$Ginv
        df <- tmp$rank
        stat.global <- t(D.df$delta) %*% inv %*% D.df$delta
        tmp <- Ginv(delta.covHWE)
        inv <- tmp$Ginv
        dfHWE <- tmp$rank
        stat.globalHWE <- t(D.df$delta) %*% inv %*% D.df$delta
    }
    D.tbl <- D.df
    D.tbl[, 1] <- a.allele[D.df[, 1]]
    D.tbl[, 2] <- b.allele[D.df[, 2]]
    D.tbl = cbind(D.tbl, var.delta = diag(delta.cov), chistat = stat.delta, 
        pval = 1 - pchisq(stat.delta, 1))
    obj <- list(D.tbl = D.tbl, chistat.global = stat.global, 
        df = df, delta.cov = delta.cov)
    oldClass(obj) <- "compositeLD"
    return(obj)
}

"print.compositeLD" <-
function (x, show.all.pairs = TRUE, digits = max(options()$digits - 
    3, 4), ...) 
{
    printBanner("Global Composite LD Test Statistic")
    cat("chi-square-stat: ", round(x$chistat.global, digits), 
        "\n")
    cat("             df: ", round(x$df, digits), "\n")
    cat("        p-value: ", round(1 - pchisq(x$chistat.global, 
        x$df), digits), "\n\n")
    if (show.all.pairs) {
        printBanner("LD Measures for all Pairs of Alleles")
        tmp <- x$D.tbl
        for (j in 3:6) {
            tmp[, j] <- round(tmp[, j], digits)
        }
        print.data.frame(tmp)
    }
}

"printBanner" <-
function (str, banner.width = 80, char.perline = 60, border = "=") 
{
    vec <- str
    new <- NULL
    onespace <- FALSE
    for (i in 1:nchar(vec)) {
        if (substring(vec, i, i) == " " && onespace == FALSE) {
            onespace <- TRUE
            new <- paste(new, substring(vec, i, i), sep = "")
        }
        else if (substring(vec, i, i) == " " && onespace == TRUE) {
            onespace <- TRUE
        }
        else {
            onespace <- FALSE
            new <- paste(new, substring(vec, i, i), sep = "")
        }
    }
    where.blank <- NULL
    indx <- 1
    for (i in 1:nchar(new)) {
        if ((substring(new, i, i) == " ")) {
            where.blank[indx] <- i
            indx <- indx + 1
        }
    }
    j <- length(where.blank) + 1
    where.blank[j] <- nchar(new)
    begin <- 1
    end <- max(where.blank[where.blank <= char.perline])
    end.ok <- is.na(end)
    if (end.ok == TRUE) {
        char.perline <- floor(banner.width/2)
        end <- max(where.blank[where.blank <= char.perline])
    }
    cat(paste(rep(border, banner.width), collapse = ""), "\n")
    repeat {
        titleline <- substring(new, begin, end)
        n <- nchar(titleline)
        if (n < banner.width) {
            n.remain <- banner.width - n
            n.left <- floor(n.remain/2)
            n.right <- n.remain - n.left
            for (i in 1:n.left) titleline <- paste(" ", titleline, 
                sep = "")
            for (i in 1:n.right) titleline <- paste(titleline, 
                " ", sep = "")
            n <- nchar(titleline)
        }
        cat(titleline, "\n")
        begin <- end + 1
        end.old <- end
        tmp <- where.blank[(end.old < where.blank) & (where.blank <= 
            end.old + char.perline + 1)]
        if (length(tmp)) 
            end <- max(tmp)
        else break
    }
    cat(paste(rep(border, banner.width), collapse = ""), "\n\n")
    invisible()
}
