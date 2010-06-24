alpha.f <- function (b, r, q, k) 
{
    k/((1 - q)^2 + 2 * q * (1 - q) * b + q^2 * r)
}

ci <- function (est, se, estH0 = 0)
{
    wald <- (est - estH0)^2/se^2
    c(est, est - 1.96 * se, est + 1.96 * se,  pchisq(wald, 1, lower.tail=FALSE))
}

cov.invlogit <- function (logit.p1, logit.p2, cov.logit) 
{
    f1 <- invlogit(logit.p1)
    f2 <- invlogit(logit.p2)
    cov.inv <- (f1 - f1^2) * (f2 - f2^2) * cov.logit
    cov.inv
}

Cox.est <- function (case, ctl, k0, initial) 
{
    fitout <- nlm(Cox.T, initial, case = case, control = ctl, 
        k = k0, hessian = TRUE)
    beta <- fitout$estimate[1]
    se.b <- sqrt(1/fitout$hessian[1, 1])
    gamma <- fitout$estimate[2]
    se.g <- sqrt(1/fitout$hessian[2, 2])
    q <- fitout$estimate[3]
    se.q <- sqrt(1/fitout$hessian[3, 3])
    alpha <- alpha.f(fitout$estimate[1], fitout$estimate[2], 
        fitout$estimate[3], k0)
    deviance <- fitout$minimum[1]
    list(alpha = alpha, beta = c(beta, se.b), gamma = c(gamma, 
        se.g), q = c(q, se.q), deviance = c(deviance, pchisq(deviance, 1, lower.tail=FALSE)))
}

Cox.T <- function (parms, case, control, k) 
{
    b <- parms[1]
    r <- parms[2]
    q <- parms[3]
    b <- pfix(b, 0.001, 19.99)
    r <- pfix(r, 0.001, 19.99)
    q <- pfix(q, 0.001, 0.999)
    a <- k/((1 - q)^2 + 2 * q * (1 - q) * b + q^2 * r)
    nt <- sum(case)
    nc <- sum(control)
    gAA <- (1 - q)^2 * a/k
    gAa <- 2 * q * (1 - q) * a * b/k
    gaa <- q^2 * a * r/k
    ecase <- c(gAA, gAa, gaa) * nt
    t <- sum((ecase - case)^2/ecase)
    gAA <- (1 - q)^2 * (1 - a)/(1 - k)
    gAa <- 2 * q * (1 - q) * (1 - a * b)/(1 - k)
    gaa <- q^2 * (1 - a * r)/(1 - k)
    econtrol <- c(gAA, gAa, gaa) * nc
    t <- t + sum((econtrol - control)^2/econtrol)
    return(t)
}

Dev <- function (parms, case, control, k) 
{
    b <- parms[1]
    r <- parms[2]
    q <- parms[3]
    b <- pfix(b, 0.001, 19.99)
    r <- pfix(r, 0.001, 19.99)
    q <- pfix(q, 0.001, 0.999)
    a <- k/((1 - q)^2 + 2 * q * (1 - q) * b + q^2 * r)
    nt <- sum(case)
    nc <- sum(control)
    gAA <- (1 - q)^2 * a/k
    gAa <- 2 * q * (1 - q) * a * b/k
    gaa <- q^2 * a * r/k
    ecase <- c(gAA, gAa, gaa) * nt
    ll <- sum(case * log(case/ecase))
    gAA <- (1 - q)^2 * (1 - a)/(1 - k)
    gAa <- 2 * q * (1 - q) * (1 - a * b)/(1 - k)
    gaa <- q^2 * (1 - a * r)/(1 - k)
    econtrol <- c(gAA, gAa, gaa) * nc
    ll <- ll + sum(control * log(control/econtrol))
    return(2 * ll)
}

Dev.est <- function (case, ctl, k0, initial) 
{
    fitout <- nlm(Dev, initial, case = case, control = ctl, k = k0, 
        hessian = TRUE)
    beta <- fitout$estimate[1]
    se.b <- sqrt(1/fitout$hessian[1, 1])
    gamma <- fitout$estimate[2]
    se.g <- sqrt(1/fitout$hessian[2, 2])
    q <- fitout$estimate[3]
    se.q <- sqrt(1/fitout$hessian[3, 3])
    alpha <- alpha.f(fitout$estimate[1], fitout$estimate[2], 
        fitout$estimate[3], k0)
    deviance <- fitout$minimum[1]
    list(alpha = alpha, beta = c(beta, se.b), gamma = c(gamma, 
        se.g), q = c(q, se.q), deviance = c(deviance, pchisq(deviance, 1, lower.tail=FALSE)))
}

DevH0dominant <- function (parms, case, control, k) 
{
    p <- parms[1]
    r <- parms[2]
    p <- pfix(p, -15, 15)
    p <- exp(p)/(1 + exp(p))
    r <- pfix(r, -10, 10)
    r <- exp(r)
    b <- r
    q <- 1 - p
    a <- k/((1 - q)^2 + 2 * q * (1 - q) * b + q^2 * r)
    nt <- sum(case)
    nc <- sum(control)
    gAA <- (1 - q)^2 * a/k
    gAa <- 2 * q * (1 - q) * a * b/k
    gaa <- q^2 * a * r/k
    ecase <- c(gAA, gAa, gaa) * nt
    if (length(ecase[ecase <= 0]) > 0) {
        ll <- 1e+10
        return(2 * ll)
    }
    ll <- sum(case * log(case/ecase))
    gAA <- (1 - q)^2 * (1 - a)/(1 - k)
    gAa <- 2 * q * (1 - q) * (1 - a * b)/(1 - k)
    gaa <- q^2 * (1 - a * r)/(1 - k)
    econtrol <- c(gAA, gAa, gaa) * nc
    if (length(econtrol[econtrol <= 0]) > 0) {
        ll <- 1e+10
        return(2 * ll)
    }
    ll <- ll + sum(control * log(control/econtrol))
    return(2 * ll)
}

DevH0dominant.est <- function (case, ctl, k0, initial) 
{
    fitout <- nlm(DevH0dominant, initial, case = case, control = ctl, 
        k = k0, hessian = TRUE)
    hessinv <- ginv(fitout$hessian)
    p <- invlogit(fitout$estimate[1])
    se.p <- se.invlogit(fitout$estimate[1], sqrt(hessinv[1, 1]))
    q <- 1 - p
    se.q <- se.p
    gamma <- exp(fitout$estimate[2])
    se.g <- se.exp(fitout$estimate[2], sqrt(hessinv[2, 2]))
    beta <- gamma
    se.b <- se.g
    alpha <- alpha.f(beta, gamma, q, k0)
    deviance <- fitout$minimum[1]
    list(alpha = alpha, beta = c(beta, se.g), gamma = c(gamma, 
        se.g), q = c(q, se.q), deviance = c(deviance, pchisq(deviance, 2, lower.tail=FALSE)))
}

DevH0recessive <- function (parms, case, control, k) 
{
    b <- 1
    p <- parms[1]
    r <- parms[2]
    p <- pfix(p, -15, 15)
    p <- exp(p)/(1 + exp(p))
    r <- pfix(r, -10, 10)
    r <- exp(r)
    q <- 1 - p
    a <- k/((1 - q)^2 + 2 * q * (1 - q) * b + q^2 * r)
    nt <- sum(case)
    nc <- sum(control)
    gAA <- (1 - q)^2 * a/k
    gAa <- 2 * q * (1 - q) * a * b/k
    gaa <- q^2 * a * r/k
    ecase <- c(gAA, gAa, gaa) * nt
    if (length(ecase[ecase <= 0]) > 0) {
        ll <- 1e+10
        return(2 * ll)
    }
    ll <- sum(case * log(case/ecase))
    gAA <- (1 - q)^2 * (1 - a)/(1 - k)
    gAa <- 2 * q * (1 - q) * (1 - a * b)/(1 - k)
    gaa <- q^2 * (1 - a * r)/(1 - k)
    econtrol <- c(gAA, gAa, gaa) * nc
    if (length(econtrol[econtrol <= 0]) > 0) {
        ll <- 1e+10
        return(2 * ll)
    }
    ll <- ll + sum(control * log(control/econtrol))
    return(2 * ll)
}

DevH0recessive.est <- function (case, ctl, k0, initial) 
{
    fitout <- nlm(DevH0recessive, initial, case = case, control = ctl, 
        k = k0, hessian = TRUE)
    hessinv <- ginv(fitout$hessian)
    beta <- 1
    p <- invlogit(fitout$estimate[1])
    se.p <- se.invlogit(fitout$estimate[1], sqrt(hessinv[1, 1]))
    q <- 1 - p
    se.q <- se.p
    gamma <- exp(fitout$estimate[2])
    se.g <- se.exp(fitout$estimate[2], sqrt(hessinv[2, 2]))
    alpha <- alpha.f(beta, gamma, q, k0)
    deviance <- fitout$minimum[1]
    list(alpha = alpha, beta = c(beta, 0), gamma = c(gamma, se.g), 
        q = c(q, se.q), deviance = c(deviance, pchisq(deviance, 2, lower.tail=FALSE)))
}

DevHaGdominant <- function (parms, case, control, k) 
{
    p0 <- parms[1]
    p1 <- parms[2]
    r <- parms[3]
    p0 <- pfix(p0, -15, 15)
    p0 <- exp(p0)/(1 + exp(p0))
    p1 <- pfix(p1, -15, 15)
    p1 <- exp(p1)/(1 + exp(p1))
    r <- pfix(r, -10, 10)
    r <- exp(r)
    b <- r
    pAA <- p0
    pAa <- p1
    paa <- 1 - p0 - p1
    a <- k/(pAA + pAa * b + paa * r)
    nt <- sum(case)
    nc <- sum(control)
    gAA <- pAA * a/k
    gAa <- pAa * a * b/k
    gaa <- paa * a * r/k
    ecase <- c(gAA, gAa, gaa) * nt
    if (length(ecase[ecase <= 0]) > 0) {
        ll <- 1e+10
        return(2 * ll)
    }
    ll <- sum(case * log(case/ecase))
    gAA <- pAA * (1 - a)/(1 - k)
    gAa <- pAa * (1 - a * b)/(1 - k)
    gaa <- paa * (1 - a * r)/(1 - k)
    econtrol <- c(gAA, gAa, gaa) * nc
    if (length(econtrol[econtrol <= 0]) > 0) {
        ll <- 1e+10
        return(2 * ll)
    }
    ll <- ll + sum(control * log(control/econtrol))
    return(2 * ll)
}

DevHaGdominant.est <- function (case, ctl, k0, initial) 
{
    fitout <- nlm(DevHaGdominant, initial, case = case, control = ctl, 
        k = k0, hessian = TRUE)
    hessinv <- ginv(fitout$hessian)
    p0 <- invlogit(fitout$estimate[1])
    se.p0 <- se.invlogit(fitout$estimate[1], sqrt(hessinv[1, 
        1]))
    p1 <- invlogit(fitout$estimate[2])
    se.p1 <- se.invlogit(fitout$estimate[2], sqrt(hessinv[2, 
        2]))
    p <- p0 + p1/2
    a12 <- cov.invlogit(fitout$estimate[1], fitout$estimate[2], 
        hessinv[1, 2])
    se.p <- sqrt(se.p0^2 + se.p1^2/4 + a12)
    q <- 1 - p
    se.q <- se.p
    gamma <- exp(fitout$estimate[3])
    se.g <- se.exp(fitout$estimate[3], sqrt(hessinv[3, 3]))
    beta <- gamma
    se.b <- se.g
    pAA <- p0
    pAa <- p1
    paa <- 1 - p0 - p1
    alpha <- k0/(pAA + pAa * beta + paa * gamma)
    deviance <- fitout$minimum[1]
    list(alpha = alpha, beta = c(beta, se.b), gamma = c(gamma, 
        se.g), q = c(q, se.q), p0 = c(p0, se.p0), p1 = c(p1, 
        se.p1), deviance = c(deviance, pchisq(deviance, 1, lower.tail=FALSE)))
}

DevHaGrecessive <- function (parms, case, control, k) 
{
    b <- 1
    p0 <- parms[1]
    p1 <- parms[2]
    r <- parms[3]
    p0 <- pfix(p0, -15, 15)
    p0 <- exp(p0)/(1 + exp(p0))
    p1 <- pfix(p1, -15, 15)
    p1 <- exp(p1)/(1 + exp(p1))
    r <- pfix(r, -10, 10)
    r <- exp(r)
    pAA <- p0
    pAa <- p1
    paa <- 1 - p0 - p1
    a <- k/(pAA + pAa * b + paa * r)
    nt <- sum(case)
    nc <- sum(control)
    gAA <- pAA * a/k
    gAa <- pAa * a * b/k
    gaa <- paa * a * r/k
    ecase <- c(gAA, gAa, gaa) * nt
    if (length(ecase[ecase <= 0]) > 0) {
        ll <- 1e+10
        return(2 * ll)
    }
    ll <- sum(case * log(case/ecase))
    gAA <- pAA * (1 - a)/(1 - k)
    gAa <- pAa * (1 - a * b)/(1 - k)
    gaa <- paa * (1 - a * r)/(1 - k)
    econtrol <- c(gAA, gAa, gaa) * nc
    if (length(econtrol[econtrol <= 0]) > 0) {
        ll <- 1e+10
        return(2 * ll)
    }
    ll <- ll + sum(control * log(control/econtrol))
    return(2 * ll)
}

DevHaGrecessive.est <- function (case, ctl, k0, initial) 
{
    fitout <- nlm(DevHaGrecessive, initial, case = case, control = ctl, 
        k = k0, hessian = TRUE)
    hessinv <- ginv(fitout$hessian)
    beta <- 1
    p0 <- invlogit(fitout$estimate[1])
    se.p0 <- se.invlogit(fitout$estimate[1], sqrt(hessinv[1, 
        1]))
    p1 <- invlogit(fitout$estimate[2])
    se.p1 <- se.invlogit(fitout$estimate[2], sqrt(hessinv[2, 
        2]))
    p <- p0 + p1/2
    a12 <- cov.invlogit(fitout$estimate[1], fitout$estimate[2], 
        hessinv[1, 2])
    se.p <- sqrt(se.p0^2 + se.p1^2/4 + a12)
    q <- 1 - p
    se.q <- se.p
    gamma <- exp(fitout$estimate[3])
    se.g <- se.exp(fitout$estimate[3], sqrt(hessinv[3, 3]))
    pAA <- p0
    pAa <- p1
    paa <- 1 - p0 - p1
    alpha <- k0/(pAA + pAa * beta + paa * gamma)
    deviance <- fitout$minimum[1]
    list(alpha = alpha, beta = c(beta, 0), gamma = c(gamma, se.g), 
        q = c(q, se.q), p0 = c(p0, se.p0), p1 = c(p1, se.p1), 
        deviance = c(deviance, pchisq(deviance, 1, lower.tail=FALSE)))
}

DevHa.recessive <- function (parms, case, control, k) 
{
    b <- 1
    p <- parms[1]
    r <- parms[2]
    f <- parms[3]
    p <- pfix(p, -15, 15)
    p <- exp(p)/(1 + exp(p))
    r <- pfix(r, -10, 10)
    r <- exp(r)
    f <- pfix(f, -15, 15)
    f <- exp(f)/(1 + exp(f))
    geno.prob <- geno.p(p, f)
    pAA <- geno.prob[1]
    pAa <- geno.prob[2]
    paa <- geno.prob[3]
    a <- k/(pAA + pAa * b + paa * r)
    nt <- sum(case)
    nc <- sum(control)
    gAA <- pAA * a/k
    gAa <- pAa * a * b/k
    gaa <- paa * a * r/k
    ecase <- c(gAA, gAa, gaa) * nt
    if (length(ecase[ecase <= 0]) > 0) {
        ll <- 1e+10
        return(2 * ll)
    }
    ll <- sum(case * log(case/ecase))
    gAA <- pAA * (1 - a)/(1 - k)
    gAa <- pAa * (1 - a * b)/(1 - k)
    gaa <- paa * (1 - a * r)/(1 - k)
    econtrol <- c(gAA, gAa, gaa) * nc
    if (length(econtrol[econtrol <= 0]) > 0) {
        ll <- 1e+10
        return(2 * ll)
    }
    ll <- ll + sum(control * log(control/econtrol))
    return(2 * ll)
}
DevHa.recessive.est <- function (case, ctl, k0, initial) 
{
    fitout <- nlm(DevHa.recessive, initial, case = case, control = ctl, 
        k = k0, hessian = TRUE)
    hessinv <- ginv(fitout$hessian)
    beta <- 1
    p <- invlogit(fitout$estimate[1])
    se.p <- se.invlogit(fitout$estimate[1], sqrt(hessinv[1, 1]))
    q <- 1 - p
    se.q <- se.p
    gamma <- exp(fitout$estimate[2])
    se.g <- se.exp(fitout$estimate[2], sqrt(hessinv[2, 2]))
    f <- invlogit(fitout$estimate[3])
    se.f <- se.invlogit(fitout$estimate[3], sqrt(hessinv[3, 3]))
    geno.prob <- geno.p(p, f)
    pAA <- geno.prob[1]
    pAa <- geno.prob[2]
    paa <- geno.prob[3]
    alpha <- k0/(pAA + pAa * beta + paa * gamma)
    deviance <- fitout$minimum[1]
    list(alpha = alpha, beta = c(beta, 0), gamma = c(gamma, se.g), 
        q = c(q, se.q), f = c(f, se.f), deviance = c(deviance, 
            pchisq(deviance, 1, lower.tail=FALSE)))
}

ExpN.recessiveH0 <- function (parms, case, control, k) 
{
    b <- 1
    q <- parms[1]
    r <- parms[2]
    p <- 1 - q
    a <- k/((1 - q)^2 + 2 * q * (1 - q) * b + q^2 * r)
    nt <- sum(case)
    nc <- sum(control)
    gAA <- (1 - q)^2 * a/k
    gAa <- 2 * q * (1 - q) * a * b/k
    gaa <- q^2 * a * r/k
    ecase <- c(gAA, gAa, gaa) * nt
    ll <- sum(case * log(case/ecase))
    gAA <- (1 - q)^2 * (1 - a)/(1 - k)
    gAa <- 2 * q * (1 - q) * (1 - a * b)/(1 - k)
    gaa <- q^2 * (1 - a * r)/(1 - k)
    econtrol <- c(gAA, gAa, gaa) * nc
    ll <- ll + sum(control * log(control/econtrol))
    list(case = case, ecase = ecase, control = control, econtrol = econtrol, 
        Dev = 2 * ll)
}

ExpN.recessiveHa <- function (parms, case, control, k) 
{
    b <- 1
    p0 <- parms[1]
    p1 <- parms[2]
    r <- parms[3]
    p2 <- 1 - p0 - p1
    a <- k/(p0 + p1 * b + p2 * r)
    nt <- sum(case)
    nc <- sum(control)
    gAA <- p0 * a/k
    gAa <- p1 * a * b/k
    gaa <- p2 * a * r/k
    ecase <- c(gAA, gAa, gaa) * nt
    ll <- sum(case * log(case/ecase))
    gAA <- p0 * (1 - a)/(1 - k)
    gAa <- p1 * (1 - a * b)/(1 - k)
    gaa <- p2 * (1 - a * r)/(1 - k)
    econtrol <- c(gAA, gAa, gaa) * nc
    ll <- ll + sum(control * log(control/econtrol))
    list(case = case, ecase = ecase, control = control, econtrol = econtrol, 
        Dev = 2 * ll)
}

geno.p <- function (p, f) 
{
    q <- 1 - p
    pAA <- p^2 + p * q * f
    pAa <- 2 * p * q * (1 - f)
    paa <- q^2 + p * q * f
    c(pAA, pAa, paa, p, f)
}

hw <- function (tt = 10, tc = 10, cc = 10) 
{
    nt <- 2 * tt + tc
    nc <- 2 * cc + tc
    pt <- nt/(nt + nc)
    pc <- nc/(nt + nc)
    n <- tt + tc + cc
    ett <- n * pt^2
    etc <- 2 * n * pt * pc
    ecc <- n * pc^2
    stat <- (ett - tt)^2/ett + (ecc - cc)^2/ecc + (etc - tc)^2/etc
    p <- pchisq(stat, 1, lower.tail=FALSE)
    list(stat = stat, p = p)
}

invlogit <- function (x = 0) exp(x)/(1 + exp(x))

logit <- function (p = 0.5) log(p/(1 - p))

pfix <- function (parm, low, high) 
{
    if (parm < low) 
        return(low)
    if (parm > high) 
        return(high)
    return(parm)
}

se.exp <- function (p, se.p) 
{
    f1 <- exp(p)
    se.delta <- sqrt((f1)^2 * se.p^2)
    se.delta
}

se.invlogit <- function (logit.p, se.logit) 
{
    f1 <- invlogit(logit.p)
    se.inv <- (f1 - f1^2) * se.logit
    se.inv
}

hwe.cc <- function(model, case, ctrl, k0, initial1, initial2)
{
model.int <- charmatch(model,c("dominant","recessive"))
if(is.na(model.int)) stop("Invalide model specification")
if(model.int==0) stop("Ambiguous model specification")

Cox1 <- Cox.est(case,ctrl,k0,initial1) 
if (model.int==1) {
   t2par <- DevH0dominant.est(case,ctrl,k0,initial2)
   p0i <- (case[1]+ctrl[1])/(sum(case)+sum(ctrl))
   p1i <- (case[2]+ctrl[2])/(sum(case)+sum(ctrl))
   initial <- c(logit(p0i), logit(p1i), log(1.0753569)) # logit(p0), logit(p1), log(gamma)
   t3par <- DevHaGdominant.est(case,ctrl,k0,initial)
}
if (model.int==2) {
   t2par <- DevH0recessive.est(case,ctrl,k0,initial2)
   p0i <- (case[1]+ctrl[1])/(sum(case)+sum(ctrl))
   p1i <- (case[2]+ctrl[2])/(sum(case)+sum(ctrl))
   initial <- c(logit(p0i), logit(p1i), log(1.0753569)) # logit(p0), logit(p1), log(gamma) 
   t3par <- DevHaGrecessive.est(case,ctrl,k0,initial)
}
lrt1.stat <- t2par$deviance[1]-t3par$deviance[1]

list(Cox=Cox1, t2par=t2par, t3par=t3par, lrt.stat=lrt1.stat, pval=pchisq(lrt1.stat,1,lower.tail=FALSE))
}
