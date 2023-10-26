get_probabilities <- function(prop, numBins=64) {
    r = range(prop, na.rm = T)
    b = seq(from=r[1], to=r[2], length.out=numBins+1 )
    yprop <- table(cut(prop, breaks=b , include.lowest=TRUE) )
    idx <- which(yprop > 0)
    yx = yprop[idx]           # remove bins with zero counts
    n = sum(yx)             # total number of counts
    p = yx/n                # empirical frequencies
    f1 = sum(yx == 1)       # number of singletons
    if (f1 == n) f1 = n-1   # avoid C=0
    C = 1 - f1/n            # estimated coverage
    pa = C*p                # coverage adjusted empirical frequencies
    la = (1-(1-pa)^n)       # probability to see a bin (species) in the sample
    
    b = b[idx]
    res <- list(probs_adj = pa, probs_bin = la, breaks = b)
    res
}

calculate_gas_phase_entropy <- function(pdata) {
    df <- read.csv(text = pdata, header=F);
    colnames(df) <- c("Area", "SHPIDX")
    nbins = 64

    SHPIDX = df$SHPIDX
    Area <- df$Area
    idx_na = which(is.na(SHPIDX))
    if (length(idx_na) > 0) {
        SHPIDX <- SHPIDX[-idx_na]
        Area <- Area[-idx_na]
    }
    rm(idx_na, df)

    shpidx_n <- NULL
    shpidx_p <- NULL
    res_p <- NULL
    res_n <- NULL
    la_p <- NULL
    la_n <- NULL
    pa_p <- NULL
    pa_n <- NULL

    # positive contrib
    idx <- which((SHPIDX < 0.99) & (SHPIDX >= 0))
    if (length(idx) > 1) {
        shpidx_p <- SHPIDX[idx]
        res_p <- get_probabilities(shpidx_p, nbins)
        la_p = res_p$probs_bin
        pa_p = res_p$probs_adj
    }


    # negative contrib
    idx <- which(SHPIDX < 0)
    if (length(idx) > 1) {
        shpidx_n <- SHPIDX[idx]
        res_n <- get_probabilities(shpidx_n, nbins)
        la_n = res_n$probs_bin
        pa_n = res_n$probs_adj
    }

    rm(shpidx_p, shpidx_n, idx)

    ntri = length(Area)

    pos_ent = numeric(ntri)
    neg_ent = numeric(ntri)

    for (j in 1:ntri) {
        
        sival = SHPIDX[j]

        if ((sival >= 0) & (sival < 1)) {
            interval_c = findInterval(sival, res_p$breaks)
            val = 0
            bval = as.double(la_p[interval_c])
            pval = as.double(pa_p[interval_c])
            xval = pval * log2(pval)/bval
            if (!is.na(xval)) {
                val = xval
            }
            pos_ent[j] = -1 * val
        } else if (sival < 0) {
            interval_c = findInterval(sival, res_n$breaks)
            val = 0
            bval = as.double(la_n[interval_c])
            pval = as.double(pa_n[interval_c])
            xval = pval * log2(pval)/bval
            if (!is.na(xval)) {
                val = xval
            }
            neg_ent[j] = -1 * val
        }
    }

    rm(res_p, res_n, la_p, pa_p, la_n, pa_n)

    x = c(101.809, 2.032308, -19.94338, 2.082685)

    tdata = as.matrix(cbind(pos_ent, neg_ent))
    gvals = tdata %*% x[c(3,4)]
    val = sum((x[2] + gvals) * Area)
    res = round(val + x[1], 4)
    res
}


# file_str <- paste(readLines("benzene.srf"), collapse="\n")
