# Shared by ascat.R and facets.R (sourced from the script's own directory; base R only):
# command-line arguments, the pieces of the both-allele purity check, the refit decision, and the arithmetic of
# the sex check and of tumour cells in the normal.

PAR_X <- c(2781479, 155701383)          # hg38 X: PAR1 end, PAR2 start
SEX_XX_MIN <- 0.3                       # X/autosome heterozygosity ratio: two X copies at or above this,
SEX_XY_MAX <- 0.1                       # one X at or below this, unclear in between

# --key=value arguments -> get(key, default). An empty value counts as not given; without a default the key is
# required.
cli_args <- function(args) {
    a <- sub("^--", "", args)
    arg <- setNames(sub("^[^=]*=", "", a), sub("=.*$", "", a))
    function(k, default = NULL) {
        if (k %in% names(arg) && nzchar(arg[[k]])) return(arg[[k]])
        if (is.null(default)) stop("missing --", k)
        default
    }
}
is_true <- function(x) tolower(x) %in% c("true", "t", "1", "yes")

# ---- sex from the normal's X heterozygosity -------------------------------------------------------------------
# ratio = (het share of X SNPs outside the PARs) / (het share of autosomal SNPs). Two X copies give roughly the
# autosomal rate (a little lower: less diversity on X), one X gives ~0 (mapping noise). A wrong sex misfits X,
# which enters the chromosome count (+X is a common HeH gain).
SEX_NOT_ASSESSED <- list(ratio = NA, data = "unknown", status = "not_assessed")
sex_call <- function(ratio, given, xx_min = SEX_XX_MIN, xy_max = SEX_XY_MAX) {
    data <- if (ratio >= xx_min) "XX" else if (ratio <= xy_max) "XY" else "unclear"
    list(ratio = ratio, data = data, status = if (data == "unclear") "unclear" else if (data == given) "ok" else "mismatch")
}

# ---- tumour cells in the normal (TIN) --------------------------------------------------------------------------
# A normal with a fraction t of tumour cells is shifted, at each germline het SNP, towards the tumour's major
# allele. For a segment with a + b copies per tumour cell:
#   normal B-allele fraction - 0.5 = delta = t(a - b) / (2(2(1 - t) + t(a + b)))
#   => t = 2 delta / ((a - b)/2 + delta (2 - a - b))
# The tumour's B-allele fraction at the same SNP says which allele is the major one, so the signed shift
# s = (normal BAF - 0.5) x sign(tumour BAF - 0.5) averages to delta (the two samples' noise is independent).
# Every clearly imbalanced segment is used: autosomal, >= 1 Mb, a - b >= 0.5, and the tumour sample's expected
# B-allele fraction >= min_dev from 0.5. Where the tumour's BAF is close to 0.5 some SNPs point the wrong way,
# which shrinks the mean by f = 1 - 2 P(wrong side) = 1 - 2 pnorm(-dev / sigma), sigma = the tumour BAF's SD at
# het SNPs in 1+1 segments; the mean is divided by f. Weighted by het SNPs; NA below min_het SNPs in total.
#   s, tb  per het SNP: the signed normal shift and the tumour's B-allele fraction
#   snps   function(i) -> the het SNPs (logical) in segment i
#   a, b   copies per tumour cell per segment (a >= b >= 0); ok: autosomal segments >= 1 Mb
tin_core <- function(s, tb, snps, a, b, p, ok, min_het = 200, min_dev = 0.12) {
    na  <- list(tin = NA, se = NA, n_het = 0)
    dev <- (p * a + 1 - p) / (p * (a + b) + 2 * (1 - p)) - 0.5     # tumour sample's BAF - 0.5, major allele
    kb  <- Reduce(`|`, lapply(which(ok & abs(a - 1) < 0.1 & abs(b - 1) < 0.1), snps), logical(length(s)))
    sigma <- if (sum(kb) >= 100) sd(tb[kb]) else 0.08
    t <- n <- v <- numeric()
    for (i in which(ok & a - b >= 0.5 & dev >= min_dev)) {
        k <- snps(i)
        if (sum(k) < 20) next
        f <- 1 - 2 * pnorm(-dev[i] / sigma)                              # sign-error shrinkage
        d <- mean(s[k]) / f; g <- (a[i] - b[i]) / 2 + d * (2 - a[i] - b[i])
        t <- c(t, 2 * d / g); n <- c(n, sum(k)); v <- c(v, var(s[k]) / f^2 * (2 / g)^2)
    }
    if (sum(n) < min_het) return(modifyList(na, list(n_het = sum(n))))
    w <- n / sum(n)
    list(tin = min(max(sum(w * t), 0), 1), se = sqrt(sum(w^2 * v / n)), n_het = sum(n))
}

# ---- the both-allele purity check ------------------------------------------------------------------------------
# The tools score a purity by how close the minor allele (ASCAT) or the EM states (FACETS) are to integers; a
# large subclonal change can then pull purity until it looks clonal. The check scores both alleles of every
# segment instead, each allele's squared distance to the nearest integer capped at cap^2, so a subclonal segment
# costs the same however far it is from a clonal state and stops steering the fit. A negative copy number is
# impossible and is not capped: an LOH segment bounds purity from below instead of setting it.
allele_off <- function(x, cap) ifelse(x < 0, x^2, pmin((x - round(x))^2, cap^2))

# A fit is feasible when every confident segment (conf) has both alleles >= -tol copies per tumour cell; a
# deletion in a share x of all cells needs purity >= ~x. The purity floor is the lowest feasible purity.
alleles_feasible <- function(major, minor, conf, tol) all(pmin(minor, major)[conf] >= -tol)

# Fraction of autosomal het SNPs (weights w, segments use) in segments with an allele > dev from 1+1.
# Purity can only be estimated from aberrant segments; below a small fraction the genome is flat.
aberrant_share <- function(a, b, w, use, dev = 0.15) {
    ab <- use & (abs(a - 1) > dev | abs(b - 1) > dev)
    if (sum(w[use]) == 0) return(0)
    sum(w[ab]) / sum(w[use])
}

# Whole-genome doubling: a doubled solution scores at least as well as its undoubled one (every state doubles,
# subclonal half-copies become whole), so the best minimum leans towards doubling. Take the lower-ploidy minimum
# (<= best / 1.6) unless the doubled one is >= doubling_ratio x better, and do not count doublings of the chosen
# fit as its runner-up. minima: data.frame(dist, purity, ploidy, ...) sorted by dist.
prefer_undoubled <- function(minima, doubling_ratio = 1.5) {
    best  <- minima[1, ]
    lower <- minima[minima$ploidy <= best$ploidy / 1.6 & minima$dist <= doubling_ratio * best$dist, , drop = FALSE]
    if (nrow(lower)) best <- lower[1, ]
    other <- minima[!(minima$purity == best$purity & minima$ploidy == best$ploidy) &
                    minima$ploidy < 1.6 * best$ploidy, , drop = FALSE]
    list(best = best, runner_up = if (nrow(other)) other[1, ] else NULL)
}

# Whether the tool's own fit is replaced by the best feasible both-allele minimum:
#   (1) the own fit is below the purity floor (infeasible), or
#   (2) it is a doubling (ploidy >= 1.6x) of the best without being >= doubling_ratio x better, or
#   (3) it scores >= 1.25x worse than the best, the best beats its runner-up by >= 1.2x (a clear answer), and
#       purity differs by >= 0.05 or ploidy by >= 0.2.
refit_decision <- function(own_purity, own_ploidy, own_dist, own_feasible, best, runner_up, doubling_ratio = 1.5) {
    if (is.null(best)) return(list(refit = FALSE, why = "no feasible minimum"))
    if (!own_feasible) return(list(refit = TRUE, why = "own fit below the purity floor"))
    if (own_ploidy >= 1.6 * best$ploidy && !(own_dist * doubling_ratio <= best$dist))
        return(list(refit = TRUE, why = "own fit is a doubling of the best"))
    vs_own    <- own_dist / best$dist
    vs_runner <- if (is.null(runner_up)) Inf else runner_up$dist / best$dist
    differs   <- abs(best$purity - own_purity) >= 0.05 || abs(best$ploidy - own_ploidy) >= 0.2
    if (vs_own >= 1.25 && vs_runner >= 1.2 && differs) return(list(refit = TRUE, why = "both-allele rule"))
    list(refit = FALSE, why = "own fit kept")
}
