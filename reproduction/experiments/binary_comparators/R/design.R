source('paper_config.R')
config <- paper_profile$config
tasks <- data.frame(id = 1:2000, setting = rep(1:4, each = 500), replication = rep(1:500, 4))
generate <- function(t, n = config$n) {
    set.seed(t$seed)
    x <- matrix(qnorm(pnorm(-2) + runif(n * 5) * (pnorm(2) - pnorm(-2))), ncol = 5)
    colnames(x) <- paste0("X", 1:5)
    z <- rbinom(n, 1, 0.5)
    g <- config$settings[t$setting, ]
    h <- drop(x %*% c(0.14, 0.12, 0.1, 0.08, 0.06))
    tr <- rbinom(n, 1, plogis(h + g$a * (z - 0.5)))
    shape <- 1 + 0.6 * sin(x[, 1] * x[, 2]) + 0.4 * sin(x[, 3] * x[, 4])
    p <- plogis(-0.4 + 1.5 * tr + drop(x %*% c(0.2, -0.25, 0.3, 0.15, -0.2)) + (2 * tr - 1) * g$b * shape * z)
    data.frame(Y = rbinom(n, 1, p), T_obs = tr, Z = z, x)
}
atomic <- function(z, p) {
    dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile("unpublished_", dirname(p))
    on.exit(if (file.exists(tmp)) unlink(tmp))
    withCallingHandlers(saveRDS(z, tmp), warning = function(w) stop(conditionMessage(w)))
    stopifnot(identical(z, readRDS(tmp)))
    stopifnot(file.rename(tmp, p), identical(z, readRDS(p)))
}
tasks$seed<-261000L+tasks$setting+100L*(tasks$replication-1L)
