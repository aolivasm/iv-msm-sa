# Main-text reproduction: required functions only.
siv_gamma <- function(r, eta) {
    exp(eta * r - eta^2/2) * siv_normalizer(0)/siv_normalizer(eta)
}

siv_gamma_inverse <- function(value, eta) {
    (log(value) + eta^2/2 - log(siv_normalizer(0)/siv_normalizer(eta)))/eta
}

siv_normalizer <- function(eta) pnorm(2 - eta) - pnorm(-2 - eta)

x5g_generate <- function(n, a, b, seed, stress = FALSE) {
    if (!stress)
        return(x5s_generate(n, a, b, seed))
    set.seed(seed)
    x <- matrix(qnorm(pnorm(-2) + runif(n * 5) * siv_normalizer(0)), ncol = 5)
    colnames(x) <- paste0("X", 1:5)
    z <- rbinom(n, 1, 0.5)
    c <- x5s_contract()
    h <- drop(x %*% c$beta)
    t <- rbinom(n, 1, plogis(h + a * (z - 0.5)))
    shape <- 1 + 0.6 * sin(x[, 1] * x[, 2]) + 0.4 * sin(x[, 3] * x[, 4])
    eta <- z * (2 * t - 1) * b * shape
    r <- eta + qnorm(pnorm(-2 - eta) + runif(n) * siv_normalizer(eta))
    loc <- drop(x %*% c$theta)
    data.frame(Y = loc + 0.75 * t + r, T_obs = t, Z = z, x, e_true = 0.5, pi0_true = plogis(h - a/2),
        pi1_true = plogis(h + a/2), gamma0_oracle = siv_gamma(loc + 0.75 * t + r - loc, -b * shape),
        gamma1_oracle = siv_gamma(loc + 0.75 * t + r - loc - 0.75, b * shape))
}

x5s_contract <- function() list(beta = c(0.14, 0.12, 0.1, 0.08, 0.06), theta = c(0.2, -0.25, 0.3, 0.15,
    -0.2), shift = 0.75, x_cols = c("X1", "X2", "X3", "X4", "X5"))

x5s_generate <- function(n, a, b, seed, contract = x5s_contract()) {
    stopifnot(n >= 1, a >= 0, b >= 0, length(contract$beta) == 5, length(contract$theta) == 5)
    set.seed(seed)
    x <- matrix(qnorm(pnorm(-2) + runif(n * 5) * siv_normalizer(0)), ncol = 5)
    colnames(x) <- contract$x_cols
    z <- rbinom(n, 1, 0.5)
    h <- drop(x %*% contract$beta)
    treatment <- rbinom(n, 1, plogis(h + a * (z - 0.5)))
    eta <- z * (2 * treatment - 1) * b
    residual <- eta + qnorm(pnorm(-2 - eta) + runif(n) * siv_normalizer(eta))
    d <- data.frame(Y = drop(x %*% contract$theta) + contract$shift * treatment + residual, T_obs = treatment,
        Z = z, x)
    nu <- x5s_oracle(d, a, b, contract)
    d$e_true <- nu$e
    d$pi0_true <- nu$pi0
    d$pi1_true <- nu$pi1
    d$gamma0_oracle <- nu$gamma0
    d$gamma1_oracle <- nu$gamma1
    d
}

