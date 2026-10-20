# Main-text reproduction dependencies; published formulas retained.
x5p_arm <- function(x, a, b, L, t, class) {
    contract <- x5s_contract()
    n <- nrow(x)
    h <- drop(as.matrix(x) %*% contract$beta)
    loc <- drop(as.matrix(x) %*% contract$theta) + 0.75 * t
    p0 <- plogis(h - a/2)
    p1 <- plogis(h + a/2)
    if (t == 0) {
        p0 <- 1 - p0
        p1 <- 1 - p1
    }
    eta <- (2 * t - 1) * b
    d <- -p0/(1 - p0)
    cl <- ((1 - p1)/L + p1)/(1 - p0)
    cu <- ((1 - p1) * L + p1)/(1 - p0)
    lower <- function(hi = 2) x5p_integral(cl, d, 1/L, eta, TRUE, hi = hi)
    upper <- function(hi = 2) x5p_integral(cu, d, L, eta, FALSE, hi = hi)
    lm <- lower()
    um <- upper()
    total <- um[, 1] - lm[, 1]
    stopifnot(all(lm[, 1] <= 1 + 1e-09), all(um[, 1] >= 1 - 1e-09), all(total > 0))
    LY <- lm[, 2] + loc * lm[, 1]
    TY <- um[, 2] - lm[, 2] + loc * total
    raw <- cbind(1, restricted_x5_features(x, class))
    trans <- solve(chol(crossprod(raw)/n))
    design <- raw %*% trans
    weight <- (1 - p0)/n
    values <- function(q) {
        r <- pmax(-2, pmin(2, q - loc))
        v <- upper(r) - lower(r)
        list(mass = v[, 1], moment = v[, 2] + loc * v[, 1])
    }
    result <- list()
    for (bound in c("lower", "upper")) {
        target <- if (bound == "lower")
            1 - lm[, 1]
        else um[, 1] - 1
        left <- rep(-2, n)
        right <- rep(2, n)
        for (k in 1:42) {
            mid <- (left + right)/2
            mass <- (upper(mid) - lower(mid))[, 1]
            right[mass >= target] <- mid[mass >= target]
            left[mass < target] <- mid[mass < target]
        }
        point <- loc + (left + right)/2
        beta <- lm.wfit(design, point, weight * total)$coefficients
        objective <- function(beta) {
            q <- drop(design %*% beta)
            v <- values(q)
            sum(weight * (q * (v$mass - target) - v$moment))
        }
        gradient <- function(beta) {
            q <- drop(design %*% beta)
            v <- values(q)
            drop(crossprod(design, weight * (v$mass - target)))
        }
        for (iter in 1:60) {
            g <- gradient(beta)
            if (max(abs(g)) < 1e-09)
                break
            q <- drop(design %*% beta)
            r <- q - loc
            gam <- siv_gamma(r, eta)
            gap <- pmin(L, cu * gam + d) - pmax(1/L, cl * gam + d)
            curvature <- ifelse(abs(r) < 2, gap * dnorm(r)/siv_normalizer(0), 0)
            H <- crossprod(design, (weight * curvature) * design)
            step <- solve(H, g)
            scale <- 1
            old <- objective(beta)
            while (scale > 1e-10 && objective(beta - scale * step) > old + 1e-14) scale <- scale/2
            beta <- beta - scale * step
        }
        grad <- max(abs(gradient(beta)))
        stopifnot(is.finite(grad), grad < 1e-08)
        endpoint <- function(q) {
            v <- values(q)
            missing <- if (bound == "lower")
                LY + v$moment + q * (1 - lm[, 1] - v$mass)
            else LY + TY - v$moment + q * (1 - um[, 1] + v$mass)
            mean(p0 * loc + (1 - p0) * missing)
        }
        q <- drop(design %*% beta)
        v <- values(q)
        norm <- if (bound == "lower")
            lm[, 1] + v$mass - 1
        else um[, 1] - v$mass - 1
        value <- endpoint(q)
        sharp <- endpoint(point)
        stopifnot(if (bound == "lower")
            value <= sharp + 1e-09
        else value >= sharp - 1e-09)
        result[[bound]] <- list(endpoint = paste0("mu", t, "_", bound), target = value, unrestricted = sharp,
            gradient = grad, iterations = iter, direct_minus_pseudo = sum(weight * q * norm), mean_abs_conditional_normalization_error = mean(abs(norm)),
            coefficients = drop(trans %*% beta))
    }
    result
}

x5p_integral <- function(c, d, k, eta, maximum, lo = -2, hi = 2) {
    n <- length(c)
    lo <- rep_len(lo, n)
    hi <- rep_len(hi, n)
    root <- siv_gamma_inverse(pmax((k - d)/c, .Machine$double.xmin), eta)
    cut <- pmax(lo, pmin(hi, root))
    piece <- function(left, right) {
        affine <- c * siv_gamma((left + right)/2, eta) + d
        choose <- if (maximum)
            affine > k
        else affine < k
        m0 <- x5p_moment(left, right, 0)
        m1 <- x5p_moment(left, right, eta)
        ifelse(choose, 1, 0) * (c * m1 + d * m0) + ifelse(choose, 0, 1) * k * m0
    }
    piece(lo, cut) + piece(cut, hi)
}

x5p_moment <- function(lo, hi, eta) {
    mass <- (pnorm(hi - eta) - pnorm(lo - eta))/siv_normalizer(eta)
    cbind(mass = mass, moment = eta * mass + (dnorm(lo - eta) - dnorm(hi - eta))/siv_normalizer(eta))
}

x5p_rule <- function(n, shift = 0) {
    stopifnot(n%%2 == 0)
    u <- vapply(first_primes(5), function(base) radical_inverse(100 + seq_len(n/2), base), numeric(n/2))
    if (shift != 0) {
        set.seed(shift)
        u <- sweep(u, 2, runif(5), "+")%%1
    }
    x <- qnorm(pnorm(-2) + rbind(u, 1 - u) * siv_normalizer(0))
    colnames(x) <- paste0("X", 1:5)
    as.data.frame(x)
}

