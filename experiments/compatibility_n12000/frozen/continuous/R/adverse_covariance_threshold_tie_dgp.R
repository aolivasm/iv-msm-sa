# Fixed observed-data laws selected for adverse five-fold covariance, with the
# default configuration using Lambda=3.5.
#
# The treated-arm upper threshold has four exact population minimizers over two
# equally likely X strata.  Their fixed-threshold influence functions have very
# different variances and asymmetric dependence with the threshold-selection
# loss.  This is the geometry that invalidates the ordinary five-fold Wald
# approximation while leaving predictable online evaluation well behaved.  The
# control outcome is fixed at zero, so the same nonregular endpoint is also the
# upper ATE endpoint: ATE_upper = mu1_upper - mu0_lower = mu1_upper.

if (!exists("make_multibranch_nonunique_threshold_spec", mode = "function")) {
  if (!file.exists("R/multibranch_nonunique_threshold_pilot.R")) {
    stop("Run this file from the repository root.")
  }
  source("R/multibranch_nonunique_threshold_pilot.R")
}
if (!exists("observed_dgp_extreme_weighted_mean", mode = "function")) {
  if (!file.exists("R/observed_data_dgp.R")) {
    stop("Run this file from the repository root.")
  }
  source("R/observed_data_dgp.R")
}

make_adverse_covariance_threshold_tie_spec <- function(
    candidate_design = c("complete_binary", "hadamard"),
    configuration_id = suppressWarnings(as.integer(Sys.getenv(
      "IVMSM_ADVERSE_COV_CONFIGURATION_ID", unset = "14487"
    )))
) {
  candidate_design <- match.arg(candidate_design)
  catalogue <- list(
    `1144` = list(
      Lambda = 3.5, strata = 3L,
      e = 0.389451395068318, pi0 = 0.501296479478478,
      pi1 = 0.710358002372019,
      atoms = c(0, 1, 2.85221378772985, 6.24160182168707),
      gamma = c(
        2.15935538802738, 1.15121155801871,
        0.419983624565659, 2.20333476743576
      ),
      fractions = c(0.000512291975821601, 0.999878572894695),
      null_fraction = 0.832705110916868
    ),
    `11655` = list(
      Lambda = 4, e = 0.635513472240418, pi0 = 0.0669972554175183,
      pi1 = 0.606139736505151,
      atoms = c(0, 0.0117977294331061, 0.223488204544330,
        1.90404303560930),
      gamma = c(1.84486313365128, 0.676467054091927,
        2.83783326332806, 1.94836059982952),
      fractions = c(0.0207724750882915, 0.999971735068110),
      null_fraction = 0.622217343167867
    ),
    `971` = list(
      Lambda = 3.5, strata = 3L,
      e = 0.768791988072917, pi0 = 0.178365446380340,
      pi1 = 0.436549376990405,
      atoms = c(0, 1, 1.50024225017987, 4.56123474405613),
      gamma = c(2.16351713785043, 2.59964473819916,
        0.784903666828766, 1.59152453826679),
      fractions = c(0.0196749510463111, 0.998724988319889),
      null_fraction = 0.472544091916643
    ),
    `20492` = list(
      Lambda = 3, e = 0.559098705779761, pi0 = 0.314138201158494,
      pi1 = 0.612588975088540,
      atoms = c(0, 0.329235053747544, 118.434109425010),
      gamma = c(1.11054825828829, 0.374135834803041,
        1.78863504558396),
      fractions = c(0.0512709565237187, 0.999993684606644),
      null_fraction = 0.702145663781557
    ),
    `30492` = list(
      Lambda = 3.5, e = 0.559098705779761, pi0 = 0.314138201158494,
      pi1 = 0.612588975088540,
      atoms = c(0, 0.329235053747544, 118.434109425010),
      gamma = c(1.11054825828829, 0.374135834803041,
        1.78863504558396),
      fractions = c(0.0512709565237187, 0.999993684606644),
      null_fraction = 0.702145663781557
    ),
    `40492` = list(
      Lambda = 4, e = 0.559098705779761, pi0 = 0.314138201158494,
      pi1 = 0.612588975088540,
      atoms = c(0, 0.329235053747544, 118.434109425010),
      gamma = c(1.11054825828829, 0.374135834803041,
        1.78863504558396),
      fractions = c(0.0512709565237187, 0.999993684606644),
      null_fraction = 0.702145663781557
    ),
    `18620` = list(
      Lambda = 4, e = 0.282695777285844, pi0 = 0.265567196602933,
      pi1 = 0.754994846737605,
      atoms = c(0, 0.674021928574480, 36.8852725474241),
      gamma = c(1.65287757317697, 0.760525178983569,
        2.00153599967629),
      fractions = c(0.0827399610867970, 0.998315706505305),
      null_fraction = 0.695812772957142
    ),
    `8667` = list(
      Lambda = 4, e = 0.174134266050532, pi0 = 0.287599009377882,
      pi1 = 0.907229816214619,
      atoms = c(0, 0.0826258499987282, 0.0883225287058485, 1),
      gamma = c(0.716714478236678, 1.04227685847079,
        0.844887190128046, 2.46969669865268),
      fractions = c(0.0010704902272184, 0.997839596841622),
      null_fraction = 0.363502326402813
    ),
    `1537` = list(
      Lambda = 3.5, e = 0.608079957943410, pi0 = 0.280753316930495,
      pi1 = 0.679673488630901,
      atoms = c(0, 1.97045687356810, 2.04104644042604, 57.6629726784295),
      gamma = c(0.859395610269433, 1.78896516423814,
        0.726687630760141, 1.94741577437939),
      fractions = c(5.27565666179695e-6, 0.999231308978455),
      null_fraction = 0.618522952666972
    ),
    `14487` = list(
      Lambda = 3.5, e = 0.800280243884772, pi0 = 0.253038700916804,
      pi1 = 0.715484388322968,
      atoms = c(0, 0.0157288646757377, 0.0299026692293297,
        20.1277539604517),
      gamma = c(0.393433416093881, 0.501535690018185,
        2.47275240608310, 1.97691004454486),
      fractions = c(4.67730548635991e-7, 0.999440143381006),
      null_fraction = 0.445036136824638
    ),
    `13375` = list(
      Lambda = 3.5, e = 0.477890119915828, pi0 = 0.404060165891424,
      pi1 = 0.810259740226338,
      atoms = c(0, 0.0164718239856142, 25.5484372716504),
      gamma = c(0.953300367684136, 0.435998935107832,
        2.24516898537874),
      fractions = c(6.20236147519471e-6, 0.999910133800562),
      null_fraction = 0.643117559654638
    ),
    `10492` = list(
      Lambda = 2.5, e = 0.559098705779761, pi0 = 0.314138201158494,
      pi1 = 0.612588975088540,
      atoms = c(0, 0.329235053747544, 118.434109425010),
      gamma = c(1.11054825828829, 0.374135834803041,
        1.78863504558396),
      fractions = c(0.0512709565237187, 0.999993684606644),
      null_fraction = 0.702145663781557
    ),
    `10389` = list(
      Lambda = 4, e = 0.538917661104351, pi0 = 0.455785786248743,
      pi1 = 0.677691581434677,
      atoms = c(0, 0.0195820510927763, 0.618220232198463),
      gamma = c(1.07987713540284, 0.373713502156947,
        1.36069360220718),
      fractions = c(0.000194600925511811, 0.995345629328569),
      null_fraction = 0.623213411073666
    ),
    `15954` = list(
      Lambda = 4, e = 0.744856553077698, pi0 = 0.351564733651467,
      pi1 = 0.463356969035806,
      atoms = c(0, 0.155231302099262, 4.43948557337045),
      gamma = c(1.64896744448929, 0.209104138019235,
        2.27422164515145),
      fractions = c(0.000922074166865379, 0.999721778345647),
      null_fraction = 0.199170866776258
    ),
    `8958` = list(
      Lambda = 3, e = 0.601130836801604, pi0 = 0.565540966694243,
      pi1 = 0.734655478531323,
      atoms = c(0, 0.233627980921659, 25.5457152893430,
        47.7748384205613),
      gamma = c(1.01221718775133, 0.629430860395054,
        0.590987608005888, 1.56652279142436),
      fractions = c(1.30746538480349e-7, 0.999875476495650),
      null_fraction = 0.223013437760528
    ),
    `9148` = list(
      Lambda = 4, e = 0.811737273614854, pi0 = 0.215852725887671,
      pi1 = 0.319004598292072,
      atoms = c(0, 822.803267862514, 827.069122840900),
      gamma = c(1.53945769440923, 1.36806521617292,
        0.782819492703632),
      fractions = c(0.0261512850845996, 0.999993786492078),
      null_fraction = 0.308993097508792,
      endpoint_direction = "lower", gap_index = 1L
    ),
    `11790` = list(
      Lambda = 4, e = 0.717287427913398, pi0 = 0.524047815464437,
      pi1 = 0.786961039246152,
      atoms = c(0, 0.340794917965102, 0.439328310357037),
      gamma = c(2.09450453205313, 0.614342907913311,
        1.55495841266777),
      fractions = c(9.92261536812348e-7, 0.899485744726691),
      null_fraction = 0.456326987245120,
      endpoint_direction = "lower", gap_index = 1L
    ),
    `12470` = list(
      Lambda = 4, e = 0.484756564460695, pi0 = 0.291715606246144,
      pi1 = 0.736452845353078,
      atoms = c(0, 773.076863125500, 775.712511014256),
      gamma = c(3.60231425508041, 3.41037830811707,
        0.320590762899988),
      fractions = c(0.0339680801900723, 0.999995432009725),
      null_fraction = 0.702055590383243,
      endpoint_direction = "lower", gap_index = 1L
    ),
    `1061` = list(
      Lambda = 2.5, e = 0.287903965292498, pi0 = 0.362189917601645,
      pi1 = 0.611924540118663,
      atoms = c(0, 599.921690757475, 668.561707463471),
      gamma = c(1.05802288975381, 0.807750166444482,
        1.80530630854555),
      fractions = c(0.0332472672731689, 0.954458578818244),
      null_fraction = 0.988462425635662,
      endpoint_direction = "lower", gap_index = 1L
    ),
    `15953` = list(
      Lambda = 3.5, e = 0.502624240526929, pi0 = 0.527683563251048,
      pi1 = 0.777845878618369,
      atoms = c(0, 10.8293261611895, 193.349095258375),
      gamma = c(0.579629747213952, 0.667792764178576,
        2.08939865995656),
      fractions = c(1.83742260642238e-6, 0.959159286462502),
      null_fraction = 0.911565135382116,
      endpoint_direction = "lower", gap_index = 1L
    ),
    # Fixed Lambda=4 law selected by the optimizer-gradient search.  Its 32
    # Hadamard thresholds are exact upper-endpoint optimizers but have distinct
    # pathwise gradients.  Configuration 1068 was preferred to the slightly
    # more adverse search winner because its minimum observed (Z,T) cell mass
    # is 0.117, which avoids making the fully estimated comparison primarily a
    # sparse-cell stress test.
    `321068` = list(
      Lambda = 4, strata = 32L, candidates = 32L,
      candidate_design = "hadamard",
      e = 0.378092418250162,
      pi0 = 0.392418883326463,
      pi1 = 0.188324497706182,
      atoms = c(0, 0.0256805883562972, 0.965597801779809),
      gamma = c(
        0.327932703358195, 1.54985943830069, 0.852734925781708
      ),
      fractions = c(0.00810424480653663, 0.998093006562091),
      null_fraction = 0.139368460448459,
      endpoint_direction = "upper", gap_index = 2L
    ),
    # Estimable version of configuration 321068.  The same fixed observed law
    # is paired with 16 orthogonal threshold branches, leaving enough data per
    # stratum for the official n=1500 nuisance and EIF regressions.
    `161068` = list(
      Lambda = 4, strata = 16L, candidates = 16L,
      candidate_design = "hadamard",
      e = 0.378092418250162,
      pi0 = 0.392418883326463,
      pi1 = 0.188324497706182,
      atoms = c(0, 0.0256805883562972, 0.965597801779809),
      gamma = c(
        0.327932703358195, 1.54985943830069, 0.852734925781708
      ),
      fractions = c(0.00810424480653663, 0.998093006562091),
      null_fraction = 0.139368460448459,
      endpoint_direction = "upper", gap_index = 2L
    )
  )
  configuration <- catalogue[[as.character(configuration_id)]]
  if (is.null(configuration)) {
    stop(
      "Unsupported adverse covariance configuration_id. Available values: ",
      paste(names(catalogue), collapse = ", "), "."
    )
  }
  atoms <- configuration$atoms / max(configuration$atoms)
  strata <- if (is.null(configuration$strata)) 2L else configuration$strata
  candidates <- if (is.null(configuration$candidates)) {
    2L
  } else {
    configuration$candidates
  }
  selected_candidate_design <- if (is.null(configuration$candidate_design)) {
    candidate_design
  } else {
    configuration$candidate_design
  }
  specification <- make_multibranch_nonunique_threshold_spec(
    sensitivity_lambda = configuration$Lambda,
    strata = strata,
    candidates = candidates,
    candidate_design = selected_candidate_design,
    endpoint_direction = if (is.null(configuration$endpoint_direction)) {
      "upper"
    } else configuration$endpoint_direction,
    treatment = 1L,
    probability_z0 = configuration$e,
    pi0 = configuration$pi0,
    pi1 = configuration$pi1,
    atoms = atoms,
    gap_index = if (is.null(configuration$gap_index)) {
      length(atoms) - 1L
    } else configuration$gap_index,
    threshold_gap_fractions = configuration$fractions,
    gamma_atoms = configuration$gamma,
    probability_null_fraction = configuration$null_fraction,
    envelope_smoothing_epsilon = 0
  )
  specification$control_outcome_scale <- 0
  specification$search_configuration_id <- as.integer(configuration_id)
  specification
}

# Overlap-adjusted fixed law selected for the layered Lambda=4 experiment.
# It retains the atoms, density ratios, threshold gap, and instrument law of
# search configuration 11655, while setting P(T=1 | Z=0) to 0.10.  The oracle
# is rebuilt after this change, so its atom probabilities, envelopes, scores,
# and exact normalization identities all correspond to the modified fixed law.
make_fixed_lambda4_covariance_threshold_tie_spec <- function(
    candidate_design = c("complete_binary", "hadamard")
) {
  candidate_design <- match.arg(candidate_design)
  specification <- make_adverse_covariance_threshold_tie_spec(
    candidate_design = candidate_design,
    configuration_id = 11655L
  )
  specification$source_search_pi0 <- specification$pi0
  specification$pi0 <- 0.10
  specification$name <-
    "fixed_lambda4_covariance_threshold_tie_cfg11655_pi0p10"
  specification$overlap_adjusted <- TRUE
  specification
}

# Fixed Lambda=3.5 law selected to separate ordinary five-fold inference from
# predictable online inference without relying on highly unequal branch
# variances.  It has three X strata, eight exactly tied upper-threshold
# branches, ordinary overlap, and a population score-SD ratio near 1.35.
make_balanced_lambda3p5_covariance_threshold_tie_spec <- function(
    candidate_design = c("complete_binary", "hadamard")
) {
  candidate_design <- match.arg(candidate_design)
  specification <- make_adverse_covariance_threshold_tie_spec(
    candidate_design = candidate_design,
    configuration_id = 1144L
  )
  specification$name <-
    "balanced_lambda3p5_covariance_threshold_tie_cfg1144"
  specification$overlap_adjusted <- FALSE
  specification
}

build_adverse_covariance_threshold_tie_oracle <- function(
    specification = make_adverse_covariance_threshold_tie_spec()
) {
  build_multibranch_nonunique_threshold_oracle(specification)
}

adverse_covariance_x_columns <- function(specification) {
  paste0("X", seq_len(specification$strata - 1L))
}

adverse_covariance_add_x_columns <- function(data, specification) {
  for (stratum in seq_len(specification$strata - 1L)) {
    data[[paste0("X", stratum)]] <- as.numeric(data$X == stratum)
  }
  data
}

generate_adverse_covariance_threshold_tie_data <- function(
    n,
    specification = make_adverse_covariance_threshold_tie_spec(),
    oracle = NULL,
    seed = NULL
) {
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 1L) {
    stop("n must be one positive integer.")
  }
  if (is.null(oracle)) {
    oracle <- build_adverse_covariance_threshold_tie_oracle(specification)
  }
  if (!is.null(seed)) set.seed(as.integer(seed))
  x <- sample.int(specification$strata, n, replace = TRUE)
  z <- stats::rbinom(n, 1L, 1 - specification$probability_z0)
  treatment_probability <- ifelse(
    z == 0L, specification$pi0, specification$pi1
  )
  treatment <- stats::rbinom(n, 1L, treatment_probability)
  probability0 <- oracle$atom_probability
  probability1_treated <- probability0 * specification$gamma_atoms
  atom <- integer(n)
  group <- interaction(treatment, z, drop = TRUE)
  for (level in levels(group)) {
    rows <- which(group == level)
    treated_z1 <- treatment[rows[1L]] == 1L && z[rows[1L]] == 1L
    atom[rows] <- sample.int(
      length(specification$atoms), length(rows), replace = TRUE,
      prob = if (treated_z1) probability1_treated else probability0
    )
  }
  outcome_scale <- ifelse(
    treatment == 1L, 1, specification$control_outcome_scale
  )
  output <- data.frame(
    Y = outcome_scale * specification$atoms[atom],
    T_obs = treatment,
    Z = z,
    X = x,
    atom = atom,
    e_true = specification$probability_z0,
    pi0_true = specification$pi0,
    pi1_true = specification$pi1,
    gamma0_oracle = 1,
    gamma1_oracle = specification$gamma_atoms[atom],
    prob_Z = 1 - specification$probability_z0,
    prob_T = treatment_probability,
    stringsAsFactors = FALSE
  )
  output <- adverse_covariance_add_x_columns(output, specification)
  attr(output, "adverse_covariance_threshold_tie_specification") <-
    specification
  output
}

compute_adverse_covariance_dual_population <- function(
    specification = make_adverse_covariance_threshold_tie_spec(),
    oracle = NULL
) {
  if (is.null(oracle)) {
    oracle <- build_adverse_covariance_threshold_tie_oracle(specification)
  }
  probability0 <- oracle$atom_probability
  probability1 <- probability0 * specification$gamma_atoms
  if (abs(sum(probability0) - 1) > 1e-10 ||
      abs(sum(probability1) - 1) > 1e-10) {
    stop("The adverse conditional outcome laws are not normalized.")
  }
  arm_results <- matrix(
    NA_real_, nrow = 2L, ncol = 2L,
    dimnames = list(treatment = c("0", "1"), bound = c("lower", "upper"))
  )
  diagnostics <- list()
  for (treatment in 0:1) {
    gamma <- if (treatment == 1L) {
      specification$gamma_atoms
    } else {
      rep(1, length(probability0))
    }
    nuisance <- list(
      e = rep(specification$probability_z0, length(probability0)),
      pi0 = rep(specification$pi0, length(probability0)),
      pi1 = rep(specification$pi1, length(probability0))
    )
    envelope <- compute_projected_envelopes(
      nuisance = nuisance,
      gamma = gamma,
      t_value = treatment,
      Lambda0 = specification$sensitivity_lambda,
      Lambda1 = specification$sensitivity_lambda,
      envelope_smoothing_epsilon = 0,
      gamma_compatibility_action = "stop"
    )
    if (any(envelope$gamma_incompatible_original) ||
        any(envelope$L > envelope$U + 1e-12)) {
      stop("The adverse law is incompatible for one treatment arm.")
    }
    outcome_values <- specification$atoms * if (
      treatment == 1L
    ) 1 else specification$control_outcome_scale
    extreme <- observed_dgp_extreme_weighted_mean(
      outcome_values,
      probability0,
      envelope$L,
      envelope$U
    )
    observed_mean <- sum(probability0 * outcome_values)
    if (treatment == 1L) {
      arm_results["1", ] <- c(
        specification$pi0 * observed_mean +
          (1 - specification$pi0) * extreme$lower,
        specification$pi0 * observed_mean +
          (1 - specification$pi0) * extreme$upper
      )
    } else {
      arm_results["0", ] <- c(
        (1 - specification$pi0) * observed_mean +
          specification$pi0 * extreme$lower,
        (1 - specification$pi0) * observed_mean +
          specification$pi0 * extreme$upper
      )
    }
    diagnostics[[treatment + 1L]] <- data.frame(
      treatment = treatment,
      minimum_envelope_width = min(envelope$U - envelope$L),
      gamma_projection_rate = mean(envelope$gamma_projected),
      extremizer_normalization_error = extreme$normalization_error,
      observed_mean_z0 = observed_mean,
      lower_target = arm_results[as.character(treatment), "lower"],
      upper_target = arm_results[as.character(treatment), "upper"],
      stringsAsFactors = FALSE
    )
  }
  endpoints <- data.frame(
    design = "adverse_covariance_complete_binary_threshold_tie",
    Lambda = specification$sensitivity_lambda,
    mu1_lower = arm_results["1", "lower"],
    mu1_upper = arm_results["1", "upper"],
    mu0_lower = arm_results["0", "lower"],
    mu0_upper = arm_results["0", "upper"],
    stringsAsFactors = FALSE
  )
  endpoints$ate_lower <- endpoints$mu1_lower - endpoints$mu0_upper
  endpoints$ate_upper <- endpoints$mu1_upper - endpoints$mu0_lower
  tied_endpoint <- paste0("mu1_", specification$endpoint_direction)
  tied_ate_endpoint <- paste0("ate_", specification$endpoint_direction)
  if (abs(endpoints[[tied_endpoint]] - oracle$target) > 1e-10) {
    stop("The finite-law and exact-tie treated targets disagree.")
  }
  if (abs(endpoints[[tied_ate_endpoint]] -
      endpoints[[tied_endpoint]]) > 1e-12) {
    stop("The adverse law must identify the tied ATE and mu1 endpoints.")
  }
  list(
    endpoints = endpoints,
    diagnostics = do.call(rbind, diagnostics),
    compatible = TRUE,
    oracle_target = oracle$target
  )
}

compute_adverse_covariance_tan_population <- function(
    specification = make_adverse_covariance_threshold_tie_spec(),
    oracle = NULL
) {
  if (is.null(oracle)) {
    oracle <- build_adverse_covariance_threshold_tie_oracle(specification)
  }
  probability0 <- oracle$atom_probability
  probability1 <- probability0 * specification$gamma_atoms
  e <- specification$probability_z0
  pi0 <- specification$pi0
  pi1 <- specification$pi1
  marginal_pi <- e * pi0 + (1 - e) * pi1
  z0_given_treated <- e * pi0 / marginal_pi
  treated_marginal <-
    z0_given_treated * probability0 +
    (1 - z0_given_treated) * probability1
  control_atoms <- specification$control_outcome_scale * specification$atoms
  x_t1 <- observed_dgp_standard_msm_arm(
    specification$atoms, treated_marginal, marginal_pi, 1L,
    specification$sensitivity_lambda
  )
  x_t0 <- observed_dgp_standard_msm_arm(
    control_atoms, probability0, marginal_pi, 0L,
    specification$sensitivity_lambda
  )
  xz_t1 <-
    e * observed_dgp_standard_msm_arm(
      specification$atoms, probability0, pi0, 1L,
      specification$sensitivity_lambda
    ) +
    (1 - e) * observed_dgp_standard_msm_arm(
      specification$atoms, probability1, pi1, 1L,
      specification$sensitivity_lambda
    )
  xz_t0 <-
    e * observed_dgp_standard_msm_arm(
      control_atoms, probability0, pi0, 0L,
      specification$sensitivity_lambda
    ) +
    (1 - e) * observed_dgp_standard_msm_arm(
      control_atoms, probability0, pi1, 0L,
      specification$sensitivity_lambda
    )
  make_row <- function(estimator, arm0, arm1) {
    output <- data.frame(
      estimator = estimator,
      Lambda = specification$sensitivity_lambda,
      mu1_lower = arm1[["lower"]],
      mu1_upper = arm1[["upper"]],
      mu0_lower = arm0[["lower"]],
      mu0_upper = arm0[["upper"]],
      stringsAsFactors = FALSE
    )
    output$ate_lower <- output$mu1_lower - output$mu0_upper
    output$ate_upper <- output$mu1_upper - output$mu0_lower
    output
  }
  rbind(
    make_row("tan_cal_x_saturated", x_t0, x_t1),
    make_row("tan_cal_xz_saturated", xz_t0, xz_t1)
  )
}

audit_adverse_covariance_threshold_tie <- function(
    specification = make_adverse_covariance_threshold_tie_spec(),
    oracle = NULL
) {
  if (is.null(oracle)) {
    oracle <- build_adverse_covariance_threshold_tie_oracle(specification)
  }
  data.frame(
    Lambda = specification$sensitivity_lambda,
    strata = specification$strata,
    candidates = specification$candidates,
    candidate_design = specification$candidate_design,
    treatment = specification$treatment,
    endpoint_direction = specification$endpoint_direction,
    control_outcome_scale = specification$control_outcome_scale,
    search_configuration_id = specification$search_configuration_id,
    source_search_pi0 = if (is.null(specification$source_search_pi0)) {
      NA_real_
    } else {
      specification$source_search_pi0
    },
    overlap_adjusted = isTRUE(specification$overlap_adjusted),
    probability_z0 = specification$probability_z0,
    pi0 = specification$pi0,
    pi1 = specification$pi1,
    minimum_xzt_probability = min(c(
      specification$probability_z0 * specification$pi0,
      specification$probability_z0 * (1 - specification$pi0),
      (1 - specification$probability_z0) * specification$pi1,
      (1 - specification$probability_z0) * (1 - specification$pi1)
    )) / specification$strata,
    minimum_atom_probability = min(oracle$atom_probability),
    maximum_atom_probability = max(oracle$atom_probability),
    target = oracle$target,
    maximum_target_range = diff(range(oracle$direct_target)),
    maximum_risk_range = diff(range(oracle$population_risk)),
    maximum_normalization_residual = max(abs(
      oracle$direct_normalization_residual
    )),
    minimum_score_standard_deviation = min(oracle$population_sd),
    maximum_score_standard_deviation = max(oracle$population_sd),
    score_standard_deviation_ratio =
      max(oracle$population_sd) / min(oracle$population_sd),
    minimum_score_correlation = min(
      oracle$population_correlation[upper.tri(
        oracle$population_correlation
      )]
    ),
    stringsAsFactors = FALSE
  )
}
