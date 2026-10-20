#!/usr/bin/env Rscript

# Prepare public-use Oregon Health Insurance Experiment samples.
#
# The primary compatibility-audit sample uses one scrambled-index adult per
# household from the administrative emergency-department cohort. This avoids
# survey-response weights and within-household dependence while preserving the
# randomized binary lottery instrument and the published Medicaid treatment.
# Only the paper's emergency-department application is prepared.

if (!requireNamespace("haven", quietly = TRUE)) {
  stop("Install the recommended haven package before preparing OHIE data.")
}

data_root <- file.path("data", "OHIE_Public_Use_Files", "OHIE_Data")
required <- file.path(data_root, c(
  "oregonhie_ed_vars.dta",
  "oregonhie_descriptive_vars.dta",
  "oregonhie_stateprograms_vars.dta"
))
if (!all(file.exists(required))) {
  stop(
    "OHIE public-use files were not found. Extract OHIE_Public_Use_Files ",
    "under data/ before running this script."
  )
}

read_numeric_dta <- function(path) {
  output <- as.data.frame(haven::read_dta(path))
  for (name in names(output)) {
    if (inherits(output[[name]], "haven_labelled")) {
      output[[name]] <- as.numeric(output[[name]])
    }
  }
  output
}

attach_by_person_id <- function(target, source, columns) {
  if (anyDuplicated(source$person_id)) {
    stop("Source data contain duplicate person_id values.")
  }
  index <- match(target$person_id, source$person_id)
  if (anyNA(index)) stop("A person_id could not be matched across OHIE files.")
  for (name in columns) target[[name]] <- source[[name]][index]
  target
}

make_analysis_covariates <- function(data, include_baseline_ed = TRUE) {
  output <- data.frame(
    age_centered = (2008 - data$birthyear_list - 45) / 10,
    female = data$female_list,
    stringsAsFactors = FALSE
  )
  if (include_baseline_ed) {
    output$baseline_ed_use <- data$any_visit_pre_ed
  }
  household_design <- stats::model.matrix(
    ~ factor(numhh_list), data = data
  )[, -1L, drop = FALSE]
  lottery_design <- stats::model.matrix(
    ~ factor(draw_lottery), data = data
  )[, -1L, drop = FALSE]
  colnames(household_design) <- paste0(
    "household_size_", sort(unique(data$numhh_list))[-1L]
  )
  colnames(lottery_design) <- paste0(
    "lottery_draw_", sort(unique(data$draw_lottery))[-1L]
  )
  output <- cbind(output, household_design, lottery_design)
  names(output) <- make.names(names(output), unique = TRUE)
  output
}

descriptive <- read_numeric_dta(required[2L])
state_programs <- read_numeric_dta(required[3L])

descriptive_columns <- c(
  "household_id", "treatment", "draw_lottery", "numhh_list",
  "birthyear_list", "female_list", "english_list", "self_list",
  "first_day_list", "have_phone_list", "pobox_list", "zip_msa_list"
)
treatment_column <- "ohp_all_ever_firstn_30sep2009"

ed <- read_numeric_dta(required[1L])
ed <- attach_by_person_id(ed, descriptive, descriptive_columns)
ed <- attach_by_person_id(ed, state_programs, treatment_column)
if (any(ed$sample_ed != 1)) stop("The public ED file should contain sample_ed=1.")

ed_covariates <- make_analysis_covariates(ed)
ed_analysis <- cbind(
  data.frame(
    person_id = ed$person_id,
    household_id = ed$household_id,
    Y = as.numeric(ed$any_visit_ed),
    T_obs = as.numeric(ed[[treatment_column]]),
    Z = as.numeric(ed$treatment),
    survey_weight = 1,
    stringsAsFactors = FALSE
  ),
  ed_covariates
)
ed_x_cols <- names(ed_covariates)
ed_published_x_cols <- c(
  "baseline_ed_use", "household_size_2", "household_size_3"
)
if (!all(ed_published_x_cols %in% ed_x_cols)) {
  stop("The published-aligned ED covariates were not constructed.")
}
if (any(!stats::complete.cases(ed_analysis))) {
  stop("Unexpected missing values in the prepared ED analysis sample.")
}

# person_id is scrambled in the public-use data. Selecting the smallest
# scrambled identifier within each household is therefore a deterministic,
# outcome-blind index-person rule.
household_order <- order(ed_analysis$household_id, ed_analysis$person_id)
household_first <- !duplicated(ed_analysis$household_id[household_order])
ed_household_index <- ed_analysis[household_order[household_first], , drop = FALSE]
rownames(ed_household_index) <- NULL

# Only 18 index households have three lottery-list members, and the public
# sample contains no treated lottery-control observation in that stratum.
# Retain the full index sample for descriptive checks, but define the primary
# overlap sample on household sizes one and two to avoid extrapolating an
# unobserved (T=1,Z=0,X) cell.
ed_household_index_overlap <- ed_household_index[
  ed_household_index$household_size_3 == 0, , drop = FALSE
]
rownames(ed_household_index_overlap) <- NULL
ed_overlap_published_x_cols <- c("baseline_ed_use", "household_size_2")
ed_overlap_expanded_x_cols <- setdiff(ed_x_cols, "household_size_3")

summarize_sample <- function(data, name) {
  cell <- as.data.frame.matrix(table(data$T_obs, data$Z))
  names(cell) <- paste0("Z", names(cell))
  cell$T_obs <- as.integer(rownames(cell))
  rownames(cell) <- NULL
  first_stage <- with(
    data,
    mean(T_obs[Z == 1]) - mean(T_obs[Z == 0])
  )
  data.frame(
    sample = name,
    observations = nrow(data),
    households = length(unique(data$household_id)),
    instrument_probability = mean(data$Z),
    treatment_probability_z0 = mean(data$T_obs[data$Z == 0]),
    treatment_probability_z1 = mean(data$T_obs[data$Z == 1]),
    first_stage = first_stage,
    outcome_mean = mean(data$Y),
    minimum_treatment_instrument_cell = min(table(data$T_obs, data$Z)),
    stringsAsFactors = FALSE
  )
}

sample_summary <- rbind(
  summarize_sample(ed_analysis, "ed_individual"),
  summarize_sample(ed_household_index, "ed_household_index"),
  summarize_sample(
    ed_household_index_overlap, "ed_household_index_overlap"
  )
)

prepared <- list(
  ed_individual = ed_analysis,
  ed_household_index = ed_household_index,
  ed_household_index_overlap = ed_household_index_overlap,
  ed_x_cols = ed_x_cols,
  ed_published_x_cols = ed_published_x_cols,
  ed_expanded_x_cols = ed_x_cols,
  ed_overlap_published_x_cols = ed_overlap_published_x_cols,
  ed_overlap_expanded_x_cols = ed_overlap_expanded_x_cols,
  sample_summary = sample_summary,
  specification = list(
    instrument = "Random selection in the Oregon Medicaid lottery",
    treatment_ed = paste(
      "Any Medicaid coverage from first lottery notification through",
      "30 September 2009"
    ),
    outcome_ed = "Any emergency-department visit during follow-up",
    household_index_rule = paste(
      "Smallest scrambled public-use person_id within household;",
      "deterministic and outcome-blind"
    ),
    overlap_rule = paste(
      "Primary ED analysis excludes household-size-three index households",
      "because the public sample has no T=1,Z=0 observation in that stratum"
    ),
    source = paste(
      "Finkelstein et al. Oregon Health Insurance Experiment public-use",
      "files, Harvard Dataverse DOI 10.7910/DVN/SJG1ED"
    )
  )
)

saveRDS(prepared, file.path("data", "ohie_application.rds"))
cat("Prepared Oregon Health Insurance Experiment samples.\n")
print(sample_summary, row.names = FALSE)
cat("Saved data/ohie_application.rds\n")
