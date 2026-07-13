## -----------------------------------------------------------------------------
## PsyMetRiC cost-effectiveness study
## Reviewer-facing implementation of the patient-level expected-value Markov model
##
## This file contains the model inputs, health-state transition logic,
## strategy effects, costing/QALY calculations, deterministic base-case analysis,
## and deterministic scenario analyses. 
##
## Version date: 13 July 2026
##
## Author: Benjamin Perry
## -----------------------------------------------------------------------------

library(tidyverse)
library(here)


# ==============================================================================
# 1. Load model inputs
# ==============================================================================

analysis_datasets <- readRDS(
  here::here(
    "overall_results", "analysis_datasets",
    "analysis_datasets_with_psymetric_thresholds.rds"
  )
)

annual_transition_grids <- readRDS(
  here::here(
    "overall_results", "incidence_models",
    "annual_transition_probability_grids",
    "primary_full_annual_transition_probability_grids_all_outcomes.rds"
  )
)

t2d_complication_grid <- readRDS(
  here::here(
    "overall_results", "incidence_models", "t2d_complications",
    "t2d_complications_annual_transition_probabilities_by_duration.rds"
  )
)

external_parameter_lookup <- readRDS(
  here::here("overall_results", "external_inputs", "external_parameter_lookup.rds")
)

all_cause_mortality_values <- readRDS(
  here::here("overall_results", "external_inputs", "all_cause_mortality_values.rds")
)

mortality_cause_allocation_values <- readRDS(
  here::here("overall_results", "external_inputs", "mortality_cause_allocation_values.rds")
)

health_state_utility_values <- readRDS(
  here::here("overall_results", "external_inputs", "health_state_utility_values.rds")
)


# Small parameter accessors used by the analysis wrappers. 
parameter_value <- function(parameters, name) {
  if (is.null(names(parameters)) || !(name %in% names(parameters))) {
    stop("Missing external parameter: ", name, call. = FALSE)
  }
  value <- parameters[[name]]
  if (length(value) == 0 || all(is.na(value))) {
    stop("Missing external parameter value: ", name, call. = FALSE)
  }
  as.numeric(value)
}

parameter_value_default <- function(parameters, name, default = NA_real_) {
  if (is.null(names(parameters)) || !(name %in% names(parameters))) {
    return(as.numeric(default))
  }
  value <- parameters[[name]]
  if (length(value) == 0 || all(is.na(value))) {
    return(as.numeric(default))
  }
  as.numeric(value)
}


# ==============================================================================
# 2. Core model
# ==============================================================================

run_psymetric_markov <- function(
    person_data,
    transition_grids,
    t2d_complication_grid,
    parameters,
    mortality_values,
    mortality_cause_allocation_values,
    utility_values,
    strategy,
    threshold_value,
    include_operational_cost = FALSE
) {
  # Fixed structural assumptions from the manuscript/base case.
  external_parameter_lookup <- parameters
  state_probability_pruning_threshold <- 1e-10
  replace_post_cvd_chronic_utility_in_event_year <- TRUE
  replace_post_cvd_chronic_cost_in_event_year <- TRUE
  apply_acute_cvd_cost_to_all_cvd_deaths <- TRUE
  enable_intervention_bmi_reclassification <- TRUE
  include_psymetric_operational_cost_scenario <- include_operational_cost

  valid_strategies <- c(
    "usual_care",
    "universal_metformin",
    "psymetric_lifestyle",
    "psymetric_metformin",
    "psymetric_lifestyle_metformin",
    "psymetric_lifestyle_glp1ra",
    "psymetric_lifestyle_glp1ra_lifelong",
    "psymetric_antipsychotic_switch"
  )
  if (length(strategy) != 1L || !strategy %in% valid_strategies) {
    stop("Unknown model strategy: ", paste(strategy, collapse = ", "), call. = FALSE)
  }
  if (length(threshold_value) != 1L || is.na(threshold_value) ||
      threshold_value <= 0 || threshold_value >= 1) {
    stop("threshold_value must be a single probability between 0 and 1.", call. = FALSE)
  }

  get_param <- function(param_name, params = external_parameter_lookup) {
    if (is.null(names(params)) || !(param_name %in% names(params))) {
      stop(
        "Missing external parameter: ",
        param_name,
        ". This name is not present in the parameter lookup. ",
        "If it is intended to be optional, call get_param_default().",
        call. = FALSE
      )
    }
    
    value <- params[[param_name]]
    
    if (is.null(value) || length(value) == 0 || all(is.na(value))) {
      stop("Missing external parameter value: ", param_name, call. = FALSE)
    }
    
    as.numeric(value)
  }

  get_param_default <- function(param_name,
                                params = external_parameter_lookup,
                                default = NA_real_) {
    if (is.null(names(params)) || !(param_name %in% names(params))) {
      return(as.numeric(default))
    }
    
    value <- params[[param_name]]
    
    if (is.null(value) || length(value) == 0 || all(is.na(value))) {
      return(as.numeric(default))
    }
    
    as.numeric(value)
  }

  clip_probability <- function(x) {
    pmin(pmax(as.numeric(x), 0), 1)
  }

  prob_to_hazard <- function(p) {
    p <- clip_probability(p)
    dplyr::case_when(
      p >= 1 ~ Inf,
      TRUE ~ -log1p(-p)
    )
  }

  hazard_to_prob <- function(h) {
    dplyr::case_when(
      is.infinite(h) ~ 1,
      TRUE ~ 1 - exp(-h)
    ) |>
      clip_probability()
  }

  apply_hazard_ratio_to_probability <- function(p, hr) {
    hazard_to_prob(prob_to_hazard(p) * hr)
  }

  normalise_model_sex <- function(x) {
    dplyr::case_when(
      as.character(x) %in% c("1", "M", "m", "Male", "male") ~ "Male",
      as.character(x) %in% c("0", "F", "f", "Female", "female") ~ "Female",
      TRUE ~ NA_character_
    )
  }

  as_binary_01_model <- function(x) {
    dplyr::case_when(
      is.na(x) ~ NA_integer_,
      as.character(x) %in% c("1", "TRUE", "True", "true", "YES", "Yes", "yes", "Y", "y") ~ 1L,
      as.character(x) %in% c("0", "FALSE", "False", "false", "NO", "No", "no", "N", "n") ~ 0L,
      TRUE ~ suppressWarnings(as.integer(as.character(x)))
    )
  }

  is_death_state <- function(state) {
    state %in% c("Death_CVD", "Death_Non_CVD")
  }

  is_living_state <- function(state) {
    !is_death_state(state)
  }

  is_t2d_state <- function(state) {
    state %in% c(
      "T2D",
      "T2D_complications",
      "Post_CVD_T2D",
      "Post_CVD_T2D_complications"
    )
  }

  is_post_cvd_state <- function(state) {
    state %in% c(
      "Post_CVD",
      "Post_CVD_Obesity",
      "Post_CVD_T2D",
      "Post_CVD_T2D_complications"
    )
  }

  get_state_mortality_hr <- function(state,
                                     params = external_parameter_lookup) {
    t2d_uncomplicated_hr <- get_param_default(
      "rr_mortality_t2d_uncomplicated_base",
      params,
      default = 1
    )
    t2d_complications_hr <- get_param_default(
      "rr_mortality_t2d_complications_base",
      params,
      default = 1
    )
    post_cvd_hr <- get_param_default(
      "rr_mortality_post_cvd_base",
      params,
      default = 1
    )
    hr_cap <- get_param_default("state_mortality_hr_cap", params, default = Inf)
    
    hr <- rep(1, length(state))
    
    hr <- hr * dplyr::case_when(
      state %in% c("T2D", "Post_CVD_T2D") ~ t2d_uncomplicated_hr,
      state %in% c("T2D_complications", "Post_CVD_T2D_complications") ~ t2d_complications_hr,
      TRUE ~ 1
    )
    
    hr <- hr * dplyr::case_when(
      is_post_cvd_state(state) ~ post_cvd_hr,
      TRUE ~ 1
    )
    
    pmin(pmax(as.numeric(hr), 0), hr_cap)
  }

  get_state_cvd_death_proportion <- function(base_proportion,
                                             state,
                                             params = external_parameter_lookup) {
    post_cvd_multiplier <- get_param_default(
      "cvd_death_allocation_post_cvd_multiplier_base",
      params,
      default = 1
    )
    allocation_cap <- get_param_default("cvd_death_allocation_cap", params, default = 1)
    
    multiplier <- dplyr::case_when(
      is_post_cvd_state(state) ~ post_cvd_multiplier,
      TRUE ~ 1
    )
    
    clip_probability(pmin(as.numeric(base_proportion) * multiplier, allocation_cap))
  }

  state_to_utility_column <- function(state) {
    dplyr::case_when(
      state == "Healthy" ~ "utility_healthy",
      state == "Obesity" ~ "utility_obesity",
      state == "T2D" ~ "utility_t2d",
      state == "T2D_complications" ~ "utility_t2d_complications",
      state == "Post_CVD" ~ "utility_post_cvd",
      state == "Post_CVD_Obesity" ~ "utility_post_cvd_obesity",
      state == "Post_CVD_T2D" ~ "utility_post_cvd_t2d",
      state == "Post_CVD_T2D_complications" ~ "utility_post_cvd_t2d_complications",
      state == "Death_CVD" ~ "utility_death_cvd",
      state == "Death_Non_CVD" ~ "utility_death_non_cvd",
      TRUE ~ NA_character_
    )
  }

  state_to_event_year_cvd_utility_column <- function(state) {
    dplyr::case_when(
      state == "Healthy" ~ "utility_cvd_event_year_from_healthy",
      state == "Obesity" ~ "utility_cvd_event_year_from_obesity",
      state == "T2D" ~ "utility_cvd_event_year_from_t2d",
      state == "T2D_complications" ~ "utility_cvd_event_year_from_t2d_complications",
      TRUE ~ NA_character_
    )
  }

  chronic_state_cost <- function(state,
                                 state_probability,
                                 params = external_parameter_lookup) {
    cost <- rep(0, length(state))
    
    cost <- cost + dplyr::if_else(
      state %in% c("Obesity", "Post_CVD_Obesity"),
      state_probability * get_param("cost_obesity_annual", params),
      0
    )
    
    cost <- cost + dplyr::if_else(
      state %in% c("T2D", "Post_CVD_T2D"),
      state_probability * get_param("cost_t2d_uncomplicated_annual", params),
      0
    )
    
    cost <- cost + dplyr::if_else(
      state %in% c("T2D_complications", "Post_CVD_T2D_complications"),
      state_probability * (
        get_param("cost_t2d_uncomplicated_annual", params) +
          get_param("cost_t2d_complications_increment_annual", params)
      ),
      0
    )
    
    cost <- cost + dplyr::if_else(
      state %in% c(
        "Post_CVD",
        "Post_CVD_Obesity",
        "Post_CVD_T2D",
        "Post_CVD_T2D_complications"
      ),
      state_probability * get_param("cost_post_cvd_annual", params),
      0
    )
    
    cost
  }

  lookup_named_utility <- function(age,
                                   sex,
                                   utility_name,
                                   utility_table) {
    lookup_data <- tibble::tibble(
      row_id_lookup = seq_along(age),
      age = as.integer(floor(age)),
      sex = normalise_model_sex(sex),
      utility_name = utility_name
    ) |>
      dplyr::mutate(
        age = pmin(pmax(age, 16L), 100L)
      )
    
    utility_long <- utility_table |>
      dplyr::select(-dplyr::any_of("psa_draw")) |>
      tidyr::pivot_longer(
        cols = dplyr::starts_with("utility_"),
        names_to = "utility_name",
        values_to = "utility_value"
      )
    
    lookup_data |>
      dplyr::left_join(
        utility_long |>
          dplyr::select(age, sex, utility_name, utility_value),
        by = c("age", "sex", "utility_name")
      ) |>
      dplyr::arrange(row_id_lookup) |>
      dplyr::pull(utility_value)
  }

  prepare_mortality_cause_allocation_by_age <- function(mortality_cause_allocation_values) {
    mortality_cause_allocation_values |>
      dplyr::mutate(
        sex = dplyr::case_when(
          sex %in% c("male", "Male", "M", "1") ~ "Male",
          sex %in% c("female", "Female", "F", "0") ~ "Female",
          TRUE ~ as.character(sex)
        )
      ) |>
      dplyr::rowwise() |>
      dplyr::mutate(
        age = list(seq.int(as.integer(age_lower), as.integer(age_upper)))
      ) |>
      tidyr::unnest(age) |>
      dplyr::ungroup() |>
      dplyr::select(age, sex, proportion_deaths_cvd)
  }

  mortality_cause_allocation_by_age <-
    prepare_mortality_cause_allocation_by_age(mortality_cause_allocation_values)

  lookup_all_cause_death_probability <- function(age,
                                                 sex,
                                                 mortality_values) {
    lookup_data <- tibble::tibble(
      row_id_lookup = seq_along(age),
      age = as.integer(floor(age)),
      sex = normalise_model_sex(sex)
    ) |>
      dplyr::mutate(
        age = pmin(pmax(age, 16L), 100L)
      )
    
    lookup_data |>
      dplyr::left_join(
        mortality_values |>
          dplyr::select(-dplyr::any_of("psa_draw")) |>
          dplyr::select(age, sex, psychosis_adjusted_qx),
        by = c("age", "sex")
      ) |>
      dplyr::arrange(row_id_lookup) |>
      dplyr::pull(psychosis_adjusted_qx) |>
      clip_probability()
  }

  lookup_cvd_death_proportion <- function(age,
                                          sex,
                                          mortality_cause_by_age) {
    lookup_data <- tibble::tibble(
      row_id_lookup = seq_along(age),
      age = as.integer(floor(age)),
      sex = normalise_model_sex(sex)
    ) |>
      dplyr::mutate(
        age = pmin(pmax(age, 16L), 100L)
      )
    
    lookup_data |>
      dplyr::left_join(
        mortality_cause_by_age,
        by = c("age", "sex")
      ) |>
      dplyr::arrange(row_id_lookup) |>
      dplyr::pull(proportion_deaths_cvd) |>
      clip_probability()
  }

  lookup_transition_probability <- function(grid,
                                            ids,
                                            start_time,
                                            missing_to_zero = FALSE) {
    lookup <- tibble::tibble(
      row_id_lookup = seq_along(ids),
      ID = ids,
      start_time = as.numeric(start_time)
    )
    
    out <- lookup |>
      dplyr::left_join(
        grid |>
          dplyr::mutate(start_time = as.numeric(start_time)) |>
          dplyr::select(ID, start_time, interval_probability),
        by = c("ID", "start_time")
      ) |>
      dplyr::arrange(row_id_lookup) |>
      dplyr::pull(interval_probability)
    
    if (missing_to_zero) {
      out <- dplyr::coalesce(out, 0)
    }
    
    clip_probability(out)
  }

  lookup_t2d_complication_probability <- function(t2d_duration_years,
                                                  complication_grid,
                                                  hr_multiplier = 1) {
    duration_lookup <- floor(as.numeric(t2d_duration_years))
    
    max_duration_available <- max(
      complication_grid$t2d_duration_start_years,
      na.rm = TRUE
    )
    
    lookup <- tibble::tibble(
      row_id_lookup = seq_along(duration_lookup),
      t2d_duration_start_years = duration_lookup
    ) |>
      dplyr::mutate(
        t2d_duration_start_years = dplyr::case_when(
          is.na(t2d_duration_start_years) ~ NA_real_,
          t2d_duration_start_years < 0 ~ NA_real_,
          t2d_duration_start_years > max_duration_available ~
            as.numeric(max_duration_available),
          TRUE ~ as.numeric(t2d_duration_start_years)
        ),
        t2d_duration_start_years = as.integer(t2d_duration_start_years)
      )
    
    p <- lookup |>
      dplyr::left_join(
        complication_grid |>
          dplyr::select(t2d_duration_start_years, interval_probability),
        by = "t2d_duration_start_years"
      ) |>
      dplyr::arrange(row_id_lookup) |>
      dplyr::pull(interval_probability)
    
    apply_hazard_ratio_to_probability(
      p = dplyr::coalesce(p, 0),
      hr = hr_multiplier
    )
  }

  make_discount_factor <- function(rate, time) {
    1 / ((1 + rate) ^ time)
  }


  mortality_cause_allocation_by_age <-
    prepare_mortality_cause_allocation_by_age(mortality_cause_allocation_values)

  get_effect_weight <- function(cycle,
                                full_effect_duration_years,
                                taper_years) {
    cycle <- as.numeric(cycle)
    full_effect_duration_years <- as.numeric(full_effect_duration_years)
    taper_years <- as.numeric(taper_years)
    
    dplyr::case_when(
      cycle < full_effect_duration_years ~ 1,
      cycle >= full_effect_duration_years + taper_years ~ 0,
      taper_years <= 0 ~ 0,
      TRUE ~ 1 - ((cycle - full_effect_duration_years + 1) / taper_years)
    ) |>
      pmin(1) |>
      pmax(0)
  }

  get_threshold_eligibility <- function(person_data,
                                        threshold_value) {
    as.integer(person_data$PsyMetRiC2_T2D >= threshold_value)
  }

  get_glp1ra_eligibility <- function(person_data,
                                     threshold_value,
                                     params) {
    threshold_positive <- get_threshold_eligibility(person_data, threshold_value)
    bmi_threshold <- get_param("glp1ra_bmi_eligibility_threshold", params)
    baseline_bmi <- as.numeric(person_data$BASELINE_BMI)
    as.integer(threshold_positive == 1L & !is.na(baseline_bmi) & baseline_bmi >= bmi_threshold)
  }

  check_antipsychotic_switch_vars <- function(person_data) {
    missing_vars <- setdiff("OLANZAPINE_SWITCH", names(person_data))
    if (length(missing_vars) > 0) {
      stop(
        "Antipsychotic-switch scenario requested, but analysis dataset is missing: ",
        paste(missing_vars, collapse = ", "),
        ". Ensure Step 1 retains OLANZAPINE_SWITCH as a clean olanzapine-only 0/1 indicator.",
        call. = FALSE
      )
    }
    invisible(TRUE)
  }

  get_antipsychotic_switch_eligibility <- function(person_data,
                                                   threshold_value) {
    threshold_positive <- get_threshold_eligibility(person_data, threshold_value)
    
    if (!"OLANZAPINE_SWITCH" %in% names(person_data)) {
      return(rep(0L, nrow(person_data)))
    }
    
    olanzapine_switch <- dplyr::coalesce(
      as_binary_01_model(person_data$OLANZAPINE_SWITCH),
      0L
    )
    as.integer(threshold_positive == 1L & olanzapine_switch == 1L)
  }

  strategy_uses_lifestyle <- function(strategy) {
    strategy %in% c(
      "psymetric_lifestyle",
      "psymetric_lifestyle_metformin",
      "psymetric_lifestyle_metformin_lifelong",
      "psymetric_lifestyle_glp1ra",
      "psymetric_lifestyle_glp1ra_lifelong"
    )
  }

  strategy_uses_metformin <- function(strategy) {
    strategy %in% c(
      "universal_metformin",
      "psymetric_metformin",
      "psymetric_lifestyle_metformin",
      "psymetric_metformin_lifelong",
      "psymetric_lifestyle_metformin_lifelong"
    )
  }

  strategy_uses_glp1ra <- function(strategy) {
    strategy %in% c(
      "psymetric_glp1ra",
      "psymetric_glp1ra_lifelong",
      "psymetric_lifestyle_glp1ra",
      "psymetric_lifestyle_glp1ra_lifelong"
    )
  }

  strategy_uses_antipsychotic_switch <- function(strategy) {
    strategy == "psymetric_antipsychotic_switch"
  }

  strategy_is_universal_metformin <- function(strategy) {
    strategy == "universal_metformin"
  }

  strategy_is_psymetric_guided <- function(strategy) {
    !(strategy %in% c("usual_care", "universal_metformin"))
  }

  get_strategy_threshold_eligibility <- function(person_data,
                                                 strategy,
                                                 threshold_value) {
    if (strategy_is_universal_metformin(strategy)) {
      return(rep(1L, nrow(person_data)))
    }
    get_threshold_eligibility(person_data, threshold_value)
  }

  get_strategy_threshold_label <- function(strategy, threshold_value) {
    dplyr::case_when(
      strategy == "usual_care" ~ "not thresholded",
      strategy_is_universal_metformin(strategy) ~ "universal",
      TRUE ~ paste0(threshold_value * 100, "%")
    )
  }

  strategy_uses_lifelong_glp1ra <- function(strategy) {
    strategy %in% c(
      "psymetric_glp1ra_lifelong",
      "psymetric_lifestyle_glp1ra_lifelong"
    )
  }

  get_metformin_effect_duration_for_strategy <- function(strategy, params = external_parameter_lookup) {
    get_param("metformin_effect_duration_years", params)
  }

  get_metformin_cost_duration_for_strategy <- function(strategy, params = external_parameter_lookup) {
    get_param("metformin_treatment_cost_duration_years", params)
  }

  get_metformin_effective_exposure_for_strategy <- function(strategy,
                                                            params = external_parameter_lookup) {
    if (strategy_is_universal_metformin(strategy)) {
      return(get_param_default(
        "effective_metformin_exposure_universal_base",
        params,
        default = get_param("effective_metformin_exposure_base", params)
      ))
    }
    get_param_default(
      "effective_metformin_exposure_psymetric_base",
      params,
      default = get_param("effective_metformin_exposure_base", params)
    )
  }

  get_metformin_direct_t2d_rr_for_strategy <- function(strategy,
                                                       params = external_parameter_lookup) {
    if (strategy_is_universal_metformin(strategy)) {
      return(get_param_default(
        "rr_t2d_metformin_if_effective_exposure_universal",
        params,
        default = get_param("rr_t2d_metformin_if_effective_exposure", params)
      ))
    }
    get_param_default(
      "rr_t2d_metformin_if_effective_exposure_psymetric",
      params,
      default = get_param("rr_t2d_metformin_if_effective_exposure", params)
    )
  }

  get_glp1ra_effect_duration_for_strategy <- function(strategy, params = external_parameter_lookup) {
    if (strategy_uses_lifelong_glp1ra(strategy)) {
      return(get_param("glp1ra_effect_duration_years_lifelong_scenario", params))
    }
    get_param("glp1ra_treatment_duration_years", params)
  }

  get_glp1ra_cost_duration_for_strategy <- function(strategy, params = external_parameter_lookup) {
    if (strategy_uses_lifelong_glp1ra(strategy)) {
      return(get_param("glp1ra_treatment_duration_years_lifelong_scenario", params))
    }
    get_param("glp1ra_treatment_duration_years", params)
  }

  get_glp1ra_cycle_cost_for_strategy <- function(strategy,
                                                 cycle,
                                                 params = external_parameter_lookup) {
    glp1ra_duration <- get_glp1ra_cost_duration_for_strategy(strategy, params)
    cycle <- as.numeric(cycle)
    if (length(cycle) != 1 || is.na(cycle)) {
      stop("cycle must be a single non-missing number in get_glp1ra_cycle_cost_for_strategy().", call. = FALSE)
    }
    if (cycle >= glp1ra_duration) {
      return(0)
    }
    dplyr::case_when(
      cycle == 0 ~ get_param("cost_glp1ra_total_year1", params),
      cycle == 1 ~ get_param("cost_glp1ra_total_year2", params),
      TRUE ~ get_param("cost_glp1ra_total_annual_after_year2", params)
    ) |>
      as.numeric()
  }

  get_expected_t2d_hr_for_strategy <- function(strategy,
                                               cycle,
                                               params = external_parameter_lookup) {
    hr <- 1
    
    if (strategy_uses_lifestyle(strategy)) {
      baseline_engaged <- get_param("baseline_effective_lifestyle_engagement_usual_care", params)
      uplift <- get_param("lifestyle_engagement_uplift_among_not_engaged_base", params)
      additional_effective_engagement <- (1 - baseline_engaged) * uplift
      
      lifestyle_full_hr <-
        1 - additional_effective_engagement *
        (1 - get_param("rr_t2d_lifestyle_if_effective_engagement", params))
      
      lifestyle_weight <- get_effect_weight(
        cycle = cycle,
        full_effect_duration_years = get_param("lifestyle_effect_duration_years", params),
        taper_years = get_param("lifestyle_effect_taper_years", params)
      )
      
      lifestyle_cycle_hr <- 1 - lifestyle_weight * (1 - lifestyle_full_hr)
      
      hr <- hr * lifestyle_cycle_hr
    }
    
    if (strategy_uses_metformin(strategy)) {
      exposure <- get_metformin_effective_exposure_for_strategy(strategy, params)
      metformin_rr <- get_metformin_direct_t2d_rr_for_strategy(strategy, params)
      
      metformin_full_hr <-
        1 - exposure *
        (1 - metformin_rr)
      
      metformin_weight <- get_effect_weight(
        cycle = cycle,
        full_effect_duration_years = get_metformin_effect_duration_for_strategy(strategy, params),
        taper_years = get_param("metformin_effect_taper_years", params)
      )
      
      metformin_cycle_hr <- 1 - metformin_weight * (1 - metformin_full_hr)
      
      hr <- hr * metformin_cycle_hr
    }
    
    if (strategy_uses_glp1ra(strategy)) {
      exposure <- get_param("effective_glp1ra_exposure_base", params)
      
      glp1ra_full_hr <-
        1 - exposure *
        (1 - get_param("rr_t2d_glp1ra_if_effective_exposure", params))
      
      glp1ra_weight <- get_effect_weight(
        cycle = cycle,
        full_effect_duration_years = get_glp1ra_effect_duration_for_strategy(strategy, params),
        taper_years = get_param("glp1ra_effect_taper_years", params)
      )
      
      glp1ra_cycle_hr <- 1 - glp1ra_weight * (1 - glp1ra_full_hr)
      
      hr <- hr * glp1ra_cycle_hr
    }
    
    if (strategy_uses_antipsychotic_switch(strategy)) {
      exposure <- get_param("effective_antipsychotic_switch_exposure", params)
      switch_full_hr <-
        1 - exposure *
        (1 - get_param_default("rr_t2d_antipsychotic_switch_if_effective", params, default = 0.80))
      
      switch_weight <- get_effect_weight(
        cycle = cycle,
        full_effect_duration_years = get_param("antipsychotic_switch_effect_duration_years", params),
        taper_years = get_param("antipsychotic_switch_effect_taper_years", params)
      )
      
      switch_cycle_hr <- 1 - switch_weight * (1 - switch_full_hr)
      
      hr <- hr * switch_cycle_hr
    }
    
    as.numeric(hr)
  }

  get_expected_cvd_hr_for_strategy <- function(strategy,
                                               cycle,
                                               params = external_parameter_lookup) {
    hr <- 1
    
    if (strategy_uses_glp1ra(strategy)) {
      exposure <- get_param("effective_glp1ra_exposure_base", params)
      
      glp1ra_full_hr <-
        1 - exposure *
        (1 - get_param("rr_cvd_glp1ra_if_effective_exposure", params))
      
      glp1ra_weight <- get_effect_weight(
        cycle = cycle,
        full_effect_duration_years = get_glp1ra_effect_duration_for_strategy(strategy, params),
        taper_years = get_param("glp1ra_effect_taper_years", params)
      )
      
      glp1ra_cycle_hr <- 1 - glp1ra_weight * (1 - glp1ra_full_hr)
      
      hr <- hr * glp1ra_cycle_hr
    }
    
    as.numeric(hr)
  }


  get_obesity_bmi_threshold <- function(ethnicity_asian) {
    ethnicity_asian <- as_binary_01_model(ethnicity_asian)
    dplyr::case_when(
      ethnicity_asian == 1L ~ 27.5,
      TRUE ~ 30
    )
  }

  get_lifestyle_bmi_multiplier_for_cycle <- function(cycle,
                                                     params = external_parameter_lookup) {
    full_multiplier <- get_param("bmi_multiplier_lifestyle_if_effective_engagement", params)
    effect_weight <- get_effect_weight(
      cycle = cycle,
      full_effect_duration_years = get_param("lifestyle_effect_duration_years", params),
      taper_years = get_param("lifestyle_effect_taper_years", params)
    )
    1 - effect_weight * (1 - full_multiplier)
  }

  get_metformin_bmi_multiplier_for_cycle <- function(strategy,
                                                     cycle,
                                                     params = external_parameter_lookup) {
    full_multiplier <- get_param("bmi_multiplier_metformin_if_effective_exposure", params)
    effect_weight <- get_effect_weight(
      cycle = cycle,
      full_effect_duration_years = get_metformin_effect_duration_for_strategy(strategy, params),
      taper_years = get_param("metformin_effect_taper_years", params)
    )
    1 - effect_weight * (1 - full_multiplier)
  }

  get_glp1ra_bmi_multiplier_for_cycle <- function(strategy,
                                                  cycle,
                                                  params = external_parameter_lookup) {
    exposed_multiplier <- get_param("bmi_multiplier_glp1ra_if_effective_exposure", params)
    after_regain_multiplier <- get_param_default(
      "bmi_multiplier_glp1ra_after_regain",
      params,
      default = 1 - (
        (1 - get_param("bmi_multiplier_glp1ra_if_effective_exposure", params)) *
          (1 - get_param("glp1ra_bmi_effect_regain_fraction_after_stopping", params))
      )
    )
    duration <- get_glp1ra_effect_duration_for_strategy(strategy, params)
    taper_years <- get_param("glp1ra_effect_taper_years", params)
    dplyr::case_when(
      cycle < duration ~ exposed_multiplier,
      taper_years <= 0 ~ after_regain_multiplier,
      cycle >= duration + taper_years ~ after_regain_multiplier,
      TRUE ~ {
        regain_fraction <- (cycle - duration + 1) / taper_years
        exposed_multiplier + regain_fraction * (after_regain_multiplier - exposed_multiplier)
      }
    ) |>
      pmin(1) |>
      pmax(0)
  }

  get_expected_bmi_reclassification_probability <- function(person_data,
                                                            strategy,
                                                            cycle,
                                                            params = external_parameter_lookup) {
    n_persons <- nrow(person_data)
    p_reclass <- rep(0, n_persons)
    
    baseline_bmi <- as.numeric(person_data$BASELINE_BMI)
    obesity_threshold <- get_obesity_bmi_threshold(person_data$ETHNICITY_ASIAN)
    threshold_eligible <- as.integer(person_data$threshold_eligible == 1L)
    glp1ra_eligible <- as.integer(person_data$glp1ra_eligible == 1L)
    
    if (strategy_uses_lifestyle(strategy) && !strategy_uses_metformin(strategy)) {
      p_lifestyle <-
        (1 - get_param("baseline_effective_lifestyle_engagement_usual_care", params)) *
        get_param("lifestyle_engagement_uplift_among_not_engaged_base", params)
      bmi_multiplier <- get_lifestyle_bmi_multiplier_for_cycle(cycle, params)
      crosses_below_threshold <- !is.na(baseline_bmi) &
        baseline_bmi >= obesity_threshold &
        baseline_bmi * bmi_multiplier < obesity_threshold
      p_reclass <- p_reclass +
        threshold_eligible * p_lifestyle * as.numeric(crosses_below_threshold)
    }
    
    if (strategy_uses_metformin(strategy) && !strategy_uses_lifestyle(strategy)) {
      p_metformin <- get_metformin_effective_exposure_for_strategy(strategy, params)
      bmi_multiplier <- get_metformin_bmi_multiplier_for_cycle(strategy, cycle, params)
      crosses_below_threshold <- !is.na(baseline_bmi) &
        baseline_bmi >= obesity_threshold &
        baseline_bmi * bmi_multiplier < obesity_threshold
      p_reclass <- p_reclass +
        threshold_eligible * p_metformin * as.numeric(crosses_below_threshold)
    }
    
    if (strategy_uses_lifestyle(strategy) && strategy_uses_metformin(strategy)) {
      p_lifestyle <-
        (1 - get_param("baseline_effective_lifestyle_engagement_usual_care", params)) *
        get_param("lifestyle_engagement_uplift_among_not_engaged_base", params)
      p_metformin <- get_metformin_effective_exposure_for_strategy(strategy, params)
      
      lifestyle_multiplier <- get_lifestyle_bmi_multiplier_for_cycle(cycle, params)
      metformin_multiplier <- get_metformin_bmi_multiplier_for_cycle(strategy, cycle, params)
      combined_cap_multiplier <- get_param_default(
        "bmi_multiplier_combined_lifestyle_metformin_if_both_effective",
        params,
        default = 1 - get_param("bmi_reduction_cap_combined_lifestyle_metformin", params)
      )
      combined_multiplier <- max(lifestyle_multiplier * metformin_multiplier, combined_cap_multiplier)
      
      crosses_lifestyle_only <- !is.na(baseline_bmi) &
        baseline_bmi >= obesity_threshold &
        baseline_bmi * lifestyle_multiplier < obesity_threshold
      crosses_metformin_only <- !is.na(baseline_bmi) &
        baseline_bmi >= obesity_threshold &
        baseline_bmi * metformin_multiplier < obesity_threshold
      crosses_both <- !is.na(baseline_bmi) &
        baseline_bmi >= obesity_threshold &
        baseline_bmi * combined_multiplier < obesity_threshold
      
      p_reclass <- p_reclass +
        threshold_eligible * (
          p_lifestyle * (1 - p_metformin) * as.numeric(crosses_lifestyle_only) +
            (1 - p_lifestyle) * p_metformin * as.numeric(crosses_metformin_only) +
            p_lifestyle * p_metformin * as.numeric(crosses_both)
        )
    }
    
    if (strategy_uses_glp1ra(strategy)) {
      p_glp1ra <- get_param("effective_glp1ra_exposure_base", params)
      bmi_multiplier <- get_glp1ra_bmi_multiplier_for_cycle(strategy, cycle, params)
      crosses_below_threshold <- !is.na(baseline_bmi) &
        baseline_bmi >= obesity_threshold &
        baseline_bmi * bmi_multiplier < obesity_threshold
      p_reclass <- p_reclass +
        glp1ra_eligible * p_glp1ra * as.numeric(crosses_below_threshold)
    }
    
    clip_probability(p_reclass)
  }

  get_expected_bmi_reduction_for_strategy <- function(person_data,
                                                      strategy,
                                                      cycle,
                                                      params = external_parameter_lookup) {
    n_persons <- nrow(person_data)
    baseline_bmi <- as.numeric(person_data$BASELINE_BMI)
    baseline_bmi[is.na(baseline_bmi)] <- 0
    expected_multiplier <- rep(1, n_persons)
    threshold_eligible <- as.integer(person_data$threshold_eligible == 1L)
    glp1ra_eligible <- as.integer(person_data$glp1ra_eligible == 1L)
    
    if (strategy_uses_lifestyle(strategy) && !strategy_uses_metformin(strategy)) {
      p_lifestyle <-
        (1 - get_param("baseline_effective_lifestyle_engagement_usual_care", params)) *
        get_param("lifestyle_engagement_uplift_among_not_engaged_base", params)
      lifestyle_multiplier <- get_lifestyle_bmi_multiplier_for_cycle(cycle, params)
      expected_multiplier <- dplyr::if_else(
        threshold_eligible == 1L,
        1 - p_lifestyle * (1 - lifestyle_multiplier),
        expected_multiplier
      )
    }
    
    if (strategy_uses_metformin(strategy) && !strategy_uses_lifestyle(strategy)) {
      p_metformin <- get_metformin_effective_exposure_for_strategy(strategy, params)
      metformin_multiplier <- get_metformin_bmi_multiplier_for_cycle(strategy, cycle, params)
      expected_multiplier <- dplyr::if_else(
        threshold_eligible == 1L,
        1 - p_metformin * (1 - metformin_multiplier),
        expected_multiplier
      )
    }
    
    if (strategy_uses_lifestyle(strategy) && strategy_uses_metformin(strategy)) {
      p_lifestyle <-
        (1 - get_param("baseline_effective_lifestyle_engagement_usual_care", params)) *
        get_param("lifestyle_engagement_uplift_among_not_engaged_base", params)
      p_metformin <- get_metformin_effective_exposure_for_strategy(strategy, params)
      lifestyle_multiplier <- get_lifestyle_bmi_multiplier_for_cycle(cycle, params)
      metformin_multiplier <- get_metformin_bmi_multiplier_for_cycle(strategy, cycle, params)
      combined_cap_multiplier <- get_param_default(
        "bmi_multiplier_combined_lifestyle_metformin_if_both_effective",
        params,
        default = 1 - get_param("bmi_reduction_cap_combined_lifestyle_metformin", params)
      )
      combined_multiplier <- max(lifestyle_multiplier * metformin_multiplier, combined_cap_multiplier)
      combined_expected_multiplier <-
        ((1 - p_lifestyle) * (1 - p_metformin) * 1) +
        (p_lifestyle * (1 - p_metformin) * lifestyle_multiplier) +
        ((1 - p_lifestyle) * p_metformin * metformin_multiplier) +
        (p_lifestyle * p_metformin * combined_multiplier)
      expected_multiplier <- dplyr::if_else(
        threshold_eligible == 1L,
        combined_expected_multiplier,
        expected_multiplier
      )
    }
    
    if (strategy_uses_glp1ra(strategy)) {
      p_glp1ra <- get_param("effective_glp1ra_exposure_base", params)
      glp1ra_multiplier <- get_glp1ra_bmi_multiplier_for_cycle(strategy, cycle, params)
      expected_multiplier <- dplyr::if_else(
        glp1ra_eligible == 1L,
        1 - p_glp1ra * (1 - glp1ra_multiplier),
        expected_multiplier
      )
    }
    
    pmax(baseline_bmi * (1 - pmin(pmax(expected_multiplier, 0), 1)), 0)
  }

  get_bmi_mediated_incidence_hr <- function(bmi_reduction,
                                            rr_per_1_bmi_unit_reduction,
                                            enabled = TRUE) {
    if (!isTRUE(enabled)) {
      return(rep(1, length(bmi_reduction)))
    }
    rr_per_1_bmi_unit_reduction <- as.numeric(rr_per_1_bmi_unit_reduction)
    if (is.na(rr_per_1_bmi_unit_reduction) || rr_per_1_bmi_unit_reduction <= 0) {
      return(rep(1, length(bmi_reduction)))
    }
    pmin(pmax(rr_per_1_bmi_unit_reduction ^ pmax(as.numeric(bmi_reduction), 0), 0), 1)
  }

  apply_intervention_bmi_reclassification <- function(transition_rows,
                                                      enabled = TRUE) {
    if (!isTRUE(enabled) || !"p_obesity_reclassification" %in% names(transition_rows)) {
      return(
        transition_rows |>
          dplyr::mutate(
            obesity_reclassified_to_non_obesity = 0,
            obesity_events_prevented_by_bmi_reclassification = 0
          )
      )
    }
    
    transition_rows <- transition_rows |>
      dplyr::mutate(
        p_obesity_reclassification = clip_probability(
          dplyr::coalesce(p_obesity_reclassification, 0)
        ),
        obesity_reclassified_to_non_obesity = 0,
        obesity_events_prevented_by_bmi_reclassification = 0
      )
    
    reclassification_candidates <- transition_rows |>
      dplyr::filter(
        to_state %in% c("Obesity", "Post_CVD_Obesity"),
        p_obesity_reclassification > 0,
        state_probability > 0
      )
    
    if (nrow(reclassification_candidates) == 0) {
      return(transition_rows)
    }
    
    non_candidates <- transition_rows |>
      dplyr::filter(
        !(to_state %in% c("Obesity", "Post_CVD_Obesity") &
            p_obesity_reclassification > 0 &
            state_probability > 0)
      )
    
    not_reclassified <- reclassification_candidates |>
      dplyr::mutate(
        state_probability = state_probability * (1 - p_obesity_reclassification)
      )
    
    reclassified <- reclassification_candidates |>
      dplyr::mutate(
        state_probability = state_probability * p_obesity_reclassification,
        to_state = dplyr::case_when(
          to_state == "Obesity" ~ "Healthy",
          to_state == "Post_CVD_Obesity" ~ "Post_CVD",
          TRUE ~ to_state
        ),
        obesity_reclassified_to_non_obesity = 1,
        obesity_events_prevented_by_bmi_reclassification = event_obesity,
        event_obesity = 0
      )
    
    dplyr::bind_rows(
      non_candidates,
      not_reclassified,
      reclassified
    ) |>
      dplyr::filter(state_probability > state_probability_pruning_threshold)
  }

  get_strategy_costs_for_cycle <- function(strategy,
                                           person_data,
                                           threshold_value,
                                           alive_start_probability,
                                           metformin_prevention_cost_probability = alive_start_probability,
                                           cycle,
                                           params = external_parameter_lookup,
                                           include_operational_cost = FALSE) {
    n_persons <- nrow(person_data)
    eligible_threshold <- get_strategy_threshold_eligibility(person_data, strategy, threshold_value)
    eligible_glp1ra <- get_glp1ra_eligibility(person_data, threshold_value, params)
    eligible_antipsychotic_switch <- get_antipsychotic_switch_eligibility(person_data, threshold_value)
    
    initial_cost <- rep(0, n_persons)
    annual_cost <- rep(0, n_persons)
    
    metformin_prevention_cost_probability <- clip_probability(
      dplyr::coalesce(
        as.numeric(metformin_prevention_cost_probability),
        as.numeric(alive_start_probability)
      )
    )
    
    if (strategy_is_psymetric_guided(strategy) && cycle == 0) {
      initial_cost <- initial_cost +
        alive_start_probability *
        get_param("cost_psymetric_risk_consultation_once", params)
    }
    
    if (isTRUE(include_operational_cost) && strategy_is_psymetric_guided(strategy)) {
      operational_denominator <- get_param_default(
        "psymetric_operational_denominator_per_trust_year",
        params,
        default = 1000
      )
      
      if (is.na(operational_denominator) || operational_denominator <= 0) {
        stop(
          "psymetric_operational_denominator_per_trust_year must be > 0 when operational costs are included.",
          call. = FALSE
        )
      }
      
      operational_age_lower <- get_param_default(
        "psymetric_operational_cost_eligibility_age_lower",
        params,
        default = 16
      )
      
      operational_age_upper_exclusive <- get_param_default(
        "psymetric_operational_cost_eligibility_age_upper_exclusive",
        params,
        default = 36
      )
      
      cycle_length <- get_param("cycle_length_years", params)
      
      attained_age_start <- as.numeric(person_data$AGE_AT_INDEX) +
        cycle * cycle_length
      attained_age_end <- attained_age_start + cycle_length
      
      psymetric_operational_eligible_person_time_fraction <-
        pmax(
          0,
          pmin(attained_age_end, operational_age_upper_exclusive) -
            pmax(attained_age_start, operational_age_lower)
        ) / cycle_length
      
      psymetric_operational_eligible_person_time_fraction <-
        pmin(pmax(psymetric_operational_eligible_person_time_fraction, 0), 1)
      
      annual_cost <- annual_cost +
        alive_start_probability *
        psymetric_operational_eligible_person_time_fraction *
        (
          get_param("cost_psymetric_operational_annual_per_trust_scenario", params) /
            operational_denominator
        )
    }
    
    if (strategy_uses_lifestyle(strategy) && cycle == 0) {
      baseline_engaged <- get_param("baseline_effective_lifestyle_engagement_usual_care", params)
      uplift <- get_param("lifestyle_engagement_uplift_among_not_engaged_base", params)
      additional_effective_engagement <- (1 - baseline_engaged) * uplift
      
      initial_cost <- initial_cost +
        alive_start_probability *
        eligible_threshold *
        additional_effective_engagement *
        get_param("cost_lifestyle_programme_per_effective_engagement", params)
    }
    
    if (strategy_uses_metformin(strategy)) {
      metformin_duration <- get_metformin_cost_duration_for_strategy(strategy, params)
      metformin_effective_exposure <- get_metformin_effective_exposure_for_strategy(strategy, params)
      
      if (cycle == 0) {
        preventive_metformin_baseline_eligible <- as.numeric(
          !is_t2d_state(as.character(person_data$baseline_markov_state))
        )
        
        initial_cost <- initial_cost +
          alive_start_probability *
          eligible_threshold *
          preventive_metformin_baseline_eligible *
          get_param_default(
            "cost_preventive_metformin_initiation_assessment_once",
            params,
            default = 0
          )
      }
      
      if (cycle < metformin_duration) {
        annual_cost <- annual_cost +
          metformin_prevention_cost_probability *
          eligible_threshold *
          (
            metformin_effective_exposure *
              get_param("cost_metformin_total_annual", params) +
              get_param_default(
                "cost_preventive_metformin_monitoring_extra_annual",
                params,
                default = 0
              )
          )
      }
    }
    
    if (strategy_uses_glp1ra(strategy)) {
      glp1ra_cost <- get_glp1ra_cycle_cost_for_strategy(strategy, cycle, params)
      
      if (glp1ra_cost > 0) {
        annual_cost <- annual_cost +
          alive_start_probability *
          eligible_glp1ra *
          get_param("effective_glp1ra_exposure_base", params) *
          glp1ra_cost
      }
    }
    
    if (strategy_uses_antipsychotic_switch(strategy)) {
      switch_exposure <- get_param("effective_antipsychotic_switch_exposure", params)
      
      if (cycle == 0) {
        initial_cost <- initial_cost +
          alive_start_probability *
          eligible_antipsychotic_switch *
          get_param("antipsychotic_switch_review_cost_once", params)
      
        annual_cost <- annual_cost +
          alive_start_probability *
          eligible_antipsychotic_switch *
          switch_exposure *
          get_param("antipsychotic_switch_excess_treatment_failure_probability", params) *
          get_param("relapse_hospitalisation_cost", params)
      }
      
      switch_duration <- get_param("antipsychotic_switch_effect_duration_years", params)
      
      if (cycle < switch_duration) {
        annual_cost <- annual_cost +
          alive_start_probability *
          eligible_antipsychotic_switch *
          switch_exposure *
          get_param("antipsychotic_switch_drug_cost_increment_annual", params)
      }
    }
    
    tibble::tibble(
      ID = person_data$ID,
      strategy_initial_cost = initial_cost,
      strategy_annual_cost = annual_cost
    )
  }

  get_strategy_qaly_adjustment_for_cycle <- function(strategy,
                                                     person_data,
                                                     threshold_value,
                                                     alive_start_probability,
                                                     cycle,
                                                     params = external_parameter_lookup) {
    if (!strategy_uses_antipsychotic_switch(strategy) || cycle != 0) {
      return(0)
    }
    
    eligible_antipsychotic_switch <- get_antipsychotic_switch_eligibility(
      person_data,
      threshold_value
    )
    
    sum(
      alive_start_probability *
        eligible_antipsychotic_switch *
        get_param("effective_antipsychotic_switch_exposure", params) *
        get_param("antipsychotic_switch_excess_treatment_failure_probability", params) *
        get_param("relapse_treatment_failure_qaly_decrement_once", params),
      na.rm = TRUE
    )
  }

  #### #### #### #### #### #### #### ####

  check_analysis_dataset <- function(data) {
    required_vars <- c(
      "ID",
      "AGE_AT_INDEX",
      "SEX",
      "BASELINE_BMI",
      "ETHNICITY_ASIAN",
      "PsyMetRiC2_T2D",
      "baseline_markov_state",
      "state_healthy_0",
      "state_obesity_0"
    )
    
    missing_vars <- setdiff(required_vars, names(data))
    
    if (length(missing_vars) > 0) {
      stop(
        "Analysis dataset is missing required variables: ",
        paste(missing_vars, collapse = ", "),
        call. = FALSE
      )
    }
    
    if (any(is.na(data$ID))) {
      stop("Analysis dataset contains missing ID values.", call. = FALSE)
    }
    
    if (anyDuplicated(data$ID) > 0) {
      stop("Analysis dataset contains duplicate ID values.", call. = FALSE)
    }
    
    if (any(is.na(data$baseline_markov_state))) {
      stop("baseline_markov_state contains missing values.", call. = FALSE)
    }
    
    if (any(is.na(normalise_model_sex(data$SEX)))) {
      stop("SEX contains values that cannot be normalised to Male/Female.", call. = FALSE)
    }
    
    invisible(TRUE)
  }

  check_transition_grids_for_imputation <- function(data,
                                                    grids_i) {
    required_outcomes <- c("t2d", "obesity", "cvd")
    missing_outcomes <- setdiff(required_outcomes, names(grids_i))
    
    if (length(missing_outcomes) > 0) {
      stop(
        "annual_transition_grids is missing: ",
        paste(missing_outcomes, collapse = ", "),
        call. = FALSE
      )
    }
    
    for (outcome in c("t2d", "cvd")) {
      missing_ids <- setdiff(data$ID, grids_i[[outcome]]$ID)
      
      if (length(missing_ids) > 0) {
        stop(
          outcome,
          " transition grid does not contain all analysis IDs. Missing n = ",
          length(missing_ids),
          ". This usually means the survival risk set differs from the Markov cohort.",
          call. = FALSE
        )
      }
    }
    
    invisible(TRUE)
  }

  #### #### #### #### #### #### #### #### #### ####
  #### Markov initial-state construction ####
  #### #### #### #### #### #### #### #### #### ####

  make_initial_state_distribution <- function(data) {
    data |>
      dplyr::transmute(
        ID,
        state = as.character(baseline_markov_state),
        t2d_duration_years = dplyr::case_when(
          state %in% c("T2D", "Post_CVD_T2D") ~ 0,
          TRUE ~ NA_real_
        ),
        state_probability = 1
      ) |>
      dplyr::filter(!is.na(state), state_probability > 0)
  }


  #### #### #### #### #### #### #### ####
  #### Transition row helper ####
  #### #### #### #### #### #### #### ####

  make_transition_rows <- function(data,
                                   to_state,
                                   probability,
                                   t2d_duration_years_to,
                                   event_t2d = 0,
                                   event_obesity = 0,
                                   event_t2d_complications = 0,
                                   event_nonfatal_cvd = 0,
                                   event_death_cvd = 0,
                                   event_death_non_cvd = 0) {
    to_state_q <- rlang::enquo(to_state)
    probability_q <- rlang::enquo(probability)
    t2d_duration_q <- rlang::enquo(t2d_duration_years_to)
    event_t2d_q <- rlang::enquo(event_t2d)
    event_obesity_q <- rlang::enquo(event_obesity)
    event_t2d_complications_q <- rlang::enquo(event_t2d_complications)
    event_nonfatal_cvd_q <- rlang::enquo(event_nonfatal_cvd)
    event_death_cvd_q <- rlang::enquo(event_death_cvd)
    event_death_non_cvd_q <- rlang::enquo(event_death_non_cvd)
    
    n_rows <- nrow(data)
    
    recycle_to_n <- function(x, name) {
      if (n_rows == 0) {
        return(x[0])
      }
      if (length(x) == 1) {
        return(rep(x, n_rows))
      }
      if (length(x) == n_rows) {
        return(x)
      }
      stop(
        "Transition helper argument '", name, "' has length ", length(x),
        " but expected length 1 or ", n_rows, ".",
        call. = FALSE
      )
    }
    
    eval_transition_arg <- function(q, name) {
      recycle_to_n(rlang::eval_tidy(q, data = data), name)
    }
    
    p_obesity_reclassification <- if ("p_obesity_reclassification" %in% names(data)) {
      data$p_obesity_reclassification
    } else {
      rep(0, n_rows)
    }
    
    tibble::tibble(
      ID = data$ID,
      from_state = data$state,
      to_state = eval_transition_arg(to_state_q, "to_state"),
      t2d_duration_years = eval_transition_arg(t2d_duration_q, "t2d_duration_years_to"),
      state_probability = eval_transition_arg(probability_q, "probability"),
      event_t2d = eval_transition_arg(event_t2d_q, "event_t2d"),
      event_obesity = eval_transition_arg(event_obesity_q, "event_obesity"),
      event_t2d_complications = eval_transition_arg(event_t2d_complications_q, "event_t2d_complications"),
      event_nonfatal_cvd = eval_transition_arg(event_nonfatal_cvd_q, "event_nonfatal_cvd"),
      event_death_cvd = eval_transition_arg(event_death_cvd_q, "event_death_cvd"),
      event_death_non_cvd = eval_transition_arg(event_death_non_cvd_q, "event_death_non_cvd"),
      p_obesity_reclassification = p_obesity_reclassification
    )
  }


  #### #### #### #### #### #### #### #### #### ####
  #### Chronic utility and cost calculation ####
  #### #### #### #### #### #### #### #### #### ####

  calculate_chronic_qalys_and_costs <- function(start_distribution,
                                                end_distribution,
                                                transition_rows,
                                                person_data,
                                                cycle,
                                                utility_values,
                                                params = external_parameter_lookup) {
    cycle_length <- get_param("cycle_length_years", params)
    age_mid <- as.numeric(person_data$AGE_AT_INDEX) + cycle + 0.5 * cycle_length
    person_lookup <- person_data |>
      dplyr::transmute(
        ID,
        SEX,
        age_mid = age_mid
      )
    
    average_state_distribution <- dplyr::bind_rows(
      start_distribution |>
        dplyr::transmute(
          ID,
          state,
          average_state_probability = 0.5 * state_probability
        ),
      end_distribution |>
        dplyr::transmute(
          ID,
          state,
          average_state_probability = 0.5 * state_probability
        )
    ) |>
      dplyr::group_by(ID, state) |>
      dplyr::summarise(
        average_state_probability = sum(average_state_probability, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::left_join(person_lookup, by = "ID") |>
      dplyr::mutate(
        utility_name = state_to_utility_column(state),
        utility = lookup_named_utility(
          age = age_mid,
          sex = SEX,
          utility_name = utility_name,
          utility_table = utility_values
        ),
        chronic_qalys = average_state_probability * utility * cycle_length,
        chronic_costs = chronic_state_cost(
          state = state,
          state_probability = average_state_probability,
          params = params
        ) * cycle_length
      )
    
    chronic_qaly_total <- sum(average_state_distribution$chronic_qalys, na.rm = TRUE)
    chronic_cost_total <- sum(average_state_distribution$chronic_costs, na.rm = TRUE)
    
    cvd_event_adjustments <- transition_rows |>
      dplyr::filter(
        event_nonfatal_cvd == 1,
        state_probability > 0
      ) |>
      dplyr::left_join(person_lookup, by = "ID") |>
      dplyr::mutate(
        source_utility_name = state_to_utility_column(from_state),
        destination_utility_name = state_to_utility_column(to_state),
        event_utility_name = state_to_event_year_cvd_utility_column(from_state),
        source_utility = lookup_named_utility(
          age = age_mid,
          sex = SEX,
          utility_name = source_utility_name,
          utility_table = utility_values
        ),
        destination_utility = lookup_named_utility(
          age = age_mid,
          sex = SEX,
          utility_name = destination_utility_name,
          utility_table = utility_values
        ),
        event_utility = lookup_named_utility(
          age = age_mid,
          sex = SEX,
          utility_name = event_utility_name,
          utility_table = utility_values
        ),
        hcc_utility_for_event_flow =
          0.5 * source_utility + 0.5 * destination_utility,
        cvd_event_qaly_adjustment = dplyr::case_when(
          isTRUE(replace_post_cvd_chronic_utility_in_event_year) ~
            state_probability *
            (event_utility - hcc_utility_for_event_flow) *
            cycle_length,
          TRUE ~ 0
        ),
        source_cost = chronic_state_cost(
          state = from_state,
          state_probability = 1,
          params = params
        ),
        destination_cost = chronic_state_cost(
          state = to_state,
          state_probability = 1,
          params = params
        ),
        hcc_cost_for_event_flow =
          0.5 * source_cost + 0.5 * destination_cost,
        cvd_event_chronic_cost_adjustment = dplyr::case_when(
          isTRUE(replace_post_cvd_chronic_cost_in_event_year) ~
            state_probability *
            (source_cost - hcc_cost_for_event_flow) *
            cycle_length,
          TRUE ~ 0
        )
      )
    
    tibble::tibble(
      chronic_qalys = chronic_qaly_total,
      chronic_costs = chronic_cost_total,
      cvd_event_qaly_adjustment = sum(
        cvd_event_adjustments$cvd_event_qaly_adjustment,
        na.rm = TRUE
      ),
      cvd_event_chronic_cost_adjustment = sum(
        cvd_event_adjustments$cvd_event_chronic_cost_adjustment,
        na.rm = TRUE
      )
    )
  }


  #### #### #### #### #### #### #### #### #### #### ####
  #### One Markov cycle for one scenario/imputation ####
  #### #### #### #### #### #### #### #### #### #### ####

  run_one_markov_cycle <- function(state_distribution,
                                   person_data,
                                   grids_i,
                                   cycle,
                                   strategy,
                                   threshold_value,
                                   mortality_values,
                                   mortality_cause_by_age,
                                   utility_values,
                                   complication_grid,
                                   params = external_parameter_lookup,
                                   return_transition_rows = FALSE,
                                   bmi_reclassification_enabled = enable_intervention_bmi_reclassification,
                                   include_operational_cost = include_psymetric_operational_cost_scenario) {
    cycle_length <- get_param("cycle_length_years", params)
    model_max_age <- get_param("model_max_age", params)
    
    person_model_data <- person_data |>
      dplyr::mutate(
        model_sex = normalise_model_sex(SEX),
        attained_age_start = as.numeric(AGE_AT_INDEX) + cycle * cycle_length,
        attained_age_mid = as.numeric(AGE_AT_INDEX) + (cycle + 0.5) * cycle_length,
        threshold_eligible = get_strategy_threshold_eligibility(person_data, strategy, threshold_value),
        glp1ra_eligible = get_glp1ra_eligibility(person_data, threshold_value, params),
        antipsychotic_switch_eligible = get_antipsychotic_switch_eligibility(person_data, threshold_value)
      )
    
    person_model_data$p_obesity_reclassification <-
      get_expected_bmi_reclassification_probability(
        person_data = person_model_data,
        strategy = strategy,
        cycle = cycle,
        params = params
      )
    
    p_death_by_person <- lookup_all_cause_death_probability(
      age = person_model_data$attained_age_start,
      sex = person_model_data$model_sex,
      mortality_values = mortality_values
    )
    
    # Force death once attained age reaches model maximum.
    p_death_by_person <- dplyr::case_when(
      person_model_data$attained_age_start >= model_max_age ~ 1,
      TRUE ~ p_death_by_person
    ) |>
      clip_probability()
    
    prop_cvd_death_by_person <- lookup_cvd_death_proportion(
      age = person_model_data$attained_age_start,
      sex = person_model_data$model_sex,
      mortality_cause_by_age = mortality_cause_by_age
    )
    
    p_t2d_by_person <- lookup_transition_probability(
      grid = grids_i$t2d,
      ids = person_model_data$ID,
      start_time = cycle,
      missing_to_zero = FALSE
    )
    
    p_obesity_by_person <- lookup_transition_probability(
      grid = grids_i$obesity,
      ids = person_model_data$ID,
      start_time = cycle,
      missing_to_zero = TRUE
    )
    
    p_cvd_by_person <- lookup_transition_probability(
      grid = grids_i$cvd,
      ids = person_model_data$ID,
      start_time = cycle,
      missing_to_zero = FALSE
    )
    
    # PSA or deterministic incidence multipliers. Deterministic runs default to 1.
    incidence_hr_t2d <- get_param_default("incidence_hr_t2d", params, default = 1)
    incidence_hr_obesity <- get_param_default("incidence_hr_obesity", params, default = 1)
    incidence_hr_cvd <- get_param_default("incidence_hr_cvd", params, default = 1)
    incidence_hr_t2d_complications <- get_param_default(
      "incidence_hr_t2d_complications",
      params,
      default = 1
    )
    
    p_t2d_by_person <- apply_hazard_ratio_to_probability(p_t2d_by_person, incidence_hr_t2d)
    p_obesity_by_person <- apply_hazard_ratio_to_probability(p_obesity_by_person, incidence_hr_obesity)
    p_cvd_by_person <- apply_hazard_ratio_to_probability(p_cvd_by_person, incidence_hr_cvd)
    
    bmi_incidence_layer_enabled <- get_param_default(
      "intervention_bmi_effects_feed_into_t2d_cvd_incidence_base",
      params,
      default = 0
    ) == 1
    expected_bmi_reduction_by_person <- get_expected_bmi_reduction_for_strategy(
      person_data = person_model_data,
      strategy = strategy,
      cycle = cycle,
      params = params
    )
    bmi_mediated_t2d_hr_by_person <- get_bmi_mediated_incidence_hr(
      bmi_reduction = expected_bmi_reduction_by_person,
      rr_per_1_bmi_unit_reduction = get_param_default(
        "rr_t2d_per_1_bmi_unit_reduction",
        params,
        default = 1
      ),
      enabled = bmi_incidence_layer_enabled
    )
    bmi_mediated_cvd_hr_by_person <- get_bmi_mediated_incidence_hr(
      bmi_reduction = expected_bmi_reduction_by_person,
      rr_per_1_bmi_unit_reduction = get_param_default(
        "rr_cvd_per_1_bmi_unit_reduction",
        params,
        default = 1
      ),
      enabled = bmi_incidence_layer_enabled
    )
    
    person_model_data <- person_model_data |>
      dplyr::mutate(
        p_death_all_cause = p_death_by_person,
        proportion_deaths_cvd = prop_cvd_death_by_person,
        p_t2d_base = p_t2d_by_person,
        p_obesity = p_obesity_by_person,
        p_cvd_base = p_cvd_by_person,
        expected_bmi_reduction_by_strategy = expected_bmi_reduction_by_person,
        bmi_mediated_t2d_hr = bmi_mediated_t2d_hr_by_person,
        bmi_mediated_cvd_hr = bmi_mediated_cvd_hr_by_person
      )
    
    cycle_data <- state_distribution |>
      dplyr::left_join(
        person_model_data |>
          dplyr::select(
            ID,
            model_sex,
            attained_age_start,
            attained_age_mid,
            threshold_eligible,
            glp1ra_eligible,
            antipsychotic_switch_eligible,
            p_death_all_cause,
            proportion_deaths_cvd,
            p_t2d_base,
            p_obesity,
            p_cvd_base,
            expected_bmi_reduction_by_strategy,
            bmi_mediated_t2d_hr,
            bmi_mediated_cvd_hr,
            p_obesity_reclassification
          ),
        by = "ID"
      ) |>
      dplyr::mutate(
        living_state = is_living_state(state),
        state_mortality_hr = get_state_mortality_hr(state, params),
        p_death_all_cause = dplyr::if_else(
          living_state,
          apply_hazard_ratio_to_probability(p_death_all_cause, state_mortality_hr),
          0
        ),
        proportion_deaths_cvd = dplyr::if_else(
          living_state,
          get_state_cvd_death_proportion(proportion_deaths_cvd, state, params),
          0
        ),
        p_death_cvd = p_death_all_cause * proportion_deaths_cvd,
        p_death_non_cvd = p_death_all_cause * (1 - proportion_deaths_cvd),
        survivor_probability_after_death =
          state_probability * pmax(1 - p_death_all_cause, 0)
      )
    
    t2d_strategy_hr <- get_expected_t2d_hr_for_strategy(
      strategy = strategy,
      cycle = cycle,
      params = params
    )
    
    cvd_strategy_hr <- get_expected_cvd_hr_for_strategy(
      strategy = strategy,
      cycle = cycle,
      params = params
    )
    
    cycle_data <- cycle_data |>
      dplyr::mutate(
        strategy_t2d_effect_eligible = dplyr::case_when(
          strategy_uses_glp1ra(strategy) ~ glp1ra_eligible == 1L,
          strategy_uses_antipsychotic_switch(strategy) ~ antipsychotic_switch_eligible == 1L,
          strategy == "usual_care" ~ FALSE,
          TRUE ~ threshold_eligible == 1L
        ),
        strategy_cvd_effect_eligible = dplyr::case_when(
          strategy_uses_glp1ra(strategy) ~ glp1ra_eligible == 1L,
          TRUE ~ FALSE
        ),
        p_t2d = dplyr::case_when(
          strategy_t2d_effect_eligible ~
            apply_hazard_ratio_to_probability(p_t2d_base, t2d_strategy_hr),
          TRUE ~ p_t2d_base
        ),
        p_t2d = apply_hazard_ratio_to_probability(p_t2d, bmi_mediated_t2d_hr),
        p_cvd = dplyr::case_when(
          state %in% c("T2D", "T2D_complications", "Post_CVD_T2D", "Post_CVD_T2D_complications") ~
            apply_hazard_ratio_to_probability(p_cvd_base, get_param("rr_cvd_after_t2d", params)),
          TRUE ~ p_cvd_base
        ),
        p_cvd = dplyr::case_when(
          strategy_cvd_effect_eligible ~
            apply_hazard_ratio_to_probability(p_cvd, cvd_strategy_hr),
          TRUE ~ p_cvd
        ),
        p_cvd = apply_hazard_ratio_to_probability(p_cvd, bmi_mediated_cvd_hr),
        p_t2d_complications = lookup_t2d_complication_probability(
          t2d_duration_years = t2d_duration_years,
          complication_grid = complication_grid,
          hr_multiplier = incidence_hr_t2d_complications
        ),
        p_t2d = clip_probability(p_t2d),
        p_obesity = clip_probability(p_obesity),
        p_cvd = clip_probability(p_cvd),
        p_t2d_complications = clip_probability(p_t2d_complications)
      )
    
    t2d_susceptible_states <- c(
      "Healthy",
      "Obesity",
      "Post_CVD",
      "Post_CVD_Obesity"
    )
    
    t2d_susceptible_data <- cycle_data |>
      dplyr::filter(
        living_state,
        state %in% t2d_susceptible_states,
        survivor_probability_after_death > 0
      )
    
    expected_t2d_events_without_strategy <- sum(
      t2d_susceptible_data$survivor_probability_after_death *
        t2d_susceptible_data$p_t2d_base,
      na.rm = TRUE
    )
    
    expected_t2d_events_with_strategy <- sum(
      t2d_susceptible_data$survivor_probability_after_death *
        t2d_susceptible_data$p_t2d,
      na.rm = TRUE
    )
    
    expected_t2d_events_prevented_by_strategy <-
      expected_t2d_events_without_strategy - expected_t2d_events_with_strategy
    
    expected_t2d_events_prevented_by_strategy_effect_eligible <- sum(
      t2d_susceptible_data$survivor_probability_after_death *
        (t2d_susceptible_data$p_t2d_base - t2d_susceptible_data$p_t2d) *
        as.numeric(t2d_susceptible_data$strategy_t2d_effect_eligible),
      na.rm = TRUE
    )
    
    n_alive_threshold_eligible_start <- sum(
      cycle_data$state_probability[
        cycle_data$living_state &
          cycle_data$threshold_eligible == 1L
      ],
      na.rm = TRUE
    )
    
    n_t2d_susceptible_threshold_eligible_start <- sum(
      t2d_susceptible_data$state_probability[
        t2d_susceptible_data$threshold_eligible == 1L
      ],
      na.rm = TRUE
    )
    
    # Death transitions first.
    death_cvd_rows <- cycle_data |>
      dplyr::filter(living_state, state_probability > 0) |>
      make_transition_rows(
        to_state = "Death_CVD",
        probability = state_probability * p_death_cvd,
        t2d_duration_years_to = NA_real_,
        event_death_cvd = 1
      )
    
    death_non_cvd_rows <- cycle_data |>
      dplyr::filter(living_state, state_probability > 0) |>
      make_transition_rows(
        to_state = "Death_Non_CVD",
        probability = state_probability * p_death_non_cvd,
        t2d_duration_years_to = NA_real_,
        event_death_non_cvd = 1
      )
    
    death_carryover_rows <- cycle_data |>
      dplyr::filter(is_death_state(state), state_probability > 0) |>
      make_transition_rows(
        to_state = state,
        probability = state_probability,
        t2d_duration_years_to = NA_real_
      )
    
    living_data <- cycle_data |>
      dplyr::filter(living_state, survivor_probability_after_death > 0)
    
    # ------------------------------------------------------------
    # Healthy source state
    # ------------------------------------------------------------
    healthy <- living_data |>
      dplyr::filter(state == "Healthy") |>
      dplyr::mutate(
        prob_after_t2d = survivor_probability_after_death * p_t2d,
        prob_no_t2d = survivor_probability_after_death * (1 - p_t2d),
        prob_no_t2d_to_obesity = prob_no_t2d * p_obesity,
        prob_no_t2d_no_obesity = prob_no_t2d * (1 - p_obesity)
      )
    
    healthy_rows <- dplyr::bind_rows(
      make_transition_rows(
        healthy,
        to_state = "Post_CVD_T2D",
        probability = prob_after_t2d * p_cvd,
        t2d_duration_years_to = 0,
        event_t2d = 1,
        event_nonfatal_cvd = 1
      ),
      make_transition_rows(
        healthy,
        to_state = "T2D",
        probability = prob_after_t2d * (1 - p_cvd),
        t2d_duration_years_to = 0,
        event_t2d = 1
      ),
      make_transition_rows(
        healthy,
        to_state = "Post_CVD_Obesity",
        probability = prob_no_t2d_to_obesity * p_cvd,
        t2d_duration_years_to = NA_real_,
        event_obesity = 1,
        event_nonfatal_cvd = 1
      ),
      make_transition_rows(
        healthy,
        to_state = "Obesity",
        probability = prob_no_t2d_to_obesity * (1 - p_cvd),
        t2d_duration_years_to = NA_real_,
        event_obesity = 1
      ),
      make_transition_rows(
        healthy,
        to_state = "Post_CVD",
        probability = prob_no_t2d_no_obesity * p_cvd,
        t2d_duration_years_to = NA_real_,
        event_nonfatal_cvd = 1
      ),
      make_transition_rows(
        healthy,
        to_state = "Healthy",
        probability = prob_no_t2d_no_obesity * (1 - p_cvd),
        t2d_duration_years_to = NA_real_
      )
    )
    
    # ------------------------------------------------------------
    # Obesity source state
    # ------------------------------------------------------------
    obesity <- living_data |>
      dplyr::filter(state == "Obesity") |>
      dplyr::mutate(
        prob_after_t2d = survivor_probability_after_death * p_t2d,
        prob_no_t2d = survivor_probability_after_death * (1 - p_t2d)
      )
    
    obesity_rows <- dplyr::bind_rows(
      make_transition_rows(
        obesity,
        to_state = "Post_CVD_T2D",
        probability = prob_after_t2d * p_cvd,
        t2d_duration_years_to = 0,
        event_t2d = 1,
        event_nonfatal_cvd = 1
      ),
      make_transition_rows(
        obesity,
        to_state = "T2D",
        probability = prob_after_t2d * (1 - p_cvd),
        t2d_duration_years_to = 0,
        event_t2d = 1
      ),
      make_transition_rows(
        obesity,
        to_state = "Post_CVD_Obesity",
        probability = prob_no_t2d * p_cvd,
        t2d_duration_years_to = NA_real_,
        event_nonfatal_cvd = 1
      ),
      make_transition_rows(
        obesity,
        to_state = "Obesity",
        probability = prob_no_t2d * (1 - p_cvd),
        t2d_duration_years_to = NA_real_
      )
    )
    
    # ------------------------------------------------------------
    # T2D source state
    # ------------------------------------------------------------
    t2d <- living_data |>
      dplyr::filter(state == "T2D") |>
      dplyr::mutate(
        prob_after_comp = survivor_probability_after_death * p_t2d_complications,
        prob_no_comp = survivor_probability_after_death * (1 - p_t2d_complications),
        next_t2d_duration = dplyr::coalesce(t2d_duration_years, 0) + cycle_length
      )
    
    t2d_rows <- dplyr::bind_rows(
      make_transition_rows(
        t2d,
        to_state = "Post_CVD_T2D_complications",
        probability = prob_after_comp * p_cvd,
        t2d_duration_years_to = NA_real_,
        event_t2d_complications = 1,
        event_nonfatal_cvd = 1
      ),
      make_transition_rows(
        t2d,
        to_state = "T2D_complications",
        probability = prob_after_comp * (1 - p_cvd),
        t2d_duration_years_to = NA_real_,
        event_t2d_complications = 1
      ),
      make_transition_rows(
        t2d,
        to_state = "Post_CVD_T2D",
        probability = prob_no_comp * p_cvd,
        t2d_duration_years_to = next_t2d_duration,
        event_nonfatal_cvd = 1
      ),
      make_transition_rows(
        t2d,
        to_state = "T2D",
        probability = prob_no_comp * (1 - p_cvd),
        t2d_duration_years_to = next_t2d_duration
      )
    )
    
    # ------------------------------------------------------------
    # T2D complications source state
    # ------------------------------------------------------------
    t2d_comp <- living_data |>
      dplyr::filter(state == "T2D_complications")
    
    t2d_comp_rows <- dplyr::bind_rows(
      make_transition_rows(
        t2d_comp,
        to_state = "Post_CVD_T2D_complications",
        probability = survivor_probability_after_death * p_cvd,
        t2d_duration_years_to = NA_real_,
        event_nonfatal_cvd = 1
      ),
      make_transition_rows(
        t2d_comp,
        to_state = "T2D_complications",
        probability = survivor_probability_after_death * (1 - p_cvd),
        t2d_duration_years_to = NA_real_
      )
    )
    
    # ------------------------------------------------------------
    # Post-CVD source state
    # ------------------------------------------------------------
    post_cvd <- living_data |>
      dplyr::filter(state == "Post_CVD") |>
      dplyr::mutate(
        prob_after_t2d = survivor_probability_after_death * p_t2d,
        prob_no_t2d = survivor_probability_after_death * (1 - p_t2d),
        prob_no_t2d_to_obesity = prob_no_t2d * p_obesity,
        prob_no_t2d_no_obesity = prob_no_t2d * (1 - p_obesity)
      )
    
    post_cvd_rows <- dplyr::bind_rows(
      make_transition_rows(
        post_cvd,
        to_state = "Post_CVD_T2D",
        probability = prob_after_t2d,
        t2d_duration_years_to = 0,
        event_t2d = 1
      ),
      make_transition_rows(
        post_cvd,
        to_state = "Post_CVD_Obesity",
        probability = prob_no_t2d_to_obesity,
        t2d_duration_years_to = NA_real_,
        event_obesity = 1
      ),
      make_transition_rows(
        post_cvd,
        to_state = "Post_CVD",
        probability = prob_no_t2d_no_obesity,
        t2d_duration_years_to = NA_real_
      )
    )
    
    # ------------------------------------------------------------
    # Post-CVD obesity source state
    # ------------------------------------------------------------
    post_cvd_obesity <- living_data |>
      dplyr::filter(state == "Post_CVD_Obesity") |>
      dplyr::mutate(
        prob_after_t2d = survivor_probability_after_death * p_t2d,
        prob_no_t2d = survivor_probability_after_death * (1 - p_t2d)
      )
    
    post_cvd_obesity_rows <- dplyr::bind_rows(
      make_transition_rows(
        post_cvd_obesity,
        to_state = "Post_CVD_T2D",
        probability = prob_after_t2d,
        t2d_duration_years_to = 0,
        event_t2d = 1
      ),
      make_transition_rows(
        post_cvd_obesity,
        to_state = "Post_CVD_Obesity",
        probability = prob_no_t2d,
        t2d_duration_years_to = NA_real_
      )
    )
    
    # ------------------------------------------------------------
    # Post-CVD T2D source state
    # ------------------------------------------------------------
    post_cvd_t2d <- living_data |>
      dplyr::filter(state == "Post_CVD_T2D") |>
      dplyr::mutate(
        prob_after_comp = survivor_probability_after_death * p_t2d_complications,
        prob_no_comp = survivor_probability_after_death * (1 - p_t2d_complications),
        next_t2d_duration = dplyr::coalesce(t2d_duration_years, 0) + cycle_length
      )
    
    post_cvd_t2d_rows <- dplyr::bind_rows(
      make_transition_rows(
        post_cvd_t2d,
        to_state = "Post_CVD_T2D_complications",
        probability = prob_after_comp,
        t2d_duration_years_to = NA_real_,
        event_t2d_complications = 1
      ),
      make_transition_rows(
        post_cvd_t2d,
        to_state = "Post_CVD_T2D",
        probability = prob_no_comp,
        t2d_duration_years_to = next_t2d_duration
      )
    )
    
    # ------------------------------------------------------------
    # Post-CVD T2D complications source state
    # ------------------------------------------------------------
    post_cvd_t2d_comp <- living_data |>
      dplyr::filter(state == "Post_CVD_T2D_complications")
    
    post_cvd_t2d_comp_rows <- make_transition_rows(
      post_cvd_t2d_comp,
      to_state = "Post_CVD_T2D_complications",
      probability = survivor_probability_after_death,
      t2d_duration_years_to = NA_real_
    )
    
    transition_rows <- dplyr::bind_rows(
      death_cvd_rows,
      death_non_cvd_rows,
      death_carryover_rows,
      healthy_rows,
      obesity_rows,
      t2d_rows,
      t2d_comp_rows,
      post_cvd_rows,
      post_cvd_obesity_rows,
      post_cvd_t2d_rows,
      post_cvd_t2d_comp_rows
    ) |>
      dplyr::mutate(
        state_probability = dplyr::coalesce(state_probability, 0)
      ) |>
      dplyr::filter(state_probability > state_probability_pruning_threshold)
    
    transition_rows <- apply_intervention_bmi_reclassification(
      transition_rows = transition_rows,
      enabled = bmi_reclassification_enabled
    )
    
    next_state_distribution <- transition_rows |>
      dplyr::group_by(ID, state = to_state, t2d_duration_years) |>
      dplyr::summarise(
        state_probability = sum(state_probability, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::filter(state_probability > state_probability_pruning_threshold)
    
    # Strategy costs are calculated at patient level and weighted by alive start
    # probability. Initial one-off costs are discounted at cycle start; annual costs
    # at cycle midpoint.
    alive_start_by_person <- state_distribution |>
      dplyr::filter(is_living_state(state)) |>
      dplyr::group_by(ID) |>
      dplyr::summarise(
        alive_start_probability = sum(state_probability, na.rm = TRUE),
        .groups = "drop"
      )
    
    # Preventive metformin costs should stop when a person transitions to T2D,
    # because annual T2D state costs already represent diabetes management after
    # diagnosis. Annual strategy costs are applied at the cycle midpoint, so use
    # the average of the T2D-free living probability at cycle start and cycle end.
    t2d_free_start_by_person <- state_distribution |>
      dplyr::filter(is_living_state(state), !is_t2d_state(state)) |>
      dplyr::group_by(ID) |>
      dplyr::summarise(
        t2d_free_start_probability = sum(state_probability, na.rm = TRUE),
        .groups = "drop"
      )
    
    t2d_free_end_by_person <- next_state_distribution |>
      dplyr::filter(is_living_state(state), !is_t2d_state(state)) |>
      dplyr::group_by(ID) |>
      dplyr::summarise(
        t2d_free_end_probability = sum(state_probability, na.rm = TRUE),
        .groups = "drop"
      )
    
    person_strategy_data <- person_data |>
      dplyr::select(ID, dplyr::everything()) |>
      dplyr::left_join(alive_start_by_person, by = "ID") |>
      dplyr::left_join(t2d_free_start_by_person, by = "ID") |>
      dplyr::left_join(t2d_free_end_by_person, by = "ID") |>
      dplyr::mutate(
        alive_start_probability = dplyr::coalesce(alive_start_probability, 0),
        t2d_free_start_probability = dplyr::coalesce(t2d_free_start_probability, 0),
        t2d_free_end_probability = dplyr::coalesce(t2d_free_end_probability, 0),
        metformin_prevention_cost_probability = 0.5 *
          (t2d_free_start_probability + t2d_free_end_probability)
      )
    
    strategy_costs <- get_strategy_costs_for_cycle(
      strategy = strategy,
      person_data = person_strategy_data,
      threshold_value = threshold_value,
      alive_start_probability = person_strategy_data$alive_start_probability,
      metformin_prevention_cost_probability = person_strategy_data$metformin_prevention_cost_probability,
      cycle = cycle,
      params = params,
      include_operational_cost = include_operational_cost
    )
    
    strategy_qaly_adjustment <- get_strategy_qaly_adjustment_for_cycle(
      strategy = strategy,
      person_data = person_strategy_data,
      threshold_value = threshold_value,
      alive_start_probability = person_strategy_data$alive_start_probability,
      cycle = cycle,
      params = params
    )
    
    chronic_qaly_costs <- calculate_chronic_qalys_and_costs(
      start_distribution = state_distribution,
      end_distribution = next_state_distribution,
      transition_rows = transition_rows,
      person_data = person_data,
      cycle = cycle,
      utility_values = utility_values,
      params = params
    )
    
    incident_t2d_events <- sum(transition_rows$state_probability * transition_rows$event_t2d, na.rm = TRUE)
    incident_obesity_events <- sum(transition_rows$state_probability * transition_rows$event_obesity, na.rm = TRUE)
    obesity_reclassified_to_non_obesity <- sum(
      transition_rows$state_probability * transition_rows$obesity_reclassified_to_non_obesity,
      na.rm = TRUE
    )
    obesity_events_prevented_by_bmi_reclassification <- sum(
      transition_rows$state_probability * transition_rows$obesity_events_prevented_by_bmi_reclassification,
      na.rm = TRUE
    )
    incident_t2d_complication_events <- sum(
      transition_rows$state_probability * transition_rows$event_t2d_complications,
      na.rm = TRUE
    )
    nonfatal_cvd_events <- sum(
      transition_rows$state_probability * transition_rows$event_nonfatal_cvd,
      na.rm = TRUE
    )
    cvd_deaths <- sum(
      transition_rows$state_probability * transition_rows$event_death_cvd,
      na.rm = TRUE
    )
    non_cvd_deaths <- sum(
      transition_rows$state_probability * transition_rows$event_death_non_cvd,
      na.rm = TRUE
    )
    first_fatal_cvd_events <- transition_rows |>
      dplyr::filter(
        event_death_cvd == 1,
        !is_post_cvd_state(from_state)
      ) |>
      dplyr::summarise(
        events = sum(state_probability, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::pull(events)
    
    if (length(first_fatal_cvd_events) == 0) {
      first_fatal_cvd_events <- 0
    }
    
    acute_cvd_death_events_for_cost <- dplyr::case_when(
      isTRUE(apply_acute_cvd_cost_to_all_cvd_deaths) ~ cvd_deaths,
      TRUE ~ first_fatal_cvd_events
    )
    
    acute_cvd_event_cost <-
      (nonfatal_cvd_events + acute_cvd_death_events_for_cost) *
      get_param("cost_cvd_event_acute", params)
    
    incident_t2d_event_cost <-
      incident_t2d_events *
      get_param("cost_t2d_diagnosis_event_once", params)
    
    undiscounted_initial_strategy_cost <-
      sum(strategy_costs$strategy_initial_cost, na.rm = TRUE)
    
    undiscounted_annual_strategy_cost <-
      sum(strategy_costs$strategy_annual_cost, na.rm = TRUE)
    
    undiscounted_chronic_cost <-
      chronic_qaly_costs$chronic_costs +
      chronic_qaly_costs$cvd_event_chronic_cost_adjustment
    
    undiscounted_event_cost <-
      acute_cvd_event_cost + incident_t2d_event_cost
    
    undiscounted_total_cost <-
      undiscounted_initial_strategy_cost +
      undiscounted_annual_strategy_cost +
      undiscounted_chronic_cost +
      undiscounted_event_cost
    
    undiscounted_strategy_qaly_adjustment <- strategy_qaly_adjustment
    
    undiscounted_total_qaly <-
      chronic_qaly_costs$chronic_qalys +
      chronic_qaly_costs$cvd_event_qaly_adjustment +
      undiscounted_strategy_qaly_adjustment
    
    discount_factor_costs_start <- make_discount_factor(
      rate = get_param("discount_rate_costs", params),
      time = cycle
    )
    
    discount_factor_costs_mid <- make_discount_factor(
      rate = get_param("discount_rate_costs", params),
      time = cycle + 0.5
    )
    
    discount_factor_qalys_mid <- make_discount_factor(
      rate = get_param("discount_rate_qalys", params),
      time = cycle + 0.5
    )
    
    discounted_initial_strategy_cost <-
      undiscounted_initial_strategy_cost * discount_factor_costs_start
    
    discounted_annual_strategy_cost <-
      undiscounted_annual_strategy_cost * discount_factor_costs_mid
    
    discounted_chronic_cost <-
      undiscounted_chronic_cost * discount_factor_costs_mid
    
    discounted_event_cost <-
      undiscounted_event_cost * discount_factor_costs_mid
    
    discounted_total_cost <-
      discounted_initial_strategy_cost +
      discounted_annual_strategy_cost +
      discounted_chronic_cost +
      discounted_event_cost
    
    discounted_total_qaly <-
      undiscounted_total_qaly * discount_factor_qalys_mid
    
    average_state_distribution_for_state_years <- dplyr::bind_rows(
      state_distribution |>
        dplyr::transmute(
          ID,
          state,
          average_state_probability = 0.5 * state_probability
        ),
      next_state_distribution |>
        dplyr::transmute(
          ID,
          state,
          average_state_probability = 0.5 * state_probability
        )
    ) |>
      dplyr::group_by(ID, state) |>
      dplyr::summarise(
        average_state_probability = sum(average_state_probability, na.rm = TRUE),
        .groups = "drop"
      )
    
    t2d_state_years <- sum(
      average_state_distribution_for_state_years$average_state_probability[
        is_t2d_state(average_state_distribution_for_state_years$state)
      ],
      na.rm = TRUE
    ) * cycle_length
    
    t2d_free_living_years <- sum(
      average_state_distribution_for_state_years$average_state_probability[
        is_living_state(average_state_distribution_for_state_years$state) &
          !is_t2d_state(average_state_distribution_for_state_years$state)
      ],
      na.rm = TRUE
    ) * cycle_length
    
    discounted_t2d_state_years <- t2d_state_years * discount_factor_qalys_mid
    discounted_t2d_free_living_years <- t2d_free_living_years * discount_factor_qalys_mid
    
    probability_mass_start <- sum(state_distribution$state_probability, na.rm = TRUE)
    probability_mass_end <- sum(next_state_distribution$state_probability, na.rm = TRUE)
    
    cycle_results <- tibble::tibble(
      cycle = cycle,
      strategy = strategy,
      threshold_value = threshold_value,
      threshold_label = get_strategy_threshold_label(strategy, threshold_value),
      n_alive_start = sum(
        state_distribution$state_probability[is_living_state(state_distribution$state)],
        na.rm = TRUE
      ),
      n_alive_end = sum(
        next_state_distribution$state_probability[is_living_state(next_state_distribution$state)],
        na.rm = TRUE
      ),
      probability_mass_start = probability_mass_start,
      probability_mass_end = probability_mass_end,
      probability_mass_error = probability_mass_end - probability_mass_start,
      n_alive_threshold_eligible_start = n_alive_threshold_eligible_start,
      n_t2d_susceptible_threshold_eligible_start = n_t2d_susceptible_threshold_eligible_start,
      t2d_strategy_hr_for_cycle = t2d_strategy_hr,
      cvd_strategy_hr_for_cycle = cvd_strategy_hr,
      expected_t2d_events_without_strategy = expected_t2d_events_without_strategy,
      expected_t2d_events_with_strategy = expected_t2d_events_with_strategy,
      expected_t2d_events_prevented_by_strategy = expected_t2d_events_prevented_by_strategy,
      expected_t2d_events_prevented_by_strategy_effect_eligible = expected_t2d_events_prevented_by_strategy_effect_eligible,
      t2d_state_years = t2d_state_years,
      t2d_free_living_years = t2d_free_living_years,
      discounted_t2d_state_years = discounted_t2d_state_years,
      discounted_t2d_free_living_years = discounted_t2d_free_living_years,
      incident_t2d = incident_t2d_events,
      incident_obesity = incident_obesity_events,
      obesity_reclassified_to_non_obesity = obesity_reclassified_to_non_obesity,
      obesity_events_prevented_by_bmi_reclassification = obesity_events_prevented_by_bmi_reclassification,
      incident_t2d_serious_non_cvd_complications = incident_t2d_complication_events,
      incident_nonfatal_cvd = nonfatal_cvd_events,
      incident_first_fatal_cvd = first_fatal_cvd_events,
      incident_first_pooled_cvd = nonfatal_cvd_events + first_fatal_cvd_events,
      deaths_cvd = cvd_deaths,
      deaths_non_cvd = non_cvd_deaths,
      undiscounted_initial_strategy_costs = undiscounted_initial_strategy_cost,
      undiscounted_annual_strategy_costs = undiscounted_annual_strategy_cost,
      undiscounted_chronic_costs = undiscounted_chronic_cost,
      undiscounted_event_costs = undiscounted_event_cost,
      undiscounted_costs = undiscounted_total_cost,
      discounted_initial_strategy_costs = discounted_initial_strategy_cost,
      discounted_annual_strategy_costs = discounted_annual_strategy_cost,
      discounted_chronic_costs = discounted_chronic_cost,
      discounted_event_costs = discounted_event_cost,
      discounted_costs = discounted_total_cost,
      undiscounted_strategy_qaly_adjustment = undiscounted_strategy_qaly_adjustment,
      undiscounted_qalys = undiscounted_total_qaly,
      discounted_qalys = discounted_total_qaly
    )
    
    list(
      next_state_distribution = next_state_distribution,
      cycle_results = cycle_results,
      transition_rows = if (isTRUE(return_transition_rows)) transition_rows else NULL
    )
  }


  # ---------------------------------------------------------------------------
  # Validate model inputs and initialise state occupancy
  # ---------------------------------------------------------------------------

  check_analysis_dataset(person_data)
  check_transition_grids_for_imputation(person_data, transition_grids)

  if (strategy_uses_antipsychotic_switch(strategy)) {
    check_antipsychotic_switch_vars(person_data)
  }

  max_cycles <- ceiling(
    get_param("model_max_age", parameters) -
      min(as.numeric(person_data$AGE_AT_INDEX), na.rm = TRUE)
  )

  state_distribution <- make_initial_state_distribution(person_data)
  cycle_results <- vector("list", max_cycles)

  for (cycle in 0:(max_cycles - 1)) {
    cycle_output <- run_one_markov_cycle(
      state_distribution = state_distribution,
      person_data = person_data,
      grids_i = transition_grids,
      cycle = cycle,
      strategy = strategy,
      threshold_value = threshold_value,
      mortality_values = mortality_values,
      mortality_cause_by_age = mortality_cause_allocation_by_age,
      utility_values = utility_values,
      complication_grid = t2d_complication_grid,
      params = parameters,
      return_transition_rows = FALSE,
      bmi_reclassification_enabled = TRUE,
      include_operational_cost = include_operational_cost
    )

    state_distribution <- cycle_output$next_state_distribution
    cycle_results[[cycle + 1]] <- cycle_output$cycle_results

    if (cycle_output$cycle_results$n_alive_end < 1e-8) {
      cycle_results <- cycle_results[seq_len(cycle + 1)]
      break
    }
  }

  cycle_results <- dplyr::bind_rows(cycle_results)

  max_mass_error <- max(abs(cycle_results$probability_mass_error), na.rm = TRUE)
  if (max_mass_error > 1e-6) {
    warning(
      "Probability mass error exceeded 1e-6; maximum absolute error = ",
      signif(max_mass_error, 4),
      call. = FALSE
    )
  }

  n_modelled <- dplyr::n_distinct(person_data$ID)

  scenario_summary <- cycle_results |>
    dplyr::summarise(
      strategy = dplyr::first(strategy),
      threshold_value = dplyr::first(threshold_value),
      threshold_label = dplyr::first(threshold_label),
      n_modelled = n_modelled,
      n_threshold_positive = sum(
        get_threshold_eligibility(person_data, threshold_value),
        na.rm = TRUE
      ),
      n_strategy_eligible = dplyr::case_when(
        strategy == "usual_care" ~ 0,
        strategy == "universal_metformin" ~ n_modelled,
        TRUE ~ sum(
          get_threshold_eligibility(person_data, threshold_value),
          na.rm = TRUE
        )
      ),
      mean_discounted_initial_strategy_costs =
        sum(discounted_initial_strategy_costs, na.rm = TRUE) / n_modelled,
      mean_discounted_annual_strategy_costs =
        sum(discounted_annual_strategy_costs, na.rm = TRUE) / n_modelled,
      mean_discounted_chronic_costs =
        sum(discounted_chronic_costs, na.rm = TRUE) / n_modelled,
      mean_discounted_event_costs =
        sum(discounted_event_costs, na.rm = TRUE) / n_modelled,
      mean_discounted_costs =
        sum(discounted_costs, na.rm = TRUE) / n_modelled,
      mean_discounted_qalys =
        sum(discounted_qalys, na.rm = TRUE) / n_modelled,
      t2d_events_per_1000 =
        1000 * sum(incident_t2d, na.rm = TRUE) / n_modelled,
      obesity_events_per_1000 =
        1000 * sum(incident_obesity, na.rm = TRUE) / n_modelled,
      t2d_complications_per_1000 =
        1000 * sum(incident_t2d_serious_non_cvd_complications, na.rm = TRUE) /
        n_modelled,
      nonfatal_cvd_events_per_1000 =
        1000 * sum(incident_nonfatal_cvd, na.rm = TRUE) / n_modelled,
      cvd_deaths_per_1000 =
        1000 * sum(deaths_cvd, na.rm = TRUE) / n_modelled,
      non_cvd_deaths_per_1000 =
        1000 * sum(deaths_non_cvd, na.rm = TRUE) / n_modelled,
      max_abs_probability_mass_error = max_mass_error,
      .groups = "drop"
    )

  final_state_distribution <- state_distribution |>
    dplyr::group_by(state) |>
    dplyr::summarise(
      final_state_probability = sum(state_probability, na.rm = TRUE),
      .groups = "drop"
    )

  list(
    summary = scenario_summary,
    cycle_results = cycle_results,
    final_state_distribution = final_state_distribution
  )
}


# ==============================================================================
# 3. Small analysis wrappers
# ==============================================================================

run_model_grid <- function(
    scenario_grid,
    analysis_datasets,
    annual_transition_grids,
    t2d_complication_grid,
    parameters,
    mortality_values,
    mortality_cause_allocation_values,
    utility_values
) {
  outputs <- vector("list", nrow(scenario_grid))

  for (i in seq_len(nrow(scenario_grid))) {
    row <- scenario_grid[i, , drop = FALSE]
    imputation <- row$imputation[[1]]

    message(
      "Running ", i, "/", nrow(scenario_grid), ": ",
      row$analysis_id[[1]], "; imputation ", imputation,
      "; ", row$strategy[[1]], "; threshold ",
      formatC(100 * row$threshold_value[[1]], format = "f", digits = 0), "%"
    )

    scenario_parameters <- parameters
    overrides <- row$parameter_overrides[[1]]
    while (is.list(overrides) && length(overrides) == 1L) {
      overrides <- overrides[[1]]
    }
    if (length(overrides) > 0) {
      if (is.null(names(overrides)) || any(names(overrides) == "")) {
        stop("Every scenario parameter override must be named.", call. = FALSE)
      }
      scenario_parameters[names(overrides)] <- as.numeric(overrides)
    }

    grids_i <- list(
      t2d = annual_transition_grids$t2d[[imputation]],
      obesity = annual_transition_grids$obesity[[imputation]],
      cvd = annual_transition_grids$cvd[[imputation]]
    )

    model_output <- run_psymetric_markov(
      person_data = analysis_datasets[[imputation]],
      transition_grids = grids_i,
      t2d_complication_grid = t2d_complication_grid,
      parameters = scenario_parameters,
      mortality_values = mortality_values,
      mortality_cause_allocation_values = mortality_cause_allocation_values,
      utility_values = utility_values,
      strategy = row$strategy[[1]],
      threshold_value = row$threshold_value[[1]],
      include_operational_cost = row$include_operational_cost[[1]]
    )

    outputs[[i]] <- model_output$summary |>
      dplyr::mutate(
        imputation = imputation,
        analysis_id = row$analysis_id[[1]],
        analysis_label = row$analysis_label[[1]],
        .before = strategy
      )
  }

  dplyr::bind_rows(outputs)
}

pool_across_imputations <- function(summary_by_imputation) {
  summary_by_imputation |>
    dplyr::group_by(
      analysis_id, analysis_label,
      strategy, threshold_value, threshold_label
    ) |>
    dplyr::summarise(
      n_imputations = dplyr::n_distinct(imputation),
      mean_n_modelled = mean(n_modelled),
      mean_n_threshold_positive = mean(n_threshold_positive),
      mean_n_strategy_eligible = mean(n_strategy_eligible),
      dplyr::across(
        dplyr::starts_with("mean_discounted_"),
        ~ mean(.x, na.rm = TRUE)
      ),
      dplyr::across(
        dplyr::ends_with("_per_1000"),
        ~ mean(.x, na.rm = TRUE)
      ),
      max_abs_probability_mass_error =
        max(max_abs_probability_mass_error, na.rm = TRUE),
      .groups = "drop"
    )
}

add_incremental_results <- function(summary_table, parameters) {
  usual_care <- summary_table |>
    dplyr::filter(strategy == "usual_care") |>
    dplyr::select(
      analysis_id,
      usual_costs = mean_discounted_costs,
      usual_qalys = mean_discounted_qalys
    )

  universal_metformin <- summary_table |>
    dplyr::filter(strategy == "universal_metformin") |>
    dplyr::select(
      analysis_id,
      universal_costs = mean_discounted_costs,
      universal_qalys = mean_discounted_qalys
    )

  strategy_3 <- summary_table |>
    dplyr::filter(strategy == "psymetric_lifestyle_metformin") |>
    dplyr::select(
      analysis_id,
      threshold_value,
      strategy_3_costs = mean_discounted_costs,
      strategy_3_qalys = mean_discounted_qalys
    )

  summary_table |>
    dplyr::left_join(usual_care, by = "analysis_id") |>
    dplyr::left_join(universal_metformin, by = "analysis_id") |>
    dplyr::left_join(strategy_3, by = c("analysis_id", "threshold_value")) |>
    dplyr::mutate(
      incremental_costs_vs_usual = mean_discounted_costs - usual_costs,
      incremental_qalys_vs_usual = mean_discounted_qalys - usual_qalys,
      icer_vs_usual = dplyr::if_else(
        abs(incremental_qalys_vs_usual) < 1e-12,
        NA_real_,
        incremental_costs_vs_usual / incremental_qalys_vs_usual
      ),
      nmb_vs_usual_25000 =
        parameter_value(parameters, "wtp_threshold_primary") *
        incremental_qalys_vs_usual - incremental_costs_vs_usual,
      nmb_vs_usual_35000 =
        parameter_value_default(
          parameters, "wtp_threshold_secondary", default = 35000
        ) * incremental_qalys_vs_usual - incremental_costs_vs_usual,
      incremental_costs_vs_universal =
        mean_discounted_costs - universal_costs,
      incremental_qalys_vs_universal =
        mean_discounted_qalys - universal_qalys,
      icer_vs_universal = dplyr::if_else(
        abs(incremental_qalys_vs_universal) < 1e-12,
        NA_real_,
        incremental_costs_vs_universal / incremental_qalys_vs_universal
      ),
      nmb_vs_universal_25000 =
        parameter_value(parameters, "wtp_threshold_primary") *
        incremental_qalys_vs_universal - incremental_costs_vs_universal,
      incremental_costs_vs_strategy_3 =
        mean_discounted_costs - strategy_3_costs,
      incremental_qalys_vs_strategy_3 =
        mean_discounted_qalys - strategy_3_qalys,
      icer_vs_strategy_3 = dplyr::if_else(
        abs(incremental_qalys_vs_strategy_3) < 1e-12,
        NA_real_,
        incremental_costs_vs_strategy_3 / incremental_qalys_vs_strategy_3
      ),
      nmb_vs_strategy_3_25000 =
        parameter_value(parameters, "wtp_threshold_primary") *
        incremental_qalys_vs_strategy_3 - incremental_costs_vs_strategy_3
    )
}

# ==============================================================================
# 4. Deterministic base case
# ==============================================================================

risk_thresholds <- c(0.05, 0.10, 0.15, 0.20)
base_case_strategies <- c(
  "usual_care",
  "universal_metformin",
  "psymetric_lifestyle",
  "psymetric_metformin",
  "psymetric_lifestyle_metformin"
)

base_case_grid <- tidyr::expand_grid(
  imputation = seq_along(analysis_datasets),
  strategy = base_case_strategies,
  threshold_value = risk_thresholds
) |>
  dplyr::filter(
    !strategy %in% c("usual_care", "universal_metformin") |
      threshold_value == min(risk_thresholds)
  ) |>
  dplyr::mutate(
    analysis_id = "base_case",
    analysis_label = "Deterministic base case",
    include_operational_cost = FALSE,
    parameter_overrides = rep(
      list(stats::setNames(numeric(0), character(0))),
      dplyr::n()
    )
  )

base_case_by_imputation <- run_model_grid(
  scenario_grid = base_case_grid,
  analysis_datasets = analysis_datasets,
  annual_transition_grids = annual_transition_grids,
  t2d_complication_grid = t2d_complication_grid,
  parameters = external_parameter_lookup,
  mortality_values = all_cause_mortality_values,
  mortality_cause_allocation_values = mortality_cause_allocation_values,
  utility_values = health_state_utility_values
)

base_case_pooled <- pool_across_imputations(base_case_by_imputation)
base_case_incremental <- add_incremental_results(
  base_case_pooled,
  external_parameter_lookup
)


# ==============================================================================
# 5. Deterministic scenario analyses described in the manuscript
# ==============================================================================

empty_overrides <- stats::setNames(numeric(0), character(0))

scenario_definitions <- tibble::tibble(
  analysis_id = c(
    "glp1ra_two_year",
    "glp1ra_lifelong_low_cost",
    "glp1ra_lifelong_base_cost",
    "glp1ra_lifelong_high_cost",
    "operational_cost_250",
    "operational_cost_500",
    "operational_cost_1000",
    "antipsychotic_switch"
  ),
  analysis_label = c(
    "Lifestyle plus GLP-1RA; two-year treatment at current cost",
    "Lifestyle plus lifelong GLP-1RA; low off-patent cost after year 2",
    "Lifestyle plus lifelong GLP-1RA; base cost after year 2",
    "Lifestyle plus lifelong GLP-1RA; high/current cost after year 2",
    "Annual PsyMetRiC operational cost; 250 eligible people per Trust-year",
    "Annual PsyMetRiC operational cost; 500 eligible people per Trust-year",
    "Annual PsyMetRiC operational cost; 1,000 eligible people per Trust-year",
    "PsyMetRiC-guided switch from olanzapine to aripiprazole"
  ),
  candidate_strategy = c(
    "psymetric_lifestyle_glp1ra",
    "psymetric_lifestyle_glp1ra_lifelong",
    "psymetric_lifestyle_glp1ra_lifelong",
    "psymetric_lifestyle_glp1ra_lifelong",
    "psymetric_lifestyle_metformin",
    "psymetric_lifestyle_metformin",
    "psymetric_lifestyle_metformin",
    "psymetric_antipsychotic_switch"
  ),
  include_operational_cost = c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE),
  parameter_overrides = list(
    empty_overrides,
    c(
      glp1ra_treatment_duration_years_lifelong_scenario = 999,
      glp1ra_effect_duration_years_lifelong_scenario = 999,
      cost_glp1ra_total_annual_after_year2 = parameter_value_default(
        external_parameter_lookup,
        "cost_glp1ra_total_annual_after_year2_low",
        default = 500
      )
    ),
    c(
      glp1ra_treatment_duration_years_lifelong_scenario = 999,
      glp1ra_effect_duration_years_lifelong_scenario = 999,
      cost_glp1ra_total_annual_after_year2 = parameter_value_default(
        external_parameter_lookup,
        "cost_glp1ra_total_annual_after_year2",
        default = 1200
      )
    ),
    c(
      glp1ra_treatment_duration_years_lifelong_scenario = 999,
      glp1ra_effect_duration_years_lifelong_scenario = 999,
      cost_glp1ra_total_annual_after_year2 = parameter_value_default(
        external_parameter_lookup,
        "cost_glp1ra_total_annual_after_year2_high",
        default = 2435
      )
    ),
    c(
      cost_psymetric_operational_annual_per_trust_scenario = 8000,
      psymetric_operational_denominator_per_trust_year = 250,
      psymetric_operational_cost_eligibility_age_lower = 16,
      psymetric_operational_cost_eligibility_age_upper_exclusive = 36
    ),
    c(
      cost_psymetric_operational_annual_per_trust_scenario = 8000,
      psymetric_operational_denominator_per_trust_year = 500,
      psymetric_operational_cost_eligibility_age_lower = 16,
      psymetric_operational_cost_eligibility_age_upper_exclusive = 36
    ),
    c(
      cost_psymetric_operational_annual_per_trust_scenario = 8000,
      psymetric_operational_denominator_per_trust_year = 1000,
      psymetric_operational_cost_eligibility_age_lower = 16,
      psymetric_operational_cost_eligibility_age_upper_exclusive = 36
    ),
    empty_overrides
  )
)

scenario_grid <- purrr::pmap_dfr(
  scenario_definitions,
  function(
      analysis_id,
      analysis_label,
      candidate_strategy,
      include_operational_cost,
      parameter_overrides
  ) {
    comparison_strategies <- unique(c(
      "usual_care",
      "universal_metformin",
      "psymetric_lifestyle_metformin",
      candidate_strategy
    ))

    tidyr::expand_grid(
      imputation = 1L,
      strategy = comparison_strategies,
      threshold_value = c(0.05, 0.10)
    ) |>
      dplyr::filter(
        !strategy %in% c("usual_care", "universal_metformin") |
          threshold_value == 0.05
      ) |>
      dplyr::mutate(
        analysis_id = .env$analysis_id,
        analysis_label = .env$analysis_label,
        include_operational_cost =
          .env$include_operational_cost &
          strategy == "psymetric_lifestyle_metformin",
        parameter_overrides = rep(list(.env$parameter_overrides), dplyr::n())
      )
  }
)

scenario_results <- run_model_grid(
  scenario_grid = scenario_grid,
  analysis_datasets = analysis_datasets,
  annual_transition_grids = annual_transition_grids,
  t2d_complication_grid = t2d_complication_grid,
  parameters = external_parameter_lookup,
  mortality_values = all_cause_mortality_values,
  mortality_cause_allocation_values = mortality_cause_allocation_values,
  utility_values = health_state_utility_values
)

scenario_incremental <- add_incremental_results(
  scenario_results,
  external_parameter_lookup
)


# ==============================================================================
# 6. Objects returned to the R session
# ==============================================================================

psymetric_markov_results <- list(
  base_case_by_imputation = base_case_by_imputation,
  base_case_pooled = base_case_pooled,
  base_case_incremental = base_case_incremental,
  scenario_results = scenario_results,
  scenario_incremental = scenario_incremental
)

psymetric_markov_results
