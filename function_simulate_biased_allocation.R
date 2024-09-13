library(coxed)
library(ggplot2)
library(tibble)
library(dplyr)
library(stringr)
library(here)
library(purrr)
library(splines)

simulate_biased_allocation <- function(seed) {
  
  set.seed(seed)
  
  # Set the baseline hazard function, simulate MELD scores according to the
  # empirical distribution of MELD in ET. Note that MELD is assumed to be
  # independent of sex. Then simulate data based on a Cox proportional hazards
  # model. Parameters for this Cox proportional hazards model were obtained
  # by estimating such a model on ET registry data.
  MELD_SCORES <- data.frame(
    meld_lab = sample(
      as.numeric(names(EMPIRICAL_DISTRIBUTION_MELD)), 
      prob=EMPIRICAL_DISTRIBUTION_MELD,
      size=TOTAL_NUMBER_OF_CANDIDATES, 
      replace=TRUE)
  )
  sim <- sim.survdata(
    N=TOTAL_NUMBER_OF_CANDIDATES, 
    `T`=TIME_HORIZON+1,
    num.data.frames = 1,
    X=MELD_SCORES-18,
    hazard.fun = function(t) {0.0010954},
    beta = exp(0.18)-1
  )
  simdata <- sim$data
  simdata$meld_lab <- simdata$meld_lab + 18
  
  # Retrieve the simulated data. Censor candidates at 10 years of waiting time.
  # Set sex for candidates randomly.
  simdata$candidate_bloodgroup <- sample(
    factor(names(BLOOD_TYPE_FREQUENCIES)), size = nrow(simdata),
    BLOOD_TYPE_FREQUENCIES, replace=TRUE
  )
  simdata$candidate_sex <- factor(
    sample(x=c('Male', 'Female'), size=nrow(simdata), prob = c(0.36, 0.64), replace=TRUE),
    levels = c('Male', 'Female')
  )
  simdata$candidate_height <- if_else(
    simdata$candidate_sex == 'Male',
    rnorm(nrow(simdata), mean = 178, sd = 7),
    rnorm(nrow(simdata), mean = 163, sd = 7)
  )
  simdata$candidate_country <- sample(factor(1:7), size = nrow(simdata), prob = P_COUNTRY_RATES, replace=TRUE)
  
  # Simulate donors. Assume that there is a donor shortage of 30%. Also assume
  # that donors are reported uniformly over the simulation period.
  donor_data = data.frame(
    bloodgroups_donor = sample(
      factor(names(BLOOD_TYPE_FREQUENCIES)),
      size = NUMBER_OF_DONORS,
      BLOOD_TYPE_FREQUENCIES,
      replace=TRUE
    ),
    donor_countries = sample(
      factor(1:7), size = NUMBER_OF_DONORS, prob = D_COUNTRY_RATES, 
      replace=TRUE)
  )
  
  # For candidates on the initial waitlist, a registration time is drawn
  # such that these candidates have not yet died at simulation start (t=0),
  # and such that they were registered before simulation start.
  # For the remaining candidates, we assume that they were registered uniformly
  # over the simulation period.
  simdata$candidate_listing_time <- if_else(
    runif(nrow(simdata)) < INITIAL_WAITLIST_SIZE/TOTAL_NUMBER_OF_CANDIDATES,
    runif(nrow(simdata), min=-simdata$y, max=0),
    round(runif(nrow(simdata), min = 0, max = TIME_HORIZON), digits=1)
  )
  simdata$candidate_death_time <- simdata$candidate_listing_time + as.numeric(simdata$y)
  simdata$transplanted <- FALSE
  
  # Make the simulation data into a data.frame, and order it descendingly
  # by MELD. Also set an ID to quickly identify needed candidates.
  setDT(simdata)
  setorder(simdata, -meld_lab, candidate_listing_time)
  simdata$id <- 1:nrow(simdata)
  
  # Prepare fixed effects for offer acceptance decisions
  X <- model.matrix(acceptance_model_formula, data=simdata)[
    , names(fixed_effects)
  ]
  simdata$lp_fixed_effects <- X %*% fixed_effects
  
  # Assume donors are reported uniformly over the simulation period.
  times_donor_reported <- round(
    runif(NUMBER_OF_DONORS, min = 0, max = TIME_HORIZON),
    digits = 2
  ) |> sort()
  
  # Allow organ offer acceptance decisions to be correlated within each donor.
  random_effects_donor_acceptance <- rnorm(
    NUMBER_OF_DONORS, sd = variance_component_donor_id
  )
  
  # Simulate the allocation of the donors.
  # Simulate the allocation of the donors.
  for (i in seq_along(times_donor_reported)) {
    
    # Get the time at which the donor was reported.
    d_reporting_time <- times_donor_reported[[i]]
    d_bloodgroup <- donor_data$bloodgroups_donor[[i]]
    d_country <- donor_data$donor_countries[[i]]
    d_random_effect <- random_effects_donor_acceptance[[i]]
    
    # Eligible candidates have to (i) have a waitlist registration,
    # (ii) not be transplanted already, (iii) listed in the same country as the donor,
    # and (iv) be blood group compatible with the donor
    elig_patients <- simdata[
      between(d_reporting_time, candidate_listing_time, candidate_death_time) &
        !transplanted &
        candidate_country == d_country &
        candidate_bloodgroup %in% BLOOD_TYPE_COMPATIBILITIES[[d_bloodgroup]]
    ]
    
    # Acceptance probabilities for patients
    elig_patients$prob_accept <- boot::inv.logit(
      elig_patients$lp_fixed_effects + d_random_effect
    )
    
    # Simulate acceptance for the candidates. The first candidate simulated
    # to accept the graft gets transplanted.
    u <- runif(n = nrow(elig_patients))
    recipient <- elig_patients$id[which.min(elig_patients$prob_accept <= u)]
    simdata[recipient, c('candidate_transplant_time', 'transplanted') := list(d_reporting_time, TRUE)]
  }
  
  # Construct censoring time
  simdata$candidate_censoring_time <- pmin(
    TIME_HORIZON,
    simdata$candidate_death_time,
    simdata$candidate_transplant_time,
    na.rm=TRUE)
  
  # Construct observed outcomes
  simdata$waitlist_death <- as.integer(!simdata$transplanted & simdata$failed)
  simdata$time_to_event <- simdata$candidate_censoring_time - simdata$candidate_listing_time
  
  # Prepare a truncated dataset
  simdata_truncated <- simdata |>
    mutate(
      transplanted_within_90d = if_else(time_to_event >= 90, 0L, transplanted),
      death_within_90d = if_else(time_to_event >= 90, 0L, waitlist_death),
      time_to_event_within_90d = pmin(time_to_event, 90)
    ) |>
    select(meld_lab, candidate_sex, transplanted_within_90d, death_within_90d, time_to_event_within_90d)
  
  
  # Prepare a truncated dataset
  return(coxph(
    Surv(time_to_event_within_90d, death_within_90d) ~ meld_lab + candidate_sex,
    data=simdata_truncated
  ) |>
    broom::tidy(conf.int=TRUE, digits=3) |>
    group_by(row_number(), term) |>
    summarize(
      estimate = glue::glue('{round(estimate, 3)} [{round(conf.low, 3)}-{round(conf.high, 3)}]'),
      p.value=p.value
    ) |>
    mutate(
      `True parameter` = recode(term, meld_lab = 0.18, candidate_sexFemale = 0), .after = term
    ) |>
    select(-1)
  )
}



