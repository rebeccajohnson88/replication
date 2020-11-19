
# Replication code for:
# What is Your Estimand? Defining the Target Quantity Connects Statistical Evidence to Theory
# Ian Lundberg, Rebecca Johnson, and Brandon Stewart
# Email: ilundberg@princeton.edu

# This file carries out the analysis for
# Concrete Estimation Example: The Family Gap in Pay

# The data for this example are available from IPUMS at
# https://cps.ipums.org/cps/
# 1. Register for an account
# 2. Put the following variables in your cart
# YEAR, SERIAL, MONTH, CPSID, ASECFLAG, ASECWTH,
# PERNUM, CPSIDP, ASECWT, RELATE, AGE, SEX, RACE,
# MARST, NCHILD, AHRSWORKT, EDUC, EARNWT, REPWTP,
# CLASSWLY, WKSWORK1, WKSWORK2, UHRSWORKLY, INCWAGE
# 3. Select the sample "IPUMS-CPS, ASEC 2019"
# 4. Select data format ".dta (Stata)"
# 5. Select cases on AGE (ages 25-44 only) and SEX (Female only)
# 6. Download the data and place it in the data subdirectory.
# You will need to rename to "cps_00040.dta" for the code below or 
# change the name of the RAW_DATA_NAME parameter

library(readstata13)
library(tidyverse)
library(reshape2)
library(foreach)
library(mgcv)


## Constants
OUTPUT_NAME= "PalWaldfogel_output.txt"

## Can either read raw data/clean or read in cleaned data
READ_RAW_DATA <- TRUE
RAW_DATA_NAME <- "cps_00040.dta"
CLEAN_DATA_NAME <- "pw_analytic.csv"
RUN_VAR_ESTIMATION_AGGREGATE <- FALSE
RUN_VAR_ESTIMATE_SUBGROUP <- FALSE


# Function to open your default browser to check when
# coding cps vars
lookup_cpsvar_inbrowser <- function(varname){
  browseURL(sprintf("https://cps.ipums.org/cps-action/variables/%s#codes_section",
                    varname))
}


## Initiate output file
sink(sprintf("output/%s", output_name))
print("Output from replication of Pal and Waldfogel")


###################################
# Prep the main data file on women 
# ages 25-44 
###################################

if(READ_RAW_DATA){
  
  ## Read in stata file
  data <- read.dta13(sprintf("data/%s", RAW_DATA_NAME),
                     convert.factors = F)
  
  # Make a data frame of the replicate weights used for variance estimation
  rep_weights <- data %>% select(pernum, starts_with("repwtp"))
  
  d_clean <- data %>%
    
    # Create different flags to filter sample
    
    ### Flags to filter to observations with person-level weight > 0 
    ### (all in this case but may differ across waves) and
    ### Filter to 25-44 (inclusive) and sex == female (2) 
    mutate(filter_isposweight = asecwt > 0,
           filter_indemrange = age >= 25 & age <= 44 & sex == 2) %>%
    
    ## apply filter before calculating income quantiles
    filter(filter_isposweight & filter_indemrange) %>%
  
    ## Mark missing incomes
    ### classwly is job codes (https://cps.ipums.org/cps-action/variables/CLASSWLY#codes_section)
    ## and >14 are wage/salary jobs
    mutate(derived_incwage = case_when(classwly > 14 ~ incwage,
                                       TRUE ~ NA_integer_), 
           ### Even for those with jobs in those categories, set to missing if has
           ### missing codes for those wages
           derived_incwage = case_when(derived_incwage != 9999999 & derived_incwage != 9999998 ~ derived_incwage,
                                       TRUE ~ NA_integer_),
           
           ## Weeks worked last year
           ### For the other years, we can just use wkswork1, which is the continuous report
           derived_wkswork1 = ifelse(wkswork1 <= 0, NA, wkswork1),
           
           ### Usual hours per week worked last year
           derived_uhrsworkly = ifelse(uhrsworkly != 999, uhrsworkly, NA),
           
           ## Create hourly wages as total income/wage divided by # of weeks x usual weekly hours
           derived_wage = derived_incwage / (derived_wkswork1 * derived_uhrsworkly),
           
           ## Truncate log wage at 1st and 99th percentile
           derived_ln_wage = log(case_when(derived_wage < quantile(derived_wage, .01, na.rm = T) ~ quantile(derived_wage, .01, na.rm = T),
                                derived_wage > quantile(derived_wage, .99, na.rm = T) ~ quantile(derived_wage, .99, na.rm = T),
                                T ~ derived_wage)),
           ## Create controls
           
           ### Code education into four levels 
           derived_educ = factor(ifelse(educ == 1 | educ == 999, NA,
                                        ## Less than high school
                                        ifelse(educ <= 60 | educ == 71, 1,
                                               ## HS degree (include diploma unclear 72)
                                               ifelse(educ == 70 | educ == 72 | educ == 73, 2,
                                                      ## Some college
                                                      ifelse(educ < 110, 3,
                                                             ## College or more
                                                             4))))),
           
           ### Family status is whether married with spouse present (TRUE),
           ### other categories (FALSE) or NA
           derived_married = ifelse(marst == 9, NA, marst == 1),
           
           ### Race is white black or other, with later including mixed race
           derived_race = factor(ifelse(race == 100, 1,
                                        ifelse(race == 200, 2, 3)),
                                 labels = c("White","Black","Other")),
           
           ## Mother is TRUE if nchild > 0, false otherwise
           derived_mother = (nchild > 0),
           
           ## For later analysis, create age^squared
           derived_agesq = age^2,
           
           ## Filters based on those covariates
           filter_nonmisseduc = !is.na(derived_educ),
           filter_nonmisswage = !is.na(derived_ln_wage)
           ) 
  
  ## Variables to include
  weights_identifiers = c("pernum", "serial", "asecwt",
                          grep("repwtp", colnames(d_clean),
                               value = TRUE)) 
  outcome = "derived_ln_wage"
  treatment = "derived_mother"
  covariates = c(sprintf("derived_%s", 
                         c("educ", "married",
                           "race", "agesq")), "age")
  filters = grep("filter", colnames(d_clean), value = TRUE)
  cols_keep = c(weights_identifiers, outcome, treatment, covariates, 
                filters)
  
  d_all <- d_clean %>%
    select(all_of(cols_keep))
  
  write.csv(d_all, sprintf("data/clean/%s", CLEAN_DATA_NAME))
} else{
  d_all <- read.csv(sprintf("data/clean/%s", CLEAN_DATA_NAME))
}

filter_vars = grep("filter", colnames(d_all), value = TRUE)
print("Sample sizes")
colSums(d_all[, filter_vars])

## Create analytic df (d) for people who pass 
## all filters (sum of true across filter vars == 
## length of filter vars)
d <- d_all[rowSums(d_all[,filter_vars]) == length(filter_vars), ]

###################################
# Note the common support problem in#
###################################

strata_vars = c("age", "derived_educ", "derived_race",
                     "derived_married")
treatment_var = "derived_mother"

support_data <- d %>%
  
  ## group by strata of covariates to check support within
  group_by(across(all_of(strata_vars))) %>%
  
  ## indicator variable checking if both levels of 
  ## treatment var have support (TRUE) otherwise false
  mutate(has_support = n_distinct(!!sym(treatment_var)) == 2) %>%
  filter(1:n() == 1) %>%
  select(all_of(c(strata_vars, "has_support"))) %>% 
  
  ## add support flags onto the main data
  right_join(d_all, by = strata_vars) %>% 
  mutate(has_support = ifelse(is.na(has_support),F,has_support))

print("Note proportion of mothers and non-mothers in the region of common support")
print(support_data %>%
        group_by(!!sym(treatment_var)) %>%
        summarize(has_support = weighted.mean(has_support, w = asecwt)))

print("Analytical sample is in the region of common support")
print(paste0("N_Total = ", 
             sum(support_data$has_support)))
print(paste0("N_Employed = ", 
             sum(support_data$has_support & !is.na(support_data$derived_ln_wage))))

print("These groups are off the region of common support")
print(data.frame(support_data %>%
                   filter(!has_support) %>%
                   select(all_of(c(treatment_var, strata_vars))) %>%
                   arrange(!!!syms(c(treatment_var, strata_vars))) %>%
                   distinct()), #rj note: i think previous was printing obs without common support; 
                                # seems like rather than obs we want all groups across the strata vars
      # The max is a very high number to tell it to print all rows
      max = 9999)

# All remaining analyses will focus on those with support only

#################
# Aggregate gap #
################

predict_finddiff <- function(one_model,
                             df_topred, treatment_varname,
                             aggregate = TRUE,
                             focal_covar = "age"){
  
  ## predict in actual group
  tmp = df_topred
  tmp[[treatment_varname]] <- TRUE
  pred_actual_focal= predict(one_model,
                             newdata = tmp)
  
  ## predict in contrast group
  tmp = df_topred
  tmp[[treatment_varname]] <- FALSE
  pred_actual_other = predict(one_model,
                              newdata = tmp)
  
  
  ## estimate is weighted difference between 
  ## predicted in focal and predicted in contrast
  ## (in our case, predicted wages for mothers with observed covar 
  ## values versus predicted wages for nonmothers with observed covar values of 
  ## mothers)
  if(aggregate){
    estimate <- weighted.mean(pred_actual_focal - pred_actual_other,
                             w = df_topred$weight)
  } else{
    estimate <- df_topred %>%
          mutate(Estimate = pred_actual_focal - pred_actual_other) %>%
          select(age, Estimate)
    
  }
  
  return(estimate)
}

fit_mods  <- function(full_data, 
                      filter_flags,
                      weight_varname,
                      treatment_varname, outcome_varname,
                      control_vec,
                      verbose = FALSE,
                      return_modfits = FALSE) {
  
  d_case <- full_data[rowSums(full_data[ , filter_flags]) == length(filter_flags),
                      ]
  
  d_case$weight <- d_case[[weight_varname]]
  
  # Normalize weights
  d_case$weight <- d_case$weight / sum(d_case$weight)
  
  # Separate vector of controls into different combinations (specific to application)
  control_vec_group = setdiff(control_vec, "derived_agesq")
  control_vec_interact = c("age", "derived_agesq")
  control_vec_nointeract = setdiff(control_vec, control_vec_interact)
  
  # Unadjusted difference in outcomes between tx and control (wages b/t mothers and nonmothers)
  unadjusted <- d_case %>%
    group_by(!!sym(treatment_varname)) %>%
    summarize(Estimate = weighted.mean(!!sym(outcome_varname),
                                       w = weight, na.rm = T), .groups= "drop") %>%
    spread(key = !!sym(treatment_varname), value = Estimate) %>%
    mutate(approach = "Unadjusted",
           Estimate = `TRUE` - `FALSE`) %>%
    select(approach, Estimate)

  # Adjusted difference in outcomes between treatment and control
  # defined by strata; weight by n treat in that strata 
  nonparametric_adjusted <- d_case %>%
    
    ## First, group by strata defined by treatment varname and control variables
    group_by(.dots = c(treatment_varname, control_vec_group)) %>%
    summarize(estimate = weighted.mean(!!sym(outcome_varname), w = weight, na.rm = T),
              weight = sum(weight), .groups = "drop") %>%
    
    ## Then, group by control variables only
    group_by(.dots = control_vec_group) %>%
    
    ## Give each stratum the total weight of treatment == 1 in this stratum
    mutate(stratum_weight = sum(weight * !!sym(treatment_varname))) %>%
    
    ## Then, group by treatment group
    group_by(!!sym(treatment_varname)) %>%
    
    ## Estimate is cell mean weighted by weight for that stratum
    summarize(Estimate = weighted.mean(estimate, w = stratum_weight),
              .groups = "drop") %>%
    
    ## Reshape to wide form to get mothers (true) - non-mothers (false)
    spread(key = !!sym(treatment_varname), value = Estimate) %>%
    mutate(approach = "Stratification",
           Estimate = `TRUE` - `FALSE`) %>%
    select(approach, Estimate)
  
  
  # Model-based approaches
  
  ## OLS all vars as additive
  ols_additive_fit <- lm(formula(sprintf("%s ~ %s + %s",
                                 outcome_varname,
                                 treatment_varname,
                                 paste(control_vec,
                                   collapse = "+"))),
                         data = d_case,
                         weights = weight)
  
  ## OLS interact treatment with age and age squared as continuous
  ols_interactive_fit <- lm(formula(sprintf("%s ~ %s + %s",
                                    outcome_varname,
                                    paste(sprintf("%s*%s", 
                                          treatment_varname, control_vec_interact),
                                          collapse = "+"),
                                    paste(control_vec_nointeract,
                                          collapse = "+"))),
                            data = d_case,
                            weights = weight)
  
  ## OLS interact with age categorical
  ols_ageFactor_fit <- lm(formula(sprintf("%s ~ %s + %s",
                                    outcome_varname,
                                    
                                    ## interact tx with factor
                                    sprintf("%s*factor(%s)", 
                                    treatment_varname, 
                                    setdiff(control_vec_interact, "derived_agesq")),
                                    
                                    ## other controls
                                    paste(control_vec_nointeract,
                                    collapse = "+"))),
                          data = d_case,
                          weights = weight)
  
  ## Gam with group-specific smoother for age
  ### Needs factor rather than logical vesion of treatment
  tmp = d_case
  tmp[[treatment_var]] = factor(tmp[[treatment_var]])
  gam_interactive_fit <- gam(formula(sprintf("%s ~ %s + s(%s, by = %s) + %s",
                                    outcome_varname,
                                    treatment_varname,
                                    setdiff(control_vec_interact, "derived_agesq"),
                                    treatment_varname,
                                    paste(control_vec_nointeract,
                                    collapse = "+"))),
                             data = tmp %>% group_by(),
                             weight = tmp$weight)
  
  # Combine models into a list and generate predictions
  all_mod = list(ols_additive_fit = ols_additive_fit,
                 ols_interactive_fit = ols_interactive_fit,
                 ols_ageFactor_fit = ols_ageFactor_fit,
                 gam_interactive_fit = gam_interactive_fit)
  
  ## filter to whatever group is focal in contrast
  df_topred = d_case %>% group_by() %>% 
    filter(!!sym(treatment_varname))
  
  ## predict
  all_pred = lapply(all_mod, predict_finddiff, 
                    df_topred = df_topred,
                    treatment_varname = treatment_varname)
  
  ## bind model-based estimates with each other
  pred_df = do.call(rbind.data.frame, all_pred)  
  colnames(pred_df) = "Estimate"
  pred_df$approach = names(all_mod)
  
  ## bind with non-model based estimates
  all_est = rbind.data.frame(unadjusted, nonparametric_adjusted,
                             pred_df)
  if(verbose){
    print("Completed one iteration of estimation")
  } 
  if(return_modfits){
    return(list(estimates = all_est, models = all_mod))
  } else{
    return(all_est)
  }
}

## run once with case weight to get point estimates
estimate_aggregate_point_all <- fit_mods(full_data = support_data,
                                    filter_flags = c("has_support"),
                                    weight_varname = "asecwt",
                                    treatment_varname = "derived_mother",
                                    outcome_varname = "derived_ln_wage",
                                    control_vec = c(strata_vars, "derived_agesq"),
                                    return_modfits = TRUE)

  
## iterate through and run w/ each of the replicate weights for later variance estimation
if(RUN_VAR_ESTIMATION_AGGREGATE){
  replicates_aggregate <- foreach(i = 1:160) %do% {
    fit_mods(weight_varname = paste0("repwtp",i),
             full_data= support_data,
             filter_flags = c("has_support"),
             treatment_varname = "derived_mother",
             outcome_varname = "derived_ln_wage",
             control_vec = c(strata_vars, "derived_agesq"),
             verbose = TRUE,
             return_modfits = TRUE)
  }
  saveRDS(replicates_aggregate, "int_objects/replicates_aggregate.RDS", compress = FALSE)
} else{
  replicates_aggregate <- readRDS("int_objects/replicates_aggregate.RDS")
}

## join point estimates onto variance estimates and
## calculate variance
estimate_aggregate_variance <- do.call(rbind.data.frame, 
                            lapply(replicates_aggregate, function(x) x$estimates)) %>%
  rename(Estimate_star = Estimate) %>%
  left_join(estimate_aggregate_point_all$estimates,
            by = "approach") %>%
  group_by(approach) %>%
  
  ## variance according to ipums methodology
  ## w/ fay's variance estimation method: https://cps.ipums.org/cps/repwt.shtml
  ## if using with different weights structure, should check on survey's
  ## specific methodology for repliace weights
  summarize(variance = 4 / 160 * sum((Estimate_star - Estimate) ^ 2), 
            .groups = "drop")

aggregate_results <- estimate_aggregate_point_all$estimates %>%
  filter(!grepl("Unadjusted", approach)) %>% # filter out unadjusted gap from estimates we present
  left_join(estimate_aggregate_variance, by = "approach") %>%
  mutate(approach = factor(case_when(grepl("Stratification", approach, ignore.case = TRUE) ~ 1,
                                     grepl("ols.*factor", approach, ignore.case = TRUE) ~ 2,
                                     grepl("gam.*interactive", approach, ignore.case = TRUE) ~ 3,
                                     grepl("ols.*interactive", approach, ignore.case = TRUE) ~ 4,
                                     grepl("ols.*additive", approach, ignore.case = TRUE) ~ 5),
                           labels = c("Stratification:\nNo estimation assumptions",
                                      "Education + Marital + Race + \n(Age indicators x Motherhood)",
                                      "Education + Marital + Race +\n(Age smooth x Motherhood)",
                                      "Education + Marital + Race +\n(Age quadratic x Motherhood)",
                                      "Education + Marital + Race +\nAge quadratic + Motherhood")))

##################################
# Estimates within age subgroups #
##################################

gen_pred_subgroup  <- function(full_data, 
                               filter_flags,
                               weight_varname,
                               treatment_varname, outcome_varname,
                               control_vec,
                               models_aggregate,
                               verbose = FALSE) {
  
  d_case <- full_data[rowSums(full_data[ , filter_flags]) == length(filter_flags), ]

  d_case$weight <- d_case[[weight_varname]]
  
  # Normalize weights
  d_case$weight <- d_case$weight / sum(d_case$weight)
  
  # Separate vector of controls into different combinations (specific to application)
  control_vec_group = setdiff(control_vec, "derived_agesq")
  control_vec_interact = c("age", "derived_agesq")
  focal_covar = "age"
  
  # Unadjusted difference in outcomes between tx and control (wages b/t mothers and nonmothers)
  unadjusted <- d_case %>%
    group_by(.dots = c(focal_covar, treatment_varname)) %>%
    summarize(Estimate = weighted.mean(!!sym(outcome_varname),
                                       w = weight, na.rm = T), .groups= "drop") %>%
    spread(key = !!sym(treatment_varname), value = Estimate) %>%
    mutate(approach = "Unadjusted",
           Estimate = `TRUE` - `FALSE`) %>%
    select(!!sym(focal_covar), approach, Estimate)
  
  # Adjusted difference in outcomes between treatment and control
  # defined by strata; weight by n treat in that strata 
  nonparametric_adjusted <- d_case %>%
    
    ## First, group by strata defined by treatment varname and control variables
    group_by(.dots = c(treatment_varname, control_vec_group)) %>%
    summarize(estimate = weighted.mean(!!sym(outcome_varname), w = weight, na.rm = T),
              weight = sum(weight), .groups = "drop") %>%
    
    ## Then, group by control variables only
    group_by(.dots = control_vec_group) %>%
    
    ## Give each stratum the total weight of treatment == 1 in this stratum
    mutate(stratum_weight = sum(weight * !!sym(treatment_varname))) %>%
    
    ## Then, group by treatment group and focal var
    group_by(.dots = c(focal_covar, treatment_varname)) %>%
    
    ## Estimate is cell mean weighted by weight for that stratum
    summarize(Estimate = weighted.mean(estimate, w = stratum_weight),
              .groups = "drop") %>%
    
    ## Reshape to wide form to get mothers (true) - non-mothers (false)
    spread(key = !!sym(treatment_varname), value = Estimate) %>%
    mutate(approach = "Stratification",
           Estimate = `TRUE` - `FALSE`) %>%
    select(!!sym(focal_covar), approach, Estimate)
  
  ## Create new prediction dataset grouped by focal variable
  # rj note: below returns id df as calling by name but
  # 0 gap predictions for a few ages, so not sure why
  df_topred = d_case %>% 
      group_by(.dots = focal_covar) %>% 
      filter(1:n() == 1) %>%
      group_by() %>%
      arrange(.dots = focal_covar)
  # df_topred <- d_case %>%
  #   group_by(age) %>%
  #   filter(1:n() == 1) %>%
  #   group_by()

  ## predict
  all_pred = lapply(models_aggregate, predict_finddiff, 
                    df_topred = df_topred,
                    treatment_varname = treatment_varname,
                    aggregate = FALSE)
  names(all_pred) = names(models_aggregate)
  
  ## bind model-based estimates with each other
  pred_df = do.call(rbind.data.frame, all_pred)  
  pred_df$approach = gsub("\\..*", "", rownames(pred_df))
  pred_df_rbind = pred_df[, c(focal_covar, "Estimate", "approach")]
  
  ## bind with non-model based estimates
  all_est = rbind.data.frame(unadjusted, nonparametric_adjusted,
                             pred_df_rbind)
  if(verbose){
    print("Completed one iteration of estimation")
  } 
  
  return(all_est)
}

estimate_subgroup_point <- gen_pred_subgroup(full_data = support_data,
                                    filter_flags = c("has_support"),
                                    weight_varname = "asecwt",
                                    treatment_varname = "derived_mother",
                                    outcome_varname = "derived_ln_wage",
                                    control_vec = c(strata_vars, "derived_agesq"),
                                    models_aggregate = estimate_aggregate_point_all$models)

## get estimated models from previous run of replicate weights
mods_aggregate_repweights = lapply(replicates_aggregate, function(x) x$models)

if(RUN_VAR_ESTIMATE_SUBGROUP){
  replicates_subgroup <- foreach(i = 1:160, .combine = "rbind") %do% {
    gen_pred_subgroup(full_data = support_data,
                      filter_flags = c("has_support"),
                      weight_varname = paste0("repwtp",i),
                      treatment_varname = "derived_mother",
                      outcome_varname = "derived_ln_wage",
                      control_vec = c(strata_vars, "derived_agesq"),
                      models_aggregate = mods_aggregate_repweights[[i]],
                      verbose = TRUE)
  }
  saveRDS(replicates_subgroup, "int_objects/replicates_subgroup.RDS", 
          compress = FALSE)
} else{
  replicates_subgroup <- readRDS("int_objects/replicates_subgroup.RDS")
}



focal_covar = "age"
estimate_subgroup_variance <- replicates_subgroup %>%
  rename(Estimate_star = Estimate) %>%
  left_join(estimate_subgroup_point,
            by = c(focal_covar,"approach")) %>%
  group_by(.dots = c(focal_covar, "approach")) %>%
  summarize(variance = 4 / 160 * sum((Estimate_star - Estimate) ^ 2))

subgroup_results <- estimate_subgroup_point %>%
  filter(approach != "Unadjusted") %>%
  left_join(estimate_subgroup_variance, by = c(focal_covar,"approach")) %>%
  mutate(approach = factor(case_when(grepl("Stratification", approach, ignore.case = TRUE) ~ 1,
                                     grepl("ols.*factor", approach, ignore.case = TRUE) ~ 2,
                                     grepl("gam.*interactive", approach, ignore.case = TRUE) ~ 3,
                                     grepl("ols.*interactive", approach, ignore.case = TRUE) ~ 4,
                                     grepl("ols.*additive", approach, ignore.case = TRUE) ~ 5),
                           labels = c("Stratification:\nNo estimation assumptions",
                                      "Education + Marital + Race + \n(Age indicators x Motherhood)",
                                      "Education + Marital + Race +\n(Age smooth x Motherhood)",
                                      "Education + Marital + Race +\n(Age quadratic x Motherhood)",
                                      "Education + Marital + Race +\nAge quadratic + Motherhood")))

##########################
# Plot the gap estimates #
##########################

for_aggregate_plot <- aggregate_results %>%
  mutate(estimand = "Aggregate\ngap",
         age = 35) %>%
  bind_rows(subgroup_results %>%
              mutate(estimand = "Age-specific\ngap")) %>%
  mutate(estimand = fct_relevel(estimand,
                                "Aggregate\ngap",
                                "Age-specific\ngap"))
for_aggregate_plot %>%
  ggplot(aes(x = age, y = Estimate,
             ymin = Estimate - qnorm(.975) * sqrt(variance),
             ymax = Estimate + qnorm(.975) * sqrt(variance))) +
  geom_hline(yintercept = 0, color = "red") +
  geom_errorbar(width = .5) +
  geom_point() +
  theme_bw() +
  #scale_y_continuous(name = "Difference in log hourly wage\n(Mother - Non-mothers)",
  #                   limits = c(-.27,.27)) +
  scale_y_continuous(name = "Controlled direct effect of motherhood\non the log hourly wage that would\nbe realized if employed",
                     limits = c(-.27,.27)) +
  xlab("Age") +
  facet_grid(estimand ~ approach) +
  theme(strip.text.y = element_text(angle = 0),
        strip.text.x = element_text(size = 8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  # Annotate within facets
  geom_text(data = for_aggregate_plot %>%
              filter(estimand == "Aggregate\ngap") %>%
              mutate(age = ifelse(approach == "Stratification:\nNo estimation assumptions",25,45),
                     Estimate = .27,
                     label = case_when(
                       approach == "Stratification:\nNo estimation assumptions" ~ "Weakest assumptions\n(most credible)",
                       approach == "Education + Marital + Race +\nAge quadratic + Motherhood" ~ "Strongest assumptions\n(least credible)"
                     )),
            aes(label = label, hjust = ifelse(approach == "Stratification:\nNo estimation assumptions",0,1)),
            vjust = 1,
            size = 3) +
  geom_text(data = for_aggregate_plot %>%
              filter(estimand == "Aggregate\ngap" &
                       approach == "Race + Marital +\n(Age smooth x Motherhood)") %>%
              mutate(Estimate = .27),
            label = "All estimation strategies\nestimate a similar\naggregate gap",
            size = 3, vjust = 1) +
  geom_text(data = for_aggregate_plot %>%
              filter(estimand == "Aggregate\ngap" & approach == "Race + Marital +\nAge quadratic + Motherhood") %>%
              mutate(age = 45,
                     Estimate = -.27),
            label = "This panel parallels the estimand and\n specification in Pal and Waldfogel (2016)",
            vjust = 0, hjust = 1,
            size = 2) +
  geom_segment(data = for_aggregate_plot %>%
                 filter(estimand == "Aggregate\ngap") %>%
                 mutate(Estimate = .205,
                        age = case_when(approach == "Race + Marital +\n(Age smooth x Motherhood)" ~ 27),
                        xend = 25),
               aes(xend = xend, yend = Estimate),
               arrow = arrow(length = unit(.15,"cm"))) +
  geom_segment(data = for_aggregate_plot %>%
                 filter(estimand == "Aggregate\ngap") %>%
                 mutate(Estimate = .205,
                        age = case_when(approach == "Race + Marital +\n(Age smooth x Motherhood)" ~ 43),
                        xend = 45),
               aes(xend = xend, yend = Estimate),
               arrow = arrow(length = unit(.15,"cm"))) +
  geom_text(data = for_aggregate_plot %>%
              filter(estimand == "Age-specific\ngap" & age == 25) %>%
              mutate(Estimate = .27,
                     label = case_when(
                       approach == "Stratification:\nNo estimation assumptions" ~ "Very uncertain\ngaps within\nsubgroups",
                       approach == "Education + Marital + Race +\n(Age smooth x Motherhood)" ~ "With this assumption,\nevidence shows a gap\nat young ages only",
                       approach == "Education + Marital + Race +\nAge quadratic + Motherhood" ~ "Additive OLS assumes\na constant effect.\nClearly this model\nis only an approximation."
                     )),
            aes(label = label),
            size = 3, hjust = 0, vjust = 1) +
  ggsave("output/all_gap_estimates.pdf",
         height = 5, width = 10)

#######################
# Cross-validated MSE #
#######################
make_cv_results  <- function(full_data, 
                             filter_flags,
                             weight_varname,
                             treatment_varname, outcome_varname,
                             control_vec,
                             nfolds = 5) {
  
  d_case <- full_data[rowSums(full_data[ , filter_flags]) == length(filter_flags), ]
  
  d_case$weight <- d_case[[weight_varname]]
  
  
  # Separate vector of controls into different combinations (specific to application)
  control_vec_group = setdiff(control_vec, "derived_agesq")
  control_vec_interact = c("age", "derived_agesq")
  control_vec_nointeract = setdiff(control_vec, control_vec_interact)
  
  ## create folded data
  folded <- d_case %>%
    group_by() %>%
    arrange(weight) %>%
    mutate(fold = rep(1:nfolds, ceiling(n() / nfolds))[1:n()])
  
  ## iterate through folds as train and test set
  foreach(f = 1:nfolds, .combine = "rbind") %do% {
    
    ## split into train and test
    train <- folded %>%
      filter(fold != f)
    test <- folded %>%
      filter(fold == f)
    
    ## estimate the models
    ols_additive_fit <- lm(formula(sprintf("%s ~ %s + %s",
                                           outcome_varname,
                                           treatment_varname,
                                           paste(control_vec,
                                                 collapse = "+"))),
                           data = train,
                           weights = weight)
    
    ols_interactive_fit <- lm(formula(sprintf("%s ~ %s + %s",
                                              outcome_varname,
                                              paste(sprintf("%s*%s", 
                                                            treatment_varname, control_vec_interact),
                                                    collapse = "+"),
                                              paste(control_vec_nointeract,
                                                    collapse = "+"))),
                              data = train,
                              weights = weight)
    
    ols_ageFactor_fit <- lm(formula(sprintf("%s ~ %s + %s",
                                            outcome_varname,
                                            
                                            ## interact tx with factor
                                            sprintf("%s*factor(%s)", 
                                                    treatment_varname, 
                                                    setdiff(control_vec_interact, "derived_agesq")),
                                            
                                            ## other controls
                                            paste(control_vec_nointeract,
                                                  collapse = "+"))),
                            data = train,
                            weights = weight)
    tmp = train
    tmp[[treatment_var]] = factor(tmp[[treatment_var]])
    gam_interactive_fit <- gam(formula(sprintf("%s ~ %s + s(%s, by = %s) + %s",
                                               outcome_varname,
                                               treatment_varname,
                                               setdiff(control_vec_interact, "derived_agesq"),
                                               treatment_varname,
                                               paste(control_vec_nointeract,
                                                     collapse = "+"))),
                               data = tmp %>% group_by(),
                               weight = tmp$weight)
    
    stratification_result <- train %>%
      group_by(.dots = c(treatment_varname, control_vec_group)) %>%
      summarize(yhat = weighted.mean(!!sym(outcome_varname), w = weight, na.rm = T),
                .groups = "drop") %>%
      inner_join(test, by = c(treatment_varname, control_vec_group)) %>%
      mutate(squared_error = (!!sym(outcome_varname) - yhat) ^ 2) %>%
      group_by() %>%
      mutate(approach = "Stratification") %>%
      select(approach, squared_error, weight)
    
    
    model_squared_errors <- test %>%
      select(all_of(c(outcome_varname, "weight"))) %>%
      mutate(ols_additive = predict(ols_additive_fit, newdata = test),
             ols_interactive = predict(ols_interactive_fit, newdata = test),
             ols_ageFactor = predict(ols_ageFactor_fit, newdata = test),
             gam_interactive = predict(gam_interactive_fit, newdata = test)) %>%
      melt(id = c(outcome_varname, "weight"), variable.name = "approach", value.name = "yhat") %>%
      mutate(squared_error = (!!sym(outcome_varname) - yhat) ^ 2)
    
    stratification_result %>%
      bind_rows(model_squared_errors)
  } %>%
    group_by(approach) %>%
    summarize(mse = weighted.mean(squared_error, w = weight, na.rm = T),
              .groups = "drop")
}

cv_point <- make_cv_results(full_data = support_data,
            filter_flags = c("has_support"),
            weight_varname = "asecwt",
            treatment_varname = "derived_mother",
            outcome_varname = "derived_ln_wage",
            control_vec = c(strata_vars, "derived_agesq"))
            
            
replicates <- foreach(i = 1:160, .combine = "rbind") %do% {
  make_cv_results(weight_varname = paste0("repwtp",i),
                  full_data = support_data,
                  filter_flags = c("has_support"),
                  treatment_varname = "derived_mother",
                  outcome_varname = "derived_ln_wage",
                  control_vec = c(strata_vars, "derived_agesq"))
}

cv_variance <- replicates %>%
  rename(mse_star = mse) %>%
  left_join(cv_point,
            by = c("approach")) %>%
  group_by(approach) %>%
  summarize(variance = 4 / 160 * sum((mse_star - mse) ^ 2))

cv_results <- cv_point %>%
  left_join(cv_variance, by = "approach") %>%
  mutate(approach = factor(case_when(approach == "Stratification" ~ 1,
                                     approach == "ols_ageFactor" ~ 2,
                                     approach == "gam_interactive" ~ 3,
                                     approach == "ols_interactive" ~ 4,
                                     approach == "ols_additive" ~ 5),
                           labels = c("Stratification:\nNo estimation assumptions",
                                      "Education + Marital + Race +\n(Age indicators x Motherhood)",
                                      "Education + Marital + Race +\n(Age smooth x Motherhood)",
                                      "Education + Marital + Race +\n(Age quadratic x Motherhood)",
                                      "Education + Marital + Race +\nAge quadratic + Motherhood")))
print(cv_results %>%
        mutate(approach = fct_rev(approach),
               best = mse == min(mse)) %>%
        as.data.frame())

sink()
