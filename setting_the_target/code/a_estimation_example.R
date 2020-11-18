
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
library(here)
library(rvest)

## Constants
OUTPUT_NAME= "PalWaldfogel_output.txt"

## Can either read raw data/clean or read in cleaned data
READ_RAW_DATA = TRUE
RAW_DATA_NAME = "cps_00040.dta"
CLEAN_DATA_NAME = "pw_analytic.csv"

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
           
           ## Filters based on those covariates
           filter_nonmisseduc = !is.na(derived_educ),
           filter_nonmisswage = !is.na(derived_ln_wage)
           ) 
  
  ## Variables to include
  weights_identifiers = c("pernum", "serial", "asecwt",
                          grep("repwtp", colnames(d_all),
                               value = TRUE)) 
  outcome = "derived_ln_wage"
  treatment = "derived_mother"
  covariates = c(sprintf("derived_%s", 
                         c("educ", "married",
                           "race")), "age")
  filters = grep("filter", colnames(d_all), value = TRUE)
  cols_keep = c(weights_identifiers, outcome, treatment, covariates, 
                filters)
  
  d_all <- d_clean %>%
    select(all_of(cols_keep))
  
  write.csv(d_all, sprintf("data/clean/%s", CLEAN_DATA_NAME))
} else{
  d_all <- read.csv(sprintf("data/clean/%s", CLEAN_DATA_NAME))
}

filter_vars = grep("filter", colnames(d), value = TRUE)
print("Sample sizes")
colSums(d[, filter_vars])

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
                   arrange(c(treatment_var, strata_vars))),
      # The max is a very high number to tell it to print all rows
      max = 9999)

# All remaining analyses will focus on those with support only

#################
# Aggregate gap #
#################

make_aggregate_results <- function(weight_name) {
  d_case <- support_data %>% filter(has_support)
  d_case$weight <- d_case[[weight_name]]
  # Normalize weights
  d_case$weight <- d_case$weight / sum(d_case$weight)
  unadjusted <- d_case %>%
    group_by(mother) %>%
    summarize(Estimate = weighted.mean(ln_wage, w = weight, na.rm = T)) %>%
    spread(key = mother, value = Estimate) %>%
    mutate(approach = "Unadjusted",
           Estimate = `TRUE` - `FALSE`) %>%
    select(approach, Estimate)
  
  nonparametric_adjusted <- d_case %>%
    group_by(mother, age, educ, race, married) %>%
    summarize(estimate = weighted.mean(ln_wage, w = weight, na.rm = T),
              weight = sum(weight)) %>%
    group_by(age, educ, race, married) %>%
    # Give this stratum the weight of mothers in this stratum
    mutate(weight = sum(weight * mother)) %>%
    group_by(mother) %>%
    summarize(Estimate = weighted.mean(estimate, w = weight)) %>%
    spread(key = mother, value = Estimate) %>%
    mutate(approach = "Stratification",
           Estimate = `TRUE` - `FALSE`) %>%
    select(approach, Estimate)
  
  ols_additive_fit <- lm(ln_wage ~ mother + age + I(age ^ 2) +
                           educ + race + married,
                         data = d_case,
                         weights = weight)
  ols_interactive_fit <- lm(ln_wage ~ mother*(age + I(age ^ 2)) +
                              educ + race + married,
                            data = d_case,
                            weights = weight)
  ols_ageFactor_fit <- lm(ln_wage ~ mother*(factor(age)) +
                            educ + race + married,
                          data = d_case,
                          weights = weight)
  gam_interactive_fit <- gam(ln_wage ~ mother + s(age, by = mother) +
                               educ + race + married,
                             data = d_case %>%
                               group_by() %>%
                               mutate(mother = factor(mother)),
                             weight = d_case$weight)
  
  to_predict <- d_case %>%
    group_by() %>%
    filter(mother)
  
  ols_additive_adjusted <- data.frame(
    approach = "OLS (additive age)",
    Estimate = weighted.mean(predict(ols_additive_fit, newdata = to_predict %>% mutate(mother = T)) -
                               predict(ols_additive_fit, newdata = to_predict %>% mutate(mother = F)),
                             w = to_predict$weight)
  )
  
  ols_interactive_adjusted <- data.frame(
    approach = "OLS (interactive age)",
    Estimate = weighted.mean(predict(ols_interactive_fit, newdata = to_predict %>% mutate(mother = T)) -
                               predict(ols_interactive_fit, newdata = to_predict %>% mutate(mother = F)),
                             w = to_predict$weight)
  )

  ols_ageFactor_adjusted <- data.frame(
    approach = "OLS (age factor)",
    Estimate = weighted.mean(predict(ols_ageFactor_fit, newdata = to_predict %>% mutate(mother = T)) -
                               predict(ols_ageFactor_fit, newdata = to_predict %>% mutate(mother = F)),
                             w = to_predict$weight)
  )
  
  gam_interactive_adjusted <- data.frame(
    approach = "GAM (interactive age)",
    Estimate = weighted.mean(predict(gam_interactive_fit, newdata = to_predict %>% mutate(mother = T)) -
                               predict(gam_interactive_fit, newdata = to_predict %>% mutate(mother = F)),
                             w = to_predict$weight)
  )
  
  return(ols_additive_adjusted %>%
           bind_rows(ols_interactive_adjusted) %>%
           bind_rows(ols_ageFactor_adjusted) %>%
           bind_rows(gam_interactive_adjusted) %>%
           bind_rows(nonparametric_adjusted))
}

estimate_aggregate_point <- make_aggregate_results("asecwt")
replicates_aggregate <- foreach(i = 1:160, .combine = "rbind") %do% {
  make_aggregate_results(paste0("repwtp",i))
}

estimate_aggregate_variance <- replicates_aggregate %>%
  rename(Estimate_star = Estimate) %>%
  left_join(estimate_aggregate_point,
            by = "approach") %>%
  group_by(approach) %>%
  summarize(variance = 4 / 160 * sum((Estimate_star - Estimate) ^ 2))

aggregate_results <- estimate_aggregate_point %>%
  left_join(estimate_aggregate_variance, by = "approach") %>%
  mutate(approach = factor(case_when(approach == "Stratification" ~ 1,
                                     approach == "OLS (age factor)" ~ 2,
                                     approach == "GAM (interactive age)" ~ 3,
                                     approach == "OLS (interactive age)" ~ 4,
                                     approach == "OLS (additive age)" ~ 5),
                           labels = c("Stratification:\nNo estimation assumptions",
                                      "Education + Marital + Race + \n(Age indicators x Motherhood)",
                                      "Education + Marital + Race +\n(Age smooth x Motherhood)",
                                      "Education + Marital + Race +\n(Age quadratic x Motherhood)",
                                      "Education + Marital + Race +\nAge quadratic + Motherhood")))

##################################
# Estimates within age subgroups #
##################################

make_subgroup_results <- function(weight_name) {
  d_case <- support_data %>% filter(has_support)
  d_case$weight <- d_case[[weight_name]]
  # Normalize weights
  d_case$weight <- d_case$weight / sum(d_case$weight)
  unadjusted <- d_case %>%
    group_by(age,mother) %>%
    summarize(Estimate = weighted.mean(ln_wage, w = weight, na.rm = T)) %>%
    spread(key = mother, value = Estimate) %>%
    mutate(approach = "Unadjusted",
           Estimate = `TRUE` - `FALSE`) %>%
    select(age, approach, Estimate)
  
  nonparametric_adjusted <- d_case %>%
    group_by(mother, age, educ, race, married) %>%
    summarize(estimate = weighted.mean(ln_wage, w = weight, na.rm = T),
              weight = sum(weight)) %>%
    group_by(age, educ, race, married) %>%
    # Give this stratum the weight of mothers in this stratum
    mutate(weight = sum(weight * mother)) %>%
    group_by(age, mother) %>%
    summarize(Estimate = weighted.mean(estimate, w = weight)) %>%
    spread(key = mother, value = Estimate) %>%
    mutate(approach = "Stratification",
           Estimate = `TRUE` - `FALSE`) %>%
    select(age, approach, Estimate)
  
  ols_additive_fit <- lm(ln_wage ~ mother + age + I(age ^ 2) +
                           educ + race + married,
                         data = d_case,
                         weights = weight)
  ols_interactive_fit <- lm(ln_wage ~ mother*(age + I(age ^ 2)) +
                              educ + race + married,
                            data = d_case,
                            weights = weight)
  ols_ageFactor_fit <- lm(ln_wage ~ mother*(factor(age)) +
                               educ + race + married,
                             data = d_case,
                             weights = weight)
  gam_interactive_fit <- gam(ln_wage ~ mother + s(age, by = mother) +
                               educ + race + married,
                             data = d_case %>%
                               group_by() %>%
                               mutate(mother = factor(mother)),
                             weight = d_case$weight)
  
  to_predict <- d_case %>%
    group_by(age) %>%
    filter(1:n() == 1) %>%
    group_by()
  
  ols_additive_adjusted <- to_predict %>%
    mutate(Estimate = predict(ols_additive_fit, newdata = to_predict %>% mutate(mother = T)) -
             predict(ols_additive_fit, newdata = to_predict %>% mutate(mother = F)),
           approach = "OLS (additive age)") %>%
    select(age, approach, Estimate)
  ols_interactive_adjusted <- to_predict %>%
    mutate(Estimate = predict(ols_interactive_fit, newdata = to_predict %>% mutate(mother = T)) -
             predict(ols_interactive_fit, newdata = to_predict %>% mutate(mother = F)),
           approach = "OLS (interactive age)") %>%
    select(age, approach, Estimate)
  ols_ageFactor_adjusted <- to_predict %>%
    mutate(Estimate = predict(ols_ageFactor_fit, newdata = to_predict %>% mutate(mother = T)) -
             predict(ols_ageFactor_fit, newdata = to_predict %>% mutate(mother = F)),
           approach = "OLS (age factor)") %>%
    select(age, approach, Estimate)
  gam_interactive_adjusted <- to_predict %>%
    mutate(Estimate = predict(gam_interactive_fit, newdata = to_predict %>% mutate(mother = T)) -
             predict(gam_interactive_fit, newdata = to_predict %>% mutate(mother = F)),
           approach = "GAM (interactive age)") %>%
    select(age, approach, Estimate)
  
  return(data.frame(unadjusted) %>%
           bind_rows(nonparametric_adjusted) %>%
           bind_rows(ols_additive_adjusted) %>%
           bind_rows(ols_interactive_adjusted) %>%
           bind_rows(ols_ageFactor_adjusted) %>%
           bind_rows(gam_interactive_adjusted))
}

estimate_subgroup_point <- make_subgroup_results("asecwt")
replicates_age <- foreach(i = 1:160, .combine = "rbind") %do% {
  make_subgroup_results(paste0("repwtp",i))
}

estimate_subgroup_variance <- replicates_age %>%
  rename(Estimate_star = Estimate) %>%
  left_join(estimate_subgroup_point,
            by = c("age","approach")) %>%
  group_by(age,approach) %>%
  summarize(variance = 4 / 160 * sum((Estimate_star - Estimate) ^ 2))

subgroup_results <- estimate_subgroup_point %>%
  filter(approach != "Unadjusted") %>%
  left_join(estimate_subgroup_variance, by = c("age","approach")) %>%
  mutate(approach = factor(case_when(approach == "Stratification" ~ 1,
                                     approach == "OLS (age factor)" ~ 2,
                                     approach == "GAM (interactive age)" ~ 3,
                                     approach == "OLS (interactive age)" ~ 4,
                                     approach == "OLS (additive age)" ~ 5),
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

make_cv_results <- function(weight_name) {
  
  with_support <- support_data %>% filter(has_support)
  with_support$weight <- with_support[[weight_name]]
  
  folded <- with_support %>%
    group_by() %>%
    arrange(weight) %>%
    mutate(fold = rep(1:5, ceiling(n() / 5))[1:n()])
  
  foreach(f = 1:5, .combine = "rbind") %do% {
    train <- folded %>%
      filter(fold != f)
    test <- folded %>%
      filter(fold == f)
    ols_additive_fit <- lm(ln_wage ~ mother + age + I(age ^ 2) +
                             educ + race + married,
                           data = train,
                           weights = weight)
    ols_interactive_fit <- lm(ln_wage ~ mother*(age + I(age ^ 2)) +
                                educ + race + married,
                              data = train,
                              weights = weight)
    ols_ageFactor_fit <- lm(ln_wage ~ mother*(factor(age)) +
                              educ + race + married,
                            data = train,
                            weights = weight)
    gam_interactive_fit <- gam(ln_wage ~ mother + s(age, by = mother) +
                                 educ + race + married,
                               data = train %>%
                                 group_by() %>%
                                 mutate(mother = factor(mother)),
                               weight = train$weight)
    
    stratification_result <- train %>%
      group_by(mother, age, educ, race, married) %>%
      summarize(yhat = weighted.mean(ln_wage, w = weight, na.rm = T)) %>%
      inner_join(test, by = c("mother","age","educ","race","married")) %>%
      mutate(squared_error = (ln_wage - yhat) ^ 2) %>%
      group_by() %>%
      mutate(approach = "Stratification") %>%
      select(approach, squared_error, weight)
    
    model_squared_errors <- test %>%
      select(ln_wage, weight) %>%
      mutate(ols_additive = predict(ols_additive_fit, newdata = test),
             ols_interactive = predict(ols_interactive_fit, newdata = test),
             ols_ageFactor = predict(ols_ageFactor_fit, newdata = test),
             gam_interactive = predict(gam_interactive_fit, newdata = test)) %>%
      melt(id = c("ln_wage","weight"), variable.name = "approach", value.name = "yhat") %>%
      mutate(squared_error = (ln_wage - yhat) ^ 2)
    
    stratification_result %>%
      bind_rows(model_squared_errors)
  } %>%
    group_by(approach) %>%
    summarize(mse = weighted.mean(squared_error, w = weight, na.rm = T))
}

cv_point <- make_cv_results("asecwt")
replicates <- foreach(i = 1:160, .combine = "rbind") %do% {
  make_cv_results(paste0("repwtp",i))
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
               best = mse == min(mse)))

sink()
