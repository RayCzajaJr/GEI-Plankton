

### get prediction skill and other summary stat values for each model used for projections


library(mgcv)
library(ggplot2)
library(gridExtra)
library(visreg)
library(DHARMa)
library(png)
library(MuMIn)

# Function to get prediction skill metrics (LOO)
calculate_metrics <- function(model, testing_data, response_var) {
  
  predicted_values <- predict.gam(model, type = "response", newdata = testing_data)
  observed_values <- testing_data[[response_var]]
  
  linear_model <- gam(observed_values ~ predicted_values, family = tw(), na.action = "na.fail")
  
  r_squared <- ((linear_model$null.deviance - linear_model$deviance) / linear_model$null.deviance) 
  correlation_coefficient <- cor(observed_values, predicted_values, use = "complete.obs", method = "spearman")
  
  return(list(
    r_squared = r_squared,
    correlation_coefficient = correlation_coefficient,
    predicted_values = predicted_values,
    observed_values = observed_values,
    linear_model = linear_model
  ))
}

# Main function with AIC and WIC from full model
calculate_prediction_skill_metrics <- function(data, response_var, predictors) {
  
  unique_years <- unique(data$YEAR)
  all_r_squared <- list()
  all_correlation_coefficients <- list()
  plots <- list()
  
  # Build the full model once on all data
  predictor_terms <- sapply(names(predictors), function(pred) {
    k_value <- predictors[[pred]]$k
    m_value <- predictors[[pred]]$m
    bs_value <- predictors[[pred]]$bs
    
    spline_call <- paste0("s(", pred)
    
    # Add optional parameters
    if (!is.null(k_value)) spline_call <- paste0(spline_call, ", k = ", k_value)
    if (!is.null(m_value)) spline_call <- paste0(spline_call, ", m = ", m_value)
    if (!is.null(bs_value)) spline_call <- paste0(spline_call, ", bs = \"", bs_value, "\"")
    
    spline_call <- paste0(spline_call, ")")
    return(spline_call)
  })
  
  formula <- as.formula(paste(response_var, "~", paste(predictor_terms, collapse = " + ")))
  
  full_model <- gam(formula, data = data, method = "REML", select = TRUE, na.action = "na.fail")
  
  # AIC and Akaike weight from full model
  full_model_aic <- AIC(full_model)
  print(paste("AIC of full model:", full_model_aic))
  
  # Akaike weight assuming it's one of multiple candidate models
  akaike_weight <- 1  # If only one model is being compared
  print(paste("Akaike weight of full model (assuming it's only model):", akaike_weight))
  
  # LOO evaluation
  for (year in unique_years) {
    
    training_data <- data[data$YEAR != year, ]
    testing_data <- data[data$YEAR == year, ]
    
    if (all(testing_data[[response_var]] == 0)) {
      print(paste("All values for year:", year, "are zero. Skipping this year."))
      next
    }
    
    model <- gam(formula, data = training_data, method = "REML", select = TRUE, na.action = "na.fail")
    metrics <- calculate_metrics(model, testing_data, response_var)
    
    if (!is.na(metrics$r_squared)) {
      print(paste("R-squared for year:", year, "is", metrics$r_squared))
      print(paste("Correlation coefficient for year:", year, "is", metrics$correlation_coefficient))
      
      all_r_squared[[as.character(year)]] <- metrics$r_squared
      all_correlation_coefficients[[as.character(year)]] <- metrics$correlation_coefficient
      
      p <- ggplot(data = data.frame(predicted = metrics$predicted_values, observed = metrics$observed_values)) +
        geom_point(aes(x = predicted, y = observed), color = "black") +
        geom_smooth(aes(x = predicted, y = observed), method = "lm", color = "black") +
        labs(x = "Predicted Values", y = "Observed Values") +
        theme_minimal() +
        theme(
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 6),
          plot.title = element_blank(),
          legend.position = "none"
        )
      
      plots[[as.character(year)]] <- p
    } else {
      print(paste("R-squared for year:", year, "is NA. Skipping this year."))
    }
  }
  
  filtered_r_squared <- unlist(all_r_squared)[!is.na(unlist(all_r_squared))]
  filtered_correlation_coefficients <- unlist(all_correlation_coefficients)[!is.na(unlist(all_correlation_coefficients))]
  
  mean_r_squared <- mean(filtered_r_squared)
  mean_correlation_coefficient <- mean(filtered_correlation_coefficients)
  
  print(paste("Mean R-squared value across non-NA years is", mean_r_squared))
  print(paste("Mean correlation coefficient across non-NA years is", mean_correlation_coefficient))
  
  grid.arrange(grobs = plots, ncol = 3)
  
  return(list(
    mean_r_squared = mean_r_squared,
    mean_correlation_coefficient = mean_correlation_coefficient,
    full_model_aic = full_model_aic,
    akaike_weight = akaike_weight,
    full_model = full_model,
    plots = plots
  ))
}

# Example usage
response_var <- "THUNNUSTHYNNUS"
predictors <- list(
  "TEMPMAX" = list(k = 3, m = NULL),
  "TEMPSURF" = list(k = 8, m = NULL),
  "CHLORSURF" = list(k = 3, m = NULL),
  "lat,lon" = list(k = NULL, m = NULL)
)
ABFTClimateResult <- calculate_prediction_skill_metrics(springplankton4, response_var, predictors)
THYNGam<- gam(THUNNUSTHYNNUS~  s(CHLORSURF, k=3)+s(TEMPSURF, k=8) +s(lon, lat)+s(TEMPMAX, k=3),
              data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
plot(THYNGam)
summary(THYNGam)
dredge(THYNGam, rank = "AIC")

### ABFT Prey 
response_var <- "THUNNUSTHYNNUS"
predictors <- list(
  "Cyclopoids" = list(k = 3, m = NULL),
  "Cladocerans" = list(k = 4, m = 1),
  "lat,lon" = list(k = NULL, m = NULL)
)

ABFTPreyResult <- calculate_prediction_skill_metrics(springplankton4, response_var, predictors)
THYNGamPrey<- gam(THUNNUSTHYNNUS~  s(Larvaceans, k=3)+s(Cladocerans, k=4, m=1) +s(lon, lat)+s(Cyclopoids, k=3),
              data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
plot(THYNGamPrey)
summary(THYNGamPrey)
dredge(THYNGamPrey, rank = "AIC")


### thunnus Climate 
response_var <- "THUNNUS"
predictors <- list(
  "TEMPSURF" = list(k = 3, m = 1),
  "TEMPMAX" = list(k = 3, m = 1),
  "SALSURF" = list(k = 3, m = 1),
  "CHLORSURF" = list(k = 3, m = NULL),
  "lat,lon" = list(k = NULL, m = NULL)
)

ThunnusClimateResult <- calculate_prediction_skill_metrics(springplankton4, response_var, predictors)
THUNGam<- gam(THUNNUS~  s(CHLORSURF, k=3)+s(TEMPSURF, k=3, m=1) +s(SALSURF, k=3, m=1) +s(lon, lat)+s(TEMPMAX, k=3),
              data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(THUNGam)
dredge(THUNGam, rank = "AIC")

### thunnus Prey 
# remove calanoids bc shrunk to zero
response_var <- "THUNNUS"
predictors <- list(
  "Larvaceans" = list(k = 3, m = 1),
  "Cyclopoids" = list(k = 3, m = NULL),
  "lat, lon" = list(k = NULL, m = NULL)
)

ThunnusPreyResult <- calculate_prediction_skill_metrics(springplankton4, response_var, predictors)
THUNGamPrey<- gam(THUNNUS~  s(Larvaceans, k=3,m=1)+s(Cyclopoids, k=3),
                  data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(THUNGamPrey)
dredge(THUNGamPrey, rank = "AIC")


### Billfish Climate 
# sal and chlor shrunk to zero so remove
response_var <- "ISTIOPHORIDAE"
predictors <- list(
  "TEMPSURF" = list(k = 3, m = 1))

BillfishClimateResult <- calculate_prediction_skill_metrics(springplankton5, response_var, predictors)
ISTGam<- gam(ISTIOPHORIDAE~s(TEMPSURF, k=3, m=1) + s(TEMPMAX, k=3) + s(lon, lat),
             data =springplankton5, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(ISTGam)
dredge(ISTGam, rank = "AIC")


### Billfish Prey 
# both copepods shrunk to zero so remove
response_var <- "ISTIOPHORIDAE"
predictors <- list(
  "Cladocerans" = list(k = 3, m = NULL),
  "lat,lon" = list(k = NULL, m = NULL)
)

BillfishPreyResult <- calculate_prediction_skill_metrics(springplankton5, response_var, predictors)
BILLGamPrey<- gam(ISTIOPHORIDAE~ s(lon, lat) +s(Cladocerans,k=3),
                 data =springplankton5, tw(),method = "REML", select = TRUE, na.action = "na.fail")
plot(BILLGamPrey)
summary(BILLGamPrey)
dredge(BILLGamPrey, rank = "AIC")


### Skipjack Climate 
# temp max shrank to zero
response_var <- "KATSUWONUSPELAMIS"
predictors <- list(
  "CHLORSURF" = list(k = 3, m = NULL),
  "SALSURF" = list(k = 3, m = 1),
  "TEMPSURF" = list(k = 3, m = 1))

SkipjackClimateResult <- calculate_prediction_skill_metrics(springplankton4, response_var, predictors)
KATGam<- gam(KATSUWONUSPELAMIS~s(TEMPSURF, k=3, m=1) + s(lon, lat)+s(SALSURF, k=3,m=1)+s(CHLORSURF, k=3),
             data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(KATGam)
dredge(KATGam, rank = "AIC")

### Skipjack Prey
# cyclopoids and cladocerans shrank to zero
response_var <- "KATSUWONUSPELAMIS"
predictors <- list(
  "Larvaceans" = list(k = NULL, m = NULL),
  "lat,lon" = list(k = NULL, m = NULL)
)

SkipjackPreyResult <- calculate_prediction_skill_metrics(springplankton4, response_var, predictors)
KATGamprey<- gam(KATSUWONUSPELAMIS~ s(lon, lat)+s(Larvaceans),
             data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(KATGamprey)
dredge(KATGamprey, rank = "AIC")


### Little Tunny Climate 
# remove sal and chlor bc shrunk to zero
response_var <- "EUTHYNNUSALLETTERATUS"
predictors <- list(
  "TEMPSURF" = list(k = 3, m = 1),
  "lat,lon" = list(k = NULL, m = NULL)
)

TunnyClimateResult <- calculate_prediction_skill_metrics(springplankton4, response_var, predictors)
EUTGam<- gam(EUTHYNNUSALLETTERATUS~ s(TEMPSURF, k=3,m=1) +s(lon, lat)+s(TEMPMAX, k=3,m=1),
             data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(EUTGam)
dredge(EUTGam, rank = "AIC")

### Little Tunny Prey
response_var <- "EUTHYNNUSALLETTERATUS"
predictors <- list(
  "Larvaceans" = list(k = 3, m = NULL),
  "lat,lon" = list(k = NULL, m = NULL)
)

TunnyPreyResult <- calculate_prediction_skill_metrics(springplankton4, response_var, predictors)
EUTGamPrey<- gam(EUTHYNNUSALLETTERATUS~  s(Calanoids, k=4, m=1) + s(Cyclopoids, k=3, m=1) + s(Larvaceans, k=3) + s(lon, lat),
                  data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(EUTGamPrey)
dredge(EUTGamPrey, rank = "AIC")


### Frigate Climate 
response_var <- "AUXIS"
predictors <- list(
  "TEMPSURF" = list(k = 3, m = 1),
  "TEMPMAX" = list(k = 3, m = NULL),
  "lat,lon" = list(k = NULL, m = NULL)
)

FrigateClimateResult <- calculate_prediction_skill_metrics(springplankton4, response_var, predictors)
AUXISGam<- gam(AUXIS~ s(TEMPSURF, k=3, m=1) +s(lon, lat)+s(TEMPMAX, k=3),
             data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(AUXISGam)
dredge(AUXISGam, rank = "AIC")


### Frigate Prey
response_var <- "AUXIS"
predictors <- list(
  "Larvaceans" = list(k = 3, m = 1),
  "Calanoids" = list(k = 3, m = 1),
  "lat,lon" = list(k = NULL, m = NULL)
)

FrigatePreyResult <- calculate_prediction_skill_metrics(springplankton4, response_var, predictors)
AUXISGamPrey<- gam(AUXIS~  s(Larvaceans, k=3)+s(Calanoids, k=3,m=1) +s(lon, lat),
                 data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(AUXISGamPrey)
dredge(AUXISGamPrey, rank = "AIC")


### Mahi mahi Climate 
response_var <- "CORYPHAENA"
predictors <- list(
  "TEMPSURF" = list(k = 3, m = 1),
  "SALSURF" = list(k = 3, m = 1),
  "CHLORSURF" = list(k = 3, m=NULL),
  "TEMPMAX" = list(k = 3, m=NULL)
)
MahiClimateResult <- calculate_prediction_skill_metrics(springplankton4, response_var, predictors)

CORYPHAENAGam<- gam(CORYPHAENA~  s(TEMPSURF, k=3,m=1) +s(SALSURF, k=3,m=1),
               data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(CORYPHAENAGam)
dredge(CORYPHAENAGam, rank = "AIC")


### Mahi mahi Prey
# calanoids shrinks to zero so remove
response_var <- "CORYPHAENA"
predictors <- list(
  "Larvaceans" = list(k = 3, m = 1),
  "Cladocerans" = list(k = 3, m = NULL),
  "Cyclopoids" = list(k = 3, m = NULL),
  "lat,lon" = list(k = NULL, m = NULL)
)
MahiPreyResult <- calculate_prediction_skill_metrics(springplankton4, response_var, predictors)
CORYPHAENAGamPrey<- gam(CORYPHAENA~ s(Larvaceans, k=3, m=1) + s(lon, lat),
                   data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(CORYPHAENAGamPrey)
dredge(CORYPHAENAGamPrey, rank = "AIC")


### Calanoid Climate 
response_var <- "Calanoids"
predictors <- list(
  "SALSURF" = list(k = 3, m = NULL),
  "TEMPMAX" = list(k = 3, m = 1),
  "CHLORSURF" = list(k = 3, m = 1),
  "TEMPSURF" = list(k = 3, m = 1),
  "TIME_EMIL" = list(k = NULL, m = NULL, bs = "cs"),
  "lat,lon" = list(k = NULL, m = NULL)
)

CalanoidClimateResult <- calculate_prediction_skill_metrics(springplankton4, response_var, predictors)
CalGam<- gam(Calanoids~  s(TEMPSURF, k=3,m=1) +s(lon, lat)+s(TEMPMAX, k=3,m=1)+s(TIME_MIL, bs="cs")+s(SALSURF,k=3,m=1)+s(CHLORSURF, k=3,m=1),
              data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(CalGam)
dredge(CalGam, rank = "AIC")


### Cladocerans Climate 
# remove chlorophyll because it shrunk to zero
response_var <- "Cladocerans"
predictors <- list(
  "SALSURF" = list(k = 3, m = NULL),
  "TEMPSURF" = list(k = 3, m = NULL),
  "TIME_EMIL" = list(k = NULL, m = NULL, bs = "cs"),
  "lat,lon" = list(k = NULL, m = NULL)
)

CladoceranClimateResult <- calculate_prediction_skill_metrics(springplankton4, response_var, predictors)
CladGam<- gam(Cladocerans~  s(TEMPSURF, k=3) +s(lon, lat)+s(SALSURF,k=3),
              data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(CladGam)
dredge(CladGam, rank = "AIC")


#### Larvaceans Climate
# remove chlorophyll because it shrunk to zero
response_var <- "Larvaceans"
predictors <- list(
  "TEMPMAX" = list(k = 3, m = 1),
  "TEMPSURF" = list(k = 3, m = 1),
  "lat,lon" = list(k = NULL, m = NULL)
)

LarvClimateResult <- calculate_prediction_skill_metrics(springplankton4, response_var, predictors)
LarvGam<- gam(Larvaceans~  s(TEMPSURF, k=3,m=1) +s(lon, lat)+s(TEMPMAX, k=3,m=1),
             data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(LarvGam)
dredge(LarvGam, rank = "AIC")


### Cyc Climate 
response_var <- "Cyclopoids"
predictors <- list(
  "TEMPSURF" = list(k = 3, m = 1),
  "TEMPMAX" = list(k = 3, m = 1),
  "CHLORSURF" = list(k = 3, m = NULL),
  "lat,lon" = list(k = NULL, m = NULL)
)
CyclopoidClimateResult <- calculate_prediction_skill_metrics(springplankton4, response_var, predictors)
CYCGam<- gam(Cyclopoids~  s(CHLORSURF, k=3)+s(TEMPSURF, k=3,m=1) +s(SALSURF, k=3) +s(lon, lat)+s(TEMPMAX, k=3,m=1),
                    data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(CYCGam)
dredge(CYCGam, rank = "AIC")





