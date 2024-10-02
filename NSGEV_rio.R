# Clear the workspace
rm(list=ls())

# Install required packages (commented out if already installed)
# install.packages(c("ROOPSD", "lubridate", "dplyr", "ggplot2", "trend", "readxl", "scales", "pals"))

# Load necessary libraries
library(ROOPSD)
library(lubridate)
library(dplyr)
library(ggplot2)
library(trend)
library(readxl)
library(scales)
library(pals)

######## Functions ########

# Function to calculate the upper bound of extreme temperature based on multiple variables 
# This is applicable if the non-stationary GEV distribution has a negative shape parameter.

# Args:
#   en3.4: (numeric) The El Niño 3.4 index, which represents sea surface temperature anomalies in the central Pacific. 
#          A higher value indicates stronger El Niño conditions.
#   gwi: (numeric) Global Warming Index, a measure of cumulative global warming over time.
#   int: (numeric) The intercept term in the linear regression between extreme temperatures and the covariates (e.g., EN3.4, GWI).
#   beta1: (numeric) Coefficient representing the effect of the en3.4 index on extreme temperatures.
#   beta2: (numeric) Coefficient representing the effect of the gwi on extreme temperatures.
#   scale: (numeric) The scale parameter of the GEV distribution, which adjusts the spread or variability of extreme temperatures.
#   shape: (numeric) The shape parameter of the GEV distribution, which defines the tail behavior. A negative value indicates a bounded distribution with an upper limit.

# Returns:
#   limit: (numeric) The calculated upper bound for the temperature extremes, according to a GEV distribution with a negative shape parameter. 
#          This value is determined based on the provided inputs, reflecting the combined effects of El Niño (en3.4), global warming (gwi), and other parameters.

upper_bound <- function(en3.4, gwi, int, beta1, beta2, scale, shape) {
  
  # Calculate the upper bound using a linear combination of EN3.4 and GWI values,
  # along with other parameters (intercept, scale, and shape).
  # The formula is: limit = intercept (int) + EN3.4 effect + GWI effect + scale/|shape|
  limit <- int + en3.4 * beta1 + gwi * beta2 + scale / abs(shape)
  
  # Return the calculated limit
  return(limit)
}


# Function to calculate the return time for extreme temperatures based on input data and GEV distribution parameters.
#
# Args:
#   data: (numeric) Observed maximum temperature for a specific year or event.
#   loc: (numeric) The location parameter of the GEV distribution, which represents the central tendency (mean) of the distribution.
#   scale: (numeric) The scale parameter of the GEV distribution, which controls the spread or variability of the extreme temperatures.
#   shape: (numeric) The shape parameter of the GEV distribution, which defines the tail behavior. A positive value indicates a heavy-tailed distribution, while a negative value suggests a bounded upper limit.
# 
# Returns:
#   Rt_F_NS_2023: (numeric) The calculated return time in years, representing the expected time interval between occurrences of the observed extreme temperature. 
#                 A higher return time indicates a rarer event.

calc_return_time <- function(data, loc, scale, shape) {
  
  # Calculate the exceedance probability for the given maximum temperature using the GEV distribution.
  # 'lower.tail = FALSE' indicates we are calculating the probability of the observed temperature exceeding the given value.
  # pgev is the GEV distribution function.
  pF_NS_2023 <- pgev(data, loc = loc, 
                     scale = scale, shape = shape, lower.tail = FALSE)
  
  # Calculate the return time (in years) by taking the inverse of the exceedance probability.
  Rt_F_NS_2023 <- 1 / pF_NS_2023
  
  # Return the calculated return time
  return(Rt_F_NS_2023)
}


# Function to perform bootstrap sampling and return time calculation for a particular year.
# This function uses bootstrap resampling to estimate the return time of extreme temperatures
# based on a non-stationary Generalized Extreme Value (GEV) distribution model.
#
# Args:
#   data_t: (data.frame) The dataset used for bootstrapping. It should have columns with names:
#          - 'date': Date of the observation
#          - 'tx': Maximum temperature
#          - 'en3.4': El Niño 3.4 index
#          - 'gwi': Global Warming Index
#   n_boot: (integer) The number of bootstrap iterations to perform.
#   EN3_4: (numeric) The El Niño 3.4 index value for the year of interest.
#   GWL: (numeric) Global warming level, where:
#        - Pre-industrial GWI is 0ºC
#        - Present GWI is approximately 1.29ºC
#        - Future GWI is projected to be around 2ºC
#   TXx: (numeric) The extreme maximum temperature value for which the return time is being calculated.
#
# Returns:
#   boot_return_times: (numeric vector) A vector of bootstrapped return times for the extreme temperature
#                      value (TXx). Each value represents the return time estimated from one bootstrap sample.

boot_return_time_calc <- function(data_t, n_boot, EN3_4, GWL, TXx) {
  
  # Initialize an empty vector to store return times from each bootstrap iteration
  boot_return_times <- numeric(n_boot)
  
  set.seed(123)  # Set seed for reproducibility
  
  for (i in 1:n_boot) {
    
    # Step 1: Create a bootstrap sample by sampling with replacement
    boot_indices <- sample(1:nrow(data_t), replace = TRUE)
    data_boot <- data_t[boot_indices, ]
    
    # Step 2: Fit a linear model with ENSO and GWI as predictors
    # The model predicts maximum temperature (tx) based on ENSO index and GWI
    model_loc_multi <- lm(tx ~ en3.4 + gwi, data = data_boot)
    loc_lm_multi <- predict(model_loc_multi)
    
    # Step 3: Fit a non-stationary GEV model to the residuals (difference between observed and predicted temperature)
    gev_NS_multi <- ROOPSD::GEV$new()$fit(data_boot$tx - loc_lm_multi)
    
    # Step 4: Extract the model coefficients and GEV parameters
    coef_BE_boot <- c(model_loc_multi$coefficients[1] + gev_NS_multi$loc, 
                      model_loc_multi$coefficients[2], 
                      model_loc_multi$coefficients[3], 
                      gev_NS_multi$scale, 
                      gev_NS_multi$shape)
    
    # Assign meaningful names to the coefficients
    names(coef_BE_boot) <- c("loc0", "loc1_enso", "loc2_gwi", "scale0", "shape0")
    
    # Step 5: Calculate the location parameter for the given GWL and EN3_4 values
    loc_boot <- coef_BE_boot["loc0"] + GWL * coef_BE_boot["loc2_gwi"] + EN3_4 * coef_BE_boot["loc1_enso"]
    
    # Create vectors for scale and shape parameters (constant across all rows)
    scale_BE_boot <- coef_BE_boot["scale0"]
    shape_BE_boot <- coef_BE_boot["shape0"]
    
    # Step 6: Calculate the return time for the extreme temperature value using the bootstrap sample
    boot_return_times[i] <- calc_return_time(TXx, loc_boot, scale_BE_boot, shape_BE_boot)
  }
  
  # Return the vector of bootstrapped return times
  return(boot_return_times)
}
# Function to perform bootstrap sampling and return time calculation for multiple return levels.
# This function estimates return times for various return levels based on a non-stationary Generalized Extreme Value (GEV) distribution model.
#
# Args:
#   data_t: (data.frame) The dataset used for bootstrapping. It should contain columns with names:
#          - 'date': Date of the observation
#          - 'tx': Maximum temperature values
#          - 'en3.4': El Niño 3.4 index values
#          - 'gwi': Global Warming Index values
#   n_boot: (integer) Number of bootstrap iterations to perform.
#   EN3_4: (numeric) El Niño 3.4 index value for the year of interest.
#   GWL: (numeric) Global warming level, where:
#        - Pre-industrial GWI is 0ºC
#        - Present GWI is approximately 1.29ºC
#        - Future GWI is projected to be around 2ºC
#   RL: (numeric vector) A vector of return levels for which return times are to be calculated. Each return level corresponds to a specific extreme temperature value.
#
# Returns:
#   boot_return_times: (matrix) A matrix of bootstrapped return times. Each row corresponds to a bootstrap iteration, and each column corresponds to a return level specified in `RL`.
#                      The matrix has dimensions (n_boot x length(RL)).

boot_return_time_calc_RL <- function(data_t, n_boot, EN3_4, GWL, RL) {
  
  # Initialize an empty matrix to store return times for each bootstrap iteration and return level
  boot_return_times <- matrix(NA, nrow = n_boot, ncol = length(RL))
  
  set.seed(123)  # Set seed for reproducibility
  
  for (i in 1:n_boot) {
    
    # Step 1: Create a bootstrap sample by sampling with replacement
    boot_indices <- sample(1:nrow(data_t), replace = TRUE)
    data_boot <- data_t[boot_indices, ]
    
    # Step 2: Fit a linear model with ENSO and GWI as predictors
    # The model predicts maximum temperature (tx) based on ENSO index and GWI
    model_loc_multi <- lm(tx ~ en3.4 + gwi, data = data_boot)
    loc_lm_multi <- predict(model_loc_multi)
    
    # Step 3: Fit a non-stationary GEV model to the residuals (observed temperature minus predicted location)
    gev_NS_multi <- ROOPSD::GEV$new()$fit(data_boot$tx - loc_lm_multi)
    
    # Step 4: Extract the model coefficients and GEV parameters
    coef_BE_boot <- c(model_loc_multi$coefficients[1] + gev_NS_multi$loc, 
                      model_loc_multi$coefficients[2], 
                      model_loc_multi$coefficients[3], 
                      gev_NS_multi$scale, 
                      gev_NS_multi$shape)
    
    # Assign meaningful names to the coefficients
    names(coef_BE_boot) <- c("loc0", "loc1_enso", "loc2_gwi", "scale0", "shape0")
    
    # Step 5: Calculate the location parameter for the given GWL and EN3_4 values
    loc_boot <- coef_BE_boot["loc0"] + GWL * coef_BE_boot["loc2_gwi"] + EN3_4 * coef_BE_boot["loc1_enso"]
    
    # Create vectors for scale and shape parameters (constant across all rows)
    scale_BE_boot <- coef_BE_boot["scale0"]
    shape_BE_boot <- coef_BE_boot["shape0"]
    
    # Step 6: Calculate return times for each return level specified in RL using the bootstrap sample
    boot_return_times[i, ] <- calc_return_time(RL, loc_boot, scale_BE_boot, shape_BE_boot)
  }
  
  # Return the matrix of bootstrapped return times
  return(boot_return_times)
}
# Function to compute the value of a linear function.
# This function calculates the result of a linear combination of two input variables plus an intercept term.
#
# Args:
#   x1: (numeric) The first predictor variable in the linear model.
#   x2: (numeric) The second predictor variable in the linear model.
#   beta1: (numeric) The coefficient (slope) for the first predictor variable (x1).
#   beta2: (numeric) The coefficient (slope) for the second predictor variable (x2).
#   beta0: (numeric) The intercept term of the linear function.
#
# Returns:
#   result: (numeric) The computed value of the linear function, given by the formula:
#           result = beta1 * x1 + beta2 * x2 + beta0

linear_function <- function(x1, x2, beta1, beta2, beta0) {
  # Compute the linear function value
  result <- beta1 * x1 + beta2 * x2 + beta0
  
  # Return the computed result
  return(result)
}


################################### Load data #############################################################


# Analysis over stations

# Daily data
setwd("C:/Users/SoledadCollazo/OneDrive - Universidad Complutense de Madrid (UCM)/Escritorio/atribucion/datos rio/")
dat <- read.csv2('tx7124_comp.csv',header=T)

code <- colnames(dat)[-1]

for (i in 1:length(code)){
  
  
  # Extract daily maximum temperature for station
  dat_station <- dat[[code[i]]]
  
  # Time
  dates <- seq(as.Date("1971/01/01"), as.Date("2024/03/20"), by = "day")
  years <- format(dates, "%Y")
  month <- format(dates, "%m")
  
  
  # Select analysis time frame
  w <- which(dates <= as.Date("2023-12-31"))
  dat_red <- dat_station[w]
  dates <- dates[w]
  
  datos <- data.frame(fecha = dates, valor = dat_red)
  
  # Extract the year from the date column
  datos <- datos %>%
    mutate(año = year(fecha))
  
  # Find the maximum value per year
  max_anual <- datos %>%
    group_by(año) %>%
    summarize(valor_maximo_anual = max(valor))
  
  # Print the results
  print(max_anual)
  
  
  # Covariables
  
  #### GWI ####
  cov <- read_xlsx("C:/Users/SoledadCollazo/OneDrive - Universidad Complutense de Madrid (UCM)/Escritorio/atribucion/atribución_igeo/data_raw/hadcrut5_wrt1850_1900.xlsx")
  
  # Apply LOESS
  loess_fit <- loess(cov$`[Celsius] annual mean of blended air_temperature_anomaly over land with sea_water_temperat` ~ cov$Time)
  
  # Predict smoothed values
  smoothed_values <- predict(loess_fit)
  
  # Plot GWI
  par(mar=c(4,4,1,1))
  
  plot(cov$`[Celsius] annual mean of blended air_temperature_anomaly over land with sea_water_temperat` ~ cov$Time, main="", col="blue", pch=16,
       ylab='GWI [ºC]', xlab='Year')
  lines(cov$Time, smoothed_values, col="red", lwd=2)
  grid()
  
  # Add smoothed GWI values to the existing covariate data
  cov <- cbind(cov, gwi = smoothed_values)
  colnames(cov)[1] <- 'año'
  
  # Merge covariate data with annual maximum values
  data <- merge(cov, max_anual)
  
  
  ### ENSO 3.4 ###
  
  cov <- read.table("C:/Users/SoledadCollazo/OneDrive - Universidad Complutense de Madrid (UCM)/Escritorio/atribucion/datos rio/enso3.4.txt",  na.strings = "-99.9")
  
  en<-as.vector(t(as.matrix(cov[2:13])))
  fec<-paste(sort(rep(cov[,1],12)),rep(1:12,nrow(cov)),sep='-')
  fecd<-data.frame(fecha=paste(years,as.numeric(month),sep='-'))
  enm<-data.frame(fecha=fec,enso=en)
  ensod<-merge(fecd,enm,sort=F)
  
  # Select analysis time frame
  dat_red <- data.frame(tx=dat_station[w],en3.4=ensod[w,2])
  
  datos <- data.frame(fecha=dates, dat_red)
  
  datos <- datos %>%
    mutate(año = year(fecha))
  
  
  # Find the maximum value per year
  max_anual_en <- datos %>%
    group_by(año) %>%
    filter(tx == max(tx))%>%  
    slice(1)  
  
  # Ver los resultados
  print(max_anual_en)
  
  data_en<-max_anual_en[,c(1:3)]
  
  
  # Data frame with the two covariables
  
  data_t<-cbind(data_en,data$gwi)
  colnames(data_t)[4]<-'gwi'
  
  
  ################################## Trend of TXx #####################################################
  
  # Fit a linear model
  model <- lm(valor_maximo_anual ~ año, data = max_anual)
  summary(model)
  
  
  # Perform Mann-Kendall test
  mk_result <- mk.test(max_anual$valor_maximo_anual)
  mk_result
  slope <- sens.slope(max_anual$valor_maximo_anual)$estimates
  
  # Create trend text
  
  if (mk_result$p.value < 0.05){
    trend_text <- paste0("Trend: ", round(slope, 3) * 10, " [ºC per decade] *")
  }else{
    trend_text <- paste0("Trend: ", round(slope, 3) * 10, " [ºC per decade]")
  }
  
  print(trend_text)
  
  # Plot the time series with customization
  ggplot(max_anual, aes(x = año, y = valor_maximo_anual)) +
    geom_line(color = "blue", linewidth = 1) +
    geom_point(color = "blue", size = 2.5) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
    labs(title = paste("Station", code[i], sep=' '),
         x = "Year",
         y = "TXx [ºC]") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) +
    annotate("text", x = 1972, y = max(max_anual$valor_maximo_anual), 
             label = trend_text, hjust = 0, vjust = 1, size = 4, color = "black")
  
  dev.print(jpeg, paste('ts_', code[i],'.jpg', sep=''), res = 700, width = 18, height = 13, units = 'cm')
  
  
  
  
  ############################## Fit GEV distribution #######################################
  
  # Stationary
  
  gev_S <- ROOPSD::GEV$new()$fit(data$valor_maximo_anual)
  gev_S$diagnostic(data$valor_maximo_anual)
  
  # Kolmogorov-Smirnov test
  gev_S$ks.test
  
  # Return period in the stationary model
  pF_S_2023 <- gev_S$sf(max(dat_station))
  Rt_F_S_2023 <- 1. / pF_S_2023
  cat(paste("Return time (2023):", round(Rt_F_S_2023, 2), "year"), end = "\n")
  
  
  # Non-stationary
  
  # 1) Covariable: GWI
  model_loc_gwi <- lm("valor_maximo_anual ~ gwi", data = data)
  loc_lm_gwi <- predict(model_loc_gwi)
  print(summary(model_loc_gwi)) 
  
  gev_NS_gwi <- ROOPSD::GEV$new()$fit(data$valor_maximo_anual - loc_lm_gwi)
  print(gev_NS_gwi$loc)
  
  gev_NS_gwi$diagnostic(data$valor_maximo_anual - loc_lm_gwi)
  
  # Kolmogorov-Smirnov test
  gev_NS_gwi$ks.test
  
  
  # 2) Covariable: EN3.4
  model_loc_enso = lm( " tx ~ en3.4" , data = data_en )
  loc_lm_enso = predict(model_loc_enso)
  print(summary(model_loc_enso)) 
  
  gev_NS_enso = ROOPSD::GEV$new()$fit( data_en$tx - loc_lm_enso )
  print(gev_NS_enso$loc)
  
  gev_NS_enso$diagnostic( data_en$tx - loc_lm_enso)
  
  # Kolmogorov-Smirnov test
  gev_NS_enso$ks.test
  
  
  # 3) Covariables: GWI and EN3.4
  model_loc_multi = lm( " tx ~ en3.4 + gwi" , data = data_t )
  loc_lm_multi = predict(model_loc_multi)
  print(summary(model_loc_multi)) 
  
  gev_NS_multi = ROOPSD::GEV$new()$fit( data_t$tx - loc_lm_multi )
  print(gev_NS_multi$loc)
  
  gev_NS_multi$diagnostic( data_en$tx - loc_lm_multi)
  
  # Kolmogorov-Smirnov test
  gev_NS_multi$ks.test
  
  
  ######################## Comparison of GEV models #############################################
  
  ## Assessing the Significance of Model Improvement
  
  # Stationary vs. Non-stationary GWI
  
  alpha <- seq(0.01, 0.3, 0.01)
  
  lll_S <- sum(gev_S$logdensity(data$valor_maximo_anual))
  lll_NS_gwi <- sum(gev_NS_gwi$logdensity(data$valor_maximo_anual - loc_lm_gwi))
  
  reject_S <- 2 * (lll_NS_gwi - lll_S) > qchisq(1 - alpha, df = 1)
  par(mfrow=c(1,1))
  plot(alpha, reject_S, yaxt = "n", ylim = c(-0.1, 1.1), col = "red")
  axis(2, at = c(0, 1), labels = c("False", "True"))
  
  
  # Stationary vs. Non-stationary EN3.4
  
  lll_NS_enso  = sum(gev_NS_enso$logdensity(data_en$tx-loc_lm_enso))
  
  reject_S = 2 * (lll_NS_enso - lll_S) > qchisq( 1 - alpha , df = 1 )
  par(mfrow=c(1,1))
  plot( alpha , reject_S , yaxt = "n" , ylim = base::c( -0.1 , 1.1 ) , col = "red" )
  axis( 2 , at = base::c(0,1) , labels = base::c("False","True") )
  
  # Non-stationary GWI vs. Non-stationary EN3.4
  
  reject_NS_enso = 2 * (lll_NS_gwi - lll_NS_enso) > qchisq( 1 - alpha , df = 1 )
  par(mfrow=c(1,1))
  plot( alpha , reject_NS_enso , yaxt = "n" , ylim = base::c( -0.1 , 1.1 ) , col = "red" )
  axis( 2 , at = base::c(0,1) , labels = base::c("False","True") )
  
  reject_NS_ensob = 2 * (lll_NS_enso - lll_NS_gwi) > qchisq( 1 - alpha , df = 1 )
  par(mfrow=c(1,1))
  plot( alpha , reject_NS_ensob , yaxt = "n" , ylim = base::c( -0.1 , 1.1 ) , col = "red" )
  axis( 2 , at = base::c(0,1) , labels = base::c("False","True") )
  
  
  # Stationary vs. Non-stationary multi
  
  lll_NS_multi  = sum(gev_NS_multi$logdensity(data_t$tx-loc_lm_multi))
  
  reject_S = 2 * (lll_NS_multi - lll_S) > qchisq( 1 - alpha , df = 1 )
  par(mfrow=c(1,1))
  plot( alpha , reject_S , yaxt = "n" , ylim = base::c( -0.1 , 1.1 ) , col = "red" )
  axis( 2 , at = base::c(0,1) , labels = base::c("False","True") )
  
  # Non-stationary multi vs. Non-stationary EN3.4
  
  reject_NS_enso = 2 * (lll_NS_multi - lll_NS_enso) > qchisq( 1 - alpha , df = 1 )
  par(mfrow=c(1,1))
  plot( alpha , reject_NS_enso , yaxt = "n" , ylim = base::c( -0.1 , 1.1 ) , col = "red" )
  axis( 2 , at = base::c(0,1) , labels = base::c("False","True") )
  
  # Non-stationary GWI vs. Non-stationary multi
  
  reject_NS_gwi = 2 * (lll_NS_multi - lll_NS_gwi) > qchisq( 1 - alpha , df = 1 )
  par(mfrow=c(1,1))
  plot( alpha , reject_NS_gwi , yaxt = "n" , ylim = base::c( -0.1 , 1.1 ) , col = "red" )
  axis( 2 , at = base::c(0,1) , labels = base::c("False","True") )
  
  
  reject_NS_gwi = 2 * (lll_NS_gwi - lll_NS_multi) > qchisq( 1 - alpha , df = 1 )
  par(mfrow=c(1,1))
  plot( alpha , reject_NS_gwi , yaxt = "n" , ylim = base::c( -0.1 , 1.1 ) , col = "red" )
  axis( 2 , at = base::c(0,1) , labels = base::c("False","True") )
  
  
  ## AIC
  
  # Log-Likelihood Values for Model Comparisons
  
  LL_station <- c(lll_S, lll_NS_enso, lll_NS_gwi, lll_NS_multi)
  
  # Number of Parameters for Each Model
  
  n_param <- c(3, 4, 4, 5)
  
  # Calculation of AIC for Each Model
  
  AIC_station <- -2 * LL_station + 2 * n_param
  
  # Create a data frame with column names
  AIC_station_df <- data.frame(
    Model = c('S', 'EN3.4', 'GWI', 'multi'),
    AIC = AIC_station
  )
  
  # Print the data frame
  print(AIC_station_df)
  
  
  ##################################### Coefficients of the NSGEV multi #########################################
  
  coef_BE = base::c( model_loc_multi$coefficients[1] + gev_NS_multi$loc , model_loc_multi$coefficients[2], model_loc_multi$coefficients[3] , gev_NS_multi$scale , gev_NS_multi$shape )
  names(coef_BE) = base::c( "loc0" , "loc1_enso" , "loc2_gwi", "scale0" , "shape0" )
  print(coef_BE)
  
  
  coef_station<-coef_BE
  
  # Function parameters
  beta1 <- coef_BE[2]  # Coefficient for en3.4
  beta2 <- coef_BE[3]  # Coefficient for gwi
  int <- coef_BE[1]    # Intercept term
  
  scale <- coef_BE[4]  # Scale parameter
  shape <- coef_BE[5]  # Shape parameter
  
  # Upper bound of the distribution
  
  if (shape < 0) {
    
    # Define the range of values for x1 (EN3.4) and x2 (GWI)
    en3.4_values <- seq(-2.5, 2.5, length.out = 1000)  # From -2.5 to 2.5 with 1000 points
    gwi_values <- seq(-0.1, 2.5, length.out = 500)  # From -0.1 to 2.5 with 500 points
    
    # Create all combinations of en3.4 and gwi
    combinations <- expand.grid(en3.4 = en3.4_values, gwi = gwi_values)
    
    # Apply the function upper_bound to each combination
    combinations$limit <- mapply(upper_bound, 
                                 en3.4 = combinations$en3.4, 
                                 gwi = combinations$gwi, 
                                 MoreArgs = list(int = int, beta1 = beta1, beta2 = beta2, scale = scale, shape = shape))
    
    # View the result
    print(combinations)
    
  }
  
  # Plot upper bound
  
  combinations$rp_cut <- cut(combinations$limit, 
                             breaks = seq(39,50,by=1), 
                             right = FALSE)
  
  # Definir una paleta de colores
  colors <- brewer.ylorrd(10)
  
  # Crear el gráfico de contornos
  # Crear el gráfico de contornos
  ggplot(combinations, aes(x = gwi, y = en3.4, z = limit)) +
    theme_bw() +
    geom_contour_filled(aes(fill = rp_cut)) +
    scale_fill_manual(values = colors,
                      breaks = levels(combinations$rp_cut)) +
    labs(title = paste("Station", code[i], sep=' '),
         x = "GWI [ºC]",
         y = "EN3.4 [ºC]",
         fill = "Upper bound [ºC]") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey40", linewidth=0.5) + 
    geom_vline(xintercept = 1.29, linetype = "dotted", color = "blue", linewidth=0.5) + 
    geom_vline(xintercept = 2, linetype = "dotted", color = "magenta", linewidth=0.5) + 
    theme(panel.ontop=TRUE,panel.grid.major = element_line(color = "gray", linetype = "dashed",linewidth=0.2), panel.grid.minor = element_blank(),
          panel.background = element_rect(color = NA, fill = NA, linetype = "dashed"),
          axis.text = element_text(size = 5),  # Tamaño del texto del eje
          axis.title = element_text(size = 5),  # Tamaño del título del eje
          axis.line = element_line(colour = "black"),
          plot.title = element_text(size=6),
          legend.text = element_text(size=4),
          legend.title = element_text(size=5),
          legend.key.height  = unit(0.3, "cm"),
          legend.key.width  = unit(0.15, "cm"),
          strip.text = element_text(size = 5),
          legend.margin=margin(0,0,0,0),
          panel.border = element_rect(color = "black", linewidth = 0.01),
          axis.ticks = element_line(linewidth = 0.05))+
    scale_x_continuous(expand = c(0,0),limits = c(-0.1,2.5)) +  # Elimina el espacio alrededor del eje x
    scale_y_continuous(expand = c(0, 0))    # Elimina el espacio alrededor del eje y
  
  dev.print(jpeg, paste(code[i],'_gev_multi_upper_bound.jpg', sep=''), res=700, height=6,width=7, units='cm')
  
  
  
  ################################## Return period ##################################################
  
  anio.max<- nrow(data_t)
  TXx2023 <- data_t$tx[anio.max]
  
  ### Pre-industrial climate
  
  loc_preind_BE  = coef_BE[1]  + data_t$en3.4[anio.max] * coef_BE[2] - 0 * coef_BE[3]
  scale_BE = rep( coef_BE[4] , nrow(data) )
  shape_BE = rep( coef_BE[5] , nrow(data) )
  
  
  Rt_F_NS_2023 = calc_return_time(TXx2023, loc_preind_BE, scale, shape)
  cat( paste( "Return time 2023 (ContraFactual - Pre-industrial):" , round(Rt_F_NS_2023,2) , "year" ) , end = "\n" )
  
  RT_station<-c()
  RT_station<-Rt_F_NS_2023
  
  # Configure number of bootstrap repetitions
  n_boot <- 1000
  
  boot_RP <- boot_return_time_calc(data_t, n_boot, EN3_4 = 2.02, GWL = 0, TXx2023) # El Niño and Pre-industrial
  
  # Calculate the 90% confidence interval
  ci_lower <- quantile(boot_RP, 0.05)
  ci_upper <- quantile(boot_RP, 0.95)
  
  
  cat(paste("Return time 2023 (ContraFactual):", round(Rt_F_NS_2023, 2), "year"), "\n")
  cat(paste("90% CI:", round(ci_lower, 2), "to", round(ci_upper, 2), "years"), "\n")
  
  # Estimating the RP for different return levels (RL)
  
  RL<-seq(34,47,by=0.1)
  
  Rt_F_NS_2023 <- calc_return_time(RL, loc_preind_BE, scale_BE, shape_BE)
  
  boot_RP_RL <- boot_return_time_calc_RL(data_t, n_boot, EN3_4 = 2.02, GWL = 0, RL) # El Niño and Pre-industrial
  
  
  # Confidence interval at 90%
  
  percentil_95 <- apply(boot_RP_RL, 2, quantile, probs = 0.95)
  percentil_5 <- apply(boot_RP_RL, 2, quantile, probs = 0.05)
  
  prein<-data.frame(ret_lev=RL,p5=percentil_5,p95=percentil_95)
  
  
  # Create an empirical cumulative distribution function (ECDF) for the residuals
  # The ECDF is based on the differences between the observed 'tx' values and the pre-industrial location parameter
  ajuste_ref <- ecdf(data_t$tx - loc_preind_BE)
  
  # Calculate the exceedance probability for each observation
  # This is done by subtracting the ECDF value from 1, which gives the quantile (i.e., the probability of exceeding the value)
  cuantil_x <- 1 - ajuste_ref(data_t$tx - loc_preind_BE)
  
  # Compute the return period (RP) based on the exceedance probability
  # The return period is the inverse of the exceedance probability, giving the expected number of years for an event of that magnitude to occur
  RP_preind_emp <- 1 / cuantil_x
  
  
  # Plot pre-industrial results
  
  par(mfrow=c(1,1))
  
  plot(Rt_F_NS_2023,RL, type='n', xlim=c(1,500), log= 'x', xlab= 'Return period (years)', xaxt='n',
       ylab="Return level (ºC)", main = paste("Station", code[i], sep=' '), sub = "Pre-industrial (black) vs. present (blue) vs. future (magenta)")
  
  axis(1,at=c(1,10,100,500),labels=c(1,10,100,500))
  
  abline(v = c(seq(1,10,1),seq(10,100,10),seq(100,500,100)), col = "gray40", lty = 3, lwd=0.5)
  
  lines(RL~Rt_F_NS_2023)
  
  aux<- data.frame(x=sort(RP_preind_emp), y=sort(data_t$tx - (data_t$gwi * coef_BE[3] + 0 * coef_BE[3] - 2 * coef_BE[2])))
  
  resul<- aggregate(y ~ x, data = aux, FUN = max)
  
  points(resul$y ~ resul$x, cex=0.5, pch=16)
  
  
  x_points <- c(percentil_95[is.finite(percentil_95)], rev(percentil_5))
  
  y_points <- c(RL[is.finite(percentil_95)],rev(RL))
  
  
  polygon(x_points, y_points, col=rgb(0.5, 0.5, 0.5, alpha=0.3), border=NA)
  
  abline(h=data$valor_maximo_anual[length(data$valor_maximo_anual)],col='darkgreen',lty=2, lwd=2)
  
  
  
  
  ###Present climate
  
  loc_present_BE  = coef_BE[1] + data$gwi * coef_BE[3] + data_t$en3.4[anio.max] * coef_BE[2]
  
  loc_pres<-loc_present_BE[anio.max]
  
  Rt_F_NS_2023 = calc_return_time(TXx2023, loc_pres, scale_BE, shape_BE)
  
  cat( paste( "Return time 2023 (Factual):" , round(Rt_F_NS_2023,2) , "year" ) , end = "\n" )
  
  RT_station <- c(RT_station, Rt_F_NS_2023)
  
  
  # Configure number of bootstrap repetitions
  n_boot <- 1000
  
  boot_RP <- boot_return_time_calc(data_t, n_boot, EN3_4 = 2.02, GWL = 1.29, TXx2023) # El Niño and Present
  
  # Calculate the 90% confidence interval
  ci_lower <- quantile(boot_RP, 0.05)
  ci_upper <- quantile(boot_RP, 0.95)
  
  
  cat(paste("Return time 2023 (Factual):", round(Rt_F_NS_2023, 2), "year"), "\n")
  cat(paste("90% CI:", round(ci_lower, 2), "to", round(ci_upper, 2), "years"), "\n")
  
  # Estimating the RP for different return levels (RL)
  
  RL<-seq(34,47,by=0.1)
  
  Rt_F_NS_2023 <- calc_return_time(RL, loc_pres, scale_BE, shape_BE)
  
  boot_RP_RL <- boot_return_time_calc_RL(data_t, n_boot, EN3_4 = 2.02, GWL = 1.29, RL) # El Niño and Present
  
  
  # Confidence interval at 90%
  
  percentil_95 <- apply(boot_RP_RL, 2, quantile, probs = 0.95)
  percentil_5 <- apply(boot_RP_RL, 2, quantile, probs = 0.05)
  
  pres<-data.frame(ret_lev=RL,p5=percentil_5,p95=percentil_95)
  
  # Create an empirical cumulative distribution function (ECDF) for the residuals
  # The ECDF is based on the differences between the observed 'tx' values and the pre-industrial location parameter
  ajuste_ref <- ecdf(data_t$tx - loc_pres)
  
  # Calculate the exceedance probability for each observation
  # This is done by subtracting the ECDF value from 1, which gives the quantile (i.e., the probability of exceeding the value)
  cuantil_x <- 1 - ajuste_ref(data_t$tx - loc_pres)
  
  # Compute the return period (RP) based on the exceedance probability
  # The return period is the inverse of the exceedance probability, giving the expected number of years for an event of that magnitude to occur
  RP_pres_emp <- 1 / cuantil_x
  
  
  # Plot present climate results
  
  
  lines(RL~Rt_F_NS_2023, col='blue')
  
  aux<- data.frame(x=sort(RP_pres_emp), y=sort(data_t$tx - (data_t$gwi * coef_BE[3] - 1.29 * coef_BE[3] - 2 * coef_BE[2])))
  
  resul<- aggregate(y ~ x, data = aux, FUN = max)
  
  points(resul$y ~ resul$x, cex=0.5, pch=16,col='blue')
  
  
  x_points <- c(percentil_95[is.finite(percentil_95)], rev(percentil_5))
  
  y_points <- c(RL[is.finite(percentil_95)],rev(RL))
  
  polygon(x_points, y_points, col=rgb(0, 0, 1, alpha=0.3), border=NA)
  
  
  
  ### Future climate
  
  
  loc_future_BE  = coef_BE[1] + 2 * coef_BE[3] + data_t$en3.4[anio.max] * coef_BE[2]
  scale_BE = coef_BE[4] 
  shape_BE = coef_BE[5]
  
  Rt_F_NS_2023 = calc_return_time(TXx2023, loc_future_BE, scale_BE, shape_BE)
  
  cat( paste( "Return time 2023 (Future):" , round(Rt_F_NS_2023,2) , "year" ) , end = "\n" )
  
  RT_station <- c(RT_station, Rt_F_NS_2023)
  
  
  # Configure number of bootstrap repetitions
  n_boot <- 1000
  
  boot_RP <- boot_return_time_calc(data_t, n_boot, EN3_4 = 2.02, GWL = 2.00, TXx2023) # El Niño and Future climate
  
  # Calculate the 90% confidence interval
  ci_lower <- quantile(boot_RP, 0.05)
  ci_upper <- quantile(boot_RP, 0.95)
  
  
  cat(paste("Return time 2023 (Future):", round(Rt_F_NS_2023, 2), "year"), "\n")
  cat(paste("90% CI:", round(ci_lower, 2), "to", round(ci_upper, 2), "years"), "\n")
  
  # Estimating the RP for different return levels (RL)
  
  RL<-seq(34,47,by=0.1)
  
  Rt_F_NS_2023 <- calc_return_time(RL, loc_future_BE, scale_BE, shape_BE)
  
  boot_RP_RL <- boot_return_time_calc_RL(data_t, n_boot, EN3_4 = 2.02, GWL = 2.00, RL) # El Niño and Future climate
  
  
  # Confidence interval at 90%
  
  percentil_95 <- apply(boot_RP_RL, 2, quantile, probs = 0.95)
  percentil_5 <- apply(boot_RP_RL, 2, quantile, probs = 0.05)
  
  fut<-data.frame(ret_lev=RL,p5=percentil_5,p95=percentil_95)
  
  # Create an empirical cumulative distribution function (ECDF) for the residuals
  # The ECDF is based on the differences between the observed 'tx' values and the pre-industrial location parameter
  ajuste_ref <- ecdf(data_t$tx - loc_future_BE)
  
  # Calculate the exceedance probability for each observation
  # This is done by subtracting the ECDF value from 1, which gives the quantile (i.e., the probability of exceeding the value)
  cuantil_x <- 1 - ajuste_ref(data_t$tx - loc_future_BE)
  
  # Compute the return period (RP) based on the exceedance probability
  # The return period is the inverse of the exceedance probability, giving the expected number of years for an event of that magnitude to occur
  RP_fut_emp <- 1 / cuantil_x
  
  
  # Plot future climate results
  
  lines(RL~Rt_F_NS_2023, col='magenta')
  
  aux<- data.frame(x=sort(RP_fut_emp), y=sort(data_t$tx - (data_t$gwi * coef_BE[3] - 2* coef_BE[3] - 2 * coef_BE[2])))
  
  resul<- aggregate(y ~ x, data = aux, FUN = max)
  
  points(resul$y ~ resul$x, cex=0.5, pch=16, col='magenta')
  
  x_points <- c(percentil_95[is.finite(percentil_95)], rev(percentil_5))
  
  y_points <- c(RL[is.finite(percentil_95)],rev(RL))
  
  polygon(x_points, y_points, col=rgb(1, 0.41, 0.71, alpha=0.3), border=NA)
  
  
  dev.print(jpeg,paste(code[i], '_gev_multi.jpg', sep=''), res=700, height=15,width=22, units='cm')
  
  
  
  
  
  # 2-dimensional Return period
  
  # Parameters of the function
  beta1 <- coef_BE[2] # Coefficient for x1: en3.4
  beta2 <- coef_BE[3] # Coefficient for x2:gwi
  beta0 <- coef_BE[1] # Independent term
  
  # Define the range of values for x1 and x2
  en3.4_values <- seq(-2.5, 2.5, length.out = 1000) # From 0 to 10 with 5 points
  gwi_values <- seq(-0.1, 2.5, length.out = 500) # From 0 to 10 with 5 points
  
  # Calculate the y-values for each combination of x1 and x2
  results <- expand.grid(x1 = en3.4_values, x2 = gwi_values)
  results$y <- mapply(linear_function, results$x1, results$x2, MoreArgs = list(beta1 = beta1, beta2 = beta2, beta0 = beta0))
  
  # Display the results
  print(results)
  
  
  prob<-c()
  
  for (j in 1:nrow(results)){
    
    pF_NS_2023 = pgev( TXx2023, loc = results[j,3] , 
                       scale = scale_BE , shape = shape_BE , lower.tail = FALSE )
    
    prob[j]<- pF_NS_2023
    
  }
  
  rp<-1/prob
  
  res_prob<-cbind(results,prob,rp)
  
  
  
  
  # Create the contour plot for probability
  
  ggplot(res_prob, aes(x = x2, y = x1, z = prob)) +
    theme_bw() +
    geom_raster(aes(fill = prob)) +  # Replace geom_contour_filled
    scale_fill_stepsn(colors = brewer.ylorrd(32),
                      breaks = seq(0,0.7,by=0.05),
                      limits = c(min(res_prob$prob), 0.7),
                      guide = guide_colorbar(barheight = unit(10, "cm"))) +  # New scale with limits
    labs(title = paste("Station", code[i], sep=' '),
         x = "GWI",
         y = "EN3.4",
         fill = "Probability") +
    theme(panel.grid.major = element_line(color = "gray", linetype = "dashed"),
          panel.grid.minor = element_blank(),
          panel.ontop = TRUE, panel.background = element_rect(color = NA, fill = NA, linetype = "dashed"),
          axis.text = element_text(size = 10),  # Tamaño del texto del eje
          axis.title = element_text(size = 12),  # Tamaño del título del eje
          axis.line = element_line(colour = "black")) +
    scale_x_continuous(expand = c(0, 0)) +  # Elimina el espacio alrededor del eje x
    scale_y_continuous(expand = c(0, 0))    # Elimina el espacio alrededor del eje y
  
  dev.print(jpeg, paste(code[i], '_gev_multi_prob.jpg', sep=''), res=700, height=15,width=18, units='cm')
  
  
  # Return period
  
  res_prob$rp_cut <- cut(res_prob$rp, 
                         breaks = c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, Inf), 
                         labels = c("[1,2)", "[2,5)", "[5,10)", "[10,20)", "[20,50)", 
                                    "[50,100)", "[100,200)", "[200,500)", "[500,1000)", "[1000,Inf)"), 
                         right = FALSE)
  
  # Define a color palette
  colors <- rev(brewer.ylorrd(10))
  
  # Create the contour plot
  ita <- ggplot(res_prob, aes(x = x2, y = x1, z = rp)) +
    theme_bw() +
    geom_contour_filled(aes(fill = rp_cut)) +
    scale_fill_manual(values = colors,
                      breaks = levels(res_prob$rp_cut)) +
    labs(title = paste("Station", code[i], sep=' '),
         x = "GWI [ºC]",
         y = "EN3.4 [ºC]",
         fill = "RP") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey40", linewidth=0.5) + 
    geom_vline(xintercept = 1.29, linetype = "dotted", color = "blue", linewidth=0.5) + 
    geom_vline(xintercept = 2, linetype = "dotted", color = "magenta", linewidth=0.5) + 
    theme(panel.ontop = TRUE, panel.grid.major = element_line(color = "gray", linetype = "dashed", linewidth = 0.2), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(color = NA, fill = NA, linetype = "dashed"),
          axis.text = element_text(size = 5),  # Axis text size
          axis.title = element_text(size = 5),  # Axis title size
          axis.line = element_line(colour = "black"),
          plot.title = element_text(size = 6),
          legend.text = element_text(size = 4),
          legend.title = element_text(size = 5),
          legend.key.height  = unit(0.3, "cm"),
          legend.key.width  = unit(0.15, "cm"),
          strip.text = element_text(size = 5),
          legend.margin = margin(0, 0, 0, 0),
          panel.border = element_rect(color = "black", linewidth = 0.01),
          axis.ticks = element_line(linewidth = 0.05)) +
    scale_x_continuous(expand = c(0, 0), limits = c(-0.1, 2.5)) +  # Removes space around the x-axis
    scale_y_continuous(expand = c(0, 0))    # Removes space around the y-axis
  
  ita
  
  dev.print(jpeg, paste(code[i], '_gev_multi_rp_bis.jpg', sep=''), res=700, height=6,width=7, units='cm')
  
  
}
