library(openair)
library(tidyverse)
library(forecast)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(mice)
library(reshape2)
library(gridExtra)
library(tseries)
library(Kendall)
library(VIM)
library(lubridate)


#Hourly Data
aurn_hourly<- importAURN(site = c('MY1','LEED','BIRR','MAN3','NOTT'),data_type='hourly',year=2019:2023,pollutant = c('no2','o3','pm10','pm2.5'),hc=FALSE,meta = FALSE)

#Removing Meterological columns from main dataset
aurn_hourly<- subset(aurn_hourly,select = -c(ws,wd,air_temp))

aurn_hourly_subset<- subset(aurn_hourly,select = -c(source,date,site,code))
aggr_plot <- aggr(aurn_hourly_subset, col=c('navyblue', 'red'), numbers=TRUE, sortVars=TRUE, labels=names(data), cex.axis=.7, gap=3, ylab=c("Histogram of missing data", "Pattern"))
summary(aggr_plot)

#Imputation for missing values using MICE technique
data_for_imputation <- aurn_hourly[, c("no2", "o3", 'pm10', "pm2.5")]
imputed_data <- mice(data_for_imputation, m = 5, method = 'pmm', maxit = 5, seed = 500)
completed_data <- complete(imputed_data, 1)
aurn_hourly[, c("imputed_no2", "imputed_o3", "imputed_pm10", "imputed_pm2.5")] <- completed_data

#Seperate dataset for imputed values
imputed_aurn_data <- subset(aurn_hourly,select = c(source,date,site,code,imputed_no2,imputed_o3,imputed_pm2.5,imputed_pm10))
summary(imputed_aurn_data)


#Histogram
hist_no2<-ggplot(imputed_aurn_data, aes(x = imputed_no2)) +
  geom_histogram(binwidth = 10, fill = "gray", color = "black", alpha = 0.7) +
  labs(title = "Histogram for NO2 concentration",
       x = "NO2 concentration",
       y = "Frequency") +
  theme_minimal()

hist_o3<-ggplot(imputed_aurn_data, aes(x = imputed_o3)) +
  geom_histogram(binwidth = 10, fill = "gray", color = "black", alpha = 0.7) +
  labs(title = "Histogram for O3 concentration",
       x = "O3 concentration",
       y = "Frequency") +
  theme_minimal()

hist_pm2.5<-ggplot(imputed_aurn_data, aes(x = imputed_pm2.5)) +
  geom_histogram(binwidth = 10, fill = "gray", color = "black", alpha = 0.7) +
  labs(title = "Histogram for PM2.5 concentration",
       x = "PM2.5 concentration",
       y = "Frequency") +
  theme_minimal()

hist_pm10<-ggplot(imputed_aurn_data, aes(x = imputed_pm10)) +
  geom_histogram(binwidth = 10, fill = "gray", color = "black", alpha = 0.7) +
  labs(title = "Histogram for PM10 concentration",
       x = "PM10 concentration",
       y = "Frequency") +
  theme_minimal()

grid.arrange(hist_no2,hist_o3,hist_pm2.5,hist_pm10,ncol=2,nrow=2)


# Calculate the percentage of missing data in the NO2 variable
missing_no2 <- sum(is.na(aurn_hourly$no2))
total_no2 <- nrow(aurn_hourly)
percentage_missing_no2 <- (missing_no2 / total_no2) * 100
percentage_missing_no2

# Calculate the percentage of missing data in the O3 variable
missing_o3 <- sum(is.na(aurn_hourly$o3))
total_o3 <- nrow(aurn_hourly)
percentage_missing_o3 <- (missing_o3 / total_o3) * 100
percentage_missing_o3

missing_pm10 <- sum(is.na(aurn_hourly$pm10))
total_pm10 <- nrow(aurn_hourly)
percentage_missing_pm10 <- (missing_pm10 / total_pm10) * 100
percentage_missing_pm10

missing_pm2.5 <- sum(is.na(aurn_hourly$pm2.5))
total_pm2.5 <- nrow(aurn_hourly)
percentage_missing_pm2.5 <- (missing_pm2.5 / total_pm2.5) * 100
percentage_missing_pm2.5

# Convert the datetime column to POSIXct
imputed_aurn_data$date <- as.POSIXct(imputed_aurn_data$date)
imputed_aurn_data$date<- as.Date(imputed_aurn_data$date)

#Extracting only year
imputed_aurn_data$year<- year(imputed_aurn_data$date)
imputed_aurn_data$year <- as.factor(imputed_aurn_data$year)

#Extracting only month
imputed_aurn_data$month <-month(imputed_aurn_data$date,label = TRUE)


#####################################
# Monthly Data
#####################################
#Summary of hourly data into monthly
monthly_data <- imputed_aurn_data%>%
 group_by(site,code,month,year,date) %>%
  summarise(
    mean_no2 = mean(imputed_no2, na.rm = TRUE),
    mean_o3 = mean(imputed_o3, na.rm = TRUE),
    mean_pm10 = mean(imputed_pm10, na.rm = TRUE),
    mean_pm2.5 = mean(imputed_pm2.5, na.rm = TRUE),
    .groups = 'drop'
  )
monthly_data$date <- as.Date(monthly_data$date)

TheilSen(imputed_aurn_data,pollutant = 'imputed_no2',deseason=TRUE,type = "site",date.format = "%m/%y",ylab='Nitrogen Dioxide (µg/m³)')
TheilSen(imputed_aurn_data,pollutant = 'imputed_o3',deseason=TRUE,type = "site",date.format = "%m/%y",,ylab='Ozone (µg/m³)')

########################################################
# Yearly Data
########################################################

# Summary of Hourly data to yearly mean values
yearly_data <- imputed_aurn_data %>%
  group_by(site,year,code) %>%
  summarise(
    mean_no2 = mean(imputed_no2, na.rm = TRUE),
    mean_o3 = mean(imputed_o3, na.rm = TRUE),
    mean_pm10 = mean(imputed_pm10, na.rm = TRUE),
    mean_pm2.5 = mean(imputed_pm2.5, na.rm = TRUE),
    .groups = 'drop'
  )

plot_correlation <- function(x, y, yearly_data, x_label, y_label) {
  ggplot(yearly_data, aes_string(x = x, y = y)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    stat_cor(method = "pearson", label.x = min(yearly_data[[x]]) * 1.1, label.y = max(yearly_data[[y]]) * 0.9) +
    labs(x = x_label, y = y_label) +
    theme_minimal()
}

p1 <- plot_correlation("mean_no2", "mean_pm2.5", yearly_data, "NO2 (µg/m³)", "PM2.5 (µg/m³)")
p2 <- plot_correlation("mean_o3", "mean_no2", yearly_data, "O3 (µg/m³)", "NO2(µg/m³)")
p3 <- plot_correlation("mean_no2", "mean_pm10", yearly_data, "NO2 (µg/m³)", "PM10 (µg/m³)")
p4 <- plot_correlation("mean_o3", "mean_pm10", yearly_data, "O3 (µg/m³)", "PM10 (µg/m³)")
ggarrange(p1, p2, p3, p4,ncol = 2, nrow = 2)

##################################################
# Line Plot for Yearly
##################################################

#Line Plot for NO2 levels at multiple sites Yearly
no2_plot<-ggline(data = yearly_data, 
       x = "year", 
       y = "mean_no2", 
       color = "site",
       ylab = "NO2 (µg/m³)",
       xlab = "Year",
       title = "NO2 Levels at Multiple Sites in UK") +
  theme_minimal()

#Line Plot for O3 levels at multiple sites Yearly
o3_plot<-ggline(data = yearly_data, 
       x = "year", 
       y = "mean_o3", 
       color = "site",
       ylab = "O3 (µg/m³)",
       xlab = "Year",
       title = "O3 Levels at Multiple Sites in UK") +
  theme_minimal()

#Line Plot for pm2.5 levels at multiple sites Yearly
pm2.5_plot<-ggline(data = yearly_data, 
       x = "year", 
       y = "mean_pm2.5", 
       color = "site",
       ylab = "pm2.5 (µg/m³)",
       xlab = "Year",
       title = "PM2.5 Levels at Multiple Sites in UK") +
  theme_minimal()

#Line Plot for pm10 levels at multiple sites Yearly
pm10_plot<-ggline(data = yearly_data, 
       x = "year", 
       y = "mean_pm10", 
       color = "site",
       ylab = "PM10 (µg/m³)",
       xlab = "Year",
       title = "PM10 Levels at Multiple Sites in UK") +
  theme_minimal()

grid.arrange(no2_plot,o3_plot,pm2.5_plot,pm10_plot)

#######################################################
# Data Visualization for the Year- 2023
#######################################################
aurn_2023<- importAURN(site = c('MY1','LEED','BIRR','MAN3','NOTT'),data_type='monthly',year=2023,pollutant = c('no2','o3','pm10','pm2.5'),hc=FALSE,meta = FALSE)
aurn_2023<- subset(aurn_2023,select = -c(o3_capture,no2_capture,pm10_capture,pm2.5_capture,uka_code))

data_2023_imputation <- aurn_2023[, c("no2", "o3", 'pm10', "pm2.5")]
imputed_data_2023 <- mice(data_2023_imputation, m = 5, method = 'pmm', maxit = 5, seed = 500)
completed_data_2023 <- complete(imputed_data_2023, 1)
aurn_2023[, c("no2", "o3", "pm10", "pm2.5")] <- completed_data_2023

aurn_2023$month <-month(aurn_2023$date,label = TRUE)
aurn_2023$year<- year(aurn_2023$date)

data_2023<- monthly_2023_data %>%
  filter(year == 2023)

data_2023$date <- as.Date(data_2023$date)

monthly_2023_data <- data_2023 %>%
  group_by(site,month,code,date) %>%
  summarise(
    mean_no2 = mean(mean_no2, na.rm = TRUE),
    mean_o3 = mean(mean_o3, na.rm = TRUE),
    mean_pm10 = mean(mean_pm10, na.rm = TRUE),
    mean_pm2.5 = mean(mean_pm2.5, na.rm = TRUE),
    .groups = 'drop'
  )

#Line plot for No2 level in year 2023
no2_2023 <- ggline(data = monthly_2023_data, 
                   x = "month", 
                   y = "mean_no2", 
                   color = "site",
                   ylab = "NO2 (µg/m³)",
                   xlab = "Month",
                   title = "NO2 Levels at Multiple Sites in UK (2023)") +
  theme_minimal()

#Line plot for O3 level in year 2023
o3_2023 <- ggline(data = monthly_2023_data, 
                  x = "month", 
                  y = "mean_o3", 
                  color = "site",
                  ylab = "O3 (µg/m³)",
                  xlab = "Month",
                  title = "O3 Levels at Multiple Sites in UK (2023)") +
  theme_minimal()

#Line plot for pm2.5 level in year 2023
pm2.5_2023 <- ggline(data = monthly_2023_data, 
                     x = "month", 
                     y = "mean_pm2.5", 
                     color = "site",
                     ylab = "pm2.5 (µg/m³)",
                     xlab = "Month",
                     title = "PM2.5 Levels at Multiple Sites in UK (2023)") +
  theme_minimal()

#Line plot for pm10 level in year 2023
pm10_2023 <- ggline(data = monthly_2023_data, 
                    x = "month", 
                    y = "mean_pm10", 
                    color = "site",
                    ylab = "pm10 (µg/m³)",
                    xlab = "Month",
                    title = "PM10 Levels at Multiple Sites in UK (2023)") +
  theme_minimal()

#Combining all the plots of year 2023 in one frame.
grid.arrange(no2_2023,o3_2023,pm2.5_2023,pm10_2023)

##################################################
# Time Series Analysis
##################################################
# Add a week column to the data frame
imputed_aurn_data <- imputed_aurn_data %>%
  mutate(week = floor_date(date, "week"),year=year(date))

# Aggregate the data by week and calculate the average pollution level
weekly_avg_pollution <- imputed_aurn_data %>%
  group_by(week,site,code,year) %>%
  summarise(no2_mean = mean(imputed_no2, na.rm = TRUE),
            .groups = 'drop')

london_week <- weekly_avg_pollution[weekly_avg_pollution$code=='MY1',]

london_week_ts<- ts(london_week$no2_mean,start = c(2018,12),end=c(2023,12),frequency = 52)

total_length <- nrow(london_week)
train_size <- floor(0.8 * total_length)

# Splitting the data into training and testing sets for mean_no2
no2_data_train <- london_week[1:train_size, ]
no2_data_test <- london_week[(train_size + 1):total_length, ]

# Creating time series objects for training and testing data
no2_ts_train <- ts(no2_data_train$no2_mean,start=c(2018,12),end=c(2022,12),frequency = 52)
no2_ts_test <- ts(no2_data_test$no2_mean,start=c(2022,12),end=c(2023,11),frequency = 52)

# Visualize the time series
plot(london_week_ts, main = "Nitrogen Dioxide Emission", ylab = "NO2 (µg/m³)", xlab = "Year")

#Trend checking
MannKendall(london_week_ts)

# Decompose the time series
no2_decompose<- stl(london_week_ts,s.window = 'periodic')
plot(no2_decompose)


# Stationarity test for original data

Acf(no2_ts_train)
Pacf(no2_ts_train)

adf_result<-adf.test(no2_ts_train)
adf_result

# first and seasonal differencing for removing trend and seasonality
combined_diff <- diff(diff(no2_ts_train, differences = 1), lag = 12)

# Plot the differenced series (to check if both trend and seasonality are removed)
plot(combined_diff, main = "Trend and Seasonality Removed", ylab = "Differenced Values", xlab = "Time")

# Stationarity test for detrended data
adf_result_1<-adf.test(combined_diff)
adf_result_1


#Acf and Pacf for no2
Acf(combined_diff) # for deciding q value
Pacf(combined_diff) # for deciding p value

# Fit ARIMA model to NO2
fit_1 <- auto.arima(combined_diff,seasonal=FALSE)
summary(fit_1)

forecast_arima<-forecast(fit_1, h=length(no2_ts_test))
plot(forecast_arima)

#SARIMA model
#fit<- auto.arima(no2_diff)
fit <- auto.arima(no2_ts_train)
summary(fit)

# Residual Analysis
checkresiduals(fit)
checkresiduals(fit_1)

# Forecasting for NO2
#forecasted_values_1 <- forecast(fit_1, h = length(no2_ts_test))
forecasted_value <- forecast(fit, h = length(no2_ts_test))
print(forecasted_value)

# Plotting the training data, test data, and forecasted values
plot(forecasted_value, main = "Forecast from SARIMA(0,1,2) (1,0,0)[52] for NO2 emission ", ylab = "NO2 (µg/m³)", xlab = "Year", xlim = c(2018, 2023))

# Adding the test data to the plot
lines(no2_ts_test, col = "blue")

# Adding the forecasted values to the plot
lines(forecasted_value$mean, col = "red")

# Adding confidence intervals for the forecast
lines(forecasted_value$lower[,2], col = "red", lty = 2)
lines(forecasted_value$upper[,2], col = "red", lty = 2)

# Adding a legend to the plot
legend("bottomleft", legend = c("Training Data", "Test Data", "Forecasted Values", "95% CI"), col = c("black", "blue", "red", "red"), lty = c(1, 1, 1,2),cex = 0.6)

#Model Evaluation
accuracy_metrics <- accuracy(forecasted_value, no2_ts_test)
print(accuracy_metrics)
