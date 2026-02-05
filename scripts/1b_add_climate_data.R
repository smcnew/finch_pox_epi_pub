# Script to add climate data to compiled bird data

library(dplyr)
library(lubridate)

# Load data
captures <- read.csv("formatted_data/combined_2022_2023_finch.csv") %>%
            mutate(year_month = as.Date(year_month)) %>%
            mutate(year_week = as.Date(year_week)) %>%
            mutate(date = as.Date(date))

climate <- read.csv("input_data/climate_puerto-ayora.csv") %>% rename(date = observation_date)

# Format data and remove columns we don't want
climate <- climate %>%
  mutate(date = as.Date(date)) %>%
  mutate(year = as.numeric(year(date))) %>%
  mutate(month = as.numeric(month(date))) %>%
  mutate(day = day(date)) %>%
  mutate(day_of_year = as.numeric(strftime(date, format = "%j"))) %>%
  mutate(year_month = floor_date(date, unit = "month")) %>%
  mutate(year_week = floor_date(date, unit = "week")) %>%
  filter(year %in% c(2021:2024)) %>%
  mutate(humidity = ifelse(humidity == 0, NA, humidity)) %>%
  select(-c(min_air_temp, max_air_temp, sea_temp, sunshine_hours, clouds))

head(climate)
# Create 2-8 week lags for temperature and humidity
# Lags need to be in days adding up to the right number of weeks

climate <- climate %>%
  mutate(
    templ2 = lag(mean_air_temp, 2*7),
    templ4 = lag(mean_air_temp, 4*7),
    templ6 = lag(mean_air_temp, 6*7),
    templ8 = lag(mean_air_temp, 8*7),
    templ10 = lag(mean_air_temp, 10*7),
    humidityl2 = lag(humidity, 2*7),
    humidityl4 = lag(humidity, 4*7),
    humidityl6 = lag(humidity, 6*7),
    humidityl8 = lag(humidity, 8*7),
    humidityl10 = lag(humidity, 10*7))

# Add weekly total for precipitation

weekly_climate_summaries <- climate %>%
    group_by(year_week) %>%
    summarize(week_precip = sum(precipitation), .groups = "drop") %>%
    mutate(precipl2 = lag(week_precip, 2),
           precipl4 = lag(week_precip, 4),
           precipl6 = lag(week_precip, 6),
           precipl8 = lag(week_precip, 8),
           precipl10 = lag(week_precip, 10)) %>%
    as.data.frame()


# add weekly averages to captures data frame

captures <- left_join(captures, weekly_climate_summaries)



# maybe cut this below...
# add daily temperature and precipitation to DF
captures <- left_join(captures, climate)

# add monthly averages to captures data frame
# captures <- left_join(captures, month_temp_summaries)


head(captures)
# write out combined dataset


write.csv(captures, "formatted_data/combined_2022_2023_finch_climate.csv", row.names = FALSE)


# Deprecated
# # aggregate average temperature, humidity, and total preciptation for each month
# month_temp_summaries <- climate %>%
#   group_by(year_month) %>%
#   summarize(mean_temp = mean(mean_air_temp),
#             mean_humid = mean(humidity),
#             total_precip = sum(precipitation), .groups = "drop") %>%
#   mutate(precip_month_lag = lag(total_precip, 1)) %>%
#   as.data.frame()





# First, aggregate average temperature, humidity, and total precip for each week
# Then, create 2-8 week lags for temp and precipitation
# weekly_climate_summaries <- climate %>%
#   group_by(year_week) %>%
#   summarize(week_temp = mean(mean_air_temp),
#             week_humidity = mean(humidity),
#             week_precip = sum(precipitation), .groups = "drop") %>%
#   mutate(templ2 = lag(week_temp, 2),
#          templ4 = lag(week_temp, 4),
#          templ6 = lag(week_temp, 6),
#          templ8 = lag(week_temp, 8),
#          templ10 = lag(week_temp, 10),
#          humidityl2 = lag(week_humidity, 2),
#          humidityl4 = lag(week_humidity, 4),
#          humidityl6 = lag(week_humidity, 6),
#          humidityl8 = lag(week_humidity, 8),
#          humidityl10 = lag(week_humidity, 10),
#          precipl2 = lag(week_precip, 2),
#          precipl4 = lag(week_precip, 4),
#          precipl6 = lag(week_precip, 6),
#          precipl8 = lag(week_precip, 8),
#          precipl10 = lag(week_precip, 10)) %>%
#   as.data.frame()



# Summary statistics of climate

aggregate(precipitation ~ year, data = climate, sum)
aggregate(mean_air_temp ~ year, data = climate, mean)

# Check for correlation of climatic variables
cor(climate$mean_air_temp, climate$humidity)^2
cor(climate$mean_air_temp, climate$precipitation)^2
cor(climate$humidity, climate$precipitation)^2


