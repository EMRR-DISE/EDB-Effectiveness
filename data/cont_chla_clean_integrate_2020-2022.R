# Emergency Drought Barrier - continuous chlorophyll data
# Clean, integrate, and calculate the daily averages and medians of the
  # continuous chlorophyll data from 2020-2022 for the stations within the
  # designated EDB regions. Data is used in figures for the 2022 EDB monitoring
  # report.
# Author: Dave Bosworth
# Contact: David.Bosworth@water.ca.gov

# Load packages
library(tidyverse)
library(readxl)
library(dataRetrieval)
library(sf)
library(here)
library(qs)
library(conflicted)

# Declare package conflict preferences 
conflicts_prefer(dplyr::filter())

# Check if we are in the correct working directory
i_am("data/cont_chla_clean_integrate_2020-2022.R")

# Set System Timezone as "Etc/GMT+8" (PST) to make it consistent with all data frames
Sys.setenv(TZ = "Etc/GMT+8")


# 1. Import Data ----------------------------------------------------------

# All 15-minute continuous chlorophyll fluorescence data collected by DWR was
  # either downloaded directly from the Water Data Library
  # (https://wdl.water.ca.gov/WaterDataLibrary/) or acquired from direct data
  # requests. Here are more specifics for each station:

  # BLP: https://wdl.water.ca.gov/WaterDataLibrary/StationDetails.aspx?Station=B9502900&source=map
  # FAL: https://wdl.water.ca.gov/WaterDataLibrary/StationDetails.aspx?Station=B9504400&source=map
  # FRK: Direct data request from DWR-DISE-CEMP. Data from the Water Quality Portal.
  # HLT: https://wdl.water.ca.gov/WaterDataLibrary/StationDetails.aspx?Station=B9545800&source=map
  # OSJ: https://wdl.water.ca.gov/WaterDataLibrary/StationDetails.aspx?Station=B9510800&source=map
  # RVB: Direct data request from DWR-DISE-CEMP. Data from the Water Quality Portal.
  # SSI: Direct data request from DWR-DISE-CEMP. Data from the Water Quality Portal.
  # TWI: Direct data request from DWR-DISE-CEMP. Data from the Water Quality Portal.

# Create a vector of all file paths for zip files containing the 15-minute
  # continuous chlorophyll data collected by DWR
fp_chla_zip <- dir(
  here("data-raw"), 
  pattern = "Cont_chla_202[0-2]\\.zip$", 
  full.names = TRUE
)

# Unzip folders containing 15-minute continuous chlorophyll data to a temporary
  # directory
walk(fp_chla_zip, ~ unzip(zipfile = .x, exdir = tempdir()))

# Create a vector of all file paths for the continuous chlorophyll data in the
  # temporary directory - csv files only
fp_chla_raw <- dir(tempdir(), full.names = TRUE, recursive = TRUE) %>% 
  str_subset("Cont_chla_202[0-2].+\\.csv$")

# Create a string to use as a regular expression pattern for the easy to import
  # files in fp_chla_raw
regexp_easy <- "BLP_Chlor_2021-01-01|FRK|RVB|SSI|TWI|HLT_Chlor_2021-01-01|OSJ_Chlor_2021-01-01"

# Import easy to import files
ls_chla_easy <- map(
  str_subset(fp_chla_raw, regexp_easy), 
  ~ read_csv(.x, col_types = cols(.default = "c"))
)

# Function to extract the three letter station codes from their file paths
extr_sta_code <- function(fp_str) {
  str_extract(fp_str, "(?<=Cont_chla_202[0-2]/)[:upper:]{3}")
}

# Subset file paths for for FAL, HLT, and OSJ collected in 2020
fp_chla_raw_fal_hlt_osj2020 <- str_subset(fp_chla_raw, "FAL_Chlor_2020|HLT_Chlor_2020|OSJ_Chlor_2020")

# Import continuous chlorophyll data collected in 2020 for FAL, HLT, and OSJ as
  # a nested dataframe
ndf_chla_fal_hlt_osj2020 <- tibble(
  Station = map_chr(fp_chla_raw_fal_hlt_osj2020, extr_sta_code),
  df_data = map(
    fp_chla_raw_fal_hlt_osj2020,
    ~ read_csv(.x, skip = 8, col_types = cols(.default = "c"))
  )
)

# Create a vector of the remaining csv file paths to be imported
fp_chla_raw_remain <- str_subset(fp_chla_raw, regexp_easy, negate = TRUE) %>% 
  str_subset("FAL_Chlor_2020|HLT_Chlor_2020|OSJ_Chlor_2020", negate = TRUE)

# Import continuous chlorophyll data for the remaining file paths as a nested
  # dataframe
ndf_chla_remain <- tibble(
  Station = map_chr(fp_chla_raw_remain, extr_sta_code),
  df_data = map(
    fp_chla_raw_remain,
    ~ read_csv(.x, skip = 1, col_types = cols(.default = "c"))
  )
)

# Import the one xlsx file containing continuous chlorophyll data
df_chla_rvb2022 <- read_excel(
  dir(tempdir(), pattern = "RVB Combined 2022\\.xlsx", full.names = TRUE, recursive = TRUE),
  na = c("", "NA")
)

# We are also using 15-minute continuous chlorophyll fluorescence data collected
  # by USGS at a few stations (MDM, SJJ). For this data, we will download and save
  # local copies of the data using the `dataRetrieval` package since some of the
  # data is provisional and may change.

# Set download_usgs to TRUE if need to download USGS continuous chlorophyll data
download_usgs <- FALSE

# Download USGS continuous chlorophyll data if necessary
if (download_usgs == TRUE) {
  # Create vectors for parameters, start and end dates
  start_date <- "2020-01-01"
  end_date <- "2022-12-31"
  params <- "32316"  # Chlorophyll concentration estimated from reference material (ug/L)
  
  # Download data for each station individually since doing it all at once
    # doesn't work correctly
  MDM2 <- as_tibble(readNWISuv("11312676", params, start_date, end_date, tz = "Etc/GMT+8"))
  SJJ <- as_tibble(readNWISuv("11337190", params, start_date, end_date, tz = "Etc/GMT+8"))
  
  # Export 15-minute data as a compressed qs file
  qsavem(MDM, SJJ, file = here("data-raw/USGS_chla_2020-2022.qs"))
  
  # Clean up
  rm(start_date, end_date, params, MDM, SJJ)
}

# Import continuous chlorophyll data collected by USGS
qload(here("data-raw/USGS_chla_2020-2022.qs"))

# Import coordinates for stations
df_coords <- read_excel(here("data/Cont_chla_station_coord.xlsx"))

# Load polygon shapefile for the EDB regions
sf_edb_reg <- read_sf(here("Spatial_data/EDB_Regions.shp")) %>%
  select(Region = Regions)


# 2. Clean and Integrate Data ---------------------------------------------

# Clean and combine data from easy to import files in ls_chla_easy
df_chla_easy <- ls_chla_easy %>% 
  # Standardize variable names
  map(
    ~ select(
      .x, 
      Station = contains("station"),
      DateTime = contains("time"),
      Chla = value,
      Qual = contains(c("Quality", "qaqc"))
    ) 
  ) %>% 
  # Combine dataframes now that all variable names are standardized
  bind_rows() %>%
  # Convert variable types making them new variables to check for parsing errors
  mutate(
    DateTime2 = parse_date_time(DateTime, orders = c("Ymd T", "mdY R"), tz = "Etc/GMT+8"),
    Chla2 = as.numeric(Chla)
  )

# Check for parsing errors in DateTime
anyNA(df_chla_easy$DateTime)
anyNA(df_chla_easy$DateTime2)
# No parsing errors identified

# Clean and combine data collected in 2020 for FAL, HLT, and OSJ
df_chla_fal_hlt_osj2020 <- ndf_chla_fal_hlt_osj2020 %>% 
  unnest(df_data) %>% 
  rename(
    DateTime = `Date Time`,
    Chla = `Chlorophyll Raw  (UT)`,
    Qual = `Quality Code`
  ) %>% 
  # Convert variable types making them new variables to check for parsing errors
  mutate(
    DateTime2 = mdy_hms(DateTime, tz = "Etc/GMT+8"),
    Chla2 = as.numeric(Chla)
  )

# Check for parsing errors in DateTime
anyNA(df_chla_fal_hlt_osj2020$DateTime)
anyNA(df_chla_fal_hlt_osj2020$DateTime2)
# No parsing errors identified

# Clean and combine data from the remaining file paths
df_chla_remain <- ndf_chla_remain %>% 
  # Standardize variable names
  mutate(
    df_data = map(
      df_data,
      ~ select(
        .x, 
        DateTime = and,
        Chla = contains("7004"),
        Qual = last_col()
      ) %>% 
      # remove first row
      slice(-1)
    )
  ) %>% 
  # Combine dataframes now that all variable names are standardized
  unnest(df_data) %>%
  # Convert variable types making them new variables to check for parsing errors
  mutate(
    DateTime2 = parse_date_time(DateTime, orders = c("mdY T", "mdY R"), tz = "Etc/GMT+8"),
    Chla2 = as.numeric(Chla)
  )

# Check for parsing errors in DateTime
anyNA(df_chla_remain$DateTime)
anyNA(df_chla_remain$DateTime2)
# No parsing errors identified

# Combine all data collected by DWR from csv files
df_chla_dwr_c1 <- bind_rows(df_chla_easy, df_chla_fal_hlt_osj2020, df_chla_remain)

# Check for parsing errors in Chla in the DWR data
identical(filter(df_chla_dwr_c1, is.na(Chla)), filter(df_chla_dwr_c1, is.na(Chla2)))
# No parsing errors identified

# Finish with general formatting and cleaning of df_chla_dwr_c1
df_chla_dwr_c2 <- df_chla_dwr_c1 %>% 
  # Clean up testing variables
  select(-c(DateTime, Chla)) %>%
  rename(
    DateTime = DateTime2,
    Chla = Chla2
  ) %>%
  # Remove NA Chla values and values less than zero
  filter(!is.na(Chla), Chla >= 0) %>%
  # Remove duplicates
  distinct()

# Inspect Qual codes in df_chla_dwr_c2
unique(df_chla_dwr_c2$Qual)
anyNA(df_chla_dwr_c2$Qual)

# Remove Chla records with Qual code of "X" (Bad data)
df_chla_dwr_c3 <- df_chla_dwr_c2 %>% 
  filter(Qual != "X") %>%
  # Remove Qual variable since all data is rated as "Good data" now
  select(-Qual)

# Clean the data set from the one xlsx file
df_chla_rvb2022_c <- df_chla_rvb2022 %>% 
  transmute(
    Station = "RVB",
    DateTime = force_tz(`Date and Time`, tzone = "Etc/GMT+8"),
    Chla = Chlorph
  ) %>% 
  drop_na(Chla)

# Clean and combine data collected by USGS
df_chla_usgs_c1 <- lst(MDM, SJJ) %>% 
  bind_rows(.id = "Station") %>% 
  select(
    Station,
    DateTime = dateTime,
    Chla =  X_32316_00000,
    Qual = X_32316_00000_cd
  )

# Inspect Qual codes in df_chla_usgs_c1
unique(df_chla_usgs_c1$Qual)
anyNA(df_chla_usgs_c1$Qual)

# Remove Qual variable in df_chla_usgs_c1 since it all appears to be "Good data"
df_chla_usgs_c2 <- df_chla_usgs_c1 %>% select(-Qual)

# Combine continuous chlorophyll data collected by DWR and USGS, and the 2022
  # RVB data and only include years 2020-2022, and months Apr-Dec
df_chla_all <- bind_rows(df_chla_dwr_c3, df_chla_rvb2022_c, df_chla_usgs_c2) %>% 
  mutate(Year = year(DateTime), .after = Station) %>% 
  filter(
    Year %in% 2020:2022,
    month(DateTime) %in% 4:12
  ) %>% 
  # Remove any duplicates
  distinct()

# Run a few quality checks on the 15-minute data before calculating daily
  # averages and medians

# Look for duplicated time stamps
df_chla_all %>%
  mutate(DateTime = round_date(DateTime, unit = "15 minute")) %>%
  count(Station, DateTime) %>%
  filter(n > 1)
# No duplicated time stamps present in data set

# Look at min and max values for each station
qc_min_max <- df_chla_all %>%
  group_by(Station, Year) %>%
  summarize(
    min_value = min(Chla),
    max_value = max(Chla)
  ) %>%
  ungroup()

View(qc_min_max)
# 1 station has some values equal to zero
# 2 stations have some values greater than 100
# A few of the stations have noisy data at times
# I tried a few outlier screening tests, but didn't perform as well as I preferred -
  # may use daily medians to visualize trends in continuous chlorophyll data
  # instead of daily means

# Look at sampling coverage to be sure we have all the data
df_chla_all %>% 
  mutate(Month = fct_drop(month(DateTime, label = TRUE))) %>% 
  count(Year, Station, Month) %>%
  complete(Year, Station, Month) %>% 
  pivot_wider(names_from = Month, values_from = n) %>% 
  print(n = 30)
# It looks like we have all the available data


# 3. Aggregate Values -----------------------------------------------------

# Assign EDB regions to each station
df_region_assign <- df_coords %>%
  select(
    Station = `CDEC code`,
    Latitude,
    Longitude
  ) %>%
  # Convert to sf object
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_join(st_make_valid(sf_edb_reg), join = st_intersects) %>%
  # Drop sf geometry column since it's no longer needed
  st_drop_geometry() %>%
  # Rename regions and assign "Outside" to stations without region assignments
  mutate(
    Region = case_match(
      Region,
      "Sacramento" ~ "Sacramento River",
      "San Joaquin" ~ "San Joaquin River",
      "Central Delta" ~ "Interior Delta",
      NA ~ "Outside"
    ),
    # Reassign FAL to the Interior Delta region
    Region = if_else(Station == "FAL", "Interior Delta", Region)
  )

# Calculate daily means and medians of continuous chlorophyll data for each station
df_chla_all_dv <- df_chla_all %>%
  mutate(Date = date(DateTime)) %>%
  summarize(
    AvgChla = mean(Chla),
    MedianChla = median(Chla),
    .by = c(Station, Year, Date)
  ) %>%
  # Add regions
  left_join(df_region_assign, join_by(Station))

# Check for NA values
df_chla_all_dv %>% summarize(across(everything(), ~ sum(is.na(.x))))
# No NA values present

# Finish cleaning the data set of daily values
cont_chla_daily_2022 <- df_chla_all_dv %>% 
  # Remove data for stations outside of the EDB regions
  filter(Region != "Outside") %>%
  # Reorder columns
  select(
    Station,
    Region,
    Year,
    Date,
    AvgChla,
    MedianChla
  ) %>% 
  arrange(Station, Year)


# 4. Save and Export Data -------------------------------------------------

# Save final data set containing continuous chlorophyll data for the EDB
  # analysis as csv file for easier diffing
cont_chla_daily_2022 %>% write_csv(here("data/cont_chla_daily_2020-2022.csv"))

# Save final data set containing continuous chlorophyll data for the EDB
  # analysis as an rds object in the GitHub repo for easier import
saveRDS(cont_chla_daily_2022, file = here("data/cont_chla_daily_2020-2022.rds"))

