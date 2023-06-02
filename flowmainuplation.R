#Data manipulation
library(tidyverse)
library(readxl)

cont = read_csv("ContinuousWQ_EDB_Eff_2021-2022_WQES.csv")
cont = mutate(cont, Date = mdy(Date))
ggplot(cont, aes(x = Date, y = AvgChla))+ geom_line()+
  facet_wrap(~Station)

ggplot(cont, aes(x = Date, y = Temperature))+ geom_line()+
  facet_wrap(~Station)

ggplot(cont, aes(x = Date, y = EC))+ geom_line()+
  facet_wrap(~Station)

dis1 = read_csv("data/2021EDBEffectivenessReport_EMPDiscreteWQ.csv") %>%
  mutate(Date = mdy(Date), Time = time(Time), DissAmmonia = as.numeric(DissAmmonia),
         DissBromide = as.numeric(DissBromide), DON = as.numeric(DON))
dis2 = read_csv("data/2021EDBEffectivenessReport_WQESDiscreteWQ.csv") %>%
  mutate(Date = mdy(Date)) %>%
  rename(LongStationName = `Long Station_Name`, ShortStationName = `Short_Station_Name`)
  

dis3 = read_excel("data/2022EDBEffectivenessReport_EMPDiscreteWQ.xlsx") %>%
  mutate(Time = time(Time)) %>%
  rename(StationCode = `Station Code`)
dis4 = read_excel("data/2022EDBEffectivenessReport_WQESDiscreteWQ.xlsx") %>%
  rename(LongStationName = `Long Station Name`, ShortStationName = `Short Station Name`, Station_Code = `Station Code`,
         Station_Number = `Station Number`, Sample_Code = `Sample Code`, Rpt_Limit = `Rpt Limit`)
str(dis1)
str(dis2)
str(dis3)
str(dis4)

EMP = bind_rows(dis1, dis3) %>%
  mutate(Program = "EMP")
WQES = bind_rows(dis2, dis4)

WQESwide = pivot_wider(WQES, id_cols = c(ShortStationName, Station_Code, Date, Month, Year), names_from = Analyte, values_from = Result,
                       values_fn = mean) %>%
  rename(Chla = `Chlorophyll a`, Pheophytin = `Pheophytin a`, DissAmmonia = `Dissolved Ammonia`, DissBromide = `Dissolved Bromide`,
         DissChloride = `Dissolved Chloride`, SpecificConductance = `Field Specific Conductance`, Temperature = `Field Water Temperature`,
         DO = `Field Dissolved Oxygen`, Turbidity = `Field Turbidity`, pH = `Field pH`, TKN = `Total Kjeldahl Nitrogen`,
         DissTKN = `Dissolved Total Kjeldahl Nitrogen`, TotPhos = `Total Phosphorus`, DissOrthophos = `Dissolved ortho-Phosphate`,
         DissNitrateNitrite = `Dissolved Nitrate + Nitrite`, DON = `Dissolved Organic Nitrogen`, DOC = `Dissolved Organic Carbon`,
         TOC = `Total Organic Carbon`, TSS = `Total Suspended Solids`, StationCode = Station_Code, Station = ShortStationName) %>%
  dplyr::select(-`*No Lab Analyses (Field Measures Only)`) %>%
  mutate(Program = "WQES")

wQ = bind_rows(EMP, WQESwide)
write.csv(wQ, "DiscreteWQ.csv")

#Data manipulation
#Flow data
DSJ = read_csv("data/DSJ_Hydrodynamics_Data_20211230_20230102.csv", 
               skip = 10, col_names = c("Time", "GageHeight", "GageHeight_QC", "Vmean", "Vmean_QC", "Discharge", "Discharge_QC", "NetFlow")) %>%
  mutate(Station = "DSJ")

FAL = read_csv("data/FAL_Hydrodynamics_Data_20211230_20230102.csv", 
               skip = 10, col_names = c("Time", "GageHeight", "GageHeight_QC", "Vmean", "Vmean_QC", "Discharge", "Discharge_QC", "NetFlow"))%>%
  mutate(Station = "FAL")

FCT = read_csv("data/FCT_Hydrodynamics_Data_20211230_20230102.csv", 
               skip = 12, col_names = c("Time", "GageHeight", "GageHeight_QC", "Vmean", "Vmean_QC", "Discharge", "Discharge_QC", "NetFlow")) %>%
  mutate(GageHeight_QC = as.character(GageHeight_QC), Vmean_QC = as.character(Vmean_QC), 
         Discharge_QC = as.character(Discharge_QC), Station = "FCT")

HOL = read_csv("data/HOL_Hydrodynamics_Data_20211230_20230102.csv", 
               skip = 10, col_names = c("Time","GageHeight", "GageHeight_QC", "Vmean", "Vmean_QC", "Discharge", "Discharge_QC", "NetFlow"))%>%
  mutate(Station = "HOL")

ORQ = read_csv("data/ORQ_Hydrodynamics_Data_20211230_20230102.csv", 
               skip = 10, col_names = c("Time","GageHeight", "GageHeight_QC", "Vmean", "Vmean_QC", "Discharge", "Discharge_QC", "NetFlow"))%>%
  mutate(Station = "ORQ")

OSJ = read_csv("data/OSJ_Hydrodynamics_Data_20211230_20230102.csv", 
               skip = 10, col_names = c("Time", "GageHeight", "GageHeight_QC", "Vmean", "Vmean_QC", "Discharge", "Discharge_QC", "NetFlow"))%>%
  mutate(Station = "OSJ")

oldflow = read_csv("Hydrodynamics_Data_2015_2020_2021.csv")

Flow = bind_rows(DSJ, FAL, FCT, HOL, ORQ, OSJ) %>%
  mutate(DateTime = mdy_hm(Time)) %>%
  select(DateTime, GageHeight, GageHeight_QC, Vmean, Vmean_QC, Discharge, Discharge_QC, NetFlow) %>%
  filter(!is.na(DateTime), DateTime > max(oldflow$DateTime))
str(Flow)

Flow2 = bind_rows(oldflow, Flow)

write.csv(Flow2, "Hydrodynamics_Data_2015_20_21_22.csv", row.names = F)
