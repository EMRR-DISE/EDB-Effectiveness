#Graph of LSZ Pseudodiaptomus copepedites


library(zooper)
library(tidyverse)
library(lubridate)
#do the pforbsei copepodites from 1995-2021

#Grab data from zooper
copeps <- Zoopsynther(Data_type = "Taxa",
                      Sources = c("EMP","STN", "FMWT"),
                      Size_class = "Meso",
                      Taxa = "Pseudodiaptomus",
                      Date_range = c("1995-01-01", "2021-12-30"))

unique(copeps$Lifestage)

#We just want the copepedites
copeps = filter(copeps, Lifestage == "Juvenile")

#filter to low salinity zone and just June-October
copepsLSZ = filter(copeps, SalSurf > 0.5, SalSurf < 6) %>%
  mutate(Month = month(Datetime), 
         Season = case_when(Month %in% c(6,7,8) ~ "Summer",
                            Month %in% c(9,10,11) ~ "Fall"),
         Season = factor(Season, levels = c("Summer", "Fall"))) %>%
  filter(Month %in% c(6,7,8,9,10))

#see what it looks like by survey
ggplot(copepsLSZ, aes(x = as.factor(Year), y = CPUE, fill = Source)) + 
  geom_boxplot()+scale_y_log10()

#Now boxplots by season
ggplot(copepsLSZ, aes(x = as.factor(Year), y = (CPUE+10), fill = Season)) + geom_boxplot()+scale_y_log10() + theme_bw()+
  xlab("Year")+ ylab("Individuals+10 per m3 (log 10 scale)")

#point and line plot
copsum = group_by(copepsLSZ, Year, Season) %>%
  summarize(med = median(CPUE), meanc = mean(CPUE), max = max(CPUE), min = min(CPUE+1), sdcop = sd(CPUE)/5, N = n())

ggplot(copsum, aes(x = Year, y = meanc, color = Season)) + geom_point(position = position_dodge(width = 1), size = 2)+
  geom_errorbar(aes(ymin = meanc-sdcop, ymax = meanc+sdcop, width = .1, group = Season), position = position_dodge(width = 1))+
  geom_line()+
  scale_y_log10() + theme_bw()+
  xlab("Year")+ ylab("Individuals+10 per m3 (log 10 scale)")

#see how many samples were collected per year
write.csv(table(copeps$Year, copeps$Source), "Samplsize.csv")

write.csv(table(copepsLSZ$Year, copepsLSZ$Source), "LSZ Samplsize.csv")

ggplot(copepsLSZ, aes(x = Year, y = (CPUE+10), color = Season)) + geom_point(alpha = 0.1)+
  geom_smooth()+
  scale_y_log10() + theme_bw()+
  xlab("Year")+ ylab("Individuals+10 per m3 (log 10 scale)")
