#TUCP zoops

library(tidyverse)
library(lubridate)
library(zooper)

#X2 from Dayflow
library(DroughtData)
View(raw_hydro_1975_2022)

#Dataframe with just the X2 values and Outflow
X2 = select(raw_hydro_1975_2022, X2, Date, Outflow)


#zooplankton dataset - just 20mm and EMP. If you want summer zooplankton you can include STN and FMWT
zoops = Zoopsynther(Data_type = 'Community', Sources = c("EMP", "20mm"), Size_class = c("Meso", "Macro"))

#Pull out the mysid data
mys = Zoopsynther(Data_type = 'Community', Sources = c("EMP"), Size_class = c("Macro"))

#remove undersampled taxa, filter to just the spring and just the most recent years
mys2 = filter(mys, Undersampled == F, month(Date) %in% c(3,4,5), Year %in% c(2000:2021))

#do the same thing with the other zooplankton
zoops2 = filter(zoops, Undersampled == F, month(Date) %in% c(3,4,5), Year %in% c(2000:2021))

#Which are the most common taxa?
taxa = group_by(zoops2, Taxlifestage) %>%
  summarize(tot = sum(CPUE))
toptax = arrange(taxa, by = -tot)[1:12,]
View(toptax)

#filter out the most common taxa and rename
zoops3 = filter(zoops2, Taxlifestage %in% toptax$Taxlifestage) %>%
  left_join(X2) %>%
  mutate(Taxa = factor(Taxlifestage, levels = sort(unique(Taxlifestage)), 
                       labels = c("Acartia sp. \nadults", "Bosmina logirostris",
                                   "Barnacle Larvae", "Other Cladocera",
                                  "Daphnia sp.", "Eurytemora affinis\nadults", "Eurrytemora affinis\ncopepedites",
                                  "Harpacticoid copepods",
                                  "Pseudodiaptomus forbesi\nadult", "Pseudodiaptomus sp.\ncopepedites", 
                                  "Sinocalanus doerrii\nadult", "Sinocalanus doerrii\ncopepedites")))

#Plot CPUE versus x2
ggplot(zoops3, aes(x = X2, y = log(CPUE+1)))+geom_point(alpha = 0.2, 
                                                        color = "seagreen")+geom_smooth(method = "lm", color = "black")+
  facet_wrap(~Taxa)+theme_bw()

#now filter to just the low salinity zone

zoopsLSZ = filter(zoops3, SalSurf <=6, SalSurf >0.5)


ggplot(zoopsLSZ, aes(x = Outflow, y = log(CPUE+1)))+geom_point(alpha = 0.2, 
                                                        color = "seagreen")+geom_smooth(method = "lm", color = "black")+
  facet_wrap(~Taxa)+theme_bw()+ xlab("Outflow - March-May")

#linear models - one for each taxon
Outmetric = zoopsLSZ  %>%
  mutate(logCPUE = log(CPUE+1))%>%
  group_by(Taxa) %>%
  summarize(intercept = coef(lm(logCPUE ~ log(Outflow)))[1],
            grad = coef(lm(logCPUE ~ log(Outflow)))[2],
            r2 = summary(lm(logCPUE ~ log(Outflow)))$r.squared,
            P =  summary(lm(logCPUE ~ log(Outflow)))$coefficients[2,4],
            Y = max(logCPUE, na.rm = T))

#add statistical significance
Outmetric =  mutate(Outmetric, sig = case_when(P >= 0.05 ~ "NS",
                         P < 0.05 & P > 0.001 ~ "*",
                         P < 0.001 ~ "**"))


#plot with the labels and formuals
ggplot(zoopsLSZ, aes(x = X2, y = log(CPUE+1)))+geom_point(alpha = 0.2, 
                                                          color = "seagreen")+geom_smooth(method = "lm", color = "black")+
  facet_wrap(~Taxa)+theme_bw()+ xlab("X2 - March-May")+
  geom_text(data = Outmetric, aes(x = 40, y = 9, 
                                  label = paste("y = x", round(grad, 3),
                                                "+", round(intercept, 3), "\n R2 =", round(r2, 4),
                                                " P = ", round(P, 4), sep = "")),
            size = 3, nudge_y = -1, inherit.aes = FALSE, hjust = 0, vjust =0)+
geom_text(data = filter(Outmetric, sig != "NS"), aes(x = 80, y = 9, 
                                label = sig),
          size = 5, color = "red", inherit.aes = FALSE, hjust = 0, vjust =1)
  


#maybe annual (spring) averages instead of all the data?
Zoopan = group_by(zoops3, Taxa, Year) %>%
  summarize(CPUE = mean(CPUE), X2 = mean(X2), Outflow = mean(Outflow), logCPUE = log(CPUE+1))

#annual averages for the low salinity zone
ZoopanLSZ = group_by(zoopsLSZ, Taxa, Year) %>%
  summarize(CPUE = mean(CPUE), X2 = mean(X2), logCPUE = log(CPUE+1), Outflow = mean(Outflow))

#plot low salinity zone zoops versus X2
ggplot(ZoopanLSZ, aes(x = X2, y = log(CPUE+1)))+geom_point(alpha = 0.2, 
                                                          color = "seagreen")+geom_smooth(method = "lm", color = "black")+
  facet_wrap(~Taxa)+theme_bw()+ xlab("Mean X2 - March-May")+ylab("log Annual mean CPUE")

#plot overall zoops versus x2
ggplot(Zoopan, aes(x = X2, y = log(CPUE+1)))+geom_point(alpha = 0.2, 
                                                           color = "seagreen")+geom_smooth(method = "lm", color = "black")+
  facet_wrap(~Taxa)+theme_bw()+ xlab("Mean X2 - March-May")+ylab("log Annual mean CPUE")

###############################################################################
#models for annual averages - using outflow instead of X2


#linear models
Outm = ZoopanLSZ  %>%
  mutate(logCPUE = log(CPUE+1))%>%
  group_by(Taxa) %>%
  summarize(intercept = coef(lm(logCPUE ~ log(Outflow)))[1],
            grad = coef(lm(logCPUE ~ log(Outflow)))[2],
            r2 = summary(lm(logCPUE ~ log(Outflow)))$r.squared,
            P =  summary(lm(logCPUE ~ log(Outflow)))$coefficients[2,4],
            Y = max(logCPUE, na.rm = T))



Outm =  mutate(Outm, sig = case_when(P >= 0.05 ~ "NS",
                                               P < 0.05 & P > 0.001 ~ "*",
                                               P < 0.001 ~ "**"))


#plot it
ggplot(ZoopanLSZ, aes(x = log(Outflow), y = log(CPUE+1)))+geom_point(alpha = 0.2, 
                                                          color = "seagreen")+
  geom_smooth(method = "lm", color = "grey", fill = "lightgrey", alpha = 0.5)+
  geom_abline(data = Outm, aes(slope = grad, intercept = intercept, color = sig))+
  scale_color_manual(values = c("blue", "grey"), guide = NULL)+
  facet_wrap(~Taxa)+theme_bw()+ xlab("Log-Transformed Mean Outflow (CFS) - March-May")+ ylab("Log Mean March-May CPUE")+
  geom_text(data = Outm, aes(x = 9, y = 9, 
                                  label = paste("y = x", round(grad, 3),
                                                "+", round(intercept, 3), "\n R2 =", round(r2, 4),
                                                " P = ", round(P, 4), sep = "")),
            size = 3, nudge_y = -1, inherit.aes = FALSE, hjust = 0, vjust =0)+
  geom_text(data = filter(Outm, sig != "NS"), aes(x = 11, y = 9, 
                                                       label = sig),
            size = 5, color = "red", inherit.aes = FALSE, hjust = 0, vjust =1)

#################################################################
#mysids

mys3 = group_by(mys2, SampleID, SalSurf, Station, Date, Year) %>%
  summarize(CPUE = mean(CPUE), logCPUE = log(CPUE+1)) %>%
  left_join(X2)

ggplot(mys3, aes(x = log(Outflow), y = logCPUE))+geom_point()+ geom_smooth()

#annual average
mysan = group_by(mys3,  Year) %>%
  summarize(CPUE = mean(CPUE), logCPUE = log(CPUE+1), X2 = mean(X2), Outflow = mean(Outflow))

ggplot(mysan,aes(x = log(Outflow), y = logCPUE))+geom_point()+ geom_smooth(method = "lm")

#restrict to low salinity zone


mys3LSZ = group_by(mys2, SampleID, SalSurf, Station, Date, Year) %>%
  summarize(CPUE = mean(CPUE), logCPUE = log(CPUE+1)) %>%
  left_join(X2)%>%
  filter(SalSurf <6, SalSurf >0.5)

ggplot(mys3LSZ, aes(x = log(Outflow), y = logCPUE))+geom_point()+ geom_smooth()

#annual average
mysanLSZ = group_by(mys3LSZ,  Year) %>%
  summarize(CPUE = mean(CPUE), logCPUE = log(CPUE+1), X2 = mean(X2), Outflow = mean(Outflow))

ggplot(mysanLSZ,aes(x = log(Outflow), y = logCPUE))+geom_point()+ geom_smooth(method = "lm")

my = lm(logCPUE ~ log(Outflow), data = mysanLSZ)
summary(my)

#################################################################################
#OK, let's do fewer taxa and do both freshwater and LSZ and maybe HSZ

Allzoops = mys3 %>%
  mutate(Taxa = "Mysids") %>%
  rbind(zoops3) %>%
  mutate(Taxa2 = case_when(Taxa %in% c("Acartia sp. \nadults")~ "Acartia",
                           Taxa %in% c("Bosmina logirostris") ~ "Bosmina",
                           Taxa %in% c("Harpacticoid copepods") ~ "Harpacticoids",
                          Taxa %in% c("Eurytemora affinis\nadults", "Eurrytemora affinis\ncopepedites") ~ "Eurytemora",
                          Taxa %in% c("Pseudodiaptomus sp.\ncopepedites",  "Pseudodiaptomus forbesi\nadult") ~ "Pseudodiaptomus",
                          Taxa %in% c("Sinocalanus doerrii\nadult", "Sinocalanus doerrii\ncopepedites") ~ 'Sinocalanus',
                          TRUE ~ Taxa)) %>%
  group_by(SampleID, SalSurf, Station, Date, Year, X2, Outflow, Taxa2) %>%
  summarize(CPUE = sum(CPUE)) %>%
  filter(!is.na(SalSurf), !Taxa2 %in% c("Barnacle Larvae", "Other Cladocera")) %>%
  mutate(Month = month(Date),
         SalCat = case_when(SalSurf <0.5 ~ "Fresh",
                            SalSurf >= 0.5 & SalSurf <= 6 ~ "Low Salinity",
                            SalSurf > 6~ "High Salinity"),
         SalCat = factor(SalCat, levels = c("Fresh", "Low Salinity", "High Salinity"), labels = c("Fresh (<0.5ppt)", "Low Salinity (0.5-6 ppt)", "High Salinity (>6 ppt)")))

#calculate annual means
AllzoopsAn = group_by(Allzoops,Year, SalCat, Taxa2) %>%
  summarise(CPUE = mean(CPUE), X2 = mean(X2), Outflow = mean(Outflow))

ggplot(AllzoopsAn, aes(x = log(Outflow), y = log(CPUE+1)))+geom_point()+geom_smooth(method = "lm")+
  facet_grid(Taxa2~SalCat, scales = "free_y")



#linear models
Outm = AllzoopsAn  %>%
  mutate(logCPUE = log(CPUE+1))%>%
  group_by(Taxa2, SalCat) %>%
  summarize(intercept = coef(lm(logCPUE ~ log(Outflow)))[1],
            grad = coef(lm(logCPUE ~ log(Outflow)))[2],
            r2 = summary(lm(logCPUE ~ log(Outflow)))$r.squared,
            P =  summary(lm(logCPUE ~ log(Outflow)))$coefficients[2,4],
            Y = max(logCPUE, na.rm = T))



Outm =  mutate(Outm, sig = case_when(P >= 0.05 ~ "NS",
                                     P < 0.05 & P > 0.001 ~ "*",
                                     P < 0.001 ~ "**"))


#I'm pretty sure this is what made it into the final TUCP documents
ggplot(AllzoopsAn, aes(x = log(Outflow), y = log(CPUE+1)))+geom_point(alpha = 0.2, 
                                                                     color = "seagreen")+
  geom_smooth(method = "lm", color = "grey", fill = "lightgrey", alpha = 0.5)+
  geom_abline(data = Outm, aes(slope = grad, intercept = intercept, color = sig))+
  scale_color_manual(values = c("blue", "grey"), guide = NULL)+
  facet_grid(Taxa2~SalCat)+theme_bw()+ xlab("Log-Transformed Mean Outflow (CFS) - March-May")+ ylab("Log Mean March-May CPUE")+
  geom_text(data = Outm, aes(x = 9, y = 9, 
                             label = paste("y = x", round(grad, 3),
                                           "+", round(intercept, 3), "\n R2 =", round(r2, 4),
                                           " P = ", round(P, 4), sep = "")),
            size = 3, nudge_y = -1, inherit.aes = FALSE, hjust = 0, vjust =0)+
  geom_text(data = filter(Outm, sig != "NS"), aes(x = 11, y = 9, 
                                                  label = sig),
            size = 5, color = "red", inherit.aes = FALSE, hjust = 0, vjust =1)

#export functions
write.csv(Outm, "ZooplanktonRegressions.csv")

baseOutflow = mean(c(16600, 10800, 7300))
TUCPOutflow = mean(c(12750, 10800, 7250))

log(baseOutflow)
log(TUCPOutflow)

##########################################################################
#quickly look at EMP's visual index

library(discretewq)

EMP = wq(Source = "EMP", Start_year = 2015, End_year = 2021)

EMP = mutate(EMP, Month = month(Date), Microcystis = factor(round(Microcystis), labels = c("Absent", "Low", "Medium", "High"))) %>%
  filter(!is.na(Microcystis))
ggplot(EMP, aes(x = Month, fill = Microcystis))+ geom_bar(position = "fill", color = "grey")+
  scale_fill_manual(values = c("antiquewhite", "yellow", "orange", "red"), name = "Microcystis\nVisual Index")+
  theme_bw()+ylab("Relative Abundace")
