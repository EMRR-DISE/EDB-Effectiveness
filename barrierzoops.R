#Zooplankton barrier analyssi

library(zooper)
library(tidyverse)
library(lubridate)
library(sf)
library(readxl)
library(RColorBrewer)
library(lme4)
library(lmerTest)
library(emmeans)
library(multcomp)
library(glmmTMB)
library(DHARMa)
library(vegan)

mypal = c(brewer.pal(8, "Dark2"), brewer.pal(8, "Set2"), "black", "white")

#pull zooplankton 2014-2021 from zooper

MyZoops <- Zoopsynther(Data_type = "Community",
                       Sources = c("EMP", "20mm","STN", "FMWT"),
                       Size_class = "Meso",
                       Date_range = c("2014-01-01", "2021-12-30"))

#subset just regions surrounding the Barrier. 
load("data/barrierregions.RData")

MyZoopssf = filter(MyZoops, !is.na(Latitude)) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, remove = F)

regions = st_transform(regions, crs = 4326)

MyZoopssf = st_join(MyZoopssf, regions) %>%
  st_drop_geometry() %>%
  filter(!is.na(Regions))

#group taxonomic names into larger categories
crosswalk = read_excel("data/zoopcrosswalk.xlsx") %>%
  distinct()

MyZoops = mutate(MyZoopssf, Taxname = str_remove(Taxname, "_UnID")) %>%
  left_join(crosswalk) %>%
  filter(!Undersampled)

#export .csv for the data package
write.csv(MyZoops, "BarrierZoops.csv")

#stacked bar plot by region and year

MyZoopssub = group_by(MyZoops, Source, SampleID, Regions, Year, Date, Station, Analy2) %>%
  summarize(CPUE = sum(CPUE))

samps = group_by(MyZoops, Regions, Year) %>%
  summarize(N = length(unique(SampleID)))

ggplot(MyZoopssub, aes(x = Regions, y = CPUE, fill = Analy2))+ geom_col()+
  geom_text(data = samps, aes(x = Regions, y = 2000000, label = N), inherit.aes = FALSE)+
  scale_fill_manual(values = mypal)+
  facet_wrap(~Year)



write.csv(MyZoopssub, "barrierzoops.csv", row.names = FALSE)


ggplot(MyZoopssub, aes(x = Regions, y = CPUE, fill = Analy2))+ geom_col()+
  geom_text(data = samps, aes(x = Regions, y = 900000, label = N), inherit.aes = FALSE)+
  facet_wrap(~Year)+ theme_bw()

#Need averages rather than totals
zoopsave = group_by(MyZoopssub, Regions, Year, Analy2) %>%
  summarize(CPUEm = mean(CPUE), sdzoop = sd(CPUE), se = sdzoop/n())

ggplot(zoopsave, aes(x = Regions, y = CPUEm, fill = Analy2))+ geom_col()+
  geom_text(data = samps, aes(x = Regions, y = 1000, label = N), inherit.aes = FALSE)+
  scale_fill_manual(values = mypal, name = "Taxon")+
  facet_wrap(~Year)+ theme_bw()


ggplot(MyZoopssub, aes(x = Station, y = CPUE, fill = Analy2))+ geom_col()+
  #geom_text(data = samps, aes(x = Regions, y = 1000, label = N), inherit.aes = FALSE)+
  scale_fill_manual(values = mypal, name = "Taxon")+
  facet_grid(Regions~Year, scales = "free_x")+ theme_bw()
#maybe we should just do June-present, since we don't have 2022

summerzoop = filter(MyZoopssub, month(Date)>5)

Szoopsave = group_by(summerzoop, Regions, Year, Analy2) %>%
  summarize(CPUEm = mean(CPUE), sdzoop = sd(CPUE), se = sdzoop/n())%>%
  mutate(Yearf = as.factor(Year))

ggplot(Szoopsave, aes(x = Yearf, y = CPUEm, fill = Analy2))+ geom_col()+
  #geom_text(data = samps, aes(x = Regions, y = 1000, label = N), inherit.aes = FALSE)+
  scale_fill_manual(values = mypal, name = "Taxon")+
  facet_wrap(~Regions)+ theme_bw()+ xlab(NULL)+ylab("Zooplankton CPUE")+
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = .5))

ggsave("plots/BarrierZoopCom.tiff", device = "tiff", width = 8, height = 7)

###############################################
#now total zooplankton catch, with error bars.
yeartypes = read_csv("data/yearassignments.csv")

EMPtots = group_by(MyZoopssub, Regions, Year, SampleID, Station) %>%
  summarize(CPUE = sum(CPUE), rCPUE = round(CPUE),Yearf = as.factor(Year)) %>%
  distinct() %>%
  left_join(yeartypes)


totsave = group_by(EMPtots, Regions, Year, ShortTerm) %>%
summarize(CPUEm = mean(CPUE), sdzoop = sd(CPUE), se = sdzoop/sqrt(n())) %>%
  mutate(Yearf = as.factor(Year))


#What's the distribution like?
ggplot(EMPtots, aes(x = CPUE)) + geom_histogram()
ggplot(EMPtots, aes(x = log(CPUE))) + geom_histogram()

#linear model based on year, region, and interaction, with station as random effect

zoop1 = lmer(log(CPUE+1) ~ Yearf*Regions  + (1|Station), data = EMPtots)

summary(zoop1)
library(car)
Anova(zoop1, type = "III")
plot(zoop1)
emmeans(zoop1, pairwise ~ Yearf|Regions)
plot(emmeans(zoop1, pairwise ~ Yearf|Regions), comparisons = T)

tuk = emmeans(zoop1, pairwise ~ Yearf|Regions)



#OK, now try it with a negative binomial instead
zoop2 = glmmTMB(rCPUE~ Yearf*Regions  + (1|Station), family = nbinom2(), data = EMPtots)
summary(zoop2)
#I like that better! Looks like what I was expecting

res = simulateResiduals(zoop2)
plot(res)

#not great
testOutliers(res, type = "bootstrap")
emmeans(zoop2, pairwise ~ Yearf|Regions)

#Beautiful. Sacramenot is lower, 2017 is differnet in the Sacramento
#no clear impact of the Barrier

library(multcomp)
mod_means <- multcomp::cld(object = tuk$emmeans,
                           Letters = letters)

ggplot(totsave, aes(x = Yearf, y = CPUEm, fill = ShortTerm)) + geom_col()+
  geom_errorbar(aes(ymin = CPUEm - se, ymax = CPUEm + se))+
  scale_fill_manual(name = "Year Type", values = c("darkorange", "tan", "skyblue"))+
  geom_text(data = mod_means, aes(x = Yearf, group = Regions, y = -100, label = .group), inherit.aes = FALSE)+
  facet_wrap(~Regions)+ theme_bw()+ ylab("Mean Zooplankton CPUE") + xlab(NULL)+
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = .5))



ggsave("plots/BarrierZooptot.tiff", device = "tiff", width = 8, height = 5)


#####################################################################
#Multivariate stats

#Create a species matrix
EMPmat = pivot_wider(MyZoopssub, id_cols = c(Regions, Year, SampleID, Station, Date), 
                     names_from = Analy2, values_from = CPUE, values_fill = 0) %>%
  left_join(EMPtots) %>%
  filter(CPUE !=0) 

EMPenv = dplyr::select(EMPmat, Regions, Year, SampleID, Station, Date) %>%
  mutate(Yearf = as.factor(Year), Regions  = as.factor(Regions))

EMPspecs = dplyr::select(ungroup(EMPmat), -Regions, -Year, -SampleID, -Station, -Date, -CPUE, -rCPUE, -Yearf, -ShortTerm, -SprNDOI,
                         -Drought, -Yr_type, -DroughtYear, -Index, -Whitepaper) 

a1 = adonis2(EMPspecs ~ Regions*Yearf, data = EMPenv)
a1
save(a1, file = "zoopadonis.RData")
#year and regions are significant, but a tiny proportion of the variance. 

#NMDS
metaMDS(EMPspecs, trymax = 200)
#Hm. Not working. Let's try relative abundance

EMPspecs2 = EMPspecs/rowSums(EMPspecs)


adonis2(EMPspecs2 ~ Regions*Yearf, data = EMPenv)
znmds = metaMDS(EMPspecs2, trymax = 200)

source("plotNMDS.R")
PlotNMDS(znmds, EMPenv, group = "Regions")
#that's weird and I dont like it. I think we can say not that much differnece bewteen regions
PlotNMDS(znmds, EMPenv, group = "Yearf")

########################################################
#do the pforbsei copepodites from 1995-2021

copeps <- Zoopsynther(Data_type = "Taxa",
                       Sources = c("EMP","STN", "FMWT"),
                       Size_class = "Meso",
                      Taxa = "Pseudodiaptomus",
                       Date_range = c("1995-01-01", "2021-12-30"))

unique(copeps$Lifestage)
copeps = filter(copeps, Lifestage == "Juvenile")
copepsLSZ = filter(copeps, SalSurf > 0.5, SalSurf < 6) %>%
  mutate(Month = month(Datetime), 
         Season = case_when(Month %in% c(6,7,8) ~ "Summer",
                            Month %in% c(9,10,11) ~ "Fall"),
         Season = factor(Season, levels = c("Summer", "Fall"))) %>%
  filter(Month %in% c(6,7,8,9,10))

ggplot(copepsLSZ, aes(x = as.factor(Year), y = CPUE, fill = Source)) + geom_boxplot()+scale_y_log10()


ggplot(copepsLSZ, aes(x = as.factor(Year), y = (CPUE+10), fill = Season)) + geom_boxplot()+scale_y_log10() + theme_bw()+
  xlab("Year")+ ylab("Individuals+10 per m3 (log 10 scale)")

copsum = group_by(copepsLSZ, Year, Season) %>%
  summarize(med = median(CPUE), meanc = mean(CPUE), max = max(CPUE), min = min(CPUE+1), sdcop = sd(CPUE)/5, N = n())

ggplot(copsum, aes(x = Year, y = meanc, color = Season)) + geom_point(position = position_dodge(width = 1), size = 2)+
  geom_errorbar(aes(ymin = meanc-sdcop, ymax = meanc+sdcop, width = .1, group = Season), position = position_dodge(width = 1))+
  geom_line()+
  scale_y_log10() + theme_bw()+
  xlab("Year")+ ylab("Individuals+10 per m3 (log 10 scale)")

write.csv(table(copeps$Year, copeps$Source), "Samplsize.csv")

write.csv(table(copepsLSZ$Year, copepsLSZ$Source), "LSZ Samplsize.csv")

ggplot(copepsLSZ, aes(x = Year, y = (CPUE+10), color = Season)) + geom_point(alpha = 0.1)+
  geom_smooth()+
  scale_y_log10() + theme_bw()+
  xlab("Year")+ ylab("Individuals+10 per m3 (log 10 scale)")
