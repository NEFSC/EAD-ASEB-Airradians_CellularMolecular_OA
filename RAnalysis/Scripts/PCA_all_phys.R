---
title: "PCA_All_Phys"
author: "Samuel Gurr"
date: "2023-08-28"
output: html_notebook
---

-   Last updates: September 13, 2023

# OA Reexposure Challenge

### May 2023: F2 Adult Bay scallops exposed to a two-week full-reciprocal pCO2 challenge

### Objective: synthesize the hemolymph flow cytometry and gill tissue lysate data 
at the individual-level resolution for principle component analysis 

## Load Libraries

```{r setup, include=TRUE}

# library(lmtest) # to receive p value from betareg model
# library(FSA) # for the Dun test post hoc for SRH non-parametric 2 way anova]
# library(emmeans)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(car)
# library(lmerTest)
library(tidyr)
# library(reshape2)
library(ggpubr)
# library(nlme)
# library(rcompanion) # to run the Schrier -Ray-Hare non parametric 2 way 
library(ggpmisc) # stat_poly for inserting equation and R2 for ggplot line 
knitr::opts_knit$set(root.dir = "C:/Users/samjg/Documents/Github_repositories/Airradians_CellularMolecular_OA/RAnalysis") # sets the working directory for the entire R markdown file - no need to reload the wd
```

## Load colorimetric datasets

* ATP assay kits
```{r load ATP assay data}
# 

# contains ATP and total protein 
ATP_master      <- read.csv(file="Output/Colorimetric_assays/ATP/ATP_Master.csv",
                          sep = ",",
                          header=TRUE)


```


* Load BCA Total Protein assay kits
```{r load BCA assay data}
# Ran on 20230829 - Plate 1 
BCA_0829_plate1      <- read.csv(file="Data/Colorimetric_assays/BCA_ATPcorrection/Runs_20230829/Plate_1/20230829_BCA_ATPcorrection_Plate1.csv",
                          sep = ",",
                          skip=2,
                          header=TRUE,
                          fileEncoding="latin1")
BCA_0829_plate1_mtx  <- as.matrix(BCA_0829_plate1[c(1:8),c(3:14)])
colnames(BCA_0829_plate1_mtx) = c("1","2","3","4","5","6","7","8","9","10","11","12")
rownames(BCA_0829_plate1_mtx) = c("A","B","C","D","E","F","G","H")
BCA_0829_plate1_table <- as.data.frame.table(BCA_0829_plate1_mtx, responseName = "value") %>% 
                            dplyr::rename(well_row=Var1, well_column=Var2, Abs_562nm=value) %>% 
                            dplyr::mutate(well=paste0(well_row,well_column),
                                          Run_date ="20230829",
                                          Plate="1") %>% 
                            dplyr::select(-c(well_row,well_column))

# Ran on 20230829 - Plate 2
BCA_0829_plate2      <- read.csv(file="Data/Colorimetric_assays/BCA_ATPcorrection/Runs_20230829/Plate_2/20230829_BCA_ATPcorrection_Plate2.csv",
                          sep = ",",
                          skip=2,
                          header=TRUE,
                          fileEncoding="latin1")
BCA_0829_plate2_mtx  <- as.matrix(BCA_0829_plate2[c(1:8),c(3:14)])
colnames(BCA_0829_plate2_mtx) = c("1","2","3","4","5","6","7","8","9","10","11","12")
rownames(BCA_0829_plate2_mtx) = c("A","B","C","D","E","F","G","H")
BCA_0829_plate2_table <- as.data.frame.table(BCA_0829_plate2_mtx, responseName = "value") %>% 
                            dplyr::rename(well_row=Var1, well_column=Var2, Abs_562nm=value) %>% 
                            dplyr::mutate(well=paste0(well_row,well_column),
                                          Run_date ="20230829",
                                          Plate="2") %>% 
                            dplyr::select(-c(well_row,well_column))

# load metadata

metadata_BCA <- read.csv(file="Data/Colorimetric_assays/BCA_ATPcorrection/Metadata_BCA_ATPcorrection.csv",
                          sep = ",",
                          header=TRUE)

# merge the rbind of all raw data with the metadat 
raw_BCA <- merge( (rbind(BCA_0829_plate1_table, 
                         BCA_0829_plate2_table)),
                  metadata_BCA) 

# plot the data to see any glaring outliers 
raw_BCA_plot <- raw_BCA %>% 
  dplyr::filter(Type %in% 'Sample') %>% 
  ggplot(aes(y = Abs_562nm, 
             x  = Scallop_ID)) +
        geom_point()
raw_BCA_plot # we see ibe liw outlier, omit the raw data of Abs_562nm < 0.5

# om dataset
raw_BCA_om <- raw_BCA %>%
                dplyr::filter(!(Type =='Sample' & Abs_562nm < 0.5))
# View(raw_BCA_om)
# write csv
write.csv(raw_BCA_om, file = "Data/Colorimetric_assays/BCA_ATPcorrection/Raw_Master_BCA_ATPcorrection.csv")

```


* Load Experiment metadata - merge with the loaded raw colorimetric data

```{r load Experimentl metadata}
Exp_metadata      <- read.csv(file="Data/Experiment_metadata.csv",sep = ",",header=TRUE)
```

# Finalize BCA data 

- Objs: 

  * (1) Merge by the common 'Scallop_ID' and retain all rows in the colorimetric data (all=T)
  - Why? We unfortunately found that that the dewar was not properly unloaded 
  on day 14 of the experiment, below we can reveal replicates that were lost (hence 'NAs' via merge below)
  
  * (2) omit NAs - empty wells (no sample or standard)
  
  * (3) separate standards from samples
  
  * (4) Run standard curve, calculate totalprotein
  
  * (5) Calculate total protien; subtract out background + correct for standard curve +  average
  
```{r Total protein data}

raw_BCA <- read.csv(file = "Data/Colorimetric_assays/BCA_ATPcorrection/Raw_Master_BCA_ATPcorrection.csv", head = T) # the raw_BCA_om file

# * (1) Merge by the common 'Scallop_ID' and retain all rows in the colorimetric data (all=T)
# - Why? We unfortunately found that that the dewar was not properly unloaded 
# on day 14 of the experiment, below we can reveal replicates that were lost (hence 'NAs' via merge below)
raw_BCA_merged    <- merge(raw_BCA, Exp_metadata, by='Scallop_ID',all=TRUE)

# View(raw_BCA_merged)
lost_gill_samples <- raw_BCA_merged %>% dplyr::filter(Abs_562nm %in% NA)
nrow(lost_gill_samples) # 29 -NOTE: NOT ALL of these were lost, some due to low yield of supernatant to run,
# the true value for lost samples 
  


# * (2) omit NAs - empty wells (no sample or standard)
raw_BCA_merged_om <- raw_BCA_merged %>% dplyr::filter(!Type %in% NA)
 
# View(raw_BCA_merged_om)

# * (3) separate standards from samples, format the standards
BCA_standards <- raw_BCA_merged_om %>% 
                    dplyr::filter(grepl('Standard', Type)) %>% # grepl for 'contains' string
                    dplyr::mutate(BCA_ug_mL = 
                                    as.numeric(gsub('.*_','',Type))) %>% # (i.e. 25 in 'Standards_25')
                    dplyr::mutate(Unique_ID = paste0('Plate_',Plate,'_', Type))



# * (4) Run standard curve, calculate totalprotein 
BCA_background_zero <- BCA_standards %>% 
                        dplyr::filter(Type %in% 'Standard_0') %>% 
                        dplyr::select(Unique_ID,Plate, BCA_ug_mL,Abs_562nm) %>% # select columns of interest
                        dplyr::group_by(Unique_ID, Plate, BCA_ug_mL) %>% # group by to get the means
                        dplyr::summarise_each(funs(mean,sd,se=sd(.)/sqrt(n()))) # get all the stats 
# Plate 1, blank to correct by is 0.08355
# Plate 2, blank to correct by is 0.08090

BCA_standards_means <- BCA_standards %>% 
                        #dplyr::filter(!Type %in% 'Standard_0') %>% 
                        dplyr::mutate(Abs_562nm_cor = 
                                     case_when(Plate == 1 ~ (Abs_562nm-0.08355),
                                               Plate == 2 ~ (Abs_562nm-0.08090) ) ) %>% 
                        dplyr::select(Unique_ID, Plate, BCA_ug_mL, Abs_562nm_cor) %>% # select columns 
                        dplyr::group_by(Unique_ID, Plate, BCA_ug_mL) %>% # group by to get the means
                        dplyr::summarise_each(funs(mean,sd,se=sd(.)/sqrt(n()))) # get all the stats 

BCA_stand_plots_quadratic <- BCA_standards_means %>% # QUADRATIC SMOOTH LINE WORKS BEST HERE (MANUFACTURERS INSTRUCTIONS)
                     #dplyr::filter(!BCA_ug_mL %in% 25) %>% # hash me out to test
                     ggplot(aes(y=mean, x=BCA_ug_mL)) + 
                        geom_point() +
                        theme_bw() +
                        labs(y= "Net Absorbance at 562nm", x = "Protein Concentration in ug/mL") +
                        #geom_line() +
                        #stat_poly_line(color='red') +
                        #geom_smooth() +
                        stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
                        stat_poly_eq(parse=T, aes(label = ..eq.label..), formula=y ~ x + I(x^2)) +
                        ggtitle('Quadratic') +
                        #stat_poly_eq(use_label(c("eq", "R2"))) +
                        facet_wrap(~Plate) 
library(polynom)
poly.calc(BCA_standards_means$BCA_ug_mL, BCA_standards_means$mean)
# 0.002129919*x - 7.316278e-06*x^2 + 3.224975e-08*x^3
f <- function(x) {
  return(0.00127*x - 1.01*10^-7*x^2)
}
ggplot(BCA_standards_means, aes(x=BCA_ug_mL, y=mean)) + 
  geom_point(size=5, col='blue') + 
  stat_function(fun = f, size=1.25, alpha=0.4)

# NOTE! I found that some of my discriminants are negative (b^2 - 4ac) so Im going to extrapolate from a linear 
# curve based on the last few standards
BCA_stand_plots_linear <- BCA_standards_means %>% # QUADRATIC SMOOTH LINE WORKS BEST HERE (MANUFACTURERS INSTRUCTIONS)
                     dplyr::filter(!BCA_ug_mL %in% c('0','25','125','250','500')) %>% # hash me out to test
                     ggplot(aes(y=mean, x=BCA_ug_mL)) + 
                        geom_point() +
                        theme_bw() +
                        labs(y= "Net Absorbance at 562nm", x = "Protein Concentration in ug/mL") +
                        #geom_line() +
                        #stat_poly_line(color='red') +
                        #geom_smooth() +
                        stat_smooth(method = "lm", formula = y ~ x, size = 1) +
                        stat_poly_eq(parse=T, aes(label = ..eq.label..), formula=y ~ x) +
                        ggtitle('Linear - for negative desriminant') +
                        #stat_poly_eq(use_label(c("eq", "R2"))) +
                        facet_wrap(~Plate) 


#print
pdf(paste("Output/Colorimetric_assays/BCA_ATPcorrection/Standard_Curve_BCA_ATPcorrection.pdf", sep =''), 
    width=10, 
    height=7)
print(ggarrange(BCA_stand_plots_quadratic,BCA_stand_plots_linear))
dev.off()

# Standard curve, Plate 1 equation y = 0.0347 + 0.00127x - 1.01x10^-7x^2 - need to solve for x!
# Standard curve, Plate 2 equation y = 0.045 + 0.00125x - 1.3x10^-7x^2 - need to solve for x!


# * (5) Calculate Total protein per sample; subtract out background + correct for standard curve +  average

# Again, Plate 1 blank to correct by is 0.08355
# Again, Plate 2 blank to correct by is 0.08090

# Standard curve, Plate 1
a1 <- -1.01*10^-7
b1 <- 0.00127
c1 <- 0.0347
# EQ: (-(b1) + sqrt( (b1^2) - (4*(((a1)-Abs_562nm_cor))*(c1)) ))/(2*a1)

# Standard curve, Plate 2
a2 <- -1.3*10^-7
b2 <- 0.00125
c2 <- 0.045
# EQ: (-(b2) + sqrt( (b2^2) - (4*a2*(c2-Abs_562nm_cor)) ) ) / (2*a2)


# linear equation plate 1 == (Abs_562nm_cor - 0.192)/0.000993
# linear equation plate 2 == (Abs_562nm_cor - 0.224)/0.000911


# IMPORTANT! we used 25 ul of the standards and 25 ul of the unknowns (samples) 
# therefore we can interpret the unknown direct to the the standard curve without having 
# to account for addition factors, fot example, if we used 5 ul unknown (sample) we would have to adjust 
# by multiplying by 5 to reach the standard curve 

V = 0.025 # 25 ul or 0.025 mL

TotalProtein_final <- raw_BCA_merged_om %>% 
                      dplyr::filter(Type %in% 'Sample') %>% # call samples
                      dplyr::mutate(Unique_ID = 
                                      paste0('Plate_',Plate,'_', Scallop_ID)) %>% # unique ID t0 group by
                      dplyr::mutate(Abs_562nm_cor = # correct the raw abs, subtract background
                                   case_when(Plate == 1 ~ (Abs_562nm-0.08355), # for plate 1
                                             Plate == 2 ~ (Abs_562nm-0.08090) ) ) %>% # for plate 2 
                      dplyr::mutate(TotalProtein_ug_mL = 
                                    case_when(
                                      # linear fr neg discrim. - luckily only two values from plate 2
                                      Scallop_ID %in% c(33, 51) ~ 
                                        ((Abs_562nm_cor - 0.224)/0.000911),
                                      # quadratic for Plate 1
                                      Plate == 1 ~ 
                                        ((-(b1) + sqrt( (b1^2) - (4*a1*(c1-Abs_562nm_cor)) ) ) / (2*a1)), 
                                      # quadratic for plate 2
                                      Plate == 2 | Scallop_ID != c(33, 51) ~ 
                                        ((-(b2) + sqrt( (b2^2) - (4*a2*(c2-Abs_562nm_cor)) ) ) / (2*a2)) ),
                                    # ug per mL concentration to ug in 25 ul sample 
                                    TotalProtein_ug = TotalProtein_ug_mL*V) %>% 
                      dplyr::group_by(Day,pCO2_history,pCO2_exposure,
                                      Unique_ID, Plate, Scallop_ID) %>% # group by to get the means
                      dplyr::summarise(mean_TotalProtein_ug = mean(TotalProtein_ug),
                                       sd_TotalProtein_ug   = sd(TotalProtein_ug),
                                       n = n(),
                                       se_TotalProtein_ug  = sd_TotalProtein_ug / sqrt(n))

# View(TotalProtein_final)
nrow(TotalProtein_final) # 62

# Sanity check Lets look at the absorbance vs. totla protein concentration data 

test <- raw_BCA_merged_om %>% 
                      dplyr::filter(Type %in% 'Sample') %>% # call samples
                      dplyr::mutate(Unique_ID = 
                                      paste0('Plate_',Plate,'_', Scallop_ID)) %>% # unique ID t0 group by
                      dplyr::mutate(Abs_562nm_cor = # correct the raw abs, subtract background
                                   case_when(Plate == 1 ~ (Abs_562nm-0.08355), # for plate 1
                                             Plate == 2 ~ (Abs_562nm-0.08090) ) ) %>% # for plate 2 
                      dplyr::mutate(TotalProtein_ug_mL = 
                                    case_when(
                                      # linear fr neg discrim. - luckily only two values from plate 2
                                      Scallop_ID %in% c(33, 51) ~ 
                                        ((Abs_562nm_cor - 0.224)/0.000911),
                                      # quadratic for Plate 1
                                      Plate == 1 ~ 
                                        ((-(b1) + sqrt( (b1^2) - (4*a1*(c1-Abs_562nm_cor)) ) ) / (2*a1)), 
                                      # quadratic for plate 2
                                      Plate == 2 | Scallop_ID != c(33, 51) ~ 
                                        ((-(b2) + sqrt( (b2^2) - (4*a2*(c2-Abs_562nm_cor)) ) ) / (2*a2)) ),
                                    # ug per mL concentration to ug in 25 ul sample 
                                    TotalProtein_ug = TotalProtein_ug_mL*V)


calc_BCA_plot <- test %>% 
  dplyr::filter(Type %in% 'Sample') %>% 
  ggplot(aes(y = Abs_562nm_cor, 
             x  = TotalProtein_ug_mL)) +
        geom_point() + 
  facet_wrap(~Plate)
calc_BCA_plot # we see ibe liw outlier, omit the raw data of Abs_562nm < 0.5

test %>% dplyr::filter(Type %in% 'Sample') %>% dplyr::filter(Plate == 2 & Abs_562nm_cor > 3)
# write csv
write.csv(TotalProtein_final, file = "Output/Colorimetric_assays/BCA_ATPcorrection/Calc_Master_BCA_ATPcorrection.csv")


```


# Finalize ATP data 

- Objs: 

  * (1) Merge by the common 'Scallop_ID' and retain all rows in the colorimetric data (all=T)
  - Why? We unfortunately found that that the dewar was not properly unloaded 
  on day 14 of the experiment, below we can reveal replicates that were lost (hence 'NAs' via merge below)
  
  * (2) omit NAs - empty wells (no sample or standard)
  
  * (3) separate standards from samples
  
  * (4) Run standard curve, calculate totalprotein
  
  * (5) Calculate ATP; subtract out background + correct for standard curve +  average +
  
  
```{r ATP data}

raw_ATP <- read.csv(file = "Data/Colorimetric_assays/ATP/Raw_Master_ATP.csv", head = T)
# * (1) Merge by the common 'Scallop_ID' and retain all rows in the colorimetric data (all=T)
# - Why? We unfortunately found that that the dewar was not properly unloaded 
# on day 14 of the experiment, below we can reveal replicates that were lost (hence 'NAs' via merge below)
raw_ATP_merged    <- merge(raw_ATP, Exp_metadata, by='Scallop_ID', all=TRUE)

lost_gill_samples <- raw_ATP_merged %>% dplyr::filter(Abs_570nm %in% NA)
nrow(lost_gill_samples) # 26 samples lost
  


# * (2) ommit NAs - empty wells (no sample or standard)
raw_ATP_merged_om <- raw_ATP_merged %>% dplyr::filter(!Type %in% NA)
 


# * (3) separate standards from samples, format the standards
ATP_standards <- raw_ATP_merged_om %>% 
                    dplyr::filter(grepl('Standard', Type)) %>% # grepl for 'contains' string
                    dplyr::mutate(ATP_nmol = 
                                    as.numeric(gsub('.*_','',Type))) %>% # (i.e. 25 in 'Standards_25')
                    dplyr::mutate(Unique_ID = paste0(Run_date,'_', Type))



# * (4) Run standard curve, calculate totalprotein 
ATP_background_zero <- ATP_standards %>% 
                        dplyr::filter(Type %in% 'Standard_0') %>% 
                        dplyr::select(Unique_ID, Run_date, ATP_nmol, Abs_570nm) %>% # select columns of interest
                        dplyr::group_by(Unique_ID, Run_date, ATP_nmol) %>% # group by to get the means
                        dplyr::summarise_each(funs(mean,sd,se=sd(.)/sqrt(n()))) # get all the stats 
# Run_date 20230828, blank to correct by is 0.04675
# Run_date 20230829, blank to correct by is 0.04680

ATP_standards_means <- ATP_standards %>% 
                        #dplyr::filter(!Type %in% 'Standard_0') %>% 
                        dplyr::mutate(Abs_570nm_cor = 
                                     case_when(Run_date == 20230828 ~ (Abs_570nm-0.04675),
                                               Run_date == 20230829 ~ (Abs_570nm-0.04680) ) ) %>% 
                        dplyr::select(Unique_ID, Run_date, ATP_nmol, Abs_570nm_cor) %>% # select columns 
                        dplyr::group_by(Unique_ID, Run_date, ATP_nmol) %>% # group by to get the means
                        dplyr::summarise_each(funs(mean,sd,se=sd(.)/sqrt(n()))) # get all the stats 

ATP_stand_plots <- ATP_standards_means %>% # LINEAR WORKS WELL
                     #dplyr::filter(!ATP_nmol %in% 0) %>%  # ommit zeros
                     ggplot(aes(y=mean, x=ATP_nmol)) + 
                        geom_point() +
                        theme_bw() +
                        #geom_line() +
                        stat_poly_line() +
                        labs(y= "OD 570nm", x = "ATP (nmol, known standards)") +
                        stat_poly_eq(use_label(c("eq", "R2"))) +
                        facet_wrap(~Run_date) 

#print
pdf(paste("Output/Colorimetric_assays/ATP/Standard_Curve_ATP.pdf", sep =''), 
    width=10, 
    height=7)
print(ATP_stand_plots)
dev.off()

# Run_date 20230828, y = 0.0385x - 0.0114 
# Run_date 20230829, y = 0.0994x - 0.0249 



# * (5) Calculate Total protein per sample; subtract out background + correct for standard curve +  average
# Background correction, Run_date 20230828, blank to correct by is 0.04675
# Background correction, Run_date 20230829, blank to correct by is 0.04680

# Standard curve, Run_date 20230828, x = ((Abs_570nm + 0.0114)/0.0385) - NOTE, below we have OD and are solving for x (ATP nmol)
# Standard curve, Run_date 20230829, x = ((Abs_570nm + 0.0249)/0.0994) - NOTE, below we have OD and are solving for x (ATP nmol)


# ATP concentration = (B/V * D) * DDF
# * B = the concentration in nmol ATP from the standard curve 
# * V = the volume added to each well 
# * D = dilution volume if sample was diluted to fit the curve (d/n apply here)
# * DDF = the deproteinization dilution factor LOOK THIS OVER WE ADDED 25 ul TO THE SAMPLE

V = 0.0005 # 50 in ul; as liters = 0.0005
# side note: ATP molecular weight = 507.18 g/mol

# IMPORTANT! we used 50 ul of the standards and 50 ul of the unknowns (samples) 
# therefore we can interpret the unknown direct to the the standard curve without having 
# to account for addition factors in order to infer nmol concentration

ATP_final <- raw_ATP_merged_om %>% 
                      dplyr::filter(Type %in% 'Sample') %>% # call samples
                      dplyr::mutate(Unique_ID = 
                                      paste0(Run_date,'_', Type, Scallop_ID)) %>% # unique ID t0 group by
                      dplyr::mutate(Abs_570nm_cor = # correct the raw abs, subtract background
                                   case_when(Run_date == 20230828 ~ (Abs_570nm-0.04675), # for plate 1
                                             Run_date == 20230829 ~ (Abs_570nm-0.04680) ) ) %>% # for plate 2 
                      dplyr::mutate(ATP_nmol = # correct the raw abs, subtract background
                                    case_when(Run_date == 20230828 ~ ((Abs_570nm + 0.0114)/0.0385), # for plate 1 - nmol to umol by div 1000
                                              Run_date == 20230829 ~ ((Abs_570nm + 0.0249)/0.0994) ) )  %>% # for plate 2 - nmol to umol by div 1000
                      dplyr::group_by(Day,pCO2_history,pCO2_exposure,
                                      Unique_ID, Run_date, Scallop_ID) %>% # group by to get the means
                      dplyr::summarise(mean_ATP_nmol = mean(ATP_nmol),
                                       sd_ATP_nmol   = sd(ATP_nmol),
                                       n = n(),
                                       se_ATP_nmol   = sd_ATP_nmol / sqrt(n))
# write csv
write.csv(ATP_final, file = "Output/Colorimetric_assays/ATP/Calc_Master_ATP.csv")

```


Finally we are here... the moment of truth! 

# MERGE ATP WITH TOTAL PROTEIN DATA 
```{r}

ATP_final <- read.csv(file = "Output/Colorimetric_assays/ATP/Calc_Master_ATP.csv", head = T)
ATP_final_prepmerge    <- ATP_final[,c(2:4,7,8)]

TotalProtein_final <- read.csv(file = "Output/Colorimetric_assays/BCA_ATPcorrection/Calc_Master_BCA_ATPcorrection.csv", head = T)
TotalProtein_prepmerge <- TotalProtein_final[,c(2:4,7,8)]

# IMPORTANT! the ATP assay used 50 ul sample whereas BCA used 25 ul of sample
# therefore to relate ATP to total protein we must divide ATP by 2

MERGED_DF <- merge(ATP_final_prepmerge, 
                   TotalProtein_prepmerge, 
                   all=T)  %>% # looks like we have 3 absent values for total protein..
                dplyr::mutate(ATP_nmol_ug_protein = (mean_ATP_nmol/2)/mean_TotalProtein_ug) # 50 ul for ATP well and 25 ul for BCA well
# View(MERGED_DF) # view the file for any NAs


# IMPORTANT! three individuals did not have enough supernatant to run total protein, listed below as NAs
# Day 1 ID 24 (low history x severe exposure)
# Day 1 ID 35 (moderate history x low exposure)
# Day 14 ID 67 (severe history x low exposure)
MERGED_DF %>% dplyr::filter(ATP_nmol_ug_protein %in% NA)

# lets see the mean_TotalProtein_ug for these treatments  
aggregate(MERGED_DF$mean_TotalProtein_ug, list(MERGED_DF$Day, 
                                               MERGED_DF$pCO2_history,
                                               MERGED_DF$pCO2_exposure
                                               ), FUN=mean, na.rm=TRUE)
# Day 1	   low	 history  severe exposure == 84.99992		(for ID 24)
# Day 1	moderate history	low	   exposure == 47.91459 (for ID 35)
# Day 14 severe	 history  low	   exposure == 40.38476  (for ID 67)

MASTER_DF <- merge(ATP_final_prepmerge, 
                   TotalProtein_prepmerge, 
                   all=T)  %>% # looks like we have 3 absent values for total protein..
                dplyr::mutate(mean_TotalProtein_ug = 
                                case_when(
                                  Scallop_ID %in% 24 ~ 84.99992, # manually based on treatment averages, no sample for protein 
                                  Scallop_ID %in% 35 ~ 47.91459, # manually based on treatment averages, no sample for protein 
                                  Scallop_ID %in% 67 ~ 40.38476, # manually based on treatment averages, no sample for protein 
                                  TRUE ~ mean_TotalProtein_ug),
                              mean_TotalProtein_ng = mean_TotalProtein_ug/1000,
                              ATP_nmol_ng_protein = (mean_ATP_nmol/2)/mean_TotalProtein_ng)

View(MASTER_DF)
# write csv
write.csv(MASTER_DF, file = "Output/Colorimetric_assays/ATP/ATP_Master.csv")

```

# FIGURES 


```{r}

MASTER_DF <- read.csv(file = "Output/Colorimetric_assays/ATP/ATP_Master.csv", head = T)

ggplot(MASTER_DF, aes(mean_ATP_nmol,mean_TotalProtein_ug)) +
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x, size = 1) +
  theme_bw() +
  facet_wrap(~pCO2_history)

ggplot(MASTER_DF, aes(mean_ATP_nmol,mean_TotalProtein_ug)) +
  geom_point()+
  stat_smooth(method = "lm", formula = y ~ x, size = 1) +
  theme_bw() +
  facet_wrap(~pCO2_exposure)

ggplot(MASTER_DF, aes(mean_ATP_nmol,mean_TotalProtein_ug)) +
  geom_point()+
  stat_smooth(method = "lm", formula = y ~ x, size = 1) +
  theme_bw() 
```


* ATP (NOT corrected)
```{r}
MASTER_DF <- read.csv(file = "Output/Colorimetric_assays/ATP/ATP_Master.csv", head = T)

Day1_ATP_RxnNorm <- MASTER_DF %>% # data %>%
  dplyr::filter(Day == 1) %>% # filter for desired date of rhte experiment
  dplyr::mutate(mean_ATP_nmol = as.numeric(mean_ATP_nmol)) %>% 
  dplyr::select(c('Day','pCO2_exposure','pCO2_history','mean_ATP_nmol')) %>% 
  na.omit() %>%   
  # we have means (from x2 assay reps per individual), now mean by treatment!
  group_by(pCO2_exposure, pCO2_history) %>% # group by columns for treatment
  dplyr::summarise( # summarise to acquire the mean and SE for plotting
    ATP_mean = mean(mean_ATP_nmol), # mean
    ATP_sd = sd(mean_ATP_nmol), # sd
    ATP_n = n(), # count
    ATP_se = ATP_sd / sqrt(ATP_n)) %>% # SE
  # plot it
  ggplot(aes(x=pCO2_exposure, y=ATP_mean, group=pCO2_history)) + # call the new mean
    geom_line(aes(group = factor(pCO2_history), 
                  linetype = pCO2_history), 
              size = 0.5, 
              position=position_dodge(.4)) +  # connect a line between variables
    scale_linetype_manual(values=c("solid", "dashed", "dotted")) +
    geom_point(aes(shape=pCO2_history, fill=pCO2_history), size = 4.5,position=position_dodge(.4)) + 
    scale_shape_manual(values=c(21, 22, 24)) + # filled circle, filled triangle, and X 
    scale_fill_manual(values=c("#009E73","#E69F00", "#CC79A7")) + # fill circles
    geom_errorbar(aes(ymin=(ATP_mean)-(ATP_se), # new means and se by treatment
                      ymax=(ATP_mean)+(ATP_se)), # new means and se by treatment
                  width=0,position=position_dodge(.4)) + # width determines the length of the end ticks
    theme_classic() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none") + 
    ggtitle("ATP (not corrected) Day 1") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 4)) + # true 0 origin
    # scale_y_continuous(expand = c(0, 0), limits = c(0, 2.5)) + # true 0 origin
    labs(y= "ATP (nmol)", x = "pCO2 Exposure")
# Day1_ATP_RxnNorm # view the plot


Day14_ATP_RxnNorm <- MASTER_DF %>% # data %>%
  dplyr::filter(Day == 14) %>% # filter for desired date of rhte experiment
  dplyr::mutate(mean_ATP_nmol = as.numeric(mean_ATP_nmol)) %>% 
  dplyr::select(c('Day','pCO2_exposure','pCO2_history','mean_ATP_nmol')) %>% 
  na.omit() %>%   
  # we have means (from x2 assay reps per individual), now mean by treatment!
  group_by(pCO2_exposure, pCO2_history) %>% # group by columns for treatment
  dplyr::summarise( # summarise to acquire the mean and SE for plotting
    ATP_mean = mean(mean_ATP_nmol), # mean
    ATP_sd = sd(mean_ATP_nmol), # sd
    ATP_n = n(), # count
    ATP_se = ATP_sd / sqrt(ATP_n)) %>% # SE
  # plot it
  ggplot(aes(x=pCO2_exposure, y=ATP_mean, group=pCO2_history)) + # call the new mean
    geom_line(aes(group = factor(pCO2_history), 
                  linetype = pCO2_history), 
              size = 0.5, 
              position=position_dodge(.4)) +  # connect a line between variables
    scale_linetype_manual(values=c("solid", "dashed", "dotted")) +
    geom_point(aes(shape=pCO2_history, fill=pCO2_history), size = 4.5,position=position_dodge(.4)) + 
    scale_shape_manual(values=c(21, 22, 24)) + # filled circle, filled triangle, and X 
    scale_fill_manual(values=c("#009E73","#E69F00", "#CC79A7")) + # fill circles
    geom_errorbar(aes(ymin=(ATP_mean)-(ATP_se), # new means and se by treatment
                      ymax=(ATP_mean)+(ATP_se)), # new means and se by treatment
                  width=0,position=position_dodge(.4)) + # width determines the length of the end ticks
    theme_classic() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none") + 
    ggtitle("ATP (not corrected) Day 14") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 4)) + # true 0 origin
    # scale_y_continuous(expand = c(0, 0), limits = c(0, 2.5)) + # true 0 origin
    labs(y= "ATP (nmol)", x = "pCO2 Exposure")
# Day14_ATP_RxnNorm # view the plot

# save plots
pdf("Output/Colorimetric_assays/ATP/ATP_uncorrected_RxnNorm.pdf", width = 7, height = 5)
ggarrange(Day1_ATP_RxnNorm,Day14_ATP_RxnNorm, ncol=2)
dev.off()
```


* ATP (corrected for total protein)
```{r}
MASTER_DF <- read.csv(file = "Output/Colorimetric_assays/ATP/ATP_Master.csv", head = T)

Day1_ATP_RxnNorm <- MASTER_DF %>% # data %>%
  dplyr::filter(Day == 1) %>% # filter for desired date of rhte experiment
  dplyr::mutate(ATP_nmol_ng_protein = as.numeric(ATP_nmol_ng_protein)) %>% 
  dplyr::select(c('Day','pCO2_exposure','pCO2_history','ATP_nmol_ng_protein')) %>% 
  na.omit() %>%   
  # we have means (from x2 assay reps per individual), now mean by treatment!
  group_by(pCO2_exposure, pCO2_history) %>% # group by columns for treatment
  dplyr::summarise( # summarise to acquire the mean and SE for plotting
    ATP_mean = mean(ATP_nmol_ng_protein), # mean
    ATP_sd = sd(ATP_nmol_ng_protein), # sd
    ATP_n = n(), # count
    ATP_se = ATP_sd / sqrt(ATP_n)) %>% # SE
  # plot it
  ggplot(aes(x=pCO2_exposure, y=ATP_mean, group=pCO2_history)) + # call the new mean
    geom_line(aes(group = factor(pCO2_history), 
                  linetype = pCO2_history), 
              size = 0.5, 
              position=position_dodge(.4)) +  # connect a line between variables
    scale_linetype_manual(values=c("solid", "dashed", "dotted")) +
    geom_point(aes(shape=pCO2_history, fill=pCO2_history), size = 4.5,position=position_dodge(.4)) + 
    scale_shape_manual(values=c(21, 22, 24)) + # filled circle, filled triangle, and X 
    scale_fill_manual(values=c("#009E73","#E69F00", "#CC79A7")) + # fill circles
    geom_errorbar(aes(ymin=(ATP_mean)-(ATP_se), # new means and se by treatment
                      ymax=(ATP_mean)+(ATP_se)), # new means and se by treatment
                  width=0,position=position_dodge(.4)) + # width determines the length of the end ticks
    theme_classic() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none") + 
    ggtitle("ATP, Day 1") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 50)) + # true 0 origin
    # scale_y_continuous(expand = c(0, 0), limits = c(0, 2.5)) + # true 0 origin
    labs(y= "ATP (nmol ng protein-1)", x = "pCO2 Exposure")
# Day1_ATP_RxnNorm # view the plot


Day14_ATP_RxnNorm <- MASTER_DF %>% # data %>%
  dplyr::filter(Day == 14) %>% # filter for desired date of rhte experiment
  dplyr::mutate(ATP_nmol_ng_protein = as.numeric(ATP_nmol_ng_protein)) %>% 
  dplyr::select(c('Day','pCO2_exposure','pCO2_history','ATP_nmol_ng_protein')) %>% 
  na.omit() %>%   
  # we have means (from x2 assay reps per individual), now mean by treatment!
  group_by(pCO2_exposure, pCO2_history) %>% # group by columns for treatment
  dplyr::summarise( # summarise to acquire the mean and SE for plotting
    ATP_mean = mean(ATP_nmol_ng_protein), # mean
    ATP_sd = sd(ATP_nmol_ng_protein), # sd
    ATP_n = n(), # count
    ATP_se = ATP_sd / sqrt(ATP_n)) %>% # SE
  # plot it
  ggplot(aes(x=pCO2_exposure, y=ATP_mean, group=pCO2_history)) + # call the new mean
    geom_line(aes(group = factor(pCO2_history), 
                  linetype = pCO2_history), 
              size = 0.5, 
              position=position_dodge(.4)) +  # connect a line between variables
    scale_linetype_manual(values=c("solid", "dashed", "dotted")) +
    geom_point(aes(shape=pCO2_history, fill=pCO2_history), size = 4.5,position=position_dodge(.4)) + 
    scale_shape_manual(values=c(21, 22, 24)) + # filled circle, filled triangle, and X 
    scale_fill_manual(values=c("#009E73","#E69F00", "#CC79A7")) + # fill circles
    geom_errorbar(aes(ymin=(ATP_mean)-(ATP_se), # new means and se by treatment
                      ymax=(ATP_mean)+(ATP_se)), # new means and se by treatment
                  width=0,position=position_dodge(.4)) + # width determines the length of the end ticks
    theme_classic() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none") + 
    ggtitle("ATP, Day 14") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 50)) + # true 0 origin
    # scale_y_continuous(expand = c(0, 0), limits = c(0, 2.5)) + # true 0 origin
    labs(y= "ATP (nmol ng protein-1)", x = "pCO2 Exposure")
# Day14_ATP_RxnNorm # view the plot

# save plots
pdf("Output/Colorimetric_assays/ATP/ATP_proteincorrect_RxnNorm.pdf", width = 7, height = 5)
ggarrange(Day1_ATP_RxnNorm,Day14_ATP_RxnNorm, ncol=2)
dev.off()

```


* Total Protein 
```{r}

MASTER_DF <- read.csv(file = "Output/Colorimetric_assays/ATP/ATP_Master.csv", head = T)


Day1_TotalProtein_RxnNorm <- MASTER_DF %>% # data %>%
  dplyr::filter(Day == 1) %>% # filter for desired date of rhte experiment
  dplyr::mutate(mean_TotalProtein_ng = as.numeric(mean_TotalProtein_ng)) %>% 
  dplyr::select(c('Day','pCO2_exposure','pCO2_history','mean_TotalProtein_ng')) %>% 
  na.omit() %>%   
  # we have means (from x2 assay reps per individual), now mean by treatment!
  group_by(pCO2_exposure, pCO2_history) %>% # group by columns for treatment
  dplyr::summarise( # summarise to acquire the mean and SE for plotting
    TP_mean = mean(mean_TotalProtein_ng), # mean
    TP_sd = sd(mean_TotalProtein_ng), # sd
    TP_n = n(), # count
    TP_se = TP_sd / sqrt(TP_n)) %>% # SE
  # plot it
  ggplot(aes(x=pCO2_exposure, y=TP_mean, group=pCO2_history)) + # call the new mean
    geom_line(aes(group = factor(pCO2_history), 
                  linetype = pCO2_history), 
              size = 0.5, 
              position=position_dodge(.4)) +  # connect a line between variables
    scale_linetype_manual(values=c("solid", "dashed", "dotted")) +
    geom_point(aes(shape=pCO2_history, fill=pCO2_history), size = 4.5,position=position_dodge(.4)) + 
    scale_shape_manual(values=c(21, 22, 24)) + # filled circle, filled triangle, and X 
    scale_fill_manual(values=c("#009E73","#E69F00", "#CC79A7")) + # fill circles
    geom_errorbar(aes(ymin=(TP_mean)-(TP_se), # new means and se by treatment
                      ymax=(TP_mean)+(TP_se)), # new means and se by treatment
                  width=0,position=position_dodge(.4)) + # width determines the length of the end ticks
    theme_classic() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none") + 
    # scale_y_continuous(expand = c(0, 0), limits = c(0, 3)) + # true 0 origin
    scale_y_continuous(expand = c(0, 0), limits = c(0, 0.15)) + # true 0 origin
    ggtitle("Total Protein, Day 1") +
    labs(y= "Total Protie (ng mL)", x = "pCO2 Exposure")
# Day1_TotalProtein_RxnNorm # view the plot


Day14_TotalProtein_RxnNorm <- MASTER_DF %>% # data %>%
  dplyr::filter(Day == 14) %>% # filter for desired date of rhte experiment
  dplyr::mutate(mean_TotalProtein_ng = as.numeric(mean_TotalProtein_ng)) %>% 
  dplyr::select(c('Day','pCO2_exposure','pCO2_history','mean_TotalProtein_ng')) %>% 
  na.omit() %>%   
  # we have means (from x2 assay reps per individual), now mean by treatment!
  group_by(pCO2_exposure, pCO2_history) %>% # group by columns for treatment
  dplyr::summarise( # summarise to acquire the mean and SE for plotting
    TP_mean = mean(mean_TotalProtein_ng), # mean
    TP_sd = sd(mean_TotalProtein_ng), # sd
    TP_n = n(), # count
    TP_se = TP_sd / sqrt(TP_n)) %>% # SE
  # plot it
  ggplot(aes(x=pCO2_exposure, y=TP_mean, group=pCO2_history)) + # call the new mean
    geom_line(aes(group = factor(pCO2_history), 
                  linetype = pCO2_history), 
              size = 0.5, 
              position=position_dodge(.4)) +  # connect a line between variables
    scale_linetype_manual(values=c("solid", "dashed", "dotted")) +
    geom_point(aes(shape=pCO2_history, fill=pCO2_history), size = 4.5,position=position_dodge(.4)) + 
    scale_shape_manual(values=c(21, 22, 24)) + # filled circle, filled triangle, and X 
    scale_fill_manual(values=c("#009E73","#E69F00", "#CC79A7")) + # fill circles
    geom_errorbar(aes(ymin=(TP_mean)-(TP_se), # new means and se by treatment
                      ymax=(TP_mean)+(TP_se)), # new means and se by treatment
                  width=0,position=position_dodge(.4)) + # width determines the length of the end ticks
    theme_classic() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none") + 
    # scale_y_continuous(expand = c(0, 0), limits = c(0, 3)) + # true 0 origin
    scale_y_continuous(expand = c(0, 0), limits = c(0, 0.15)) + # true 0 origin
    ggtitle("Total Protein, Day 14") +
    labs(y= "Total Protien (ng mL)", x = "pCO2 Exposure")
# Day14_TotalProtein_RxnNorm # view the plot

# save plots
pdf("Output/Colorimetric_assays/BCA_ATPcorrection/TotalProtein_RxnNorm.pdf", width = 7, height = 5)
ggarrange(Day1_TotalProtein_RxnNorm,Day14_TotalProtein_RxnNorm, ncol=2)
dev.off()


```



```{r stats - for loop for pCO2 history}

MASTER_DF <- read.csv(file = "Data/Colorimetric_assays/ATP_Master.csv", head = T)

# for loop d
Day  <- as.data.frame(list(unique(as.character(MASTER_DF$Day))))

# for loop i
pCO2history_treatments  <- as.data.frame(list(unique(paste(MASTER_DF$pCO2_history))))

# for loop m
Colorimetric_assays          <- as.data.frame(list(c('ATP_nmol_ng_protein',
                                                'mean_TotalProtein_ng')))

# prepare data for analysis
data <- MASTER_DF %>% 
  dplyr::select(c("Day","pCO2_history","pCO2_exposure",
                  "ATP_nmol_ng_protein","mean_TotalProtein_ng")) 

# Call the cumulative dataframe that we will write to in the for loop below
df_total              <- data.frame() # start dataframe 
stats.table           <- data.frame(matrix(nrow = 1, ncol = 10)) # create dataframe to save cumunalitively during for loop
colnames(stats.table) <- c('Day', 
                          'pCO2_history',
                          'Colorimetric_assays', 
                          'Test', 
                          'Shapiro', 
                          'Levenes',
                          'DF',
                          'Fvalue',
                          'chisquared',
                          'Pvalue') # names for comuns in the for loop


# RUN TE FOR LOOP!

for (d in 1:nrow(Day)) {
  
  x <- (Day[d,])
  data_by_Day <- data %>% 
        dplyr::filter(Day %in% x)
  
  for (i in 1:nrow(pCO2history_treatments)) {
    
    # filter the dataset 
    data_by_history <- data_by_Day %>% 
      dplyr::filter(pCO2_history %in% pCO2history_treatments[i,1]) 
      #dplyr::filter(pCO2_exposure %in% gsub('.* ', '', treatments[i,1])) 
    
      for (m in 1:nrow(Colorimetric_assays)) {
        
        # if loop for probe type-speicif testing based on the naure of the data..
          
        colNum   <- which( colnames(data_by_history)==Colorimetric_assays[m,])
        data_by_history[,colNum] = as.numeric(data_by_history[,colNum])
        modAOV   <- aov(lm(data_by_history[,colNum]~pCO2_exposure, data = data_by_history))
        shapTEST <- shapiro.test(resid(modAOV))
        levTEST  <- leveneTest(modAOV)
        
            # within continuous data.. run either ANOVA or Kruskal Wallis based on normality testing..
            if(shapTEST$p.value> 0.05 & levTEST$`Pr(>F)`[1] > 0.05) { # if normal - run anova
              
              modSUMMARY     <- summary(modAOV) # create anova table summary
              modSUMMARYDF   <- signif(modSUMMARY[[1]][["Df"]][1], 2) # cal p value for anova table - 2 sig figs
              modSUMMARYFval <- signif(modSUMMARY[[1]][["F value"]][1], 4) # cal p value for anova table - 2 sig figs
              modSUMMARYchi  <- NA  
              modSUMMARYPval <- signif(modSUMMARY[[1]][["Pr(>F)"]][1], 4) # cal p value for anova table - 2 sig figs
              testRun        <- "One-way ANOVA"
              
            } else { # if non-normal (either shapiro OR levenes defied) - run the Kruskal Wallis
              
              modSUMMARY     <- kruskal.test(data_by_history[,colNum] ~ pCO2_exposure, data = data_by_history) 
              modSUMMARYDF   <- signif(modSUMMARY[[2]], 2) # cal p value for Kruskal wallis test table - 2 sig figs
              modSUMMARYFval <- NA
              modSUMMARYchi  <- signif(modSUMMARY[[1]], 4) # cal p value for Kruskal wallis test table - 2 sig figs
              modSUMMARYPval <- signif(modSUMMARY[[3]], 4) # cal p value for Kruskal wallis test table - 2 sig figs
              testRun        <- "Kruskal Wallis"
            } # close the if else loop for probe type-specific testing (we are still in the the 'm' for loop!)
        
        stats.table$Day           <- data_by_history$Day[1]
        stats.table$pCO2_history  <- pCO2history_treatments[i,]
        stats.table$Colorimetric_assays   <- Colorimetric_assays[m,]
        stats.table$Test          <- testRun
        if (shapTEST[1] == 'NA') {
          stats.table$Shapiro = 'NA'} else {stats.table$Shapiro <- shapTEST$p.value}
        if (levTEST[1,1] == 'NA') {
          stats.table$Levenes = 'NA'} else {stats.table$Levenes <- levTEST$`Pr(>F)`[1]}
        stats.table$DF            <- modSUMMARYDF
        stats.table$Fvalue        <- modSUMMARYFval
        stats.table$chisquared    <- modSUMMARYchi
        stats.table$Pvalue        <- modSUMMARYPval
        
        df       <- data.frame(stats.table) # name dataframe for this single row
        df_total <- rbind(df_total,df) #bind to a cumulative list dataframe
        # print(df_total) # print to monitor progress      
        
          } # close 'm' for loop - through each flow cy probe
    
    } # close 'i' for loop - through each pCO2 history (low, moderate, severe)
  
} # close 'd' for loop - through each sampling data (5/2/2023 and 5/16/2023)

# View(df_total)
df_total
  # summary(aov(lm(SYBR_PI_count_dead~pCO2_`exposure,data=low)))
```



```{r stats run Two-Way tests}
library(rcompanion) # to run the Schrier -Ray-Hare non parametric 2 way 

# Day 1 data
MASTER_DF_d1 <- MASTER_DF %>% dplyr::filter(Day %in% 1)

# test ATP
Day1_ATP_AOV <- aov(lm(ATP_nmol_ng_protein ~ pCO2_history*pCO2_exposure, data=MASTER_DF_d1))
shapiro.test(resid(Day1_ATP_AOV)) # 0.0001624- non normal
leveneTest(Day1_ATP_AOV) # 0.5982 PASS variance test
Day1_ATP_AOV_SRH <- scheirerRayHare(ATP_nmol_ng_protein ~ pCO2_history*pCO2_exposure, data=MASTER_DF_d1)
Day1_ATP_AOV_SRH
#                            Df Sum Sq      H p.value
# pCO2_history                2  229.7 1.3318 0.51381
# pCO2_exposure               2  252.1 1.4616 0.48151
# pCO2_history:pCO2_exposure  4 1293.3 7.4976 0.11182
# Residuals                  36 5814.8


# Day 14 data
MASTER_DF_d14 <- MASTER_DF %>% dplyr::filter(Day %in% 14)

# test ATP
Day14_ATP_AOV <- aov(lm(ATP_nmol_ng_protein ~ pCO2_history*pCO2_exposure, data=MASTER_DF_d14))
shapiro.test(resid(Day14_ATP_AOV)) # 0.9328 normal
leveneTest(Day14_ATP_AOV) # 9.574e-07 *** NO PASS variance test
Day14_ATP_AOV_SRH <- scheirerRayHare(ATP_nmol_ng_protein ~ pCO2_history*pCO2_exposure, data=MASTER_DF_d14)
Day14_ATP_AOV_SRH
#                            Df Sum Sq      H p.value
# pCO2_history                2  49.80 1.4229 0.49094
# pCO2_exposure               2  26.25 0.7500 0.68729
# pCO2_history:pCO2_exposure  4 142.20 4.0629 0.39757
# Residuals                  11 446.75


```

# one way anovas within history, exposure effect
```{r, one way ANOVA effects of exposure}
library(rcompanion) # to run the Schrier -Ray-Hare non parametric 2 way 

# Day 1 data
MASTER_DF_d1 <- MASTER_DF %>% dplyr::filter(Day %in% 1)

MASTER_DF_d1_LOW <- MASTER_DF_d1 %>% dplyr::filter(pCO2_history %in% 'low')
Day1_ATP_AOV_LOW <- aov(lm(ATP_nmol_ng_protein ~ pCO2_exposure, data=MASTER_DF_d1_LOW))
shapiro.test(resid(Day1_ATP_AOV_LOW)) # 0.7476 normal
leveneTest(Day1_ATP_AOV_LOW) #  0.6173  PASS variance test
summary(Day1_ATP_AOV_LOW)
#               Df Sum Sq Mean Sq F value Pr(>F)  
# pCO2_exposure  2  308.7  154.35   3.846 0.0512 .
# Residuals     12  481.5   40.13

MASTER_DF_d1_MOD <- MASTER_DF_d1 %>% dplyr::filter(pCO2_history %in% 'moderate')
Day1_ATP_AOV_MOD <- aov(lm(ATP_nmol_ng_protein ~ pCO2_exposure, data=MASTER_DF_d1_MOD))
shapiro.test(resid(Day1_ATP_AOV_MOD)) # 0.09701 normal
leveneTest(Day1_ATP_AOV_MOD) #  0.5746  PASS variance test
summary(Day1_ATP_AOV_MOD)
#               Df Sum Sq Mean Sq F value Pr(>F)
# pCO2_exposure  2   71.2    35.6   0.273  0.766
# Residuals     12 1565.0   130.4

MASTER_DF_d1_SEV <- MASTER_DF_d1 %>% dplyr::filter(pCO2_history %in% 'severe')
Day1_ATP_AOV_SEV <- aov(lm(ATP_nmol_ng_protein ~ pCO2_exposure, data=MASTER_DF_d1_SEV))
shapiro.test(resid(Day1_ATP_AOV_SEV)) # 0.005418 NON normal
leveneTest(Day1_ATP_AOV_SEV) #  0.3667 PASS variance test
Day1_ATP_SEV_KW <- kruskal.test(ATP_nmol_ng_protein ~ pCO2_exposure, data=MASTER_DF_d1_SEV)
Day1_ATP_SEV_KW
# Kruskal-Wallis chi-squared = 0.96, df = 2, p-value = 0.6188






# Day 14 data
MASTER_DF_d14 <- MASTER_DF %>% dplyr::filter(Day %in% 14)

MASTER_DF_d14_LOW <- MASTER_DF_d14 %>% dplyr::filter(pCO2_history %in% 'low')
Day14_ATP_AOV_LOW <- aov(lm(ATP_nmol_ng_protein ~ pCO2_exposure, data=MASTER_DF_d14_LOW))
shapiro.test(resid(Day14_ATP_AOV_LOW)) # 0.5464 normal
leveneTest(Day14_ATP_AOV_LOW) #   < 2.2e-16 *** NO PASS variance test
Day14_ATP_LOW_KW <- kruskal.test(ATP_nmol_ng_protein ~ pCO2_exposure, data=MASTER_DF_d14_LOW)
Day14_ATP_LOW_KW
# Kruskal-Wallis chi-squared = 0.4, df = 2, p-value = 0.8187

MASTER_DF_d14_MOD <- MASTER_DF_d14 %>% dplyr::filter(pCO2_history %in% 'moderate')
Day14_ATP_AOV_MOD <- aov(lm(ATP_nmol_ng_protein ~ pCO2_exposure, data=MASTER_DF_d14_MOD))
shapiro.test(resid(Day14_ATP_AOV_MOD)) # 0.1213 normal
leveneTest(Day14_ATP_AOV_MOD) #  < 2.2e-16 ***  NO PASS variance test
Day14_ATP_MOD_KW <- kruskal.test(ATP_nmol_ng_protein ~ pCO2_exposure, data=MASTER_DF_d14_MOD)
Day14_ATP_MOD_KW
# Kruskal-Wallis chi-squared = 3.6, df = 2, p-value = 0.1653

MASTER_DF_d14_SEV <- MASTER_DF_d14 %>% dplyr::filter(pCO2_history %in% 'severe')
Day14_ATP_AOV_SEV <- aov(lm(ATP_nmol_ng_protein ~ pCO2_exposure, data=MASTER_DF_d14_SEV))
shapiro.test(resid(Day14_ATP_AOV_SEV)) # 0.792 normal
leveneTest(Day14_ATP_AOV_SEV) #  0.01186 *  PASS variance test
Day14_ATP_SEV_KW <- kruskal.test(ATP_nmol_ng_protein ~ pCO2_exposure, data=MASTER_DF_d14_SEV)
Day14_ATP_SEV_KW
# Kruskal-Wallis chi-squared = 2.9182, df = 2, p-value = 0.2324

```

# one way anovas within exposure, history effect
```{r, one way ANOVA effects of history}
library(rcompanion) # to run the Schrier -Ray-Hare non parametric 2 way 

# Day 1 data
MASTER_DF_d1 <- MASTER_DF %>% dplyr::filter(Day %in% 1)

MASTER_DF_d1_LOW <- MASTER_DF_d1 %>% dplyr::filter(pCO2_exposure %in% 'low')
Day1_ATP_AOV_LOW <- aov(lm(ATP_nmol_ng_protein ~ pCO2_history, data=MASTER_DF_d1_LOW))
shapiro.test(resid(Day1_ATP_AOV_LOW)) # 0.8536 normal
leveneTest(Day1_ATP_AOV_LOW) #  0.5721  PASS variance test
summary(Day1_ATP_AOV_LOW)
#              Df Sum Sq Mean Sq F value Pr(>F)  
# pCO2_history  2  341.0  170.48   4.706  0.031 *
# Residuals    12  434.7   36.23
TukeyHSD(Day1_ATP_AOV_LOW)
#                      diff        lwr       upr     p adj
# moderate-low    11.313717   1.157924 21.469509 0.0291432
# severe-low       8.164200  -1.991592 18.319993 0.1222530
# severe-moderate -3.149516 -13.305309  7.006276 0.6938877


MASTER_DF_d1_MOD <- MASTER_DF_d1 %>% dplyr::filter(pCO2_exposure %in% 'moderate')
Day1_ATP_AOV_MOD <- aov(lm(ATP_nmol_ng_protein ~ pCO2_history, data=MASTER_DF_d1_MOD))
shapiro.test(resid(Day1_ATP_AOV_MOD)) # 0.009256 NON normal
leveneTest(Day1_ATP_AOV_MOD) #  0.567  PASS variance test
Day1_ATP_MOD_KW <- kruskal.test(ATP_nmol_ng_protein ~ pCO2_history, data=MASTER_DF_d1_MOD)
Day1_ATP_MOD_KW
# Kruskal-Wallis chi-squared = 0.98, df = 2, p-value = 0.6126

MASTER_DF_d1_SEV <- MASTER_DF_d1 %>% dplyr::filter(pCO2_exposure %in% 'severe')
Day1_ATP_AOV_SEV <- aov(lm(ATP_nmol_ng_protein ~ pCO2_history, data=MASTER_DF_d1_SEV))
shapiro.test(resid(Day1_ATP_AOV_SEV)) # 0.258 NON normal
leveneTest(Day1_ATP_AOV_SEV) #  0.4152 PASS variance test
summary(Day1_ATP_AOV_SEV)
#              Df Sum Sq Mean Sq F value Pr(>F)
# pCO2_history  2   20.2   10.11    0.15  0.863
# Residuals    12  810.0   67.50 





# Day 14 data
MASTER_DF_d14 <- MASTER_DF %>% dplyr::filter(Day %in% 14)

MASTER_DF_d14_LOW <- MASTER_DF_d14 %>% dplyr::filter(pCO2_exposure %in% 'low')
Day14_ATP_AOV_LOW <- aov(lm(ATP_nmol_ng_protein ~ pCO2_history, data=MASTER_DF_d14_LOW))
shapiro.test(resid(Day14_ATP_AOV_LOW)) # 0.8218 normal
leveneTest(Day14_ATP_AOV_LOW) #    0.001459 ** NO PASS variance test
Day14_ATP_LOW_KW <- kruskal.test(ATP_nmol_ng_protein ~ pCO2_history, data=MASTER_DF_d14_LOW)
Day14_ATP_LOW_KW
# Kruskal-Wallis chi-squared = 0.5, df = 2, p-value = 0.7788

MASTER_DF_d14_MOD <- MASTER_DF_d14 %>% dplyr::filter(pCO2_exposure %in% 'moderate')
Day14_ATP_AOV_MOD <- aov(lm(ATP_nmol_ng_protein ~ pCO2_history, data=MASTER_DF_d14_MOD))
shapiro.test(resid(Day14_ATP_AOV_MOD)) # 0.683 normal
leveneTest(Day14_ATP_AOV_MOD) #   < 2.2e-16 *** NO PASS variance test
Day14_ATP_MOD_KW <- kruskal.test(ATP_nmol_ng_protein ~ pCO2_history, data=MASTER_DF_d14_MOD)
Day14_ATP_MOD_KW
# Kruskal-Wallis chi-squared = 1.8, df = 2, p-value = 0.4066

MASTER_DF_d14_SEV <- MASTER_DF_d14 %>% dplyr::filter(pCO2_exposure %in% 'severe')
Day14_ATP_AOV_SEV <- aov(lm(ATP_nmol_ng_protein ~ pCO2_history, data=MASTER_DF_d14_SEV))
shapiro.test(resid(Day14_ATP_AOV_SEV)) # 0.3469 normal
leveneTest(Day14_ATP_AOV_SEV) #  0.0001975 *** NO  PASS variance test
Day14_ATP_SEV_KW <- kruskal.test(ATP_nmol_ng_protein ~ pCO2_history, data=MASTER_DF_d14_SEV)
Day14_ATP_SEV_KW
# Kruskal-Wallis chi-squared = 2.8333, df = 2, p-value = 0.2425

```

