library(tidyverse)
library(lme4)
library(dplyr)
library(MASS)
library(lmerTest)
library(performance)
library(detectseparation)

#ULT Crass/Specific Only###############################################################################
#ULT Crass/ Only data for Specific (1C intervals) Assays

#Uploading/subsetting
specific_crass<-read.csv(file="C:/Users/jessi/OneDrive/Move to personal/Ambrosia 2024-2025/Upper Lethal Limit Trials/Spreadsheets_for_analysis/GLM_Specific_Data.csv", header=TRUE, sep=",")
str(specific_crass)
specific_crass$Species<-as.factor(specific_crass$Species)
specific_crass$Host_species<-as.factor(specific_crass$Host_species)
specific_crass$Lifestage<-as.factor(specific_crass$Lifestage)
specific_crass$Rep.<-as.factor(specific_crass$Rep.)
str(specific_crass)
summary(specific_crass)

#fit model
#Overall fit
specific_crass_subset<-specific_crass %>% 
  filter(Species =="X. crassiusculus")

fit_specific_crass_glmer <- glmer(Survival_24hr ~ Temp+Host_species+Lifestage + (1|Rep.),
                                  data=specific_crass_subset,family=binomial(link="logit"))
fit_specific_crass_glm <- glm(Survival_24hr ~ Temp+Host_species+Lifestage,
                              data=specific_crass_subset,family=binomial(link="logit"))

anova(fit_specific_crass_glmer, fit_specific_crass_glm)
#glm is better so use that summary
summary(fit_specific_crass_glm)

#check for perfect separation
fit_specific_crass_separation <- glm(Survival_24hr ~ Temp+Host_species+Lifestage,
                                     data=specific_crass_subset,family=binomial(link="logit"), method="detect_separation")
print(fit_specific_crass_separation)

#Models for individual graphs
#Adults
specific_crass_adult_subset<-specific_crass_subset %>% 
  filter(Lifestage =="Adult")

pdf(file="C:/Users/jessi/OneDrive/Move to personal/Ambrosia 2024-2025/ISPM15_Ambrosia_Beetles/Graph_Crass_Adult_Specific.pdf", width=6.5, height=4.25)
ggplot(specific_crass_adult_subset, aes(x=Temp, y=Survival_24hr, color=Host_species))+
  geom_jitter(height=0)+
  xlab("Temperature (C)")+
  ylab("Survival (Dead/Alive)")+
  labs(color="Host Species")+
  geom_smooth(method="glm", method.args=list(family="binomial"), se=FALSE)+
  scale_color_manual(values=c("Beech"="royalblue2","Red Maple"="red2"))+
  geom_hline(data=data.frame(c(0.50)), aes(yintercept=c(0.50)), color="darkgray", linetype="dashed")+theme_classic()+
  theme(axis.text=element_text(size=15))+theme(axis.title=element_text(size=14))+
  theme(legend.text = element_text(size=13))+ theme(legend.title = element_text(size=13))
dev.off()

pdf(file="C:/Users/jessi/OneDrive/Move to personal/Ambrosia 2024-2025/ISPM15_Ambrosia_Beetles/Graph_Crass_Adult_Specific_black.pdf", width=6.5, height=4.25)
ggplot(specific_crass_adult_subset, aes(x=Temp, y=Survival_24hr, color=Host_species))+
  geom_jitter(height=0)+
  xlab("Temperature (C)")+
  ylab("Survival (Dead/Alive)")+
  labs(color="Host Species")+
  geom_smooth(aes(linetype=Host_species), method="glm", linewidth=1.5,method.args=list(family="binomial"), se=FALSE)+
  scale_color_manual(values=c("Beech"="darkgrey","Red Maple"="black"))+
  geom_hline(data=data.frame(c(0.50)), aes(yintercept=c(0.50)), color="darkgray", linetype="dashed")+theme_classic()+
  theme(axis.text=element_text(size=15))+theme(axis.title=element_text(size=14))+
  theme(legend.text = element_text(size=13))+ theme(legend.title = element_text(size=13))
dev.off()

#Larvae
specific_crass_larvae_subset<-specific_crass_subset %>% 
  filter(Lifestage =="Larvae")

pdf(file="C:/Users/jessi/OneDrive/Move to personal/Ambrosia 2024-2025/ISPM15_Ambrosia_Beetles/Graph_Crass_Larvae_Specific.pdf", width=6.5, height=4.25)
ggplot(specific_crass_larvae_subset, aes(x=Temp, y=Survival_24hr, color=Host_species))+
  geom_jitter(height=0)+
  xlab("Temperature (C)")+
  ylab("Survival (Dead/Alive)")+
  labs(color="Host Species")+
  geom_smooth(method="glm", method.args=list(family="binomial"), se=FALSE)+
  scale_color_manual(values=c("Beech"="royalblue2","Red Maple"="red2"))+
  geom_hline(data=data.frame(c(0.50)), aes(yintercept=c(0.50)), color="darkgray", linetype="dashed")+theme_classic()+
  theme(axis.text=element_text(size=15))+theme(axis.title=element_text(size=14))+
  theme(legend.text = element_text(size=13))+ theme(legend.title = element_text(size=13))
dev.off()

pdf(file="C:/Users/jessi/OneDrive/Move to personal/Ambrosia 2024-2025/ISPM15_Ambrosia_Beetles/Graph_Crass_Larvae_Specific_black.pdf", width=6.5, height=4.25)
ggplot(specific_crass_larvae_subset, aes(x=Temp, y=Survival_24hr, color=Host_species))+
  geom_jitter(height=0)+
  xlab("Temperature (C)")+
  ylab("Survival (Dead/Alive)")+
  labs(color="Host Species")+
  geom_smooth(aes(linetype=Host_species), method="glm", linewidth=1.5, method.args=list(family="binomial"), se=FALSE)+
  scale_color_manual(values=c("Beech"="darkgrey","Red Maple"="black"))+
  geom_hline(data=data.frame(c(0.50)), aes(yintercept=c(0.50)), color="darkgray", linetype="dashed")+theme_classic()+
  theme(axis.text=element_text(size=15))+theme(axis.title=element_text(size=14))+
  theme(legend.text = element_text(size=13))+ theme(legend.title = element_text(size=13))
dev.off()

#Determine LD99

#Adult crass on beech
LD_crass_adult_beech_subset<-specific_crass_adult_subset %>% 
  filter(Host_species =="Beech")
LD_crass_adult_beech <- glm(Survival_24hr ~ Temp, 
                            data=LD_crass_adult_beech_subset,family=binomial(link="logit"))

ld99_crass_adult_beech <- dose.p(cf = 1:2, LD_crass_adult_beech, p = 0.99)
ld99_crass_adult_beech

#Adult crass on maple
LD_crass_adult_maple_subset<-specific_crass_adult_subset %>% 
  filter(Host_species =="Red Maple")
LD_crass_adult_maple <- glm(Survival_24hr ~ Temp, 
                            data=LD_crass_adult_maple_subset,family=binomial(link="logit"))

ld99_crass_adult_maple <- dose.p(cf = 1:2, LD_crass_adult_maple, p = 0.99)
ld99_crass_adult_maple

#Larval crass on beech
LD_crass_larvae_beech_subset<-specific_crass_larvae_subset %>% 
  filter(Host_species =="Beech")
LD_crass_larvae_beech <- glm(Survival_24hr ~ Temp, 
                             data=LD_crass_larvae_beech_subset,family=binomial(link="logit"))

ld99_crass_larvae_beech <- dose.p(cf = 1:2, LD_crass_larvae_beech, p = 0.99)
ld99_crass_larvae_beech

#Larval crass on maple
LD_crass_larvae_maple_subset<-specific_crass_larvae_subset %>% 
  filter(Host_species =="Red Maple")
LD_crass_larvae_maple <- glm(Survival_24hr ~ Temp, 
                             data=LD_crass_larvae_maple_subset,family=binomial(link="logit"))

ld99_crass_larvae_maple <- dose.p(cf = 1:2, LD_crass_larvae_maple, p = 0.99)
ld99_crass_larvae_maple


#ULT germanus/Specific Only##############################################################################
#ULT germanus/ Only data for Specific (1C intervals) Assays

#Uploading/subsetting
specific_germanus<-read.csv(file="C:/Users/jessi/OneDrive/Move to personal/Ambrosia 2024-2025/Upper Lethal Limit Trials/Spreadsheets_for_analysis/GLM_Specific_data.csv", header=TRUE, sep=",")
str(specific_germanus)
specific_germanus$Species<-as.factor(specific_germanus$Species)
specific_germanus$Host_species<-as.factor(specific_germanus$Host_species)
specific_germanus$Lifestage<-as.factor(specific_germanus$Lifestage)
specific_germanus$Rep.<-as.factor(specific_germanus$Rep.)
str(specific_germanus)
summary(specific_germanus)

#fit model
#Overall fit
specific_germanus_subset<-specific_germanus %>% 
  filter(Species =="X. germanus")

fit_specific_germanus_glmer <- glmer(Survival_24hr ~ Temp+Host_species+Lifestage + (1|Rep.),
                                     data=specific_germanus_subset,family=binomial(link="logit"))
fit_specific_germanus_glm <- glm(Survival_24hr ~ Temp+Host_species+Lifestage,
                                 data=specific_germanus_subset,family=binomial(link="logit"))

anova(fit_specific_germanus_glmer, fit_specific_germanus_glm)
#glm is better so use that summary
summary(fit_specific_germanus_glm)

#check for perfect separation
fit_specific_germanus_separation <- glm(Survival_24hr ~ Temp+Host_species+Lifestage,
                                        data=specific_germanus_subset,family=binomial(link="logit"), method="detect_separation")
print(fit_specific_germanus_separation)

#Models for individual graphs
#Adults
specific_germanus_adult_subset<-specific_germanus_subset %>% 
  filter(Lifestage =="Adult")

pdf(file="C:/Users/jessi/OneDrive/Move to personal/Ambrosia 2024-2025/ISPM15_Ambrosia_Beetles/Graph_Germanus_Adult_Specific.pdf", width=6.5, height=4.25)
ggplot(specific_germanus_adult_subset, aes(x=Temp, y=Survival_24hr, color=Host_species))+
  geom_jitter(height=0)+
  xlab("Temperature (C)")+
  ylab("Survival (Dead/Alive)")+
  labs(color="Host Species")+
  geom_smooth(method="glm", method.args=list(family="binomial"), se=FALSE)+
  scale_color_manual(values=c("Beech"="royalblue2","Red Maple"="red2"))+
  geom_hline(data=data.frame(c(0.50)), aes(yintercept=c(0.50)), color="darkgray", linetype="dashed")+theme_classic()+
  theme(axis.text=element_text(size=15))+theme(axis.title=element_text(size=14))+
  theme(legend.text = element_text(size=13))+ theme(legend.title = element_text(size=13))
dev.off()

pdf(file="C:/Users/jessi/OneDrive/Move to personal/Ambrosia 2024-2025/ISPM15_Ambrosia_Beetles/Graph_Germanus_Adult_Specific_black.pdf", width=6.5, height=4.25)
ggplot(specific_germanus_adult_subset, aes(x=Temp, y=Survival_24hr, color=Host_species))+
  geom_jitter(height=0)+
  xlab("Temperature (C)")+
  ylab("Survival (Dead/Alive)")+
  labs(color="Host Species")+
  geom_smooth(aes(linetype=Host_species), method="glm", linewidth=1.5, method.args=list(family="binomial"), se=FALSE)+
  scale_color_manual(values=c("Beech"="darkgrey","Red Maple"="black"))+
  geom_hline(data=data.frame(c(0.50)), aes(yintercept=c(0.50)), color="darkgray", linetype="dashed")+theme_classic()+
  theme(axis.text=element_text(size=15))+theme(axis.title=element_text(size=14))+
  theme(legend.text = element_text(size=13))+ theme(legend.title = element_text(size=13))
dev.off()

#Larvae
specific_germanus_larvae_subset<-specific_germanus_subset %>% 
  filter(Lifestage =="Larvae")

pdf(file="C:/Users/jessi/OneDrive/Move to personal/Ambrosia 2024-2025/ISPM15_Ambrosia_Beetles/Graph_Germanus_Larvae_Specific.pdf", width=6.5, height=4.25)
ggplot(specific_germanus_larvae_subset, aes(x=Temp, y=Survival_24hr, color=Host_species))+
  geom_jitter(height=0)+
  xlab("Temperature (C)")+
  ylab("Survival (Dead/Alive)")+
  labs(color="Host Species")+
  geom_smooth(method="glm", method.args=list(family="binomial"), se=FALSE)+
  scale_color_manual(values=c("Beech"="royalblue2","Red Maple"="red2"))+
  geom_hline(data=data.frame(c(0.50)), aes(yintercept=c(0.50)), color="darkgray", linetype="dashed")+theme_classic()+
  theme(axis.text=element_text(size=15))+theme(axis.title=element_text(size=14))+
  theme(legend.text = element_text(size=13))+ theme(legend.title = element_text(size=13))
dev.off()

pdf(file="C:/Users/jessi/OneDrive/Move to personal/Ambrosia 2024-2025/ISPM15_Ambrosia_Beetles/Graph_Germanus_Larvae_Specific_black.pdf", width=6.5, height=4.25)
ggplot(specific_germanus_larvae_subset, aes(x=Temp, y=Survival_24hr, color=Host_species))+
  geom_jitter(height=0)+
  xlab("Temperature (C)")+
  ylab("Survival (Dead/Alive)")+
  labs(color="Host Species")+
  geom_smooth(aes(linetype=Host_species), method="glm", linewidth=1.5, method.args=list(family="binomial"), se=FALSE)+
  scale_color_manual(values=c("Beech"="darkgrey","Red Maple"="black"))+
  geom_hline(data=data.frame(c(0.50)), aes(yintercept=c(0.50)), color="darkgray", linetype="dashed")+theme_classic()+
  theme(axis.text=element_text(size=15))+theme(axis.title=element_text(size=14))+
  theme(legend.text = element_text(size=13))+ theme(legend.title = element_text(size=13))
dev.off()

#Determine LD99

#Adult germanus on beech
LD_germanus_adult_beech_subset<-specific_germanus_adult_subset %>% 
  filter(Host_species =="Beech")
LD_germanus_adult_beech <- glm(Survival_24hr ~ Temp, 
                               data=LD_germanus_adult_beech_subset,family=binomial(link="logit"))

ld99_germanus_adult_beech <- dose.p(cf = 1:2, LD_germanus_adult_beech, p = 0.99)
ld99_germanus_adult_beech

#Adult germanus on maple
LD_germanus_adult_maple_subset<-specific_germanus_adult_subset %>% 
  filter(Host_species =="Red Maple")
LD_germanus_adult_maple <- glm(Survival_24hr ~ Temp, 
                               data=LD_germanus_adult_maple_subset,family=binomial(link="logit"))

ld99_germanus_adult_maple <- dose.p(cf = 1:2, LD_germanus_adult_maple, p = 0.99)
ld99_germanus_adult_maple

#Larval germanus on beech
LD_germanus_larvae_beech_subset<-specific_germanus_larvae_subset %>% 
  filter(Host_species =="Beech")
LD_germanus_larvae_beech <- glm(Survival_24hr ~ Temp, 
                                data=LD_germanus_larvae_beech_subset,family=binomial(link="logit"))

ld99_germanus_larvae_beech <- dose.p(cf = 1:2, LD_germanus_larvae_beech, p = 0.99)
ld99_germanus_larvae_beech

#Larval germanus on maple
LD_germanus_larvae_maple_subset<-specific_germanus_larvae_subset %>% 
  filter(Host_species =="Red Maple")
LD_germanus_larvae_maple <- glm(Survival_24hr ~ Temp, 
                                data=LD_germanus_larvae_maple_subset,family=binomial(link="logit"))

ld99_germanus_larvae_maple <- dose.p(cf = 1:2, LD_germanus_larvae_maple, p = 0.99)
ld99_germanus_larvae_maple

#ULT maiche/Specific Only##############################################################################
#ULT maiche/ Only data for Specific (1C intervals) Assays

#Uploading/subsetting
specific_maiche<-read.csv(file="C:/Users/jessi/OneDrive/Move to personal/Ambrosia 2024-2025/Upper Lethal Limit Trials/Spreadsheets_for_analysis/GLM_Specific_data.csv", header=TRUE, sep=",")
str(specific_maiche)
specific_maiche$Species<-as.factor(specific_maiche$Species)
specific_maiche$Host_species<-as.factor(specific_maiche$Host_species)
specific_maiche$Lifestage<-as.factor(specific_maiche$Lifestage)
specific_maiche$Rep.<-as.factor(specific_maiche$Rep.)
str(specific_maiche)
summary(specific_maiche)

#fit model
#Overall fit
specific_maiche_subset<-specific_maiche %>% 
  filter(Species =="A. maiche")

fit_specific_maiche_glmer <- glmer(Survival_24hr ~ Temp+Host_species+Lifestage + (1|Rep.),
                                   data=specific_maiche_subset,family=binomial(link="logit"))
fit_specific_maiche_glm <- glm(Survival_24hr ~ Temp+Host_species+Lifestage,
                               data=specific_maiche_subset,family=binomial(link="logit"))

anova(fit_specific_maiche_glmer, fit_specific_maiche_glm)
#glm is better so use that summary
summary(fit_specific_maiche_glm)

#check for perfect separation
fit_specific_maiche_separation <- glm(Survival_24hr ~ Temp+Host_species+Lifestage,
                                      data=specific_maiche_subset,family=binomial(link="logit"), method="detect_separation")
print(fit_specific_maiche_separation)

#Models for individual graphs
#Adults
specific_maiche_adult_subset<-specific_maiche_subset %>% 
  filter(Lifestage =="Adult")

pdf(file="C:/Users/jessi/OneDrive/Move to personal/Ambrosia 2024-2025/ISPM15_Ambrosia_Beetles/Graph_Maiche_Adult_Specific.pdf", width=6.5, height=4.25)
ggplot(specific_maiche_adult_subset, aes(x=Temp, y=Survival_24hr, color=Host_species))+
  geom_jitter(height=0)+
  xlab("Temperature (C)")+
  ylab("Survival (Dead/Alive)")+
  labs(color="Host Species")+
  geom_smooth(method="glm", method.args=list(family="binomial"), se=FALSE)+
  scale_color_manual(values=c("Beech"="royalblue2","Red Maple"="red2"))+
  geom_hline(data=data.frame(c(0.50)), aes(yintercept=c(0.50)), color="darkgray", linetype="dashed")+theme_classic()+
  theme(axis.text=element_text(size=15))+theme(axis.title=element_text(size=14))+
  theme(legend.text = element_text(size=13))+ theme(legend.title = element_text(size=13))
dev.off()

pdf(file="C:/Users/jessi/OneDrive/Move to personal/Ambrosia 2024-2025/ISPM15_Ambrosia_Beetles/Graph_Maiche_Adult_Specific_black.pdf", width=6.5, height=4.25)
ggplot(specific_maiche_adult_subset, aes(x=Temp, y=Survival_24hr, color=Host_species))+
  geom_jitter(height=0)+
  xlab("Temperature (C)")+
  ylab("Survival (Dead/Alive)")+
  labs(color="Host Species")+
  geom_smooth(aes(linetype=Host_species), method="glm", linewidth=1.5, method.args=list(family="binomial"), se=FALSE)+
  scale_color_manual(values=c("Beech"="darkgrey","Red Maple"="black"))+
  geom_hline(data=data.frame(c(0.50)), aes(yintercept=c(0.50)), color="darkgray", linetype="dashed")+theme_classic()+
  theme(axis.text=element_text(size=15))+theme(axis.title=element_text(size=14))+
  theme(legend.text = element_text(size=13))+ theme(legend.title = element_text(size=13))
dev.off()

#Larvae
specific_maiche_larvae_subset<-specific_maiche_subset %>% 
  filter(Lifestage =="Larvae")

pdf(file="C:/Users/jessi/OneDrive/Move to personal/Ambrosia 2024-2025/ISPM15_Ambrosia_Beetles/Graph_Maiche_Larvae_Specific.pdf", width=6.5, height=4.25)
ggplot(specific_maiche_larvae_subset, aes(x=Temp, y=Survival_24hr, color=Host_species))+
  geom_jitter(height=0)+
  xlab("Temperature (C)")+
  ylab("Survival (Dead/Alive)")+
  labs(color="Host Species")+
  geom_smooth(method="glm", method.args=list(family="binomial"), se=FALSE)+
  scale_color_manual(values=c("Beech"="royalblue2","Red Maple"="red2"))+
  geom_hline(data=data.frame(c(0.50)), aes(yintercept=c(0.50)), color="darkgray", linetype="dashed")+theme_classic()+
  theme(axis.text=element_text(size=15))+theme(axis.title=element_text(size=14))+
  theme(legend.text = element_text(size=13))+ theme(legend.title = element_text(size=13))
dev.off()

pdf(file="C:/Users/jessi/OneDrive/Move to personal/Ambrosia 2024-2025/ISPM15_Ambrosia_Beetles/Graph_Maiche_Larvae_Specific_black.pdf", width=6.5, height=4.25)
ggplot(specific_maiche_larvae_subset, aes(x=Temp, y=Survival_24hr, color=Host_species))+
  geom_jitter(height=0)+
  xlab("Temperature (C)")+
  ylab("Survival (Dead/Alive)")+
  labs(color="Host Species")+
  geom_smooth(aes(linetype=Host_species), method="glm", linewidth=1.5, method.args=list(family="binomial"), se=FALSE)+
  scale_color_manual(values=c("Beech"="darkgrey","Red Maple"="black"))+
  geom_hline(data=data.frame(c(0.50)), aes(yintercept=c(0.50)), color="darkgray", linetype="dashed")+theme_classic()+
  theme(axis.text=element_text(size=15))+theme(axis.title=element_text(size=14))+
  theme(legend.text = element_text(size=13))+ theme(legend.title = element_text(size=13))
dev.off()

#Determine LD99

#Adult maiche on beech
LD_maiche_adult_beech_subset<-specific_maiche_adult_subset %>% 
  filter(Host_species =="Beech")
LD_maiche_adult_beech <- glm(Survival_24hr ~ Temp, 
                             data=LD_maiche_adult_beech_subset,family=binomial(link="logit"))

ld99_maiche_adult_beech <- dose.p(cf = 1:2, LD_maiche_adult_beech, p = 0.99)
ld99_maiche_adult_beech

#Adult maiche on maple
LD_maiche_adult_maple_subset<-specific_maiche_adult_subset %>% 
  filter(Host_species =="Red Maple")
LD_maiche_adult_maple <- glm(Survival_24hr ~ Temp, 
                             data=LD_maiche_adult_maple_subset,family=binomial(link="logit"))

ld99_maiche_adult_maple <- dose.p(cf = 1:2, LD_maiche_adult_maple, p = 0.99)
ld99_maiche_adult_maple

#Larval maiche on beech
LD_maiche_larvae_beech_subset<-specific_maiche_larvae_subset %>% 
  filter(Host_species =="Beech")
LD_maiche_larvae_beech <- glm(Survival_24hr ~ Temp, 
                              data=LD_maiche_larvae_beech_subset,family=binomial(link="logit"))

ld99_maiche_larvae_beech <- dose.p(cf = 1:2, LD_maiche_larvae_beech, p = 0.99)
ld99_maiche_larvae_beech

#ULT crass/Nonspecific Only##############################################################################
#ULT crass/ Only data for Nonspecific (4C intervals) Assays

#Uploading/subsetting
nonspecific_crass<-read.csv(file="C:/Users/jessi/OneDrive/Move to personal/Ambrosia 2024-2025/Upper Lethal Limit Trials/Spreadsheets_for_analysis/GLM_NonSpecific_Data.csv", header=TRUE, sep=",")
str(nonspecific_crass)
summary(nonspecific_crass)
nonspecific_crass$Species<-as.factor(nonspecific_crass$Species)
nonspecific_crass$Host_species<-as.factor(nonspecific_crass$Host_species)
nonspecific_crass$Lifestage<-as.factor(nonspecific_crass$Lifestage)
str(nonspecific_crass)


#fit model
#Overall fit
nonspecific_crass_subset<-nonspecific_crass %>% 
  filter(Species =="X. crassiusculus")
fit_nonspecific_crass <- glm(Survival_24hr ~ Temp + Host_species + Lifestage,
                             data=nonspecific_crass_subset,family=binomial(link=logit))
summary(fit_nonspecific_crass)


nonspecific_crass_larvae_subset<-nonspecific_crass_subset %>% 
  filter(Lifestage =="Larvae")
nonspecific_crass_adult_subset<-nonspecific_crass_subset %>% 
  filter(Lifestage =="Adult")

#Models for individual graphs
pdf(file="C:/Users/jessi/OneDrive/Move to personal/Ambrosia 2024-2025/ISPM15_Ambrosia_Beetles/Graph_Crass_NonSpecific_Larvae.pdf", width=6.5, height=4.25)
ggplot(nonspecific_crass_larvae_subset, aes(x=Temp, y=Survival_24hr, color=Host_species))+
  geom_jitter(height=0)+
  xlab("Temperature (C)")+
  ylab("Survival (Dead/Alive)")+
  labs(color="Host Species")+
  geom_smooth(method="glm", method.args=list(family="binomial"), se=FALSE)+
  scale_color_manual(values=c("Beech"="royalblue2","Red Maple"="red2"))+
  geom_hline(data=data.frame(c(0.50)), aes(yintercept=c(0.50)), color="darkgray", linetype="dashed")+theme_classic()+
  theme(axis.text=element_text(size=15))+theme(axis.title=element_text(size=14))+
  theme(legend.text = element_text(size=13))+ theme(legend.title = element_text(size=13))
dev.off()

#ULT germanus/Nonspecific Only##############################################################################
#ULT germanus/ Only data for Nonspecific (4C intervals) Assays

#Uploading/subsetting
nonspecific_germanus<-read.csv(file="C:/Users/jessi/OneDrive/Move to personal/Ambrosia 2024-2025/Upper Lethal Limit Trials/Spreadsheets_for_analysis/GLM_NonSpecific_Data.csv", header=TRUE, sep=",")
str(nonspecific_germanus)
summary(nonspecific_germanus)
nonspecific_germanus$Species<-as.factor(nonspecific_germanus$Species)
nonspecific_germanus$Host_species<-as.factor(nonspecific_germanus$Host_species)
nonspecific_germanus$Lifestage<-as.factor(nonspecific_germanus$Lifestage)
str(nonspecific_germanus)


#fit model
#Overall fit
nonspecific_germanus_subset<-nonspecific_germanus %>% 
  filter(Species =="X. germanus")
fit_nonspecific_germanus <- glm(Survival_24hr ~ Temp + Host_species + Lifestage,
                                data=nonspecific_germanus_subset,family=binomial(link=logit))
summary(fit_nonspecific_germanus)


nonspecific_germanus_larvae_subset<-nonspecific_germanus_subset %>% 
  filter(Lifestage =="Larvae")
nonspecific_germanus_adult_subset<-nonspecific_germanus_subset %>% 
  filter(Lifestage =="Adult")

#Models for individual graphs
pdf(file="C:/Users/jessi/OneDrive/Move to personal/Ambrosia 2024-2025/ISPM15_Ambrosia_Beetles/Graph_Germanus_NonSpecific_Larvae.pdf", width=6.5, height=4.25)
ggplot(nonspecific_germanus_larvae_subset, aes(x=Temp, y=Survival_24hr, color=Host_species))+
  geom_jitter(height=0)+
  xlab("Temperature (C)")+
  ylab("Survival (Dead/Alive)")+
  labs(color="Host Species")+
  geom_smooth(method="glm", method.args=list(family="binomial"), se=FALSE)+
  scale_color_manual(values=c("Beech"="royalblue2","Red Maple"="red2"))+
  geom_hline(data=data.frame(c(0.50)), aes(yintercept=c(0.50)), color="darkgray", linetype="dashed")+theme_classic()+
  theme(axis.text=element_text(size=15))+theme(axis.title=element_text(size=14))+
  theme(legend.text = element_text(size=13))+ theme(legend.title = element_text(size=13))
dev.off()

#ULT maiche/Nonspecific Only##############################################################################
#ULT maiche/ Only data for Nonspecific (4C intervals) Assays

#Uploading/subsetting
nonspecific_maiche<-read.csv(file="C:/Users/jessi/OneDrive/Move to personal/Ambrosia 2024-2025/Upper Lethal Limit Trials/Spreadsheets_for_analysis/GLM_NonSpecific_Data.csv", header=TRUE, sep=",")
str(nonspecific_maiche)
summary(nonspecific_maiche)
nonspecific_maiche$Species<-as.factor(nonspecific_maiche$Species)
nonspecific_maiche$Host_species<-as.factor(nonspecific_maiche$Host_species)
nonspecific_maiche$Lifestage<-as.factor(nonspecific_maiche$Lifestage)
str(nonspecific_maiche)


#fit model
#Overall fit
nonspecific_maiche_subset<-nonspecific_maiche %>% 
  filter(Species =="A. maiche")
fit_nonspecific_maiche <- glm(Survival_24hr ~ Temp + Host_species + Lifestage,
                              data=nonspecific_maiche_subset,family=binomial(link=logit))
summary(fit_nonspecific_maiche)


nonspecific_maiche_larvae_subset<-nonspecific_maiche_subset %>% 
  filter(Lifestage =="Larvae")
nonspecific_maiche_adult_subset<-nonspecific_maiche_subset %>% 
  filter(Lifestage =="Adult")

#Models for individual graphs
pdf(file="C:/Users/jessi/OneDrive/Move to personal/Ambrosia 2024-2025/ISPM15_Ambrosia_Beetles/Graph_Maiche_NonSpecific_Adult.pdf", width=6.5, height=4.25)
ggplot(nonspecific_maiche_adult_subset, aes(x=Temp, y=Survival_24hr, color=Host_species))+
  geom_jitter(height=0)+
  xlab("Temperature (C)")+
  ylab("Survival (Dead/Alive)")+
  labs(color="Host Species")+
  geom_smooth(method="glm", method.args=list(family="binomial"), se=FALSE)+
  scale_color_manual(values=c("Beech"="royalblue2","Red Maple"="red2"))+
  geom_hline(data=data.frame(c(0.50)), aes(yintercept=c(0.50)), color="darkgray", linetype="dashed")+theme_classic()+
  theme(axis.text=element_text(size=15))+theme(axis.title=element_text(size=14))+
  theme(legend.text = element_text(size=13))+ theme(legend.title = element_text(size=13))
dev.off()

