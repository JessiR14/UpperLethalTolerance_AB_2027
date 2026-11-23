
#No error bars because separation issues are causing problems with its calculations

library(ggpattern)

ld99<-read.csv(file="C:/Users/jessi/OneDrive/Move to personal/Ambrosia 2024-2025/Upper Lethal Limit Trials/Spreadsheets_for_analysis/LD datasheet.csv", header=TRUE, sep=",")
summary(ld99)
ld99$Group<-as.factor(ld99$Group)
ld99$Species<-as.factor(ld99$Species)
str(ld99)
##"X. crassiusculus"="crosshatch", "X. germanus"="circle")

crass_subset<-ld99 %>% 
  filter(Species =="X. crassiusculus")
germanus_subset<-ld99 %>% 
  filter(Species =="X. germanus")
maiche_subset<-ld99 %>% 
  filter(Species =="A. maiche")

pdf(file="C:/Users/jessi/OneDrive/Move to personal/Ambrosia 2024-2025/ISPM15_Ambrosia_Beetles/ld99_graph_black_maiche.pdf", width=4, height=4.25)
ggplot(maiche_subset, aes(x=Group, y=LD99, pattern=Species))+
  geom_col_pattern(pattern_density=0.1, pattern_fill="black", fill="lightgrey", pattern_color="black", width=0.6, colour="black")+
  scale_pattern_manual(values=c("A. maiche"="stripe"))+
  coord_cartesian(ylim=c(0,60))+
  xlab("Group")+
  ylab("Lethal Temp. at 99%")+
  scale_x_discrete(labels=c("LV-RM","LV-AB","AD-RM","AD-AB"),limits=rev)+
  geom_hline(data=data.frame(c(56)), aes(yintercept=c(56)), color="darkgray", linetype="dashed")+theme_classic()+
  theme(axis.text=element_text(size=16))+theme(axis.title=element_text(size=16), axis.text.x=element_text(angle=45,vjust=1,hjust=1, size=16))
dev.off()

pdf(file="C:/Users/jessi/OneDrive/Move to personal/Ambrosia 2024-2025/ISPM15_Ambrosia_Beetles/ld99_graph_black_germanus.pdf", width=4, height=4.25)
ggplot(germanus_subset, aes(x=Group, y=LD99, pattern=Species))+
  geom_col_pattern(pattern_density=0.1, pattern_fill="black", fill="lightgrey", pattern_color="black", width=0.6, colour="black")+
  scale_pattern_manual(values=c("X. germanus"="circle"))+
  coord_cartesian(ylim=c(0,60))+
  xlab("Group")+
  ylab("Lethal Temp. at 99%")+
  scale_x_discrete(labels=c("LV-RM","LV-AB","AD-RM","AD-AB"),limits=rev)+
  geom_hline(data=data.frame(c(56)), aes(yintercept=c(56)), color="darkgray", linetype="dashed")+theme_classic()+
  theme(axis.text=element_text(size=16))+theme(axis.title=element_text(size=16), axis.text.x=element_text(angle=45,vjust=1,hjust=1, size=16))
dev.off()

pdf(file="C:/Users/jessi/OneDrive/Move to personal/Ambrosia 2024-2025/ISPM15_Ambrosia_Beetles/ld99_graph_black_crass.pdf", width=4, height=4.25)
ggplot(crass_subset, aes(x=Group, y=LD99, pattern=Species))+
  geom_col_pattern(pattern_density=0.1, pattern_fill="black", fill="lightgrey", pattern_color="black", width=0.6, colour="black")+
  scale_pattern_manual(values=c("X. crassiusculus"="crosshatch"))+
  coord_cartesian(ylim=c(0,60))+
  xlab("Group")+
  ylab("Lethal Temp. at 99%")+
  scale_x_discrete(labels=c("LV-RM","LV-AB","AD-RM","AD-AB"),limits=rev)+
  geom_hline(data=data.frame(c(56)), aes(yintercept=c(56)), color="darkgray", linetype="dashed")+theme_classic()+
  theme(axis.text=element_text(size=16))+theme(axis.title=element_text(size=16), axis.text.x=element_text(angle=45,vjust=1,hjust=1, size=16))
dev.off()
