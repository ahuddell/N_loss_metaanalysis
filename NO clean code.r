#Nitric oxide data analysis
setwd("~/N Loss Metaanalysis/analysis code and data")

library(lme4)
library(car)
library(MuMIn)
library(lmerTest)
library(MASS)
library(ggplot2)
#library(stringr)
library(dplyr)

#standardizing variables 2x sd per Gelman reccomendation 
z.trans<-function(x) {(x - mean(x, na.rm=T))/(2*sd(x, na.rm=T))}


#theme for plots
theme_default <- function(axis_text_size = 13) {
  theme(text=element_text(size=16, colour="black"),
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.text.x=element_text(size=axis_text_size, colour="black"),
        axis.text.y=element_text(size=axis_text_size, colour="black"),
        axis.title.y=element_text(angle=90, vjust=-0.5),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.background=element_rect(fill="white", color="white"),
        legend.position="none") 
}
# NO ----------------------------------------------------------------------
no<-read.csv("no_database.csv")

#Q1/Q2 How do NO3-, N2O and NO losses from tropical agricultural systems respond to N inputs? (2) Is the relationship of N losses relative to N inputs different between tropical and temperate regions? --------
# set up dataframe for eq. 1
no_full<-no
NO_model_full=data.frame(row.names = 1:nrow(no_full))
NO_model_full$site_lat <-no_full$lat_decimal #no trans. needed
NO_model_full$tropical<-(no_full$tropical) #no trans. needed
NO_model_full$log_no_season<-log10(no_full$NO_kg_seas_zero_rm) #logging response variable
NO_model_full$N_in_kg_untrans<-(no_full$N_in_kg) #z-transforming
NO_model_full$N_in_kg<-z.trans(no_full$N_in_kg) #z-transforming
NO_model_full$duration_d<-z.trans(no_full$duration_d) #z-transforming

#remove incomplete cases for lmer function
dim(NO_model_full) 
NO_model_full<-na.omit(NO_model_full)
dim(NO_model_full) #removed imcomplete cases

#########eq. 1
meQ1=lmer( (log_no_season)~  duration_d  +  N_in_kg*tropical+ (1|site_lat), data=NO_model_full, na.action=NULL)
summary(meQ1) 
r.squaredGLMM(meQ1)

#population prediction intervals to calculate CIs from p.257 Bolker 2008
x<- 1:2000
cis <- 0.95
lowb <- 0.5-(cis/2)
upb <- 0.5+(cis/2)

vmat<-mvrnorm(1000,mu=fixef(meQ1), Sigma=vcov(meQ1))
a<-vmat[,3] #coef for N fertilizer (if temperate)
b<- vmat[,2] #coef for study duration
c<-vmat[,1] #global intercept 
f<-vmat[,4] #coef for tropical=1
g<-vmat[,5] #coef for interaction N in if tropical

n_bar<-mean(no_full$N_in_kg, na.rm=T)
n_sd<-sd(no_full$N_in_kg, na.rm=T)

#function for CI curves
yQ2_tropical_CIs=function(a,b,c,f,g,x) {10^(c)*10^(f)*10^((a+g)*((x-n_bar)/(2*n_sd)))}
yQ2_temperate_CIs=function(a,b,c,f,g,x) {10^(c)*10^((a)*((x-n_bar)/(2*n_sd)))}

#calculating CIs for tropical curve
dist=array(dim=c(1000,length(x)))
for (i in 1:1000) {
  dist[i,] = yQ2_tropical_CIs(a=a[i],b=b[i],c=c[i],f=f[i],g=g[i],x=x)
}
civec_yQ2_trop <- array(dim=c(2,length(x)))
for(j in 1:length(x)){
  civec_yQ2_trop [,j] <- quantile(dist[,j],c(lowb,upb),na.rm=TRUE)
}

#calculating CIs for temperate curve
dist=array(dim=c(1000,length(x)))
for (i in 1:1000) {
  dist[i,] = yQ2_temperate_CIs(a=a[i],b=b[i],c=c[i],f=f[i],g=g[i],x=x)
}
civec_yQ2_temp <- array(dim=c(2,length(x)))
for(j in 1:length(x)){
  civec_yQ2_temp [,j] <- quantile(dist[,j],c(lowb,upb),na.rm=TRUE)
}

#eq. 1 fitted model outputs for plotting
a<-fixef(meQ1)[3] #coef for N fertilizer (if temperate)
b<- fixef(meQ1)[2]#coef for study duration
c<-fixef(meQ1)[1] #global intercept 
f<- fixef(meQ1)[4]#coef for tropical=1
g<- fixef(meQ1)[5] #coef for interaction N in if tropical

##########plotting on linear scale
#functions to calculate tropical and temperate model (eq. 1) fit curves
yQ2_tropical=function(x) {10^(c)*10^(f)*10^((a +g)*((x-n_bar)/(2*n_sd)))}
yQ2_temperate=function(x) {10^(c)*10^(a*((x-n_bar)/(2*n_sd)))}

# reordering tropical and temperate levels
no_full$tropical<-factor(as.factor(no_full$tropical), rev(levels(as.factor(no_full$tropical))))
levels(as.factor(no_full$tropical)) #tropical now comes first

#calculating curves and CIs
ytrop<-yQ2_tropical(x)
ytemp<- yQ2_temperate(x)
CItroplb<-civec_yQ2_trop [1,]
CItropub<-civec_yQ2_trop [2,]
CItemplb<-civec_yQ2_temp [1,]
CItempub<-civec_yQ2_temp [2,]
df<-data.frame(x,ytrop,ytemp,CItroplb, CItropub, CItemplb, CItempub)

#plot
nop<-ggplot(no_full, aes(x=N_in_kg,y=NO_kg_seas_zero_rm)) +
  geom_point(aes(colour=(tropical), shape=(tropical)), alpha=.55,
             position = "jitter", size=2) +
   ylab(expression(paste('Nitric oxide (kg NO-N ',~ha^-1, season^-1,')'))) +
  xlab(expression(paste('N inputs (kg N ha'^'-1', 
                        paste('season'^'-1'), ')'))) +
  ggtitle('Nitric oxide losses') +
  theme_default() +
  ylim(0,15) + xlim (0,300)+ 
  geom_segment(aes(x=0,xend=0,y=0,yend=15), colour="black") +
  geom_segment(aes(x=0,xend=300,y=0,yend=0),colour="black") +
  scale_shape_manual(values=c(16,17)) +
  geom_line(data=df, aes(x,ytrop), col='red3', size=1)+
  geom_line(data=df, aes(x,ytemp), col='dodgerblue', size=1) +
  geom_ribbon(data = df, aes(x=x, y=ytrop,
                             ymin = CItroplb, ymax = CItropub), 
              fill = "red3", alpha=.12) + 
  geom_ribbon(data = df, aes(x=x, y=ytemp,
                             ymin = CItemplb, ymax = CItempub), 
              fill = "dodgerblue", alpha=.12) +
  theme(legend.position="right", legend.key=element_blank()) +
  scale_color_manual(name= 'Region',
                     labels= c('Tropical', 'Temperate'),
                     values = c('red3', 'dodgerblue')) +
  scale_shape_manual(name= 'Region',
                     labels= c('Tropical', 'Temperate'),
                     values = c(16,17))
nop

#plotting full range of Fig 2 for SFig. 1
#calculating max of x axis range for tropical and temperate
N_in_trop<-  no_full %>% filter(tropical==1) %>% dplyr::select(N_in_kg)
N_in_temp <-no_full %>% filter(tropical==0) %>% dplyr::select(N_in_kg)
x_trop<-1:max(N_in_trop)
x_temp<-1:max(N_in_temp)
ytrop<-yQ2_tropical(x_trop)
ytemp<- yQ2_temperate(x_temp)
CItroplb<-civec_yQ2_trop [1,(1:max(N_in_trop))]
CItropub<-civec_yQ2_trop [2,(1:max(N_in_trop))]
CItemplb<-civec_yQ2_temp [1,(1:max(N_in_temp))]
CItempub<-civec_yQ2_temp [2,(1:max(N_in_temp))]
df_trop<-data.frame(x_trop,ytrop,CItroplb, CItropub)
df_temp<-data.frame(x_temp,ytemp, CItemplb, CItempub)

#plotting full range of Fig 2 for SFig. 1
nop_full<-ggplot(no_full, aes(x=N_in_kg,y=NO_kg_seas_zero_rm)) +
  geom_point(aes(colour=factor(tropical), shape=factor(tropical)), alpha=.6,
             position = "jitter", size=2) +
    ylab(expression(paste('Nitric oxide (kg NO-N ',~ha^-1, season^-1,')'))) +
  xlab(expression(paste('N inputs (kg N ha'^'-1', 
                        paste('season'^'-1'), ')'))) +
  ggtitle('Nitric oxide losses') +
  theme_default() +
  geom_segment(aes(x=0,xend=0,y=0,yend=40), colour="black") +
  geom_segment(aes(x=0,xend=900,y=0,yend=0),colour="black") +
  geom_line(data=df_trop, aes(x_trop,ytrop), col='red3')+
  geom_line(data=df_temp, aes(x_temp,ytemp), col='dodgerblue') +
  geom_ribbon(data = df_trop, aes(x=x_trop, y=ytrop,
                                  ymin = CItroplb, ymax = CItropub),
              fill = "red3", alpha=.12) +
  geom_ribbon(data = df_temp, aes(x=x_temp, y=ytemp,
                                  ymin = CItemplb, ymax = CItempub),
              fill = "dodgerblue", alpha=.12) +
  theme(legend.position="right", legend.key=element_blank()) +
  scale_color_manual(name= 'Region',
                     labels= c('Tropical', 'Temperate'),
                     values = c('red3', 'dodgerblue')) +
  scale_shape_manual(name= 'Region',
                     labels= c('Tropical', 'Temperate'),
                     values = c(16,17)) 
  
nop_full

##################
#library(ggpubr)

#only re-run to update plot
#Q2plot<-ggarrange(no3p,n2op,nop, ncol=3, nrow=1, 
#                   labels=c("a","b","c"), 
#                   common.legend = TRUE, legend="right")
#ggsave('Fig 2 Q2plot.png', plot = Q2plot, device = NULL,
#        path = NULL, scale = 3, width = 8.7, height = 4.5, 
#        units = 'cm', dpi = 500, limitsize = TRUE)

#Tropical N losses at 50 and 150 kg N
y50<-yQ2_tropical(50) #no3 losses at 50 kg N inputs
y150<-yQ2_tropical(150) #no3 losses at 150 kg N inputs
(y150-y50)/y50*100 #% increase

#calculation number of observations in Figs 1, 2, S1
unique(no$lat_decimal)
no %>% filter (N_in_kg<=300, NO_kg_seas_zero_rm<=15) %>% summarize (n())
no %>% summarize (n())


#############################################################################################
# Q3 Which environmental and management control the magnitude of N losses in tropical agroecosystems?   -------------------------------------------
#############################################################################################
no_trop<-no%>% filter(tropical > 0)

#set up dataframe for eq. 2 model
NO_model_tropical<-data.frame(row.names = 1:length(no_trop$NO_kg_seas_zero_rm))
NO_model_tropical$site_lat<-no_trop$lat_decimal # site ID; no transformation needed
NO_model_tropical$crop_type<-no_trop$crop_type # crop type; no transformation needed
NO_model_tropical$log_g_NO_day<-log10(no_trop$NO_kg_seas_zero_rm*1000/no_trop$duration_d) # crop type; no transformation needed
NO_model_tropical$split_app<-no_trop$split_app #split=1 (no trans. needed)
NO_model_tropical$org_Nfix_in<-no_trop$org_Nfix_in #organic/Nfix in=1  (no trans. needed)
NO_model_tropical$N_in_kg<-z.trans(no_trop$N_in_kg) #z-transforming 
NO_model_tropical$duration_d<-z.trans(no_trop$duration_d) #z-transforming 
NO_model_tropical$precip<-z.trans(no_trop$MAP_mm_combined) #z-transforming 

#remove incomplete cases for lmer function
dim(NO_model_tropical) 
NO_model_tropical<-na.omit(NO_model_tropical)
dim(NO_model_tropical) #removed incomplete cases 
#removed sand because it reduces the sample size so much that it causes rank deficiencies

#test for collinearity
vif_test=lmer(log_g_NO_day~ N_in_kg  + precip + split_app + org_Nfix_in +
                (1|site_lat), data=NO_model_tropical, na.action=NULL) 
vif(vif_test) #all covariates below 3

meq2=lmer(log_g_NO_day~ N_in_kg  + precip + split_app + org_Nfix_in +  
            crop_type+(1|site_lat), data=NO_model_tropical, 
            na.action=NULL) 
summary(meq2) 
r.squaredGLMM(meq2)

#plotting key relationships from eq2
#continous variables plot

#precipitation/irrigation input plot
no_trop$log_no_g_day<-log10(no_trop$NO_kg_seas_zero_rm*1000/no_trop$duration_d) #y
x=800:4200
precip<- mean(no_trop$MAP_mm_combined,na.rm=T) #to unstandardize sand
precip_sd<-sd(no_trop$MAP_mm_combined,na.rm=T) #to unstandardize sand
precip_fit_SE<- sqrt(diag(vcov(meq2)))[3] #SE amounts
yprecip<-function(x){fixef(meq2)[1]+fixef(meq2)[3]*((x-precip)/(2*precip_sd))}
yprecip_hat<-yprecip(x=x)
yprecip_hat_lb<-yprecip_hat-rep(precip_fit_SE, length(yprecip_hat))
yprecip_hat_ub<-yprecip_hat+rep(precip_fit_SE, length(yprecip_hat))
df<-data.frame(x,yprecip_hat,yprecip_hat_lb,yprecip_hat_ub)

#precip plot
no_precip<- ggplot(no_trop, aes(x=MAP_mm_combined,y=log_no_g_day)) +
  geom_point(aes(alpha=.6)) +
  ylab(expression(paste('nitric oxide (g NO-N ',~ha^-1,~d^-1, ')'))) +
  xlab(expression(paste('mean annual precipitation (mm)'))) +
  theme_default() +
  ylim(0,2) + xlim (800,4200)+
  geom_segment(aes(x=800,xend=800,y=0,yend=2), colour="black") +
  geom_segment(aes(x=800,xend=4200,y=0,yend=0),colour="black") +
  geom_line(data=df, aes(x,yprecip_hat), col='red3')+
  geom_ribbon(data=df,aes(x=x, y=yprecip_hat, ymin = yprecip_hat_lb, ymax = yprecip_hat_ub),
              fill = "red3", alpha=.12) 
no_precip<-no_precip + scale_y_continuous(breaks=c(0,1,2),
                        labels=c(0,10,100))
no_precip

#count of plot
length(na.omit(no_trop$MAP_mm_combined))

##only run to update Figure 3
 library(ggpubr)
 Q3_contin<-ggarrange(no3_precip, no3_sand, n2o_precip,
                     n2o_sand, no_precip,ncol=2, nrow=3,
                     labels=c("a","b","c","d","e"))
 ggsave('Fig 3.jpeg', plot = Q3_contin, device = NULL,
       path = NULL, scale = 1, width = 9, height = 11,
       dpi = 500, limitsize = TRUE, units = "in")


#Fig 4 binary/dummy crops plot
coef <- data.frame(data = 
                   c(fixef(meq2)[c('(Intercept)','(Intercept)',
                   'split_app','(Intercept)','org_Nfix_in','crop_typefallow',
                   'crop_typepasture','crop_typerice','crop_typetree crop', 
                   'crop_typelegume')]) +
                     c(rep(0,2),fixef(meq2)['(Intercept)'], 0, 
                       rep(fixef(meq2)['(Intercept)'],6))
)
#adding the intercept (cereal) to all model fits
meq2.se <- sqrt(diag(vcov(meq2))) #standard errors from meq2
meq2.names <- row.names(vcov(meq2)) #variable names of meq2
meq2.se.comb<-as.data.frame(cbind(meq2.names,meq2.se)) #data frame of parameter SEs
meq2.se.comb<-data.frame(meq2.names,meq2.se) #data frame of parameter SEs

SE <- c(meq2.se.comb[meq2.names=='(Intercept)',2],
        meq2.se.comb[meq2.names=='(Intercept)',2],
        meq2.se.comb[meq2.names=='split_app',2],
        meq2.se.comb[meq2.names=='(Intercept)',2],
        meq2.se.comb[meq2.names=='org_Nfix_in',2],
        meq2.se.comb[meq2.names=='crop_typefallow',2],
        meq2.se.comb[meq2.names=='crop_typepasture',2],
        meq2.se.comb[meq2.names=='crop_typerice',2],
        meq2.se.comb[meq2.names=='crop_typetree crop',2],
        meq2.se.comb[meq2.names=='crop_typelegume',2])

meq2_outputs<- cbind(coef, SE)
colnames(meq2_outputs) <- c('coef', 'SE')
meq2_outputs$names<-c('cereal','sa0','sa1','org0','org1', 'fallow',
                      'pasture', 'rice', 'tree crop', 'legume')
meq2_outputs


#new plot with y data and model fit

#creating individual dfs for binary variables to join to crop types
discrete_var<-no_trop %>%
  dplyr::select(org_Nfix_in ,split_app, log_no_g_day) 
org1<-discrete_var%>% filter(org_Nfix_in==1) %>% 
  mutate(crop_type= rep("org1",3))
org0<-discrete_var%>% filter(org_Nfix_in==0)%>% 
  mutate(crop_type= rep("org0",46))
sa1<-discrete_var%>% filter(split_app==1)%>% 
  mutate(crop_type= rep("sa1",15))
sa0<-discrete_var%>% filter(split_app==0)%>% 
  mutate(crop_type= rep("sa0",34))


#joining the data frames
discrete_var2<-no_trop %>%
  dplyr::select(org_Nfix_in ,split_app, crop_type, log_no_g_day) 

discrete_var3<-rbind(discrete_var2,org0,org1,sa1,sa0)

#reordering the factors based on the model outputs 
discrete_var3$crop_type<-factor(discrete_var3$crop_type,
                                  levels=c('cereal','sa0','sa1','org0','org1',
                                           'fallow','pasture', 'rice',
                                           'tree crop', 'legume'))

#x-axis names
x_axis<-c('cereal/intercept', 'split application=0','split application=1',
          'organic/N-fix in=0', 'organic/N-fix in=1', 'fallow',  
          'pasture','rice', 'tree crop', 'legume')

#first removing the 2 "other crop types" which are not estimated by the model
discrete_var3<-discrete_var3%>% filter(crop_type!='other')        


#fig 4--we removed the violin plots from these because the outline covered data points
no_Q3_disc<-  ggplot() +
  geom_violin(data=discrete_var3,
              aes(x=crop_type,y=log_no_g_day), fill="light gray", 
                  colour = "light gray" )+
  geom_jitter(data=discrete_var3, aes(x=crop_type,y=log_no_g_day), 
              width=.05, alpha=.3) +
  ylab(expression(paste('nitric oxide (g NO-N',~~ha^-1,~d^-1, ')'))) +
  geom_pointrange(data = meq2_outputs,aes(reorder(names, coef), coef, 
                  ymin = coef-SE, ymax =  coef+SE), col='red') +
  scale_colour_manual(values = coef) +
  scale_x_discrete('' ,labels=x_axis) +
  geom_hline(yintercept=0, color = "black") +
  theme(axis.text.x = element_text(angle=30, hjust=1)) +
  geom_segment(aes(x=0,xend=0,y=0,yend=2), colour="black") +
  theme_default(axis_text_size = 13) +
  geom_vline(xintercept=c(1.5,3.5, 5.5), linetype='longdash')

no_Q3_disc<-no_Q3_disc + scale_y_continuous(breaks= c(0,1,2),
                                            labels=c(0,10,100))
no_Q3_disc

#count of plot
length(na.omit(no_trop$crop_type))

####only run to update the plot
library(ggpubr)
# Fig4<-ggarrange(no3_Q3_disc,n2o_Q3_disc,no_Q3_disc,
#                   ncol=1, nrow=3, labels=c("a","b","c"))
# ggsave('Fig 4.png', plot = Fig4, device = NULL,
#                        width=3, height = 6, units = 'in', path = NULL,
#                       scale=2.2, dpi = 500, limitsize = TRUE)
# 

#tables--------------------------------------------

rm(list = ls()) #clear workspace

#standardizing variables 2x sd per Gelman reccomendation 
z.trans<-function(x) {(x - mean(x, na.rm=T))/(2*sd(x, na.rm=T))}

#theme for plots
theme_default <- function(axis_text_size = 13) {
  theme(text=element_text(size=16, colour="black"),
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.text.x=element_text(size=axis_text_size, colour="black"),
        axis.text.y=element_text(size=axis_text_size, colour="black"),
        axis.title.y=element_text(angle=90, vjust=-0.5),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.background=element_rect(fill="white", color="white"),
        legend.position="none")
}

no_full<-read.csv("no_database.csv")

#number of papers and observations
no_trop<-no_full %>% filter(tropical > 0)
levels(no_trop$title) #100 unique titles
dim(no_trop) #616 observations


#regional breakdown for table 1
 region_breakdown_tropical<-no_trop %>%
   group_by(region) %>%
   summarise(count = n())
 
 region_breakdown_tropical$pct<-region_breakdown_tropical$count/
   sum(region_breakdown_tropical$count)*100
 region_breakdown_tropical
 
 #percent ox/ult
 perc_oxult1=(sum(no_trop$ox_ult==1, na.rm=T)/length(no_trop$ox_ult)*100)
 perc_oxult1
 
#crop type breakdown for table S1
crop_breakdown_tropical<-no_trop %>%
  group_by(crop_type) %>%
  summarise(count = n())

crop_breakdown_tropical$pct<-crop_breakdown_tropical$count/
  sum(crop_breakdown_tropical$count)*100
crop_breakdown_tropical

#irrigation, split applicaiton, and organic inputs breakdown
irrig_trop<- no_trop %>% filter(irrig_1_0==1) %>%
  summarise(count = n())
irrig_pct<-irrig_trop$count/
  nrow(no_trop)*100

split_app_trop<- no_trop %>%  filter(split_app==1) %>%
  summarise(count = n())
split_app_pct<-split_app_trop$count/
  nrow(no_trop)*100

org_Nfix_trop<- no_trop %>%  filter(org_Nfix_in==1) %>%
  summarise(count = n())
org_Nfix_in_pct<-org_Nfix_trop$count/
  nrow(no_trop)*100


# supplemental figures ----------------------------------------------------
# Q2plot<-ggarrange(no3p_full,n2op_full,nop_full, ncol=3, nrow=1,
#                   labels=c("a","b","c"),
#                   common.legend = T, legend="right")
# ggsave('Fig S1.png', plot = Q2plot, device = NULL,
#        path = NULL, scale = 4, width = 8.7, height = 4,
#        units = 'cm', dpi = 500, limitsize = TRUE)

no_model_trop_plot=data.frame(row.names = 1:49)
no_model_trop_plot$N_in_kg<-(no_trop$N_in_kg)
no_model_trop_plot$duration_d<-(no_trop$duration_d)
no_model_trop_plot$MAP_mm_combined<-no_trop$MAP_mm_combined
no_model_trop_plot$irrigation<-no_trop$irrig_1_0 
no_model_trop_plot$split_app<-no_trop$split_app 
no_model_trop_plot$org_Nfix_in<-no_trop$org_Nfix_in 
no_model_trop_plot$reg<-no_trop$reg

#removing NA's and using complete cases only
dim(no_model_trop_plot) # 170, 5
no_model_trop_plot<-na.omit(no_model_trop_plot)

#N inputs boxplot
p1<-ggplot(no_model_trop_plot, aes(x=" " , y=N_in_kg)) + 
  geom_boxplot() +
  ylab(expression(paste('N inputs (kg N '~ha^-1,~y^-1, ')'))) +
  xlab(" ") +
  theme(text=element_text(size=16, colour="black"),
        plot.title = element_text(size = 16, hjust = 0.5), 
        axis.text.x=element_text(size=16, colour="black"),
        axis.text.y=element_text(size=16, colour="black")) +
  theme(panel.background=element_rect(fill="white", color="white"),
        legend.position="none") +
  geom_segment(aes(x=0,xend=0,y=0,yend=400), colour="black") +
  geom_segment(aes(x=0,xend=2,y=0,yend=0),colour="black") 
p1

#count of plot
length(na.omit(no_model_trop_plot$N_in_kg))

#MAP boxplot
p2<-ggplot(no_model_trop_plot, aes(x=" " , y=MAP_mm_combined)) + 
  geom_boxplot() +
  ylab(expression(paste('mean annual precipitation (mm '~~~y^-1, ')'))) +
  xlab(" ") +
  theme(text=element_text(size=16, colour="black"),
        plot.title = element_text(size = 16, hjust = 0.5), 
        axis.text.x=element_text(size=16, colour="black"),
        axis.text.y=element_text(size=16, colour="black")) +
  theme(panel.background=element_rect(fill="white", color="white"),
        legend.position="none") +
  geom_segment(aes(x=0,xend=0,y=0,yend=4000), colour="black") +
  geom_segment(aes(x=0,xend=2,y=0,yend=0),colour="black") 
p2

#duration_d boxplot
p4<-ggplot(no_model_trop_plot, aes(x=" " , y=duration_d)) + 
  geom_boxplot() +
  ylab(expression(paste('study duration (days)'))) +
  xlab(" ") +
  theme(text=element_text(size=16, colour="black"),
        plot.title = element_text(size = 16, hjust = 0.5), 
        axis.text.x=element_text(size=16, colour="black"),
        axis.text.y=element_text(size=16, colour="black")) +
  theme(panel.background=element_rect(fill="white", color="white"),
        legend.position="none") +
  geom_segment(aes(x=0,xend=0,y=0,yend=500), colour="black") +
  geom_segment(aes(x=0,xend=2,y=0,yend=0),colour="black") 
p4

################# publishing plots
library(ggpubr)
#only run to update the plot
SI_Q3_no_inputs<-ggarrange(p1,p2,p4,ncol=3, nrow=1, labels=c("a","b","c")) 
#SI_Q3_no_inputs<-annotate_figure(SI_Q3_no_inputs, fig.lab=("Summary of covariates for Nitrous oxide Model for Equation 2 \n    "), fig.lab.size = 16)
SI_Q3_no_inputs
ggsave('SI_Q3_no_inputs2.png', plot = SI_Q3_no_inputs, device = NULL, path = NULL, scale = 1, width = NA, height = NA, dpi = 300, limitsize = TRUE)

# emission factors --------------------------------------------------------

#calculing actual EF's for available data
EF_dat <- NULL 
for (i in unique(no_full$lat_decimal)){
  tempdf<-no_full %>% 
    filter(lat_decimal==i)  #subset each site by latitude
  tempdf2<-tempdf[tempdf$N_in_kg==0,] #isolate observations at N=0
  y0<-mean(tempdf2$NO_kg_seas_zero_rm)
  tempdf<-tempdf[tempdf$N_in_kg>0,] #removing x=0 points
  tempdf$EF<-((tempdf$NO_kg_seas_zero_rm-y0)/tempdf$N_in_kg)*100
  EF_dat<-rbind(EF_dat,tempdf)
}
EF_dat<-as.data.frame(EF_dat)
EF_dat<-dplyr::select(EF_dat,N_in_kg, NO_kg_seas_zero_rm, duration_d, 
                      EF,tropical,lat_decimal)

#finding minimum positive EF value to replace negative EFs for log-transformation
pos_EF<-EF_dat %>% filter(EF>0) #subsetting positive EFs
min(pos_EF$EF)
# minimum positive EF is 0.0327381; will replace non-positive EFs with half of this minimum for log-transformation
EF_dat_alt<-EF_dat$EF #making a vector of just EFs
EF_dat_alt[EF_dat_alt<=0] <- 0.0327381*.5 # replacing nonpositive EFs with half of minimum EF
EF_dat$EF_alt<-EF_dat_alt

#remove NAs to run model
dim(EF_dat)
EF_dat<-na.omit(EF_dat) #remove NAs
dim(EF_dat)

#z-transforming continuous covariates
EF_dat$N_in_kg.z<-z.trans(EF_dat$N_in_kg)
EF_dat$duration_d.z<-z.trans(EF_dat$duration_d)

#######################################################
######################Log-transformation approach:
########################################################

#finding minimum positive EF value to replace negative EFs for log-transformation
pos_EF<-EF_dat %>% filter(EF>0) #subsetting positive EFs
min(pos_EF$EF)
# minimum positive EF is 0.004; will replace non-positive EFs with half of this minimum for log-transformation
EF_dat_alt<-EF_dat$EF #making a vector of just EFs
EF_dat_alt[EF_dat_alt<=0] <- 0.004*.5 # replacing nonpositive EFs with half of minimum EF
EF_dat$EF_alt<-EF_dat_alt

#supplementary eq. 2 modeling EFs similar to as in eq. 1
mEF_seq2=lmer((log10(EF_alt))~  duration_d.z + tropical*N_in_kg.z + 
                (1|lat_decimal), data=EF_dat, na.action=NULL)
summary(mEF_seq2) 

##population prediction intervals to calculate CIs from p.257 Bolker 2008
x<- 1:1000
cis <- 0.95
lowb <- 0.5-(cis/2)
upb <- 0.5+(cis/2)
vmat<-mvrnorm(1000,mu=fixef(mEF_seq2), Sigma=vcov(mEF_seq2))

a<-vmat[,4] #coef for N fertilizer (if temperate)
b<-vmat[,2] #coef for study duration
c<-vmat[,1] #global intercept 
f<-vmat[,3] #coef for tropical=1
g<-vmat[,5] #coef for interaction N in if tropical

n_bar<-mean(EF_dat$N_in_kg, na.rm=T) #original mean for N in to unstandardize
n_sd<-sd(EF_dat$N_in_kg, na.rm=T) #original sd for N in to unstandardize

d_bar<-mean(EF_dat$duration_d, na.rm=T) #original mean for study duration to unstandardize
d_sd<-sd(EF_dat$duration_d, na.rm=T) #original sd for study duration to unstandardize

#function for CI curves
mEF_tropical_CIs=function(a,b,c,f,g,x){10^(c)*10^(f)*10^((a+g)*((x-n_bar)/(2*n_sd)))*b*((365-d_bar)/(2*d_sd))}
mEF_temperate_CIs=function(a,b,c,f,x) {10^(c)*10^((a)*((x-n_bar)/(2*n_sd)))*b*((365-d_bar)/(2*d_sd))}

#calculating CI values for tropical curve
dist=array(dim=c(1000,length(x)))
for (i in 1:1000) {
  dist[i,] = mEF_tropical_CIs(a=a[i],b=b[i],c=c[i],f=f[i],g=g[i],x=x)
}
civec_mEF_trop <- array(dim=c(2,length(x)))
for(j in 1:length(x)){
  civec_mEF_trop [,j] <- quantile(dist[,j],c(lowb,upb),na.rm=TRUE)
}

#calculating CI values for temperate curve
dist=array(dim=c(1000,length(x)))
for (i in 1:1000) {
  dist[i,] = mEF_temperate_CIs(a=a[i],b=b[i],c=c[i],f=f[i],x=x)
}
civec_mEF_temp <- array(dim=c(2,length(x)))
for(j in 1:length(x)){
  civec_mEF_temp [,j] <- quantile(dist[,j],c(lowb,upb),na.rm=TRUE)
}

 #supplemental eq. 2 fitted model outputs for plotting
a<-fixef(mEF_seq2)[4] #coef for N fertilizer (if temperate)
b<-fixef(mEF_seq2)[2] #coef for study duration
c<-fixef(mEF_seq2)[1] #global intercept 
f<-fixef(mEF_seq2)[3] #coef for tropical=1
g<-fixef(mEF_seq2)[5] #coef for interaction N in if tropical

##########plotting on linear scale
#functions to calculate tropical and temperate model (eq. 1) fit curves
mEF_tropical=function(x) {10^(c)*10^(f)*10^((a+g)*((x-n_bar)/(2*n_sd)))*b*((365-d_bar)/(2*d_sd))}
mEF_temperate=function(x) {10^(c)*10^((a)*((x-n_bar)/(2*n_sd)))*b*((365-d_bar)/(2*d_sd))}

#calculating curves and CIs
ytrop<-mEF_tropical(x)
ytemp<- mEF_temperate(x)
CItroplb<-civec_mEF_trop [1,]
CItropub<-civec_mEF_trop [2,]
CItemplb<-civec_mEF_temp [1,]
CItempub<-civec_mEF_temp [2,]
df<-data.frame(x,ytrop,ytemp,CItroplb, CItropub, CItemplb, CItempub)

#plot restricted range  fig. S5 "emission factor" ((y-y intercept)/x*100) ~x
x=1:1000
no_ef<-ggplot() +
  geom_point(data=EF_dat, aes(x=N_in_kg, y=EF,colour=factor(tropical)), alpha=.6,
             position = "jitter", size=2) +
  ylab("emission factor nitric oxide-N (%)") +
  xlab(expression(paste('N inputs (kg N ha'^'-1',
                        paste('season'^'-1'), ')'))) +
  theme_default() +
  geom_line(data=df, aes(x,ytrop), col='red3')+
  geom_line(data=df, aes(x,ytemp), col='dodgerblue') +
  geom_segment(aes(x=0,xend=0,y=-2,yend=25), colour="black") +
  ylim(-2,25) + xlim (0,300)+
  geom_segment(aes(x=0,xend=300,y=0,yend=0),colour="black") +
  scale_shape_manual(values=c(1,2)) +
  geom_ribbon(data = df, aes(x=x, y=ytrop,
                             ymin = CItroplb, ymax = CItropub),
              fill = "red3", alpha=.12) +
  geom_ribbon(data = df, aes(x=x, y=ytemp,
                             ymin = CItemplb, ymax = CItempub),
              fill = "dodgerblue", alpha=.12)

no_ef

#only run to update plots
#S5plot<-ggarrange(no3_ef,n2o_ef,no_ef, ncol=3, nrow=1, labels=c("a","b","c"))
#ggsave('S5plot.png', plot = S5plot, device = NULL, path = NULL, scale = 1.3, width = 10, height = 4, units = 'in', dpi = 300, limitsize = TRUE)


#exploring EFs calculated from fitted actual EFs as in SFig 5 specific points
x=c(5,10,25,50,100,150,200)
seq2_trop_EF<-mEF_tropical(x)
seq2_trop_EF

seq2_temp_EF<- mEF_temperate(x)
seq2_temp_EF

##############linear approaches

#linear regression approach based on the calculated EF's for each site
mEF_seq3=lm(EF_alt~ duration_d.z + tropical*N_in_kg.z, data=EF_dat, na.action=NULL)
summary(mEF_seq3)

#supplementary eq. 3 fitted model outputs for plotting
a<-coef(mEF_seq3)[4] #coef for N fertilizer (if temperate)
b<-coef(mEF_seq3)[2] #coef for duration
c<-coef(mEF_seq3)[1] #global intercept 
f<-coef(mEF_seq3)[3] #coef for tropical=1
g<-coef(mEF_seq3)[5] #coef for interaction N in if tropical

##########plotting on linear scale
#functions to calculate tropical and temperate model (eq. 1) fit curves; estimated for 365 days
mEF_seq3_tropical=function(x) {c+f+(a+g)*((x-n_bar)/(2*n_sd))+b*((365-d_bar)/(2*d_sd))}
mEF_seq3_temperate=function(x) {c+a*((x-n_bar)/(2*n_sd))+b*((365-d_bar)/(2*d_sd))}

#exploring EFs calculated from fitted actual EFs as in SFig 6 specific points
x=c(5,10,25,50,100,150,200)
mEF_seq3_trop_EF<-mEF_seq3_tropical(x)
mEF_seq3_trop_EF

mEF_seq3_temp_EF<-mEF_seq3_temperate(x)
mEF_seq3_temp_EF

plot(EF_dat$EF_alt~EF_dat$N_in_kg.z, ylim=c(0,200), 
     col=as.factor(EF_dat$tropical))
abline(c,a)
abline((c+f), (a+g), col='red')

###########################
#supplementary eq. 4--same approach as seq3 but without separating tropical and temperate data
#linear regression approach based on the calculated EF's for each site
mEF_seq4=lm(EF_alt~ duration_d.z + N_in_kg.z, data=EF_dat, na.action=NULL)
summary(mEF_seq4)


#supplemental eq 4 fitted model outputs for plotting
a<-coef(mEF_seq4)[3] #coef for N fertilizer (if temperate)
b<-coef(mEF_seq4)[2] #coef for duration
c<-coef(mEF_seq4)[1] #global intercept 

##########plotting on linear scale
plot(EF_dat$EF_alt~EF_dat$N_in_kg.z, ylim=c(0,200))
abline(c,a)

#functions to calculate tropical and temperate model (eq. 1) fit curves; estimated for 365 days
mEF_seq4=function(x) {c+a*((x-n_bar)/(2*n_sd))+b*((365-d_bar)/(2*d_sd))}

#exploring EFs calculated from fitted actual EFs as in SFig 5 specific points
x=c(5,10,25,50,100,150,200)
mEF_seq4_EF<-mEF_seq4(x)
mEF_seq4_EF

################################
#supplementary eq. 5 linear regression approach based on the overall N loss~ N input data
mEF_seq5=lm(NO_kg_seas_zero_rm~ tropical*N_in_kg + duration_d, data=no_full, na.action=NULL)
summary(mEF_seq5)

plot(no_full$NO_kg_seas_zero_rm~no_full$N_in_kg, ylim=c(0,200), 
     col=as.factor(no_full$tropical))
abline(coef(mEF_seq5)[1],coef(mEF_seq5)[3])
abline((coef(mEF_seq5)[1]+coef(mEF_seq5)[2]), 
       (coef(mEF_seq5)[3]+coef(mEF_seq5)[5]), col='red')

#################
#supplementary eq. 6 linear regression on the overall N loss~ N input data WITHOUT TEMP/TROP split
mEF_seq6=lm( NO_kg_seas_zero_rm~N_in_kg + duration_d, data=no_full, na.action=NULL)
summary(mEF_seq6)


mEF_seq6_pred=function(x) {c+a*x+b*365}

###########
#######this was copied from above########### is redundant?
#######
#plot FULL RANGE fig. S6 "emission factor" alternative calculation ((y-y intercept)/x*100) ~x
x=1:800
no_ef<-ggplot() +
  geom_point(data=EF_dat, aes(x=N_in_kg, y=EF,colour=factor(tropical)), alpha=.6,
             position = "jitter", size=2) +
  ylab("emission factor nitric oxide-N (%)") +
  xlab(expression(paste('N inputs (kg N ha'^'-1',
                        paste('season'^'-1'), ')'))) +
  theme_default() +
  geom_line(data=df, aes(x,(ytrop-ytrop[1])/x*100), col='red3')+
  geom_line(data=df, aes(x,(ytemp-ytemp[1])/x*100), col='dodgerblue') +
  geom_segment(aes(x=0,xend=0,y=-2,yend=15), colour="black") +
  geom_segment(aes(x=0,xend=800,y=0,yend=0),colour="black") +
  ylim(-2,15) + xlim (0,800)+
  scale_shape_manual(values=c(1,2)) 

no_ef
##############
##############

####################### Fig S6
#plot full range  fig. S6 based on seq6 with a shaded 95% confidence interval
x=1:900
d = data.frame(no_full, predict(mEF_seq6, interval="prediction"))

no_seq6<-ggplot(d,aes(x=N_in_kg,y=NO_kg_seas_zero_rm)) +
  geom_point(aes(colour=factor(tropical)), alpha=.6,position = "jitter", size=2) +
  geom_smooth( method="lm", colour="black") +
  ylab(expression(paste('Nitric oxide (kg NO-N',~ha^-1, season^-1,')    '))) +
  
  xlab(expression(paste('N inputs (kg N ha'^'-1',
                        paste('season'^'-1'), ')'))) +
  theme_default()+
  geom_segment(aes(x=0,xend=0,y=0,yend=20), colour="black") +
  geom_segment(aes(x=0,xend=400,y=0,yend=0),colour="black") +
  xlim(0,400)+ ylim(0,20)+
  scale_shape_manual(values=c(1,2))  

no_seq6

##################
#only re-run to update plot
#SF6plot<-ggarrange(no3_seq6,n2o_seq6,no_seq6, ncol=3, nrow=1, labels=c("a","b","c"))
#ggsave('SFig 6 plot.png', plot = SF7plot, device = NULL, path = NULL, scale = 1.3,
#       width = 10, height = 4, units = 'in', dpi = 300, limitsize = TRUE)

