#ammonia data analysis
setwd("~/N Loss Metaanalysis/analysis code and data")

library(lme4)
library(car)
library(MuMIn)
library(lmerTest)
library(MASS)
library(ggplot2)
library(dplyr)

#standardizing variables 2x sd per Gelman reccomendation
z.trans <- function(x) {
  (x - mean(x, na.rm = T)) / (2 * sd(x, na.rm = T))
}

#theme for plots
theme_default <- function(axis_text_size = 13) {
  theme(
    text = element_text(size = 16, colour = "black"),
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.text.x = element_text(size = axis_text_size, colour = "black"),
    axis.text.y = element_text(size = axis_text_size, colour = "black"),
    axis.title.y = element_text(angle = 90, vjust = -0.5),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    panel.background = element_rect(fill = "white", color = "white"),
    legend.position = "none"
  )
}
# eq. 1 ---------------------------------------------------------------------
nh3 <- read.csv("nh3_database.csv")

#set up dataframe for eq. 1
nh3_full <- nh3
nh3_model_full = data.frame(row.names = 1:nrow(nh3_full))
nh3_model_full$site_lat <-
  nh3_full$lat_decimal #no transformation needed
nh3_model_full$log_nh3_kg_season <-
  log10(nh3_full$nh3_kg_seas_zero.rm) #no transformation needed
nh3_model_full$tropical <-
  (nh3_full$tropical) #no transformation needed
nh3_model_full$N_in_kg <- z.trans(nh3_full$N_in_kg) #z-transforming
nh3_model_full$duration_d <-
  z.trans(nh3_full$duration_d) #z-transforming


#########eq. 1
meQ1 = lmer(
  (log_nh3_kg_season) ~  duration_d  + tropical * N_in_kg + (1 |
                                                               site_lat),
  data = nh3_model_full,
  na.action = NULL
)
summary(meQ1)
r.squaredGLMM(meQ1)

#population prediction intervals to calculate CIs from p.257 Bolker 2008
x <- 1:2000
cis <- 0.95
lowb <- 0.5 - (cis / 2)
upb <- 0.5 + (cis / 2)

vmat <- mvrnorm(1000, mu = fixef(meQ1), Sigma = vcov(meQ1))
a <- vmat[, 4] #coef for N fertilizer (if temperate)
b <- vmat[, 2] #coef for study duration
c <- vmat[, 1] #global intercept
f <- vmat[, 3] #coef for tropical=1
g <- vmat[, 5] #coef for interaction N in if tropical

n_bar <- mean(nh3_full$N_in_kg, na.rm = T)
n_sd <- sd(nh3_full$N_in_kg, na.rm = T)

#function for CI curves
yQ2_tropical_CIs = function(a, c, f, g, x) {
  10 ^ (c) * 10 ^ (f) * 10 ^ ((a + g) * ((x - n_bar) / (2 * n_sd)))
}
yQ2_temperate_CIs = function(a, c, f, g, x) {
  10 ^ (c) * 10 ^ ((a) * ((x - n_bar) / (2 * n_sd)))
}

#calculating CI values for tropical curve
dist = array(dim = c(1000, length(x)))
for (i in 1:1000) {
  dist[i, ] = yQ2_tropical_CIs(
    a = a[i],
    c = c[i],
    f = f[i],
    g = g[i],
    x = x
  )
}
civec_yQ2_trop <- array(dim = c(2, length(x)))
for (j in 1:length(x)) {
  civec_yQ2_trop [, j] <- quantile(dist[, j], c(lowb, upb), na.rm = TRUE)
}

#calculating CIs for temperate curve
dist = array(dim = c(1000, length(x)))
for (i in 1:1000) {
  dist[i, ] = yQ2_temperate_CIs(
    a = a[i],
    c = c[i],
    f = f[i],
    g = g[i],
    x = x
  )
}
civec_yQ2_temp <- array(dim = c(2, length(x)))
for (j in 1:length(x)) {
  civec_yQ2_temp [, j] <- quantile(dist[, j], c(lowb, upb), na.rm = TRUE)
}

#eq. 1 fitted model outputs for plotting
a <- fixef(meQ1)[4] #coef for N fertilizer (if temperate)
b <- fixef(meQ1)[2]#coef for study duration
c <- fixef(meQ1)[1] #global intercept
f <- fixef(meQ1)[3]#coef for tropical=1
g <- fixef(meQ1)[5] #coef for interaction N in if tropical

##########plotting on linear scale
#functions to calculate tropical and temperate model (eq. 1) fit curves
yQ2_tropical = function(x) {
  10 ^ (c) * 10 ^ (f) * 10 ^ ((a + g) * ((x - n_bar) / (2 * n_sd)))
}
yQ2_temperate = function(x) {
  10 ^ (c) * 10 ^ ((a * ((x - n_bar) / (2 * n_sd))))
}

# reordering tropical and temperate levels
nh3_full$tropical <-
  factor(as.factor(nh3_full$tropical), rev(levels(as.factor(nh3_full$tropical))))
levels(as.factor(nh3_full$tropical)) #tropical now comes first

#calculating curves and CIs
ytrop <- yQ2_tropical(x)
ytemp <- yQ2_temperate(x)
CItroplb <- civec_yQ2_trop [1, ]
CItropub <- civec_yQ2_trop [2, ]
CItemplb <- civec_yQ2_temp [1, ]
CItempub <- civec_yQ2_temp [2, ]
df <- data.frame(x, ytrop, ytemp, CItroplb, CItropub, CItemplb, CItempub)

#Figure 2
nh3p <- ggplot(nh3_full, aes(x = N_in_kg, y = nh3_kg_seas_zero.rm)) +
  geom_point(
    aes(colour = factor(tropical), shape = factor(tropical)),
    alpha = .55,
    position = "jitter",
    size = 2
  ) +
  ylab(expression(paste(
    'Ammonia (kg ', ' ', NH[3], '-N',  ~ ha ^ -1, season ^ -1, ')    '
  ))) +
  xlab(expression(paste(
    'N inputs (kg N ha' ^ '-1',
    paste('season' ^ '-1'), ')'
  ))) +
  ggtitle('Ammonia losses') +
  theme_default() +
  ylim(0, 40) + xlim (0, 300) +
  geom_segment(aes(
    x = 0,
    xend = 0,
    y = 0,
    yend = 40
  ), colour = "black") +
  geom_segment(aes(
    x = 0,
    xend = 300,
    y = 0,
    yend = 0
  ), colour = "black") +
  scale_shape_manual(values = c(16, 17)) +
  geom_line(data = df,
            aes(x, ytrop),
            col = 'red3',
            size = 1) +
  geom_line(data = df,
            aes(x, ytemp),
            col = 'dodgerblue',
            size = 1) +
  geom_ribbon(
    data = df,
    aes(
      x = x,
      y = ytrop,
      ymin = CItroplb,
      ymax = CItropub
    ),
    fill = "red3",
    alpha = .12
  ) +
  geom_ribbon(
    data = df,
    aes(
      x = x,
      y = ytemp,
      ymin = CItemplb,
      ymax = CItempub
    ),
    fill = "dodgerblue",
    alpha = .12
  ) +
  theme(legend.position = "right", legend.key = element_blank()) +
  scale_color_manual(
    name = 'Region',
    labels = c('Tropical', 'Temperate'),
    values = c('red3', 'dodgerblue')
  ) +
  scale_shape_manual(
    name = 'Region',
    labels = c('Tropical', 'Temperate'),
    values = c(16, 17)
  )
nh3p


# eq. 2 -------------------------------------------------------------------

nh3_trop <- nh3 %>% filter(tropical > 0)

#removing empty crop type factor
nh3_trop <- subset(nh3_trop, crop_type != "")
table((nh3_trop$crop_type))
nh3_trop$crop_type <- droplevels(nh3_trop$crop_type)
table((nh3_trop$crop_type))

#set up dataframe for eq. 2 model
nh3_trop$log_g_nh3_day <-
  log10(nh3_trop$nh3_kg_seas_zero.rm * 1000 / nh3_trop$duration_d)
nh3_model_tropical <-
  data.frame(row.names = 1:length(nh3_trop$nh3_kg_seas_zero.rm))
nh3_model_tropical$site_lat <-
  nh3_trop$lat_decimal # site ID; no transformation needed
nh3_model_tropical$log_g_nh3_day <-
  nh3_trop$log_g_nh3_day #no transformation needed
nh3_model_tropical$crop_type <-
  nh3_trop$crop_type #no transformation needed
nh3_model_tropical$irrigation <-
  nh3_trop$irrig_1_0 #irrigation=1 (no trans. needed)
nh3_model_tropical$split_app <-
  nh3_trop$split_app #split=1 (no trans. needed)
nh3_model_tropical$org_Nfix_in <-
  nh3_trop$org_Nfix_in #organic/Nfix in=1  (no trans. needed)
nh3_model_tropical$N_in_kg <- z.trans(nh3_trop$N_in_kg)
nh3_model_tropical$sand <-
  z.trans(nh3_trop$pct_sand_combined) #z-transforming
nh3_model_tropical$precip <-
  z.trans(nh3_trop$MAP_mm_combined) #z-transforming
nh3_model_tropical$duration_d_temp <-
  z.trans(nh3_trop$duration_d) #z-transforming
nh3_model_tropical$pH <- z.trans(nh3_trop$pH) #z-transforming

#remove incomplete cases for lmer function
dim(nh3_model_tropical) #193 10
nh3_model_tropical <- na.omit(nh3_model_tropical)
dim(nh3_model_tropical) #removed incomplete cases

#test for collinearity
vif_test = lmer(
  log_g_nh3_day ~ N_in_kg + sand +  precip + irrigation + split_app +
    org_Nfix_in + pH + (1 |
                          site_lat),
  data = nh3_model_tropical,
  na.action = NULL
)
vif(vif_test) #all covariates below 3

meq2 = lmer(
  log_g_nh3_day ~ N_in_kg + sand + precip + irrigation + split_app +
    org_Nfix_in  + pH + crop_type + (1 |
                                       site_lat),
  data = nh3_model_tropical,
  na.action = NULL
)

summary(meq2)
r.squaredGLMM(meq2)

#Fig. 3
x = 1:100
sand <-
  mean(nh3_trop$pct_sand_combined, na.rm = T) #to unstandardize sand
sand_sd <-
  sd(nh3_trop$pct_sand_combined, na.rm = T) #to unstandardize sand
sand_fit_SE <- sqrt(diag(vcov(meq2)))[3] #SE amounts
ysand <-
  function(x) {
    fixef(meq2)[1] + fixef(meq2)[3] * ((x - sand) / (2 * sand_sd))
  }
ysand_hat <- ysand(x = x)
ysand_hat_lb <- ysand_hat - rep(sand_fit_SE, length(ysand_hat))
ysand_hat_ub <- ysand_hat + rep(sand_fit_SE, length(ysand_hat))
df <- data.frame(x, ysand_hat, ysand_hat_lb, ysand_hat_ub)

#sand plot
nh3_sand <-
  ggplot(nh3_trop, aes(x = pct_sand_combined, y = log_g_nh3_day)) +
  geom_point(aes(alpha = .6)) +
  ylab('') +
  xlab(expression(paste('soil texture (% sand)'))) +
  theme_default() +
  xlim (0, 100) +
  geom_hline(yintercept = -2, color = "black") +
  geom_vline(xintercept = 0, color = "black") +
  geom_line(data = df, aes(x, ysand_hat), col = 'red3') +
  geom_ribbon(
    data = df,
    aes(
      x = x,
      y = ysand_hat,
      ymin = ysand_hat_lb,
      ymax = ysand_hat_ub
    ),
    fill = "red3",
    alpha = .12
  ) +
  scale_y_continuous(
    labels = c(0.0001, 0.01, 1, 100, '10,000', 10 ^ 6, 10 ^ 8, 10 ^ 10),
    breaks = c(-4, -2, 0, 2, 4, 6, 8, 10),
    limits = c(-2, 5)
  )

nh3_sand

#precipitation/irrigation input plot
x = 700:3300
precip <- mean(nh3_trop$MAP_mm_combined, na.rm = T) #to unstandardize
precip_sd <- sd(nh3_trop$MAP_mm_combined, na.rm = T) #to unstandardize
precip_fit_SE <- sqrt(diag(vcov(meq2)))[4] #SE amounts
yprecip <-
  function(x) {
    fixef(meq2)[1] + fixef(meq2)[4] * ((x - precip) / (2 * precip_sd))
  }
yprecip_hat <- yprecip(x = x)
yprecip_hat_lb <- yprecip_hat - rep(precip_fit_SE, length(yprecip_hat))
yprecip_hat_ub <- yprecip_hat + rep(precip_fit_SE, length(yprecip_hat))
df <- data.frame(x, yprecip_hat, yprecip_hat_lb, yprecip_hat_ub)

#precip plot
nh3_precip <-
  ggplot(nh3_trop, aes(x = MAP_mm_combined, y = log_g_nh3_day)) +
  geom_point(aes(alpha = .6)) +
  ylab(expression(paste('ammonia (g ',  ~ NH[3], ' - N',  ~ ha ^ -1,  ~
                          d ^ -1, ')'))) +
  xlab(expression(paste('mean annual precipitation (mm)'))) +
  theme_default() +
  geom_line(data = df, aes(x, yprecip_hat), col = 'red3') +
  geom_ribbon(
    data = df,
    aes(
      x = x,
      y = yprecip_hat,
      ymin = yprecip_hat_lb,
      ymax = yprecip_hat_ub
    ),
    fill = "red3",
    alpha = .12
  ) +
  xlim (700, 4600) +
  geom_hline(yintercept = -2, color = "black") +
  geom_vline(xintercept = 700, color = "black") +
  scale_y_continuous(
    labels = c(0.0001, 0.01, 1, 100, '10,000', 10 ^ 6, 10 ^ 8, 10 ^ 10),
    breaks = c(-4, -2, 0, 2, 4, 6, 8, 10),
    limits = c(-2, 5)
  )
nh3_precip

#pH plot
x = c(4, 5, 6, 7.35)
pH <- mean(nh3_trop$pH, na.rm = T) #to unstandardize
pH_sd <- sd(nh3_trop$pH, na.rm = T) #to unstandardize
pH_fit_SE <- sqrt(diag(vcov(meq2)))[8] #SE amounts
ypH <- function(x) {
  fixef(meq2)[1] + fixef(meq2)[8] * ((x - pH) / (2 * pH_sd))
}
ypH_hat <- ypH(x = x)
ypH_hat_lb <- ypH_hat - rep(pH_fit_SE, length(ypH_hat))
ypH_hat_ub <- ypH_hat + rep(pH_fit_SE, length(ypH_hat))
df <- data.frame(x, ypH_hat, ypH_hat_lb, ypH_hat_ub)

#pH plot
nh3_pH <- ggplot(nh3_trop, aes(x = pH, y = log_g_nh3_day)) +
  geom_point(aes(alpha = .6)) +
  ylab(" ") +
  xlab(expression(paste('pH'))) +
  theme_default() +
  geom_hline(yintercept = -2, color = "black") +
  geom_vline(xintercept = 4, color = "black") +
  geom_line(data = df, aes(x, ypH_hat), col = 'red3') +
  geom_ribbon(
    data = df,
    aes(
      x = x,
      y = ypH_hat,
      ymin = ypH_hat_lb,
      ymax = ypH_hat_ub
    ),
    fill = "red3",
    alpha = .12
  ) +
  scale_y_continuous(
    labels = c(0.0001, 0.01, 1, 100, '10,000', 10 ^ 6, 10 ^ 8, 10 ^ 10),
    breaks = c(-4, -2, 0, 2, 4, 6, 8, 10),
    limits = c(-2, 5)
  ) +
  scale_x_continuous(breaks = c(4, 5, 6, 7, 8), limits = c(4, 8.5))

nh3_pH


#Fig 4 binary/dummy crops plot
coef <-
  data.frame(data = c(fixef(meq2)[c(
    '(Intercept)',
    '(Intercept)',
    'split_app',
    '(Intercept)',
    'org_Nfix_in',
    '(Intercept)',
    'irrigation',
    'crop_typeother',
    'crop_typeCF rice'
  )]) +
    c(
      rep(0, 2),
      fixef(meq2)['(Intercept)'],
      0,
      fixef(meq2)['(Intercept)'],
      0,
      fixef(meq2)['(Intercept)'],
      fixef(meq2)['(Intercept)'],
      fixef(meq2)['(Intercept)']
    ))
#adding the intercept (cereal) to all model fits

meq2.se <- sqrt(diag(vcov(meq2))) #standard errors from meq2
meq2.names <- row.names(vcov(meq2)) #variable names of meq2
meq2.se.comb <-
  as.data.frame(cbind(meq2.names, meq2.se)) #data frame of parameter SEs
meq2.se.comb <-
  data.frame(meq2.names, meq2.se) #data frame of parameter SEs

SE <- c(
  meq2.se.comb[meq2.names == '(Intercept)', 2],
  meq2.se.comb[meq2.names == '(Intercept)', 2],
  meq2.se.comb[meq2.names == 'split_app', 2],
  meq2.se.comb[meq2.names == '(Intercept)', 2],
  meq2.se.comb[meq2.names == 'org_Nfix_in', 2],
  meq2.se.comb[meq2.names == '(Intercept)', 2],
  meq2.se.comb[meq2.names == 'irrigation', 2],
  meq2.se.comb[meq2.names == 'crop_typeother', 2, ],
  meq2.se.comb[meq2.names == 'crop_typeCF rice', 2, ]
)

meq2_outputs <- cbind(coef, SE)
colnames(meq2_outputs) <- c('coef', 'SE')
meq2_outputs$names <-
  c('cereal',
    'sa0',
    'sa1',
    'org0',
    'org1',
    'irr0',
    'irr1',
    'other',
    'CF rice')
meq2_outputs

#new plot with y data and model fit
#creating individual dfs for binary variables to join to crop types
discrete_var <- nh3_trop %>%
  dplyr::select(org_Nfix_in , split_app, irrig_1_0, log_g_nh3_day)
org1 <- discrete_var %>% filter(org_Nfix_in == 1) %>%
  mutate(crop_type = rep("org1", 83))
org0 <- discrete_var %>% filter(org_Nfix_in == 0) %>%
  mutate(crop_type = rep("org0", 110))
sa1 <- discrete_var %>% filter(split_app == 1) %>%
  mutate(crop_type = rep("sa1", 49))
sa0 <- discrete_var %>% filter(split_app == 0) %>%
  mutate(crop_type = rep("sa0", 144))
irr1 <- discrete_var %>% filter(irrig_1_0 == 1) %>%
  mutate(crop_type = rep("irr1", 68))
irr0 <- discrete_var %>% filter(irrig_1_0 == 0) %>%
  mutate(crop_type = rep("irr0", 125))

#joining the data frames
discrete_var2 <- nh3_trop %>%
  dplyr::select(org_Nfix_in , split_app, irrig_1_0, crop_type, log_g_nh3_day)

discrete_var3 <- rbind(discrete_var2, org0, org1, sa1, sa0, irr0, irr1)


#reordering the factors based on the model outputs
levels(discrete_var3$crop_type)
discrete_var3$crop_type <-
  factor(
    discrete_var3$crop_type,
    levels = c(
      'cereal',
      'sa0',
      'sa1',
      'org0',
      'org1',
      'irr0',
      'irr1',
      'fallow',
      'other',
      'CF rice'
    )
  )


#x-axis names
x_axis <-
  c(
    'cereal/intercept',
    'split application=0',
    'split application=1',
    'organic/N-fix in=0',
    'organic/N-fix in=1',
    'irrigation=0',
    'irrigation=1',
    'other',
    'continuously flooded rice'
  )

#fig. 4
nh3_Q3_disc <- ggplot() +
  geom_violin(
    data = discrete_var3,
    aes(x = crop_type, y = log_g_nh3_day),
    fill = "light gray",
    colour = "light gray"
  ) +
  geom_jitter(
    data = discrete_var3,
    aes(x = crop_type, y = log_g_nh3_day),
    width = .05,
    alpha = .3
  ) +
  ylab(expression(paste('ammonia (g ',  ~ NH[3], ' - N',  ~ ha ^ -1,  ~
                          d ^ -1, ')'))) +
  geom_pointrange(data = meq2_outputs,
                  aes(
                    reorder(names, coef),
                    coef,
                    ymin = coef - SE,
                    ymax =  coef + SE
                  ),
                  col = 'red') +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  scale_x_discrete('' , labels = x_axis) +
  geom_segment(aes(
    x = 0,
    xend = 0,
    y = -2,
    yend = 4
  ), colour = "black") +
  geom_hline(yintercept = -2, color = "black") +
  theme_default(axis_text_size = 13) +
  annotate(
    "text",
    x = 7,
    y = 4.7,
    size = 8,
    label = "**"
  ) +
  geom_vline(xintercept = c(1.5, 3.5, 5.5, 7.5), linetype = 'longdash')



nh3_Q3_disc <-
  nh3_Q3_disc + scale_y_continuous(breaks = c(-2, 0, 2, 4),
                                   labels = c(0.01, 1, 100, '10,000'))
nh3_Q3_disc
