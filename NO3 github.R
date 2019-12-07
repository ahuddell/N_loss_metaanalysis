#Nitrate leaching data analysis
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
no3 <- read.csv("no3_database.csv")

#set up dataframe for eq. 1
no3_full <- no3
NO3_model_full = data.frame(row.names = 1:nrow(no3_full))
NO3_model_full$site_lat <-
  no3_full$lat_decimal #no transformation needed
NO3_model_full$log_NO3_kg_season <-
  log10(no3_full$NO3_kg_seas_zero_rm) #no transformation needed
NO3_model_full$N_in_kg <- (no3_full$N_in_kg)
NO3_model_full$duration_d <- (no3_full$duration_d)
NO3_model_full$sample.depth <- (no3_full$sample_cm)
NO3_model_full$tropical <- (no3_full$tropical)
NO3_model_full$N_in_kg <-
  z.trans(NO3_model_full$N_in_kg) #z-transforming
NO3_model_full$duration_d <-
  z.trans(NO3_model_full$duration_d) #z-transforming
NO3_model_full$sample.depth <-
  z.trans(no3_full$sample_cm) #z-transforming

#########eq. 1
meq1 = lmer((log_NO3_kg_season) ~  duration_d + tropical * N_in_kg + sample.depth + (1 |
                                                                                       site_lat),
            data = NO3_model_full,
            na.action = NULL
)
summary(meq1)
r.squaredGLMM(meq1)

#population prediction intervals to calculate CIs from p.257 Bolker 2008
x <- 0:1000
cis <- 0.95
lowb <- 0.5 - (cis / 2)
upb <- 0.5 + (cis / 2)

vmat <- mvrnorm(1000, mu = fixef(meq1), Sigma = vcov(meq1))
a <- vmat[, 4] #coef for N fertilizer (if temperate)
c <- vmat[, 1] #global intercept
f <- vmat[, 3] #coef for tropical=1
g <- vmat[, 6] #coef for interaction N in if tropical

n_bar <-
  mean(no3_full$N_in_kg, na.rm = T) #original mean for N in to unstandardize
n_sd <-
  sd(no3_full$N_in_kg, na.rm = T) #original sd for N in to unstandardize

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

#calculating CI values for temperate curve
dist = array(dim = c(1000, length(x)))
for (i in 1:1000) {
  dist[i, ] = yQ2_temperate_CIs(a = a[i],
                                c = c[i],
                                f = f[i],
                                x = x)
}
civec_yQ2_temp <- array(dim = c(2, length(x)))
for (j in 1:length(x)) {
  civec_yQ2_temp [, j] <- quantile(dist[, j], c(lowb, upb), na.rm = TRUE)
}

#eq. 1 fitted model outputs for plotting
a <- fixef(meq1)[4] #coef for N fertilizer (if temperate)
b <- fixef(meq1)[2] #coef for study duration
c <- fixef(meq1)[1] #global intercept
f <- fixef(meq1)[3] #coef for tropical=1
g <- fixef(meq1)[6] #coef for interaction N in if tropical
h <- fixef(meq1)[5] #coef for sample depth

##########plotting on linear scale
#functions to calculate tropical and temperate model (eq. 1) fit curves
yQ2_tropical = function(x) {
  10 ^ (c) * 10 ^ (f) * 10 ^ ((a + g) * ((x - n_bar) / (2 * n_sd)))
}
yQ2_temperate = function(x) {
  10 ^ (c) * 10 ^ ((a) * ((x - n_bar) / (2 * n_sd)))
}

# reordering tropical and temperate levels
no3_full$tropical <-
  factor(as.factor(no3_full$tropical), rev(levels(as.factor(no3_full$tropical))))
levels(as.factor(no3_full$tropical)) #tropical now comes first

#calculating curves and CIs
ytrop <- yQ2_tropical(x)
ytemp <- yQ2_temperate(x)
CItroplb <- civec_yQ2_trop [1, ]
CItropub <- civec_yQ2_trop [2, ]
CItemplb <- civec_yQ2_temp [1, ]
CItempub <- civec_yQ2_temp [2, ]
df <- data.frame(x, ytrop, ytemp, CItroplb, CItropub, CItemplb, CItempub)

#Figure 2
no3p <- ggplot(no3_full, aes(x = N_in_kg, y = NO3_kg_seas_zero_rm)) +
  geom_point(
    aes(colour = factor(tropical), shape = factor(tropical)),
    alpha = .6,
    position = "jitter",
    size = 2
  ) +
  ylab(expression(paste(
    "Nitrate (kg N ", O[3], '' ^ '-',
    ' - N',  ~ ha ^ -1,  ~ season ^ -1, ')'
  ))) +
  xlab(expression(paste(
    'N inputs (kg N ha' ^ '-1',
    paste('season' ^ '-1'), ')'
  ))) +
  ggtitle('Nitrate leaching losses') +
  theme_default() +
  ylim(0, 60) + xlim (0, 300) +
  geom_segment(aes(
    x = 0,
    xend = 0,
    y = 0,
    yend = 60
  ), colour = "black") +
  geom_segment(aes(
    x = 0,
    xend = 300,
    y = 0,
    yend = 0
  ), colour = "black") +
  geom_line(data = df, aes(x, ytrop), col = 'red3') +
  geom_line(data = df, aes(x, ytemp), col = 'dodgerblue') +
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
  scale_color_manual(
    name = 'Region',
    labels = c('Tropical', 'Temperate'),
    values = c('red3', 'dodgerblue')
  ) +
  scale_shape_manual(
    name = 'Region',
    labels = c('Tropical', 'Temperate'),
    values = c(16, 17)
  ) +
  theme(legend.position = "right", legend.key = element_blank())
no3p

# eq. 2 -------------------------------------------------------------------
no3_trop <-
  no3 %>% filter(tropical > 0) #removing temperate comparison data
levels(no3_trop$crop_type)

#removing empty crop type factor
no3_trop <- subset(no3_trop, crop_type != "")
table((no3_trop$crop_type))
no3_trop$crop_type <- droplevels(no3_trop$crop_type)
table((no3_trop$crop_type))

#set up dataframe for eq. 2 model
NO3_model_tropical = data.frame(row.names = 1:nrow(no3_trop))
NO3_model_tropical$site_lat <-
  as.factor(no3_trop$lat_decimal) #no transformation needed
NO3_model_tropical$crop_type <-
  no3_trop$crop_type #no transformation needed
NO3_model_tropical$log_NO3_g_day <-
  log10(no3_trop$NO3_kg_seas_zero_rm * 1000 / no3_trop$duration_d) #no transformation needed
NO3_model_tropical$log_NO3_kg_season <-
  log10(no3_trop$NO3_kg_seas_zero_rm)
NO3_model_tropical$irrigation <-
  no3_trop$irrig_1_0 #irrigation=1 (no trans. needed)
NO3_model_tropical$split_app <-
  no3_trop$split_app #split=1 (no trans. needed)
NO3_model_tropical$pH <- no3_trop$pH #pH
NO3_model_tropical$org_Nfix_in <-
  no3_trop$org_Nfix_in #organic/Nfix in=1  (no trans. needed)
NO3_model_tropical$N_in_kg <-
  z.trans(no3_trop$N_in_kg) #z-transforming
NO3_model_tropical$duration_d <-
  z.trans(no3_trop$duration_d) #z-transforming
NO3_model_tropical$sand <-
  z.trans(no3_trop$pct_sand_combined) #z-transforming
NO3_model_tropical$water_input_mm <-
  z.trans(no3_trop$water_input_mm) #z-transforming
NO3_model_tropical$pH <- z.trans(no3_trop$pH) #z-transforming
NO3_model_tropical$sample_depth <-
  z.trans(no3_trop$sample_cm) #z-transforming

#remove incomplete cases for lmer function
dim(NO3_model_tropical)
NO3_model_tropical <- na.omit(NO3_model_tropical)
dim(NO3_model_tropical) #removed imcomplete cases

#test for collinearity
vif_test = lmer((log_NO3_g_day) ~ N_in_kg + sand + water_input_mm + irrigation
                + split_app + org_Nfix_in + pH + sample_depth + (1 |
                                                                   site_lat),
                data = NO3_model_tropical,
                na.action = NULL
)

vif(vif_test) #all covariates below 3

# removed irrigation  + added sample depth
meq2 = lmer((log_NO3_g_day) ~ N_in_kg + sand + water_input_mm  + sample_depth + split_app
            + org_Nfix_in + pH + crop_type + (1 |
                                                site_lat) ,
            data = NO3_model_tropical,
            na.action = NULL
)

summary(meq2)
r.squaredGLMM(meq2)

#fig. 3
#inputs for graph
x = 1:100 #x
no3_trop$log_NO3_g_day <-
  log10(no3_trop$NO3_kg_seas_zero_rm * 1000 / no3_trop$duration_d) #y
sand <-
  mean(no3_trop$pct_sand_combined, na.rm = T) #to unstandardize sand
sand_sd <-
  sd(no3_trop$pct_sand_combined, na.rm = T) #to unstandardize sand
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
no3_sand <-
  ggplot(no3_trop, aes(x = pct_sand_combined, y = log_NO3_g_day)) +
  geom_point(aes(alpha = .6)) +
  ylab(" ") +
  xlab(expression(paste('soil texture (% sand)'))) +
  annotate(
    "text",
    x = 85,
    y = 3.8,
    size = 14,
    label = "**"
  ) +
  theme_default() +
  ylim(0, 4) + xlim (0, 100) +
  geom_segment(aes(
    x = 0,
    xend = 0,
    y = 0,
    yend = 4
  ), colour = "black") +
  geom_segment(aes(
    x = 0,
    xend = 100,
    y = 0,
    yend = 0
  ), colour = "black") +
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
    labels = c(1, 10, 100, '1,000', '10,000'),
    breaks = c(0, 1, 2, 3, 4),
    limits = c(0, 4)
  )
no3_sand

#precipitation/irrigation input plot
water_input_mm <- fixef(meq2)[4]
water_input_mm_SE <- sqrt(diag(vcov(meq2)))[4]
x = 0:3000
precip <-
  mean(no3_trop$water_input_mm, na.rm = T) #to unstandardize sand
precip_sd <-
  sd(no3_trop$water_input_mm, na.rm = T) #to unstandardize sand
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
no3_precip <-
  ggplot(no3_trop, aes(x = water_input_mm, y = log_NO3_g_day)) +
  geom_point(aes(alpha = .6)) +
  ylab(expression(paste(
    'nitrate (g N', O[3], '' ^ '-',
    ' - N',  ~ ha ^ -1,  ~ d ^ -1, ')'
  ))) +
  xlab(expression(paste('precipitation and irrigation inputs (mm)'))) +
  theme_default() +
  annotate(
    "text",
    x = 2600,
    y = 3.8,
    size = 14,
    label = "**"
  ) +
  xlim (0, 4600) +
  geom_hline(yintercept = 0, color = "black") +
  geom_vline(xintercept = 0, color = "black") +
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
  scale_y_continuous(
    labels = c(1, 10, 100, '1,000', '10,000'),
    breaks = c(0, 1, 2, 3, 4),
    limits = c(0, 4)
  )
no3_precip

#pH plot
x = c(4, 5, 6, 7.45)
pH <- mean(no3_trop$pH, na.rm = T) #to unstandardize
pH_sd <- sd(no3_trop$pH, na.rm = T) #to unstandardize
pH_fit_SE <- sqrt(diag(vcov(meq2)))[8] #SE amounts
ypH <- function(x) {
  fixef(meq2)[1] + fixef(meq2)[8] * ((x - pH) / (2 * pH_sd))
}
ypH_hat <- ypH(x = x)
ypH_hat_lb <- ypH_hat - rep(pH_fit_SE, length(ypH_hat))
ypH_hat_ub <- ypH_hat + rep(pH_fit_SE, length(ypH_hat))
df <- data.frame(x, ypH_hat, ypH_hat_lb, ypH_hat_ub)

#pH plot
no3_pH <- ggplot(no3_trop, aes(x = pH, y = log_NO3_g_day)) +
  geom_point(aes(alpha = .6)) +
  ylab(" ") +
  xlab(expression(paste('pH'))) +
  theme_default() +
  geom_hline(yintercept = 0, color = "black") +
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
    labels = c(1, 10, 100, '1,000', '10,000'),
    breaks = c(0, 1, 2, 3, 4),
    limits = c(0, 4)
  ) +
  scale_x_continuous(breaks = c(4, 5, 6, 7, 8), limits = c(4, 8.5))
no3_pH

#Fig 4 binary/dummy crops plot
coef <-
  data.frame(data = c(fixef(meq2)[c(
    '(Intercept)',
    '(Intercept)',
    'split_app',
    '(Intercept)',
    'org_Nfix_in',
    'crop_typefallow',
    'crop_typeother',
    'crop_typepasture',
    'crop_typerice',
    'crop_typetree crop'
  )]) +
    c(rep(0, 2), fixef(meq2)['(Intercept)'], 0,
      rep(fixef(meq2)['(Intercept)'], 6)))
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
  meq2.se.comb[meq2.names == 'crop_typefallow', 2],
  meq2.se.comb[meq2.names == 'crop_typeother', 2],
  meq2.se.comb[meq2.names == 'crop_typepasture', 2],
  meq2.se.comb[meq2.names == 'crop_typerice', 2],
  meq2.se.comb[meq2.names == 'crop_typetree crop', 2]
)

meq2_outputs <- cbind(coef, SE)
colnames(meq2_outputs) <- c('coef', 'SE')
meq2_outputs$names <- c(
  'cereal',
  'sa0',
  'sa1',
  'org0',
  'org1',
  'fallow',
  'other',
  'pasture',
  'rice',
  'tree crop'
)
meq2_outputs

#creating individual dfs for binary variables to join to crop types
discrete_var <- no3_trop %>%
  dplyr::select(org_Nfix_in , split_app, log_NO3_g_day)
org1 <- discrete_var %>% filter(org_Nfix_in == 1) %>%
  mutate(crop_type = rep("org1", 65))
org0 <- discrete_var %>% filter(org_Nfix_in == 0) %>%
  mutate(crop_type = rep("org0", 128))
sa1 <- discrete_var %>% filter(split_app == 1) %>%
  mutate(crop_type = rep("sa1", 110))
sa0 <- discrete_var %>% filter(split_app == 0) %>%
  mutate(crop_type = rep("sa0", 83))

#joining the data frames
discrete_var2 <- no3_trop %>%
  dplyr::select(org_Nfix_in , split_app, crop_type, log_NO3_g_day)

discrete_var3 <- rbind(discrete_var2, org0, org1, sa1, sa0)

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
      'fallow',
      'other',
      'pasture',
      'rice',
      'tree crop'
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
    'fallow',
    'other',
    'pasture',
    'rice',
    'tree crop'
  )

#fig 4
no3_Q3_disc <- ggplot() +
  geom_violin(
    data = discrete_var3,
    aes(x = crop_type, y = log_NO3_g_day),
    fill = 'light gray',
    colour = 'light gray'
  ) +
  geom_jitter(
    data = discrete_var3,
    aes(x = crop_type, y = log_NO3_g_day),
    width = .05,
    alpha = .3
  ) +
  ylab(expression(paste(
    ' nitrate (g N ', O[3], '' ^ '-',
    ' - N',  ~ ha ^ -1,  ~ d ^ -1, ')'
  ))) +
  geom_pointrange(data = meq2_outputs,
                  aes(
                    reorder(names, coef),
                    coef,
                    ymin = coef - SE,
                    ymax =  coef + SE
                  ),
                  col = 'red') +
  scale_x_discrete(" ", labels = x_axis) +
  geom_hline(yintercept = 0, color = "black") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  geom_segment(aes(
    x = 0,
    xend = 0,
    y = 0,
    yend = 4
  ), colour = "black") +
  theme_default(axis_text_size = 13) +
  geom_vline(xintercept = c(1.5, 3.5, 5.5), linetype = 'longdash')

no3_Q3_disc <-
  no3_Q3_disc +  scale_y_continuous(labels = c('0', '10', '100', '1,000', '10,000'))
no3_Q3_disc
