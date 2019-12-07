#Nitric oxide data analysis
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
# eq. 1 ----------------------------------------------------------------------
no <- read.csv("no_database.csv")

# set up dataframe for eq. 1
no_full <- no
NO_model_full = data.frame(row.names = 1:nrow(no_full))
NO_model_full$site_lat <- no_full$lat_decimal #no trans. needed
NO_model_full$tropical <- (no_full$tropical) #no trans. needed
NO_model_full$log_no_season <-
  log10(no_full$NO_kg_seas_zero_rm) #logging response variable
NO_model_full$N_in_kg_untrans <- (no_full$N_in_kg) #z-transforming
NO_model_full$N_in_kg <- z.trans(no_full$N_in_kg) #z-transforming
NO_model_full$duration_d <-
  z.trans(no_full$duration_d) #z-transforming

#########eq. 1
meQ1 = lmer(
  (log_no_season) ~  duration_d  +  N_in_kg * tropical + (1 |
                                                            site_lat),
  data = NO_model_full,
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
a <- vmat[, 3] #coef for N fertilizer (if temperate)
b <- vmat[, 2] #coef for study duration
c <- vmat[, 1] #global intercept
f <- vmat[, 4] #coef for tropical=1
g <- vmat[, 5] #coef for interaction N in if tropical

n_bar <- mean(no_full$N_in_kg, na.rm = T)
n_sd <- sd(no_full$N_in_kg, na.rm = T)

#function for CI curves
yQ2_tropical_CIs = function(a, b, c, f, g, x) {
  10 ^ (c) * 10 ^ (f) * 10 ^ ((a + g) * ((x - n_bar) / (2 * n_sd)))
}
yQ2_temperate_CIs = function(a, b, c, f, g, x) {
  10 ^ (c) * 10 ^ ((a) * ((x - n_bar) / (2 * n_sd)))
}

#calculating CIs for tropical curve
dist = array(dim = c(1000, length(x)))
for (i in 1:1000) {
  dist[i, ] = yQ2_tropical_CIs(
    a = a[i],
    b = b[i],
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
    b = b[i],
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
a <- fixef(meQ1)[3] #coef for N fertilizer (if temperate)
b <- fixef(meQ1)[2]#coef for study duration
c <- fixef(meQ1)[1] #global intercept
f <- fixef(meQ1)[4]#coef for tropical=1
g <- fixef(meQ1)[5] #coef for interaction N in if tropical

##########plotting on linear scale
#functions to calculate tropical and temperate model (eq. 1) fit curves
yQ2_tropical = function(x) {
  10 ^ (c) * 10 ^ (f) * 10 ^ ((a + g) * ((x - n_bar) / (2 * n_sd)))
}
yQ2_temperate = function(x) {
  10 ^ (c) * 10 ^ (a * ((x - n_bar) / (2 * n_sd)))
}

# reordering tropical and temperate levels
no_full$tropical <-
  factor(as.factor(no_full$tropical), rev(levels(as.factor(no_full$tropical))))
levels(as.factor(no_full$tropical)) #tropical now comes first

#calculating curves and CIs
ytrop <- yQ2_tropical(x)
ytemp <- yQ2_temperate(x)
CItroplb <- civec_yQ2_trop [1, ]
CItropub <- civec_yQ2_trop [2, ]
CItemplb <- civec_yQ2_temp [1, ]
CItempub <- civec_yQ2_temp [2, ]
df <- data.frame(x, ytrop, ytemp, CItroplb, CItropub, CItemplb, CItempub)

#Figure 2
nop <- ggplot(no_full, aes(x = N_in_kg, y = NO_kg_seas_zero_rm)) +
  geom_point(
    aes(colour = (tropical), shape = (tropical)),
    alpha = .55,
    position = "jitter",
    size = 2
  ) +
  ylab(expression(paste(
    'Nitric oxide (kg NO-N ',  ~ ha ^ -1, season ^ -1, ')'
  ))) +
  xlab(expression(paste(
    'N inputs (kg N ha' ^ '-1',
    paste('season' ^ '-1'), ')'
  ))) +
  ggtitle('Nitric oxide losses') +
  theme_default() +
  ylim(0, 15) + xlim (0, 300) +
  geom_segment(aes(
    x = 0,
    xend = 0,
    y = 0,
    yend = 15
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
nop


# eq. 2 --------------------------------------------------------------------

no_trop <- no %>% filter(tropical > 0)

#set up dataframe for eq. 2 model
NO_model_tropical <-
  data.frame(row.names = 1:length(no_trop$NO_kg_seas_zero_rm))
NO_model_tropical$site_lat <-
  no_trop$lat_decimal # site ID; no transformation needed
NO_model_tropical$crop_type <-
  no_trop$crop_type # crop type; no transformation needed
NO_model_tropical$log_g_NO_day <-
  log10(no_trop$NO_kg_seas_zero_rm * 1000 / no_trop$duration_d) # crop type; no transformation needed
NO_model_tropical$split_app <-
  no_trop$split_app #split=1 (no trans. needed)
NO_model_tropical$org_Nfix_in <-
  no_trop$org_Nfix_in #organic/Nfix in=1  (no trans. needed)
NO_model_tropical$pH <-
  no_trop$pH #organic/Nfix in=1  (no trans. needed)
NO_model_tropical$N_in_kg <- z.trans(no_trop$N_in_kg) #z-transforming
NO_model_tropical$duration_d <-
  z.trans(no_trop$duration_d) #z-transforming
NO_model_tropical$precip <-
  z.trans(no_trop$MAP_mm_combined) #z-transforming
NO_model_tropical$pH <- z.trans(no_trop$pH) #z-transforming

#remove incomplete cases for lmer function
dim(NO_model_tropical)
NO_model_tropical <- na.omit(NO_model_tropical)
dim(NO_model_tropical) #removed incomplete cases
#removed sand because it reduces the sample size so much that it causes rank deficiencies

#test for collinearity
vif_test = lmer(
  log_g_NO_day ~ N_in_kg  + precip + split_app + org_Nfix_in + pH +
    (1 |
       site_lat),
  data = NO_model_tropical,
  na.action = NULL
)
vif(vif_test) #all covariates below 3

meq2 = lmer(
  log_g_NO_day ~ N_in_kg  + precip + split_app + org_Nfix_in + pH +
    crop_type + (1 | site_lat),
  data = NO_model_tropical,
  na.action = NULL
)
summary(meq2)
r.squaredGLMM(meq2)


#fig. 4

#precipitation/irrigation input plot
no_trop$log_no_g_day <-
  log10(no_trop$NO_kg_seas_zero_rm * 1000 / no_trop$duration_d) #y
x = 800:4200
precip <-
  mean(no_trop$MAP_mm_combined, na.rm = T) #to unstandardize sand
precip_sd <-
  sd(no_trop$MAP_mm_combined, na.rm = T) #to unstandardize sand
precip_fit_SE <- sqrt(diag(vcov(meq2)))[3] #SE amounts
yprecip <-
  function(x) {
    fixef(meq2)[1] + fixef(meq2)[3] * ((x - precip) / (2 * precip_sd))
  }
yprecip_hat <- yprecip(x = x)
yprecip_hat_lb <- yprecip_hat - rep(precip_fit_SE, length(yprecip_hat))
yprecip_hat_ub <- yprecip_hat + rep(precip_fit_SE, length(yprecip_hat))
df <- data.frame(x, yprecip_hat, yprecip_hat_lb, yprecip_hat_ub)

#precip plot
no_precip <- ggplot(no_trop, aes(x = MAP_mm_combined, y = log_no_g_day)) +
  geom_point(aes(alpha = .6)) +
  ylab(expression(paste('nitric oxide (g NO-N ',  ~ ha ^ -1,  ~ d ^ -1, ')'))) +
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
  geom_hline(yintercept = 0, color = "black") +
  geom_vline(xintercept = 700, color = "black") +
  xlim (700, 4600) +
  scale_y_continuous(
    labels = c(.1, 1, 10, 100),
    breaks = c(-1, 0, 1, 2),
    limits = c(-2, 2)
  )

no_precip

#pH plot
x = c(4, 5, 6, 7, 8.5)
pH <- mean(no_trop$pH, na.rm = T) #to unstandardize
pH_sd <- sd(no_trop$pH, na.rm = T) #to unstandardize
pH_fit_SE <- sqrt(diag(vcov(meq2)))[8] #SE amounts
ypH <- function(x) {
  fixef(meq2)[1] + fixef(meq2)[8] * ((x - pH) / (2 * pH_sd))
}
ypH_hat <- ypH(x = x)
ypH_hat_lb <- ypH_hat - rep(pH_fit_SE, length(ypH_hat))
ypH_hat_ub <- ypH_hat + rep(pH_fit_SE, length(ypH_hat))
df <- data.frame(x, ypH_hat, ypH_hat_lb, ypH_hat_ub)

#pH plot
no_pH <- ggplot(no_trop, aes(x = pH, y = log_no_g_day)) +
  geom_point(aes(alpha = .6)) +
  ylab("") +
  xlab(expression(paste('pH'))) +
  theme_default() +
  ylim(-2, 4) + xlim (4, 8.5) +
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
  geom_hline(yintercept = 0, color = "black") +
  geom_vline(xintercept = 4, color = "black") +
  scale_y_continuous(
    labels = c(0.1, 1, 10, 100),
    breaks = c(-1, 0, 1, 2),
    limits = c(-1, 2)
  )
no_pH


#Fig 4 binary/dummy crops plot
coef <- data.frame(data =
                     c(fixef(meq2)[c(
                       '(Intercept)',
                       '(Intercept)',
                       'split_app',
                       '(Intercept)',
                       'org_Nfix_in',
                       'crop_typefallow',
                       'crop_typepasture',
                       'crop_typelegume'
                     )]) +
                     c(rep(0, 2), fixef(meq2)['(Intercept)'], 0,
                       rep(fixef(meq2)['(Intercept)'], 4)))
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
  meq2.se.comb[meq2.names == 'crop_typepasture', 2],
  meq2.se.comb[meq2.names == 'crop_typelegume', 2]
)

meq2_outputs <- cbind(coef, SE)
colnames(meq2_outputs) <- c('coef', 'SE')
meq2_outputs$names <- c('cereal',
                        'sa0',
                        'sa1',
                        'org0',
                        'org1',
                        'fallow',
                        'pasture',
                        'legume')
meq2_outputs

#creating individual dfs for binary variables to join to crop types
discrete_var <- no_trop %>%
  dplyr::select(org_Nfix_in , split_app, log_no_g_day)
org1 <- discrete_var %>% filter(org_Nfix_in == 1) %>%
  mutate(crop_type = rep("org1", 3))
org0 <- discrete_var %>% filter(org_Nfix_in == 0) %>%
  mutate(crop_type = rep("org0", 46))
sa1 <- discrete_var %>% filter(split_app == 1) %>%
  mutate(crop_type = rep("sa1", 15))
sa0 <- discrete_var %>% filter(split_app == 0) %>%
  mutate(crop_type = rep("sa0", 34))


#joining the data frames
discrete_var2 <- no_trop %>%
  dplyr::select(org_Nfix_in , split_app, crop_type, log_no_g_day)

discrete_var3 <- rbind(discrete_var2, org0, org1, sa1, sa0)

#reordering the factors based on the model outputs
discrete_var3$crop_type <- factor(
  discrete_var3$crop_type,
  levels = c(
    'cereal',
    'sa0',
    'sa1',
    'org0',
    'org1',
    'fallow',
    'pasture',
    'legume'
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
    'pasture',
    'legume'
  )

#first removing the 2 "other crop types" which are not estimated by the model
discrete_var3 <- discrete_var3 %>% filter(crop_type != 'other')


#Fig 4--we removed the violin plots from these because the outline covered data points
no_Q3_disc <-  ggplot() +
  geom_violin(
    data = discrete_var3,
    aes(x = crop_type, y = log_no_g_day),
    fill = "light gray",
    colour = "light gray"
  ) +
  geom_jitter(
    data = discrete_var3,
    aes(x = crop_type, y = log_no_g_day),
    width = .05,
    alpha = .3
  ) +
  ylab(expression(paste('nitric oxide (g NO-N',  ~  ~ ha ^ -1,  ~ d ^ -1, ')'))) +
  geom_pointrange(data = meq2_outputs,
                  aes(
                    reorder(names, coef),
                    coef,
                    ymin = coef - SE,
                    ymax =  coef + SE
                  ),
                  col = 'red') +
  scale_colour_manual(values = coef) +
  scale_x_discrete('' , labels = x_axis) +
  geom_hline(yintercept = 0, color = "black") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  geom_segment(aes(
    x = 0,
    xend = 0,
    y = 0,
    yend = 2
  ), colour = "black") +
  theme_default(axis_text_size = 13) +
  annotate(
    "text",
    x = 7,
    y = 2,
    size = 8,
    label = "*"
  ) +
  geom_vline(xintercept = c(1.5, 3.5, 5.5), linetype = 'longdash')

no_Q3_disc <- no_Q3_disc + scale_y_continuous(breaks = c(0, 1, 2),
                                              labels = c(0, 10, 100))
no_Q3_disc

