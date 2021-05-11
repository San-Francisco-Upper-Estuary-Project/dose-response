################## Mixture Analysis Diazinon and Malathion
################## Eric Lawrence 9/17/2020

################## Data from Laetz et al. 2009 and Laetz et al. 2013
################## Some code copied from MixtureAnalysis_EJL.R

library(drc)
library(dplyr)
library(plotly)

setwd("C:/Users/ej/Dropbox/New Risk/Laetz_Data/Laetz data/Lawrence_analysis_2020")

############################################################################# Single Chemical Analysis
##########################################################
##########################################################

################## Diazinon #####################################

diaz <- read.csv("Diazinon_R.csv")

# Using only 2005 data
diaz <- diaz %>%
  filter(Year == "2005")

summary(diaz)
plot(diaz$ache~diaz$conc)

############################## Model Construction

diaz.mod <- drm(diaz$ache~diaz$conc, fct=W2.4())

plot(diaz.mod)
qqnorm(resid(diaz.mod))
qqline(resid(diaz.mod))

plot(resid(diaz.mod)~predict(diaz.mod))
abline(h=0)

############################## Model Selection

mselect(diaz.mod, list(LL.3(), LL.4(), LL.5(), EXD.2(), EXD.3(), W1.3(), W1.4(), W2.3(), W2.4()) )
#### W2.4 selected for best fit

############################## Prediction interval

# Sequence of doses to predict
pre.seq <- exp(seq(log(0.00001), log(1000), length = 100))


pre.d <- expand.grid(conc = pre.seq)

summary(pre.seq)

# Make predictions
pre.p <- predict(diaz.mod, newdata = pre.d, interval = "prediction")

# Add predictions to dataframe
pre.d$p <- pre.p[,1]
pre.d$p.l <- pre.p[,2]
pre.d$p.u <- pre.p[,3]
pre.d

summary(pre.d)

# Make models for upper and lower predictions

pre.l.m <- drm(pre.d$p.l~pre.d$conc, fct = W2.4())
pre.u.m <- drm(pre.d$p.u ~ pre.d$conc, fct = W2.4())


############################## Figure for Diazinon

plot(diaz.mod, type = "confidence", ylab = "AchE (Percent Control)",
     xlab = "Concentration (ug/L)",
     ylim = c(0,165),
     main = "Diazinon, W2.4, 95% CI")
plot(diaz.mod, type = "all", add = TRUE)

#Prediction Interval

plot(pre.l.m, add = TRUE, type = "none", col = "blue")
plot(pre.u.m, add=TRUE, type = "none", col = "blue")

#Legend

legend("topright", c("Prediction Interval", "95% Confidence Interval"),
       lty=c(1,1),
       lwd=c(2.5,2.5), 
       col=c('blue', 'gray'),
       bty = "n")



############################## Model Parameters

diaz.mod

# Call:
#   drm(formula = diaz$ache ~ diaz$conc, fct = W2.4())
# 
#   Coefficients:
#   b:(Intercept)  c:(Intercept)  d:(Intercept)  e:(Intercept)  
#        -0.9116        24.7515        99.9396        35.6787   

summary(diaz.mod)

################################ Diazinon EC50

diaz.ec <- ED(diaz.mod, c(5, 10, 20, 50), interval = "delta")
diaz.ec

#       Estimate Std. Error   Lower   Upper
# e:1:5   10.7078     4.2001  2.3426 19.0729
# e:1:10  14.2912     4.5578  5.2135 23.3689
# e:1:20  21.1686     4.9518 11.3062 31.0311
# e:1:50  53.3358    12.2935 28.8511 77.8204

ec50.d <- 53.3358

################################ Diazinon Toxic Units

diaz.s.tu <- diaz$conc / ec50.d
plot(diaz$ache ~ diaz.s.tu)

diaz.s.tu.m <- drm(diaz$ache ~ diaz.s.tu, fct=W2.4())
plot(diaz.s.tu.m)


############################################################ Malathion
############################

mala <- read.csv("Malathion_R.csv")
### I selected records from the Laetz et al. 2009 data that allowed me to reconstruct the EC50 value for malathion 
### recorded in Laetz et al. 2009 all other malathion records removed.


summary(mala)
plot(mala$ache~mala$conc)



plot(mala$ache~mala$conc)

mala.mod <- drm(mala$ache~mala$conc, fct=LL.3())

mselect(mala.mod, list(LL.3(), LL.4(), LL.5(), EXD.2(), EXD.3(), W1.3(), W1.4(), W2.3(), W2.4()) )

#### LL.3 selected for best fit

mala.mod

############################## Prediction interval

# Sequence of doses to predict
pre.seq <- exp(seq(log(0.00001), log(1000), length = 100))


pre.d <- expand.grid(conc = pre.seq)

summary(pre.seq)

# Make predictions
pre.p <- predict(mala.mod, newdata = pre.d, interval = "prediction")

# Add predictions to dataframe
pre.d$p <- pre.p[,1]
pre.d$p.l <- pre.p[,2]
pre.d$p.u <- pre.p[,3]
pre.d

summary(pre.d)

# Make models for upper and lower predictions

pre.l.m <- drm(pre.d$p.l~pre.d$conc, fct = LL.3())
pre.u.m <- drm(pre.d$p.u ~ pre.d$conc, fct = LL.3())


######################### Malathion Figure

plot(mala.mod)

plot(mala.mod, type = "confidence", ylab = "AchE (Percent Control)",
     xlab = "ug/L",
     ylim = c(0,165),
     main = "Malathion, LL3, 95% CI")
plot(mala.mod, type = "all", add = TRUE)

#Prediction Interval

plot(pre.l.m, add = TRUE, type = "none", col = "blue")
plot(pre.u.m, add=TRUE, type = "none", col = "blue")

#Legend

legend("topright", c("Prediction Interval", "95% Confidence Interval"),
       lty=c(1,1),
       lwd=c(2.5,2.5), 
       col=c('blue', 'gray'),
       bty = "n")



summary(mala.mod)

############################### Malathion EC50

mala.ec <- ED(mala.mod, c(5, 10, 20, 50), interval = "delta")
mala.ec

#        Estimate Std. Error    Lower    Upper
# e:1:5  18.38975   3.532989 11.35467 25.42483
# e:1:10 26.98071   4.109170 18.79830 35.16311
# e:1:20 40.90029   4.600384 31.73975 50.06083
# e:1:50 83.28927   5.939810 71.46160 95.11695

ec50.m <- 83.29


############################### Malathion Toxic Units


mala.s.tu <- mala$conc / ec50.m
plot(mala$ache ~ mala.s.tu)

mala.s.tu.mod <- drm(mala$ache ~ mala.s.tu, fct=LL.3())
plot(diaz.s.tu.m)


##########################################################################################################
################################################## TU approach for Diaz + Mala

mdmix <- read.csv("donmon_R.csv")

summary(mdmix)
plot(mdmix$pct_control~mdmix$moles)

# Convert concentrations to toxic units

tu.m <- mdmix$meas_mon / ec50.m
tu.m

tu.d <- mdmix$meas_diaz / ec50.d
tu.d

plot(tu.d ~ tu.m)

plot(tu.d ~ tu.m,
     xlab = "Malathion TU",
     ylab = "Diazinon TU",
     xlim = c(0,1),
     ylim = c(0,1),
     main = "Isobole Plot, Toxic Units for Malathion and Diazinon"
)


#################### Model Construction


tu.sum <- tu.m + tu.d
tu.sum

plot(mdmix$pct_control ~ tu.sum)

tu.mix.m <- drm(mdmix$pct_control ~ tu.sum, fct = LL.3())
plot(tu.mix.m)
tu.mix.m

#   Coefficients:
#   b:(Intercept)  d:(Intercept)  e:(Intercept)  
#       10.18484       99.24952        0.05025  

mselect(tu.mix.m, list(LL.3(), LL.4(), LL.5(), EXD.2(), EXD.3(), W1.3(), W1.4(), W2.3(), W2.4()))

## LL.3 best fit


summary(mdmix)

summary(tu.sum)

################### TU EC50

tu.mix.ED <- ED(tu.mix.m, c(5, 10, 20, 50), interval = "delta")
tu.mix.ED

#         Estimate  Std. Error      Lower      Upper
# e:1:5  0.03763569 0.002717258 0.03223117 0.04304020
# e:1:10 0.04050065 0.002229490 0.03606629 0.04493502
# e:1:20 0.04385722 0.001643798 0.04058777 0.04712667
# e:1:50 0.05025212 0.001026211 0.04821103 0.05229322

######### Prediction Interval

pre.seq <- expand.grid(seq(from = 0, to = 2.5, by = 0.01))
pre.d <- expand.grid(tu = pre.seq)

summary(pre.seq)

pre.d
pre.p <- predict(tu.mix.m, newdata = pre.d, interval = "prediction")

pre.d$p <- pre.p[,1]
pre.d$p.l <- pre.p[,2]
pre.d$p.u <- pre.p[,3]
pre.d

summary(pre.d)

plot(pre.d$p.l~pre.d$tu, add=TRUE)

pre.l.m <- drm(pre.d$p.l~pre.d$Var1, fct = LL.3())

pre.u.m <- drm(pre.d$p.u ~ pre.d$Var1, fct = LL.5())

######################################### TU Figure

plot(tu.mix.m, type = "confidence", ylab = "AchE (Percent Control)",
     xlab = "Toxic Units (Summed)",
     xlim = c(0, 1),
     ylim = c(0, 150),
     main = "Malathion + Diazinon Mixture, LL3",
     add = FALSE)
plot(tu.mix.m, type = "all", add = TRUE)

# Prediction interval

plot(pre.l.m, add = TRUE, type = "none", col = "blue")
plot(pre.u.m, add=TRUE, type = "none", col = "blue")

legend("topright", c("Prediction Interval", "95% Confidence Interval"),
       lty=c(1,1),
       lwd=c(2.5,2.5), 
       col=c('blue', 'gray'),
       bty = "n")

######################################################################################
####################################### 3D Scatterplot
############### ####################### Using plotly

donmon <- read.csv("DONMON_2008_resp_surf.csv")

## Set the the axes
axx <- list(nticks = 8, range = c(1.4, 0), title = "Malathion (ug/L)")
axy <- list(nticks = 8, range = c(0, 3), title = "Diazinon (ug/L")
axz <- list(nticks = 5, range = c(0,200), title = "AChE (pct control)")

## Make total concentration factor to use in grouping and legend
total_conc <- as.factor(donmon$meas_mon + donmon$meas_diaz)

## Construct Figure
fig <- plot_ly(x = donmon$meas_mon, y = donmon$meas_diaz,
               z = donmon$pct_control, color = ~total_conc)
fig <- fig %>% layout(scene= list(xaxis = axx, yaxis=axy,
                                  zaxis = axz),
                      legend = list(title = list(text = "Concentration Sum (ug/L)")))

fig
