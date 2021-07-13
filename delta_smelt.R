################## Delta Smelt Dose-Response
################## Eric Lawrence 7/5/2021

################## 
################## 

library(drc)
library(tidyverse)

smelt <- read.csv("smelt_range_finding_survival_data.csv")

tibble(smelt)



########### bifenthrin #########

###########################################################################################
############################################################################################

############################## Model Construction

bif <- smelt %>%
  filter(chemical == "bifenthrin") %>%
  mutate(percent_m = dead / total)

for (i in 1:length(bif$percent_m)) {
  if (bif$percent_m[i] == 0) {
    bif$percent_m[i] == 0.000001
  }
}

plot(bif$percent_m ~ bif$dose.adj)

bif.m <- drm(bif$percent_m ~ bif$dose.adj, fct = LL.3())

plot(bif.m)
qqnorm(resid(bif.m))
qqline(resid(bif.m))

plot(resid(bif.m)~predict(bif.m))
abline(h=0)

############################## Model Selection

mselect(bif.m, list(LL.3(), LL.4(), LL.5(), EXD.2(), EXD.3(), W1.3(), W1.4(), W2.3(), W2.4()) )
#### LL.4 selected for best fit

############################## Prediction interval

# Sequence of doses to predict
pre.seq <- exp(seq(log(0.000001), log(0.01), length = 50))


pre.d <- expand.grid(conc = pre.seq)

summary(pre.seq)

# Make predictions
pre.p <- predict(bif.m, newdata = pre.d, interval = "prediction")
pre.p

pre.p <- predict(bif.m, interval = "prediction", level = 0.5)


# Add predictions to dataframe
pre.d$p <- pre.p[,1]
pre.d$p.l <- pre.p[,2]
pre.d$p.u <- pre.p[,3]
pre.d

summary(pre.d)

# Make models for upper and lower predictions

pre.l.m <- drm(pre.d$p.l~pre.d$conc, fct=LL.3())
pre.u.m <- drm(pre.d$p.u ~ pre.d$conc, fct=LL.3())


############################## Figure for Bifenthrin

plot(bif.m, type = "confidence", ylab = "Percent Mortality",
     xlab = "Concentration (mg/L)",
     #ylim = c(0,165),
     main = "Bifenthrin, LL.4")
plot(bif.m, type = "all", add = TRUE)
plot(pre.d$p.l~pre.d$conc, add = T)

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