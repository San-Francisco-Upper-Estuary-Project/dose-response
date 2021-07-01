####### You first need to install the drc package ######

install.packages("drc")

#############################################################################

library("drc") #load package 

#read in Striped Bass and Chinook Salmon toxicity data
chinook = read.csv("Delta_Conc_Response/chinooktoxdataM.csv")
bass = read.csv("Delta_Conc_Response/basstoxdataM.csv")

#########################################################################

########### create dose-response models for striped bass data ###########

#########################################################################

## Copper ##
bassCu_drm = drm(percentmortality~conc, data = bass[1:76, 1:2], fct = AR.2())
x1 = seq(0,100, by = 0.1)
new1 = data.frame(dose = x1)
bassCu_pred_ints = predict(bassCu_drm, interval = "prediction", 
                                     newdata = new1) #calculates 95% prediction intervals for curve

# Selecting model: best is AR.2
mselect(bassCu_drm, list(LL.5(), LL.4(), LL.3(), LL.2(), W1.4(), W1.3(), AR.2(), AR.3(), LN.2(), LN.3(), 
                         LN.4(), EXD.2(), EXD.3(), W2.3(), W2.4()))

plot(bassCu_drm, log = '', xlab = expression(paste("Copper Concentration (",mu,"M)")), 
     ylab = "Mortality (%)",
     type = "confidence", #adds 95% confidence intervals to plot
     main = "Copper Toxicity to Striped Bass",
     lwd = 2, cex.lab = 1.3, cex.axis = 1.4, xlim = c(0,100), ylim = c(0,100))

points(percentmortality~conc, data = bass[1:76, 1:2], pch = 1, lwd = 1.5) #plot data 
lines(bassCu_pred_ints[,2]~x1, lty = 2) #plot lower prediction bounds
lines(bassCu_pred_ints[,3]~x1, lty = 2) #plot upper prediction bounds

legend("bottomright", c(expression("0.143 μM LC"[20]), expression("0.445 μM LC"[50]), expression("2.96 μM LC"[99]),
                        "95% Prediction Interval", "95% Confidence Interval"), #change after running ED code below
       lty = c(NA, NA, NA, 2, NA), fill = c(NA, NA, NA, NA,"gray65"), border = "white", bty = 'n', 
       merge = TRUE, cex = 1.3)


## Cadmium ##
bassCd_drm = drm(percentmortality~conc, data = bass[77:232, 1:2], fct = W2.4())
x2 = seq(0,250, by = 1)
new2 = data.frame(dose = x2)
bassCd_pred_ints = predict(bassCd_drm, interval = "prediction", newdata = new2)

# Selecting model: best is W2.4
mselect(bassCd_drm, list(LL.5(), LL.4(), LL.3(), LL.2(), W1.4(), W1.3(), AR.2(), AR.3(), LN.2(), LN.3(), 
                         LN.4(), EXD.2(), EXD.3(), W2.3(), W2.4()))

plot(bassCd_drm, log = '', xlab = expression(paste("Cadmium Concentration (",mu,"M)")), 
     ylab = "Mortality (%)",
     type = "confidence", #adds 95% confidence intervals to plot
     main = "Cadmium Toxicity to Striped Bass",
     lwd = 2, cex.lab = 1.3, cex.axis = 1.4, xlim = c(0,250), ylim = c(0,100))

points(percentmortality~conc, data = bass[77:232, 1:2], pch = 1, lwd = 1.5) #plot data 
lines(bassCd_pred_ints[,2]~x2, lty = 2) #plot lower prediction bounds
lines(bassCd_pred_ints[,3]~x2, lty = 2) #plot upper prediction bounds

legend("right", c(expression("0.0193 μM LC"[20]), expression("0.0598 μM LC"[50]), expression("0.398 μM LC"[99]),
                        "95% Prediction Interval", "95% Confidence Interval"), #change after running ED code below
       lty = c(NA, NA, NA, 2, NA), fill = c(NA, NA, NA, NA,"gray65"), border = "white", bty = 'n', 
       merge = TRUE, cex = 1.3)


## Molinate ##
bassMol_drm = drm(percentmortality~conc, data = bass[233:280,1:2], fct = AR.3())
x3 = seq(0,225, by = 0.1)
new3 = data.frame(dose = x3)
bassMol_pred_ints = predict(bassMol_drm, interval = "prediction", 
                            newdata = new3) #calculates 95% prediction intervals for curve

# Selecting model: best is AR.3
mselect(bassMol_drm, list(LL.5(), LL.4(), LL.3(), LL.2(), W1.4(), W1.3(), AR.2(), AR.3(), LN.2(), LN.3(), 
                         LN.4(), EXD.2(), EXD.3(), W2.3(), W2.4()))

plot(bassMol_drm, log = '', xlab = expression(paste("Molinate Concentration (",mu,"M)")), 
     ylab = "Mortality (%)",
     type = "confidence", #adds 95% confidence intervals to plot
     main = "Molinate Toxicity to Striped Bass",
     lwd = 2, cex.lab = 1.3, cex.axis = 1.4, xlim = c(0,225), ylim = c(0,100))

points(percentmortality~conc, data = bass[233:280,1:2], pch = 1, lwd = 1.5) #plot data 
lines(bassMol_pred_ints[,2]~x3, lty = 2) #plot lower prediction bounds
lines(bassMol_pred_ints[,3]~x3, lty = 2) #plot upper prediction bounds

legend("bottomright", c(expression("15.39 μM LC"[20]), expression("47.81 μM LC"[50]), expression("317.66 μM LC"[99]),
                        "95% Prediction Interval", "95% Confidence Interval"), #change after running ED code below
       lty = c(NA, NA, NA, 2, NA), fill = c(NA, NA, NA, NA,"gray65"), border = "white", bty = 'n', 
       merge = TRUE, cex = 1.3)

## Thiobencarb ##
bassThio_drm = drm(percentmortality~conc, data = bass[281:339,1:2], fct = LN.4())
x4 = seq(0,15, by = 0.01)
new4 = data.frame(dose = x4)
bassThio_pred_ints = predict(bassThio_drm, interval = "prediction", 
                             newdata = new4) #calculates 95% prediction intervals for curve

# Selecting model: best is LN.4
mselect(bassThio_drm, list(LL.5(), LL.4(), LL.3(), LL.2(), W1.4(), W1.3(), AR.2(), AR.3(), LN.2(), LN.3(), 
                         LN.4(), EXD.2(), EXD.3(), W2.3(), W2.4()))

plot(bassThio_drm, log = '', xlab = expression(paste("Thiobencarb Concentration (",mu,"M)")), 
     ylab = "Mortality (%)",
     type = "confidence", #adds 95% confidence intervals to plot
     main = "Thiobencarb Toxicity to Striped Bass",
     lwd = 2, cex.lab = 1.3, cex.axis = 1.4, xlim = c(0,15), ylim = c(0,100))

points(percentmortality~conc, data = bass[281:339,1:2], pch = 1, lwd = 1.5) #plot data 
lines(bassThio_pred_ints[,2]~x4, lty = 2) #plot lower prediction bounds
lines(bassThio_pred_ints[,3]~x4, lty = 2) #plot upper prediction bounds

legend("bottomright", c(expression("0.685 μM LC"[20]), expression("2.13 μM LC"[50]), expression("14.14 μM LC"[99]),
                        "95% Prediction Interval", "95% Confidence Interval"), #change after running ED code below
       lty = c(NA, NA, NA, 2, NA), fill = c(NA, NA, NA, NA,"gray65"), border = "white", bty = 'n', 
       merge = TRUE, cex = 1.3)

## Getting all EC20, 50, and 99 values
ALLbass_drm <- drm(percentmortality~conc, toxicant, data = bass, fct = c(AR.2(), W2.4(), AR.3(), LN.4()))
ED(ALLbass_drm, c(20,50,99), interval="delta", level = 0.95)


#########################################################################

########## create dose-response models for chinook salmon data ##########

#########################################################################

## Esfenvalerate ##
chinEsf_drm = drm(percentmortality~conc, data = chinook[1:4, 1:2], fct = W2.3())
x5 = seq(0,0.5, by = 0.00001)
new5 = data.frame(dose = x5)
chinEsf_pred_ints = predict(chinEsf_drm, interval = "prediction", 
                           newdata = new5) #calculates 95% prediction intervals for curve

# Selecting model: best is W2.3
mselect(chinEsf_drm, list(LL.5(), LL.4(), LL.3(), LL.2(), W1.4(), W1.3(), AR.2(), AR.3(), LN.2(), LN.3(), 
                         LN.4(), EXD.2(), EXD.3(), W2.3(), W2.4()))

plot(chinEsf_drm, log = '', xlab = expression(paste("Esfenvalerate Concentration (",mu,"M)")), 
     ylab = "Mortality (%)",
     type = "confidence", #adds 95% confidence intervals to plot
     main = "Esfenvalerate Toxicity to Chinook Salmon",
     lwd = 2, cex.lab = 1.3, cex.axis = 1.4, xlim = c(0,0.5), ylim = c(0,100))

points(percentmortality~conc, data = chinook[1:4, 1:2], pch = 1, lwd = 1.5) #plot data 
lines(chinEsf_pred_ints[,2]~x5, lty = 2) #plot lower prediction bounds
lines(chinEsf_pred_ints[,3]~x5, lty = 2) #plot upper prediction bounds

legend("bottomright", c(expression("0.0144 μM LC"[20]), expression("0.0159 μM LC"[50]), expression("NA μM LC"[99]),
                        "95% Prediction Interval", "95% Confidence Interval"), #change after running ED code below
       lty = c(NA, NA, NA, 2, NA), fill = c(NA, NA, NA, NA,"gray65"), border = NA, bty = 'n', 
       merge = TRUE, cex = 1.3)

ED(chinEsf_drm, c(20,50,99), interval = "delta")


## Triclopyr Butoxyethyl Ester ##
chinTri_drm = drm(percentmortality~conc, data = chinook[5:19, 1:2], fct = LL.3())
x6 = seq(0,110, by = 0.001)
new6 = data.frame(dose = x6)
chinTri_pred_ints = predict(chinTri_drm, interval = "prediction", 
                            newdata = new6) #calculates 95% prediction intervals for curve

# Selecting model: best is LL.3
mselect(chinTri_drm, list(LL.5(), LL.4(), LL.3(), LL.2(), W1.4(), W1.3(), AR.2(), AR.3(), LN.2(), LN.3(), 
                          LN.4(), EXD.2(), EXD.3(), W2.3(), W2.4()))

plot(chinTri_drm, log = '', xlab = expression(paste("Triclopyr Butoxyethyl Ester Concentration (",mu,"M)")), 
     ylab = "Mortality (%)",
     type = "confidence", #adds 95% confidence intervals to plot
     main = "Triclopyr Toxicity to Chinook Salmon",
     lwd = 2, cex.lab = 1.3, cex.axis = 1.4, xlim = c(0,110), ylim = c(0,100))

points(percentmortality~conc, data = chinook[5:19, 1:2], pch = 1, lwd = 1.5) #plot data 
lines(chinTri_pred_ints[,2]~x6, lty = 2) #plot lower prediction bounds
lines(chinTri_pred_ints[,3]~x6, lty = 2) #plot upper prediction bounds

legend("top", c(expression("4.44 μM LC"[20]), expression("4.55 μM LC"[50]), expression("NA μM LC"[99]),
                        "95% Prediction Interval", "95% Confidence Interval"), #change after running ED code below
       lty = c(NA, NA, NA, 2, NA), fill = c(NA, NA, NA, NA,"gray65"), border = NA, bty = 'n', 
       merge = TRUE, cex = 1.3)

ED(chinTri_drm, c(20,50,99), interval = "delta")


## Diazinon ##
chinDia_drm = drm(percentmortality~conc, data = chinook[20:23, 1:2], fct = W1.3())
x7 = seq(0,200, by = 0.1)
new7 = data.frame(dose = x7)
chinDia_pred_ints = predict(chinDia_drm, interval = "prediction", 
                            newdata = new7) #calculates 95% prediction intervals for curve

# Selecting model: best is W1.3
mselect(chinDia_drm, list(LL.5(), LL.4(), LL.3(), LL.2(), W1.4(), W1.3(), AR.2(), AR.3(), LN.2(), LN.3(), 
                          LN.4(), EXD.2(), EXD.3(), W2.3(), W2.4()))

plot(chinDia_drm, log = '', xlab = expression(paste("Diazinon Concentration (",mu,"M)")), 
     ylab = "Mortality (%)",
     type = "confidence", #adds 95% confidence intervals to plot
     main = "Diazinon Toxicity to Chinook Salmon",
     lwd = 2, cex.lab = 1.3, cex.axis = 1.4, xlim = c(0,200), ylim = c(0,100))

points(percentmortality~conc, data = chinook[20:23, 1:2], pch = 1, lwd = 1.5) #plot data 
lines(chinDia_pred_ints[,2]~x7, lty = 2) #plot lower prediction bounds
lines(chinDia_pred_ints[,3]~x7, lty = 2) #plot upper prediction bounds

legend("bottomright", c(expression("68.36 μM LC"[20]), expression("78.49 μM LC"[50]), expression("157.16 μM LC"[99]),
                        "95% Prediction Interval", "95% Confidence Interval"), #change after running ED code below
       lty = c(NA, NA, NA, 2, NA), fill = c(NA, NA, NA, NA,"gray65"), border = NA, bty = 'n', 
       merge = TRUE, cex = 1.3)

ED(chinDia_drm, c(20,50,99), interval = "delta")


## Dinoseb ##
chinDino_drm = drm(percentmortality~conc, data = chinook[24:31, 1:2], fct = AR.2())
x8 = seq(0,5, by = 0.001)
new8 = data.frame(dose = x8)
chinDino_pred_ints = predict(chinDino_drm, interval = "prediction", 
                            newdata = new8) #calculates 95% prediction intervals for curve

# Selecting model: best is AR.2
mselect(chinDino_drm, list(LL.5(), LL.4(), LL.3(), LL.2(), W1.4(), W1.3(), AR.2(), AR.3(), LN.2(), LN.3(), 
                          LN.4(), EXD.2(), EXD.3(), W2.3(), W2.4()))

plot(chinDino_drm, log = '', xlab = expression(paste("Dinoseb Concentration (",mu,"M)")), 
     ylab = "Mortality (%)",
     type = "confidence", #adds 95% confidence intervals to plot
     main = "Dinoseb Toxicity to Chinook Salmon",
     lwd = 2, cex.lab = 1.3, cex.axis = 1.4, xlim = c(0,5), ylim = c(0,100))

points(percentmortality~conc, data = chinook[24:31, 1:2], pch = 1, lwd = 1.5) #plot data 
lines(chinDino_pred_ints[,2]~x8, lty = 2) #plot lower prediction bounds
lines(chinDino_pred_ints[,3]~x8, lty = 2) #plot upper prediction bounds

legend("bottomright", c(expression("0.286 μM LC"[20]), expression("0.888 μM LC"[50]), expression("5.90 μM LC"[99]),
                        "95% Prediction Interval", "95% Confidence Interval"), #change after running ED code below
       lty = c(NA, NA, NA, 2, NA), fill = c(NA, NA, NA, NA,"gray65"), border = NA, bty = 'n', 
       merge = TRUE, cex = 1.3)

ED(chinDino_drm, c(20,50,99), interval = "delta")


## Copper ##
chinCu_drm = drm(percentmortality~conc, data = chinook[32:37, 1:2], fct = W1.4())
x9 = seq(0,1, by = 0.0001)
new9 = data.frame(dose = x9)
chinCu_pred_ints = predict(chinCu_drm, interval = "prediction", 
                           newdata = new9) #calculates 95% prediction intervals for curve

# Selecting model: best is W1.4
mselect(chinCu_drm, list(LL.5(), LL.4(), LL.3(), LL.2(), W1.4(), W1.3(), AR.2(), AR.3(), LN.2(), LN.3(), 
                         LN.4(), EXD.2(), EXD.3(), W2.3(), W2.4()))

plot(chinCu_drm, log = '', xlab = expression(paste("Copper Concentration (",mu,"M)")), 
     ylab = "Mortality (%)",
     type = "confidence", #adds 95% confidence intervals to plot
     main = "Copper Toxicity to Chinook Salmon",
     lwd = 2, cex.lab = 1.3, cex.axis = 1.4, xlim = c(0,1), ylim = c(0,100))

points(percentmortality~conc, data = chinook[32:37, 1:2], pch = 1, lwd = 1.5) #plot data 
lines(chinCu_pred_ints[,2]~x9, lty = 2) #plot lower prediction bounds
lines(chinCu_pred_ints[,3]~x9, lty = 2) #plot upper prediction bounds

legend("bottomright", c(expression("0.263 μM LC"[20]), expression("0.406 μM LC"[50]), expression("3.59 μM LC"[99]),
                        "95% Prediction Interval", "95% Confidence Interval"), #change after running ED code below
       lty = c(NA, NA, NA, 2, NA), fill = c(NA, NA, NA, NA,"gray65"), border = NA, bty = 'n', 
       merge = TRUE, cex = 1.3)

ED(chinCu_drm, c(20,50,99), interval = "delta")


## Zinc ##
chinZn_drm = drm(percentmortality~conc, data = chinook[38:43, 1:2], fct = AR.3())
x10 = seq(0,10, by = 0.01)
new10 = data.frame(dose = x10)
chinZn_pred_ints = predict(chinZn_drm, interval = "prediction", 
                           newdata = new10) #calculates 95% prediction intervals for curve

# Selecting model: best is AR.3
mselect(chinZn_drm, list(LL.5(), LL.4(), LL.3(), LL.2(), W1.4(), W1.3(), AR.2(), AR.3(), LN.2(), LN.3(), 
                         LN.4(), EXD.2(), EXD.3(), W2.3(), W2.4()))

plot(chinZn_drm, log = '', xlab = expression(paste("Zinc Concentration (",mu,"M)")), 
     ylab = "Mortality (%)",
     type = "confidence", #adds 95% confidence intervals to plot
     main = "Zinc Toxicity to Chinook Salmon",
     lwd = 2, cex.lab = 1.3, cex.axis = 1.4, xlim = c(0,10), ylim = c(0,100))

points(percentmortality~conc, data = chinook[38:43, 1:2], pch = 1, lwd = 1.5) #plot data 
lines(chinZn_pred_ints[,2]~x10, lty = 2) #plot lower prediction bounds
lines(chinZn_pred_ints[,3]~x10, lty = 2) #plot upper prediction bounds

legend("topleft", c(expression("NA μM LC"[20]), expression("NA μM LC"[50]), expression("NA μM LC"[99]),
                        "95% Prediction Interval", "95% Confidence Interval"), #change after running ED code below
       lty = c(NA, NA, NA, 2, NA), fill = c(NA, NA, NA, NA,"gray65"), border = "white", bty = 'n', 
       merge = TRUE, cex = 1.3)

ED(chinZn_drm, c(20,50,99), interval = "delta")


## Cadmium ##
chinCd_drm = drm(percentmortality~conc, data = chinook[44:49, 1:2], fct = W1.4())
x11 = seq(0,0.5, by = 0.0001)
new11 = data.frame(dose = x11)
chinCd_pred_ints = predict(chinCd_drm, interval = "prediction", 
                           newdata = new11) #calculates 95% prediction intervals for curve

# Selecting model: best is W1.4
mselect(chinCd_drm, list(LL.5(), LL.4(), LL.3(), LL.2(), W1.4(), W1.3(), AR.2(), AR.3(), LN.2(), LN.3(), 
                         LN.4(), EXD.2(), EXD.3(), W2.3(), W2.4()))

plot(chinCd_drm, log = '', xlab = expression(paste("Cadmium Concentration (",mu,"M)")), 
     ylab = "Mortality (%)",
     type = "confidence", #adds 95% confidence intervals to plot
     main = "Cadmium Toxicity to Chinook Salmon",
     lwd = 2, cex.lab = 1.3, cex.axis = 1.4, xlim = c(0,0.5), ylim = c(0,100))

points(percentmortality~conc, data = chinook[44:49, 1:2], pch = 1, lwd = 1.5) #plot data 
lines(chinCd_pred_ints[,2]~x11, lty = 2) #plot lower prediction bounds
lines(chinCd_pred_ints[,3]~x11, lty = 2) #plot upper prediction bounds

legend("topleft", c(expression("0.0105 μM LC"[20]), expression("NA μM LC"[50]), expression("NA μM LC"[99]),
                        "95% Prediction Interval", "95% Confidence Interval"), #change after running ED code below
       lty = c(NA, NA, NA, 2, NA), fill = c(NA, NA, NA, NA,"gray65"), border = "white", bty = 'n', 
       merge = TRUE, cex = 1.3)

ED(chinCd_drm, c(20,50,99), interval = "delta")


## DDT ##
chinDDT_drm = drm(percentmortality~conc, data = chinook[50:54, 1:2], fct = LN.4())
x12 = seq(0,1250, by = 1)
new12 = data.frame(dose = x12)
chinDDT_pred_ints = predict(chinDDT_drm, interval = "prediction", 
                           newdata = new12) #calculates 95% prediction intervals for curve

# Selecting model: best is LN.4
mselect(chinDDT_drm, list(LL.5(), LL.4(), LL.3(), LL.2(), W1.4(), W1.3(), AR.2(), AR.3(), LN.2(), LN.3(), 
                         LN.4(), EXD.2(), EXD.3(), W2.3(), W2.4()))

plot(chinDDT_drm, log = '', xlab = expression(paste("DDT Concentration (",mu,"M)")), 
     ylab = "Mortality (%)",
     type = "confidence", #adds 95% confidence intervals to plot
     main = "DDT Toxicity to Chinook Salmon",
     lwd = 2, cex.lab = 1.3, cex.axis = 1.4, xlim = c(0,1250), ylim = c(0,100))

points(percentmortality~conc, data = chinook[50:54, 1:2], pch = 1, lwd = 1.5) #plot data 
lines(chinDDT_pred_ints[,2]~x12, lty = 2) #plot lower prediction bounds
lines(chinDDT_pred_ints[,3]~x12, lty = 2) #plot upper prediction bounds

legend("bottomright", c(expression("153.4 μM LC"[20]), expression("238.6 μM LC"[50]), expression("808.4 μM LC"[99]),
                    "95% Prediction Interval", "95% Confidence Interval"), #change after running ED code below
       lty = c(NA, NA, NA, 2, NA), fill = c(NA, NA, NA, NA,"gray65"), border = "white", bty = 'n', 
       merge = TRUE, cex = 1.3)

ED(chinDDT_drm, c(20,50,99), interval = "delta")


## 2,4-D Butoxyethanol Ester ##
chin24d_drm = drm(percentmortality~conc, data = chinook[55:60, 1:2], fct = LN.4())
x13 = seq(0,1, by = 0.001)
new13 = data.frame(dose = x13)
chin24d_pred_ints = predict(chin24d_drm, interval = "prediction", 
                            newdata = new13) #calculates 95% prediction intervals for curve

# Selecting model: best is LN.4
mselect(chin24d_drm, list(LL.5(), LL.4(), LL.3(), LL.2(), W1.4(), W1.3(), AR.2(), AR.3(), LN.2(), LN.3(), 
                          LN.4(), EXD.2(), EXD.3(), W2.3(), W2.4()))

plot(chin24d_drm, log = '', xlab = expression(paste("2,4-D Butoxyethanol Ester Concentration (",mu,"M)")), 
     ylab = "Mortality (%)",
     type = "confidence", #adds 95% confidence intervals to plot
     main = "2,4-D Toxicity to Chinook Salmon",
     lwd = 2, cex.lab = 1.3, cex.axis = 1.4, xlim = c(0,1), ylim = c(0,100))

points(percentmortality~conc, data = chinook[55:60, 1:2], pch = 1, lwd = 1.5) #plot data 
lines(chin24d_pred_ints[,2]~x13, lty = 2) #plot lower prediction bounds
lines(chin24d_pred_ints[,3]~x13, lty = 2) #plot upper prediction bounds

legend("bottomright", c(expression("0.212 μM LC"[20]), expression("0.245 μM LC"[50]), expression("NA μM LC"[99]),
                        "95% Prediction Interval", "95% Confidence Interval"), #change after running ED code below
       lty = c(NA, NA, NA, 2, NA), fill = c(NA, NA, NA, NA,"gray65"), border = NA, bty = 'n', 
       merge = TRUE, cex = 1.3)

ED(chin24d_drm, c(20,50,99), interval = "delta")


## Malathion ##
chinMal_drm = drm(percentmortality~conc, data = chinook[61:96, 1:2], fct = LL.3())
x14 = seq(0,1, by = 0.001)
new14 = data.frame(dose = x14)
chinMal_pred_ints = predict(chinMal_drm, interval = "prediction", 
                            newdata = new14) #calculates 95% prediction intervals for curve

# Selecting model: best is LL.3
mselect(chinMal_drm, list(LL.5(), LL.4(), LL.3(), LL.2(), W1.4(), W1.3(), AR.2(), AR.3(), LN.2(), LN.3(), 
                          LN.4(), EXD.2(), EXD.3(), W2.3(), W2.4()))

plot(chinMal_drm, log = '', xlab = expression(paste("Malathion Concentration (",mu,"M)")), 
     ylab = "Mortality (%)",
     type = "confidence", #adds 95% confidence intervals to plot
     main = "Malathion Toxicity to Chinook Salmon",
     lwd = 2, cex.lab = 1.3, cex.axis = 1.4, xlim = c(0,1), ylim = c(0,100))

points(percentmortality~conc, data = chinook[61:96, 1:2], pch = 1, lwd = 1.5) #plot data 
lines(chinMal_pred_ints[,2]~x14, lty = 2) #plot lower prediction bounds
lines(chinMal_pred_ints[,3]~x14, lty = 2) #plot upper prediction bounds

legend("bottomright", c(expression("0.318 μM LC"[20]), expression("0.397 μM LC"[50]), expression("0.832 μM LC"[99]),
                        "95% Prediction Interval", "95% Confidence Interval"), #change after running ED code below
       lty = c(NA, NA, NA, 2, NA), fill = c(NA, NA, NA, NA,"gray65"), border = "white", bty = 'n', 
       merge = TRUE, cex = 1.3)

ED(chinMal_drm, c(20,50,99), interval = "delta")


## Getting all EC20, 50, and 99 values (might be inaccurate)
ALLchin_drm <- drm(percentmortality~conc, toxicant, data = chinook, fct = c(W2.3(), LL.3(), W1.3(), AR.2(), W1.4(),
                                                                            AR.3(), W1.4(), LN.4(), LN.4(), LL.3()))
ED(ALLchin_drm, c(20,50,99), interval="delta", level = 0.95)
