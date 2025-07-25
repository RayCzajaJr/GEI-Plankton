

### make partial effects plot for each relationship from all final GAMs
### get partial deviance explained values for each relationship from all final GAMs
### get DHARMa residual plots for each model


library(mgcv)
library(visreg)
library(ggplot2)
library(gridExtra)
library(DHARMa)
library(visreg)

#Auxis gam 1
AUXGam1<- gam(AUXIS~ s(Calanoids, k=3,m=1) + s(Larvaceans, k=3)+
               s(lon, lat),
             data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(AUXGam1)

res1 <- simulateResiduals(AUXGam1)
plot(res1)

common_params1 <- list(
  partial = FALSE,
  rug = FALSE,
  band = TRUE,
  line = list(col = "darkgoldenrod"),
  fill = list(col = rgb(240/255, 230/255, 140/255, alpha = 0.6)),
  type = c("contrast")
)
setwd("~/Desktop")

png("AuxisGamPlot.png", width = 2400, height = 1600, res = 200)  

par(mfrow = c(1, 2),
    oma = c(5, 4, 1, 1),
    mar = c(4, 4, 1, 1),
    cex.lab = 1.3,
    mgp = c(2, 0.5, 0)) 

do.call(visreg, c(list(AUXGam1, "Calanoids"), common_params1,
                  list(xlab = "Calanoid Abundance", ylab = expression("Frigate Tuna Partial Effects"))))
do.call(visreg, c(list(AUXGam1, "Larvaceans"), common_params1,
                  list(xlab = "Larvacean Abundance", ylab = expression("Frigate Tuna Partial Effects"))))
print(plot) 
dev.off() 

#Auxis gam 2
AUXGam2<- gam(AUXIS~s(TEMPSURF, k=3, m=1)+s(TEMPMAX, k=3)+
               s(lon, lat),
             data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(AUXGam2)

res2 <- simulateResiduals(AUXGam2)
plot(res2)

common_params2 <- list(
  partial = FALSE,
  rug = FALSE,
  band = TRUE,
  line = list(col = "royalblue4"),
  fill = list(col = "lightblue1"),
  type = c("contrast")
)
setwd("~/Desktop")

png("AuxisGamPlot2.png", width = 2400, height = 1600, res = 200)  

par(mfrow = c(1, 2),
    oma = c(5, 4, 1, 1),
    mar = c(4, 4, 1, 1),
    cex.lab = 1.3,
    mgp = c(2, 0.5, 0)) 

do.call(visreg, c(list(AUXGam2, "TEMPSURF"), common_params2,
                  list(xlab = "Surface Temperature (°C)", ylab = expression("Frigate Tuna Partial Effects"))))
do.call(visreg, c(list(AUXGam2, "TEMPMAX"), common_params2,
                  list(xlab = "Temperature at Depth (°C)", ylab = expression("Frigate Tuna Partial Effects"))))
print(plot) 
dev.off()  



#Mahi mahi 

#mahi mahi gam
CORYGam1<- gam(CORYPHAENA~s(Larvaceans, k=3,m=1)+ s(lon, lat),
              data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(CORYGam1)
res1 <- simulateResiduals(CORYGam1)
plot(res1)

setwd("~/Desktop")

png("CoryGamPlot1.png", width = 2400, height = 1600, res = 200)  

par(
  oma = c(5, 4, 1, 1),
  mar = c(4, 4, 1, 1),
  cex.lab = 1.3,
  mgp = c(2, 0.5, 0)) 

do.call(visreg, c(list(CORYGam1, "Larvaceans"), common_params1,
                  list(xlab = "Larvacean Abundance", ylab = expression("Mahi-Mahi Partial Effects"))))
print(plot) 
dev.off()  


#mahi mahi gam
CORYGam2<- gam(CORYPHAENA~s(TEMPSURF, k=3, m=1) + s(SALSURF, k=3, m=1),
             data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(CORYGam2)
res2 <- simulateResiduals(CORYGam2)
plot(res2)

setwd("~/Desktop")

png("CoryGamPlot2.png", width = 2400, height = 1600, res = 200)  

par(mfrow = c(1, 2),
    oma = c(5, 4, 1, 1),
    mar = c(4, 4, 1, 1),
    cex.lab = 1.3,
    mgp = c(2, 0.5, 0)) 

do.call(visreg, c(list(CORYGam2, "TEMPSURF"), common_params2,
                  list(xlab = "Surface Temperature (°C)", ylab = expression("Mahi-Mahi Partial Effects"))))
do.call(visreg, c(list(CORYGam2, "SALSURF"), common_params2,
                  list(xlab = "Surface Salinity", ylab = expression("Mahi-Mahi Partial Effects"))))
print(plot) 
dev.off()  

# Little tunny

# Fit the Gam
EUTGam1 <- gam(EUTHYNNUSALLETTERATUS ~ s(Calanoids, k=4, m=1) + s(Cyclopoids, k=3, m=1) + s(Larvaceans, k=3) + s(lon, lat),
              data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(EUTGam1)

res1 <- simulateResiduals(EUTGam1)
plot(res1)

setwd("~/Desktop")

png("TunnyGamPlot1.png", width = 2400, height = 1600, res = 200)  

par(mfrow = c(2, 2),
    oma = c(5, 4, 1, 1),
    mar = c(4, 4, 1, 1),
    cex.lab = 1.3,
    mgp = c(2, 0.5, 0)) 

do.call(visreg, c(list(EUTGam1, "Calanoids"), common_params1,
                  list(xlab = "Calanoid Abundance", ylab = expression("Little Tunny Partial Effects"))))
do.call(visreg, c(list(EUTGam1, "Cyclopoids"), common_params1,
                  list(xlab = "Cyclopoid Abundance", ylab = expression("Little Tunny Partial Effects"))))
do.call(visreg, c(list(EUTGam1, "Larvaceans"), common_params1,
                  list(xlab = "Larvacean Abundance", ylab = expression("Little Tunny Partial Effects"))))
print(plot) 
dev.off()  

# Fit the Gam
EUTGam2 <- gam(EUTHYNNUSALLETTERATUS ~s(TEMPSURF, k=3, m=1)+s(TEMPMAX, k=3,m=1)+
                 s(lon, lat),
              data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(EUTGam2)

res2 <- simulateResiduals(EUTGam2)
plot(res2)

setwd("~/Desktop")

png("TunnyGamPlot2.png", width = 2400, height = 1600, res = 200)  

par(mfrow = c(1, 2),
    oma = c(5, 4, 1, 1),
    mar = c(4, 4, 1, 1),
    cex.lab = 1.3,
    mgp = c(2, 0.5, 0)) 

do.call(visreg, c(list(EUTGam2, "TEMPSURF"), common_params2,
                  list(xlab = "Surface Temperature (°C)", ylab = expression("Little Tunny Partial Effects"))))
do.call(visreg, c(list(EUTGam2, "TEMPMAX"), common_params2,
                  list(xlab = "Temperature at Depth (°C)",  ylab = expression("Little Tunny Partial Effects"))))
print(plot) 
dev.off()  

#MBillfish 

#bilfish gam
springplankton5 <- springplankton4 %>%
  filter(!(ISTIOPHORIDAE > 0 & TEMPSURF < 25))

#fit climate gam
BILLGam2<- gam(ISTIOPHORIDAE~s(TEMPSURF, k=3, m=1)  + s(TEMPMAX, k=3) + s(lon, lat),
               data =springplankton5, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(BILLGam2)

res2 <- simulateResiduals(BILLGam2)
plot(res2)

setwd("~/Desktop")

png("BillGamPlot2.png", width = 2400, height = 1600, res = 200)  

par(mfrow = c(1,2),
    oma = c(5,4,1,1),
    mar = c(4,4,1,1),
    cex.lab = 1.3,
    mgp = c(2, 0.5, 0)) 

do.call(visreg, c(list(BILLGam2, "TEMPSURF"), common_params2,
                  list(xlab = "Surface Temperature (°C)", ylab = expression("Billfish Partial Effects"))))
do.call(visreg, c(list(BILLGam2, "TEMPMAX"), common_params2,
                  list(xlab = "Temperature at Depth (°C)", ylab = expression("Billfish Partial Effects"))))
print(plot) 
dev.off()  

#fit prey gam
BILLGam1<- gam(ISTIOPHORIDAE~ s(Cladocerans, k=3) + s(lon, lat),
              data =springplankton5, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(BILLGam1)

res1 <- simulateResiduals(BILLGam1)
plot(res1)

setwd("~/Desktop")

png("BillGamPlot1.png", width = 2400, height = 1600, res = 200)  

par(mfrow = c(1,1),
    oma = c(5,4,1,1),
    mar = c(4,4,1,1),
    cex.lab = 1.3,
    mgp = c(2, 0.5, 0)) 

do.call(visreg, c(list(BILLGam1, "Cladocerans"), common_params1,
                  list(xlab = "Cladoceran Abundance", ylab = expression("Billfish Partial Effects"))))
print(plot) 
dev.off()  

#Skipjack tuna

#fit the prey gam
KATSUGam1<- gam(KATSUWONUSPELAMIS~ s(Larvaceans) + s(lon, lat),
              data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(KATSUGam1)

res1 <- simulateResiduals(KATSUGam1)
plot(res1)

setwd("~/Desktop")

png("SkipjackGamPlot1.png", width = 2400, height = 1600, res = 200)  

par(
  oma = c(5, 4, 1, 1),
  mar = c(4, 4, 1, 1),
  cex.lab = 1.3,
  mgp = c(2, 0.5, 0)) 

do.call(visreg, c(list(KATSUGam1, "Larvaceans"), common_params1,
                  list(xlab = "Larvacean Abundance", ylab = expression("Skipjack Tuna Partial Effects"))))
print(plot) 
dev.off() 

#fit the climate gam
KATSUGam2<- gam(KATSUWONUSPELAMIS~s(CHLORSURF, k=3) + s(SALSURF, k=3, m=1) + s(TEMPSURF, k=3, m=1),
                data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(KATSUGam2)

res2 <- simulateResiduals(KATSUGam2)
plot(res2)

setwd("~/Desktop")

png("SkipjackGamPlot2.png", width = 2400, height = 1600, res = 200)  

par(mfrow = c(2, 2),
    oma = c(5, 4, 1, 1),
    mar = c(4, 4, 1, 1),
    cex.lab = 1.3,
    mgp = c(2, 0.5, 0)) 

do.call(visreg, c(list(KATSUGam2, "TEMPSURF"), common_params2,
                  list(xlab = "Surface Temperature (°C)", ylab = expression("Skipjack Tuna Partial Effects"))))
do.call(visreg, c(list(KATSUGam2, "SALSURF"), common_params2,
                  list(xlab = "Surface Salinity", ylab = expression("Skipjack Tuna Partial Effects"))))
do.call(visreg, c(list(KATSUGam2, "CHLORSURF"), common_params2,
                  list(xlab = "Surface Chlorophyll (μg/L)", ylab = expression("Skipjack Tuna Partial Effects"))))
print(plot) 
dev.off()  



#other true tuna

# Fit the prey gam
THUNNUSGam1 <- gam(THUNNUS ~s(Larvaceans, k=3, m=1) +s(Cyclopoids, k=3),
                   data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(THUNNUSGam1)

res1 <- simulateResiduals(THUNNUSGam1)
plot(res1)

setwd("~/Desktop")

png("thunnusGamPlot1.png", width = 2400, height = 1600, res = 200)  

par(mfrow = c(1, 2),
    oma = c(5, 4, 1, 1),
    mar = c(4, 4, 1, 1),
    cex.lab = 1.3,
    mgp = c(2, 0.5, 0)) 

do.call(visreg, c(list(THUNNUSGam1, "Larvaceans"), common_params1,
                  list(xlab = "Larvacean Abundance", ylab = expression("Other True Tuna Partial Effects"))))
do.call(visreg, c(list(THUNNUSGam1, "Cyclopoids"), common_params1,
                  list(xlab = "Cyclopoid Abundance", ylab = expression("Other True Tuna Partial Effects"))))
print(plot) 
dev.off() 

# Fit the climate gam
THUNNUSGam2 <- gam(THUNNUS ~ s(TEMPSURF, k=3, m=1) + s(SALSURF, k=3, m=1) + s(TEMPMAX, k=3) +s(CHLORSURF, k=3)+ s(lon, lat),
                  data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(THUNNUSGam2)

res2 <- simulateResiduals(THUNNUSGam2)
plot(res2)

setwd("~/Desktop")

png("thunnusGamPlot2.png", width = 2400, height = 1600, res = 200)  

par(mfrow = c(2, 2),
    oma = c(5, 4, 1, 1),
    mar = c(4, 4, 1, 1),
    cex.lab = 1.3,
    mgp = c(2, 0.5, 0)) 

do.call(visreg, c(list(THUNNUSGam2, "TEMPSURF"), common_params2,
                  list(xlab = "Surface Temperature (°C)", ylab = expression("Other True Tuna Partial Effects"))))
do.call(visreg, c(list(THUNNUSGam2, "SALSURF"), common_params2,
                  list(xlab = "Surface Salinity", ylab = expression("Other True Tuna Partial Effects"))))
do.call(visreg, c(list(THUNNUSGam2, "TEMPMAX"), common_params2,
                  list(xlab = "Temperature at Depth (°C)", ylab = expression("Other True Tuna Partial Effects"))))
do.call(visreg, c(list(THUNNUSGam2, "CHLORSURF"), common_params2,
                  list(xlab = "Surface Chlorophyll (μg/L)", ylab = expression("Other True Tuna Partial Effects"))))
print(plot) 
dev.off() 


#ABFT
#fit the prey  gam
ABFTGam1 <- gam(THUNNUSTHYNNUS~  s(Larvaceans, k=3) +s(Cladocerans, k=4, m=1) +s(lon, lat)+s(Cyclopoids, k=3),
                data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(ABFTGam1)

res1 <- simulateResiduals(ABFTGam1)
plot(res1)

setwd("~/Desktop")

png("abftGamPlot1.png", width = 2400, height = 1600, res = 200)  

par(mfrow = c(2,2),
    oma = c(5,4,1,1),
    mar = c(4,4,1,1),
    cex.lab = 1.3,
    mgp = c(2, 0.5, 0)) 

do.call(visreg, c(list(ABFTGam1, "Larvaceans"), common_params1,
                  list(xlab = "Larvacean Abundance", ylab = expression("Atlantic Bluefin Tuna Partial Effects"))))
do.call(visreg, c(list(ABFTGam1, "Cladocerans"), common_params1,
                  list(xlab = "Cladoceran Abundance", ylab = expression("Atlantic Bluefin Tuna Partial Effects"))))
do.call(visreg, c(list(ABFTGam1, "Cyclopoids"), common_params1,
                  list(xlab = "Cyclopoid Abundance", ylab = expression("Atlantic Bluefin Tuna Partial Effects"))))
print(plot) 
dev.off() 

#fit the climate  gam
ABFTGam2 <- gam(THUNNUSTHYNNUS~  s(TEMPSURF, k=8) +s(CHLORSURF, k=3, m=1) +s(TEMPMAX,k=3,m=1)+s(lon, lat),
               data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(ABFTGam2)

res2 <- simulateResiduals(ABFTGam2)
plot(res2)

setwd("~/Desktop")

png("abftGamPlot2.png", width = 2400, height = 1600, res = 200)  

par(mfrow = c(2,2),
    oma = c(5,4,1,1),
    mar = c(4,4,1,1),
    cex.lab = 1.3,
    mgp = c(2, 0.5, 0)) 

do.call(visreg, c(list(ABFTGam2, "TEMPSURF"), common_params2,
                  list(xlab = "Surface Temperature (°C)", ylab = expression("Atlantic Bluefin Tuna Partial Effects"))))
do.call(visreg, c(list(ABFTGam2, "TEMPMAX"), common_params2,
                  list(xlab = "Temperature at Depth (°C)", ylab = expression("Atlantic Bluefin Tuna Partial Effects"))))
do.call(visreg, c(list(ABFTGam2, "CHLORSURF"), common_params2,
                  list(xlab = "Surface Chlorophyll (μg/L)", ylab = expression("Atlantic Bluefin Tuna Partial Effects"))))
print(plot) 
dev.off() 


#calanoids
CALAGam <- gam(Calanoids ~ s(TEMPSURF, k=3, m=1) + s(SALSURF, k=3) + s(TEMPMAX, k=3, m=1) +s(CHLORSURF, k=3, m=1)+ s(lon, lat) +
                 s(TIME_EMIL, bs = "cc"), data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
plot(CALAGam)

res2 <- simulateResiduals(CALAGam)
plot(CALAGam)

setwd("~/Desktop")

png("CalaGamPlot.png", width = 2400, height = 1600, res = 200)  

par(mfrow = c(2, 2),
    oma = c(5, 4, 1, 1),
    mar = c(4, 4, 1, 1),
    cex.lab = 1.3,
    mgp = c(2, 0.5, 0)) 

do.call(visreg, c(list(CALAGam, "TEMPSURF"), common_params2,
                  list(xlab = "Surface Temperature (°C)", ylab = expression("Calanoid Partial Effects"))))
do.call(visreg, c(list(CALAGam, "SALSURF"), common_params2,
                  list(xlab = "Surface Salinity", ylab = expression("Calanoid Partial Effects"))))
do.call(visreg, c(list(CALAGam, "TEMPMAX"), common_params2,
                  list(xlab = "Temperature at Depth (°C)", ylab = expression("Calanoid Partial Effects"))))
do.call(visreg, c(list(CALAGam, "CHLORSURF"), common_params2,
                  list(xlab = "Surface Chlorophyll (μg/L)", ylab = expression("Calanoid Partial Effects"))))
print(plot) 
dev.off() 

#cyclopoids
CYCGam <- gam(Cyclopoids ~ s(TEMPSURF, k=3, m=1) + s(SALSURF, k=3) + s(TEMPMAX, k=3, m=1) +s(CHLORSURF, k=3)+ s(lon, lat)+s(TIME_EMIL, bs="cs"),
               data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
res2 <- simulateResiduals(CYCGam)
plot(CYCGam)

setwd("~/Desktop")

png("CycGamPlot.png", width = 2400, height = 1600, res = 200)  

par(mfrow = c(2, 2),
    oma = c(5, 4, 1, 1),
    mar = c(4, 4, 1, 1),
    cex.lab = 1.3,
    mgp = c(2, 0.5, 0)) 

do.call(visreg, c(list(CYCGam, "TEMPSURF"), common_params2,
                  list(xlab = "Surface Temperature (°C)", ylab = expression("Cyclopoid Partial Effects"))))
do.call(visreg, c(list(CYCGam, "SALSURF"), common_params2,
                  list(xlab = "Surface Salinity", ylab = expression("Cyclopoid Partial Effects"))))
do.call(visreg, c(list(CYCGam, "TEMPMAX"), common_params2,
                  list(xlab = "Temperature at Depth (°C)", ylab = expression("Cyclopoid Partial Effects"))))
do.call(visreg, c(list(CYCGam, "CHLORSURF"), common_params2,
                  list(xlab = "Surface Chlorophyll (μg/L)", ylab = expression("Cyclopoid Partial Effects"))))
print(plot) 
dev.off() 

#larvaceans
LARVGam <- gam(Larvaceans~ s(TEMPSURF, k=3, m=1) + s(TEMPMAX, k=3, m=1) + s(lon, lat),
              data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
res2 <- simulateResiduals(LARVGam)
plot(res2)

setwd("~/Desktop")

png("LarvGamPlot.png", width = 2400, height = 1600, res = 200)  

par(mfrow = c(1, 2),
    oma = c(5, 4, 1, 1),
    mar = c(4, 4, 1, 1),
    cex.lab = 1.3,
    mgp = c(2, 0.5, 0)) 

do.call(visreg, c(list(LARVGam, "TEMPSURF"), common_params2,
                  list(xlab = "Surface Temperature (°C)", ylab = expression("Larvacean Partial Effects"))))
do.call(visreg, c(list(LARVGam, "TEMPMAX"), common_params2,
                  list(xlab = "Temperature at Depth (°C)", ylab = expression("Larvacean Partial Effects"))))
print(plot) 
dev.off() 

#cladocerans
CLADGam <- gam(Cladocerans~ s(TEMPSURF, k=3) + s(SALSURF, k=3)  + s(lon, lat),
               data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
res2 <- simulateResiduals(CLADGam)
plot(res2)

setwd("~/Desktop")

png("CLADGamPlot.png", width = 2400, height = 1600, res = 200)  

par(mfrow = c(1, 2),
    oma = c(5, 4, 1, 1),
    mar = c(4, 4, 1, 1),
    cex.lab = 1.3,
    mgp = c(2, 0.5, 0)) 

do.call(visreg, c(list(CLADGam, "TEMPSURF"), common_params2,
                  list(xlab = "Surface Temperature (°C)", ylab = expression("Cladoceran Partial Effects"))))
do.call(visreg, c(list(CLADGam, "SALSURF"), common_params2,
                  list(xlab = "Surface Salinity", ylab = expression("Cladoceran Partial Effects"))))
print(plot) 
dev.off() 



### get partial DE values

calculate_gam_contributions <- function(gam_model, data) {
  library(mgcv)
  
  # Check if the model is a GAM
  if (!inherits(gam_model, "gam")) {
    stop("The input model must be a GAM (fitted with mgcv::gam).")
  }
  
  # Extract the response variable and predictors from the full model
  response <- as.character(attr(terms(gam_model), "variables"))[2]
  full_formula <- formula(gam_model)
  
  # Deviance of the null model (intercept-only model)
  null_model <- gam(as.formula(paste(response, "~ 1")), 
                    data = data, 
                    family = gam_model$family, 
                    method = gam_model$method)
  null_deviance <- deviance(null_model)
  
  # Deviance of the full model
  full_deviance <- deviance(gam_model)
  
  # Extract smooth terms
  smooth_terms <- attr(terms(full_formula), "term.labels")
  
  # Retrieve summary statistics for the full model
  model_summary <- summary(gam_model)
  edf_values <- model_summary$s.table[, "edf"]
  p_values <- model_summary$s.table[, "p-value"]
  
  # Calculate partial deviance explained for each term
  partial_deviances <- sapply(smooth_terms, function(term) {
    # Create a reduced formula without the current term
    reduced_terms <- setdiff(smooth_terms, term)
    reduced_formula <- as.formula(paste(response, "~", paste(reduced_terms, collapse = " + ")))
    
    # Fit the reduced model
    reduced_model <- gam(reduced_formula, 
                         data = data, 
                         family = gam_model$family, 
                         method = gam_model$method)
    
    # Deviance of the reduced model
    reduced_deviance <- deviance(reduced_model)
    
    # Partial deviance explained
    (reduced_deviance - full_deviance) / null_deviance
  })
  
  # Convert PDE to percentages
  partial_deviances <- partial_deviances * 100
  
  # Combine results into a data frame
  results <- data.frame(
    Predictor = smooth_terms,
    EDF = edf_values,
    P_Value = p_values,
    Partial_Deviance_Explained = partial_deviances
  )
  
  # Format the results for clarity
  results <- results[order(-results$Partial_Deviance_Explained), ] # Order by PDE
  rownames(results) <- NULL  # Clean up row names
  
  return(results)
}

#Auxis 
results <- calculate_gam_contributions(AUXGam, data = springplankton4)
print(results)

results <- calculate_gam_contributions(AUXGam2, data = springplankton4)
print(results)

#Mahi Mahi
results <- calculate_gam_contributions(CORYGam, data = springplankton4)
print(results)

results <- calculate_gam_contributions(CORYGam2, data = springplankton4)
print(results)

#Thunnus
results <- calculate_gam_contributions(THUNGam, data = springplankton4)
print(results)

results <- calculate_gam_contributions(THUNGam2, data = springplankton4)
print(results)

#Skipjack
results <- calculate_gam_contributions(KATSUGam, data = springplankton4)
print(results)

results <- calculate_gam_contributions(KATSUGam2, data = springplankton4)
print(results)

#Tunny
results <- calculate_gam_contributions(EUTHGam, data = springplankton4)
print(results)

results <- calculate_gam_contributions(EUTHGam2, data = springplankton4)
print(results)

#ABFT
results <- calculate_gam_contributions(THYNGam, data = springplankton4)
print(results)

results <- calculate_gam_contributions(THYNGam2, data = springplankton4)
print(results)

#billfish
results <- calculate_gam_contributions(ISTGam, data = springplankton4)
print(results)

results <- calculate_gam_contributions(ISTGam2, data = springplankton4)
print(results)

#calanoids
results <- calculate_gam_contributions(CalGam, data = springplankton4)
print(results)

#cyclopoids
results <- calculate_gam_contributions(CycGam, data = springplankton4)
print(results)

# the function is not working well for cladocerans and larvaceans, so lets get those values manually
#cladocerans
results <- calculate_gam_contributions(CladGam, data = springplankton4)
print(results)
CladGam <- gam(Cladocerans~ s(TEMPSURF, k=3) + s(SALSURF, k=3)+ s(lon, lat),
               data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(CladGam)

CladGamNoTemp <- gam(Cladocerans~ s(SALSURF, k=3)+ s(lon, lat),
               data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(CladGamNoTemp) #42.8 - 40.4 = 2.4

CladGamNoSal <- gam(Cladocerans~ s(TEMPSURF, k=3)+ s(lon, lat),
                     data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(CladGamNoSal) #42.8 - 37.1 = 5.7

CladGamNoGPS <- gam(Cladocerans~ s(TEMPSURF, k=3) + s(SALSURF, k=3),
               data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(CladGamNoGPS) #42.8 - 24.2 = 18.6


#larvaceans
results <- calculate_gam_contributions(LarvGam, data = springplankton4)
print(results)
LarvGam<- gam(Larvaceans~ s(TEMPSURF, k=3, m=1) + s(TEMPMAX, k=3, m=1) + s(lon, lat),
              data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(LarvGam)
LarvGamNoSurf<- gam(Larvaceans~ s(TEMPMAX, k=3, m=1) + s(lon, lat),
              data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(LarvGamNoSurf) # 21.2 - 15.1 = 6.1
LarvGamNoMax<- gam(Larvaceans~ s(TEMPSURF, k=3, m=1) + s(lon, lat),
                   data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(LarvGamNoMax) # 21.2 - 11.2 = 10.0
LarvGamNoGPS<- gam(Larvaceans~ s(TEMPSURF, k=3, m=1) + s(TEMPMAX, k=3, m=1),
              data =springplankton4, tw(),method = "REML", select = TRUE, na.action = "na.fail")
summary(LarvGamNoGPS) # 21.2 - 13.9 = 7.3


########################################



#figure out times of day the predict mean values 

# Function to find the predictor value corresponding to the average response
find_time_for_avg_response <- function(model, data, time_var, other_vars) {
  # Arguments:
  # model      : Fitted GAM model (e.g., from mgcv::gam())
  # data       : Data frame used for the model fitting
  # time_var   : The name of the predictor variable of interest (e.g., "TIME_EMIL")
  # other_vars : Vector of names of other predictors to hold constant
  
  # Create a sequence for the time variable
  time_seq <- seq(min(data[[time_var]], na.rm = TRUE), 
                  max(data[[time_var]], na.rm = TRUE), 
                  length.out = 100)
  
  # Calculate mean values for the other predictors
  means <- sapply(other_vars, function(var) mean(data[[var]], na.rm = TRUE))
  
  # Create a prediction data frame
  pred_data <- data.frame(matrix(ncol = length(other_vars) + 1, nrow = length(time_seq)))
  colnames(pred_data) <- c(other_vars, time_var)
  
  # Populate the data frame with mean values for other predictors and time sequence
  for (var in other_vars) {
    pred_data[[var]] <- means[var]
  }
  pred_data[[time_var]] <- time_seq
  
  # Predict the response variable
  predicted <- predict(model, newdata = pred_data, type = "response")
  
  # Find the value of the time variable corresponding to the average response
  avg_response <- mean(predicted)
  closest_time <- time_seq[which.min(abs(predicted - avg_response))]
  
  # Return the closest time value
  return(closest_time)
}

calatime_var <- "TIME_EMIL"
calaother_vars <- c("TEMPSURF", "SALSURF", "TEMPMAX","CHLORSURF", "lon", "lat")
calaavg_time <- find_time_for_avg_response(CALAGam, springplankton4, calatime_var, calaother_vars)
print(calaavg_time)
