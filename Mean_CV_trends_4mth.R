#######################################################################################################
#####            This code is for the figure entitled "Mean and CV trends before death            ##### 
#####                                for the 4-month variable list"                               #####
#######################################################################################################
library(RColorBrewer)

rm(list = ls())

source("Functions_HD_variability_paper.R")

#### Load data ####
clean.HD$date_visit <- as.Date(clean.HD$date_visit)
clean.HD$date_death <- as.Date(clean.HD$date_death)


## variable list
selected.biomarkersHD <- c("wbc", "hb", "hct", "mch", "mchc", "mcv", "pltlt", "rdw", "ca", "rbc", "creat", 
                           "Glucose", "k", "Phosphates.sériques", "sodium", "alb")
nms.4months = c("WBC", "Hemoglobin", "Hematocrit", "MCH", "MCHC", "MCV", "Platelets",
                "RDW", "Calcium", "RBC", "Creatinine", "Glucose", "Potassium", 
                "Phosphate","Sodium", "Albumin")

# Select complete observations
clean.HD2 <- clean.HD[complete.cases(clean.HD[ , selected.biomarkersHD]), ]

# Transform variables not normally distributed
if ("pltlt" %in% selected.biomarkersHD) clean.HD2$pltlt <- sqrt(clean.HD2$pltlt)
clean.HD2[ , which(colnames(clean.HD2) %in% c("wbc", "rdw", "Glucose"))] <- 
  log(clean.HD2[ , which(colnames(clean.HD2) %in% c("wbc", "rdw", "Glucose"))])

# Add date of last contact
dat <- ddply(clean.HD2, .(id_no), lastvis)

### Plot
tiff("Mean_CV_trends_4mth.tiff", width = 270, height = 180, "mm", res = 300, compression = "lzw")
opar <- par(mfrow = c(2, 3), mar = c(3, 3, 0.5, 0.5), oma = c(0.2, 0.2, 0.2, 0.2), 
            mgp = c(1.5, 0.4, 0), tcl = -0.3, cex = 1.1)

# Set color palette 
cols <- colorRampPalette(brewer.pal(12, "Paired"))
myPal <- rev(cols(length(selected.biomarkersHD)))

###### Panels A and B (means and CVs)
# Calculate biomarker means per year
dat0 <- dat
dat0[ , selected.biomarkersHD] <- scale(dat0[ , selected.biomarkersHD])
dat.means <- ddply(dat0, .(id_no), cv.year, vars = selected.biomarkersHD, cv = F, min.vis.nb = 2, 
                   covars = c("id_no", "sex", "date_death", "date_birth", "date_last", "diabetes", "fu_length"))

# Add "mean" to column names
colnames(dat.means)[which(colnames(dat.means) %in% selected.biomarkersHD)] <- 
  unlist(lapply(colnames(dat.means)[which(colnames(dat.means) %in% selected.biomarkersHD)], 
                function (x) paste(x, "mean", sep = "_")))

# Calculate biomarker CVs per year
dat.cvs <- ddply(dat, .(id_no), cv.year, vars = selected.biomarkersHD, cv = T, min.vis.nb = 2, 
                 covars = c("id_no", "sex", "date_death", "date_birth", "date_last", "diabetes", "fu_length"))

# Apply log transformation on CVs
dat.cvs[ , selected.biomarkersHD] <- apply(dat.cvs[ , selected.biomarkersHD], 2, log_zeros)

# Some variables have no variation, leading to outliers
dat.cvs[ , selected.biomarkersHD] <- apply(dat.cvs[ , selected.biomarkersHD], 2, uniform_val)

# Add "CV" to column names
colnames(dat.cvs)[which(colnames(dat.cvs) %in% selected.biomarkersHD)] <- 
  unlist(lapply(colnames(dat.cvs)[which(colnames(dat.cvs) %in% selected.biomarkersHD)], 
                function (x) paste(x, "CV", sep = "_")))

# Combine means with CVs
dat.both <- merge(dat.means, 
                  dat.cvs[ , -which(colnames(dat.cvs) %in% c("sex","date_death", "date_birth", "date_last", 
                                                             "diabetes", "fu_length", "age_visit", "nb_vis", 
                                                             "vis_sd", "status"))],
                  by = c("id_no", "years"))


# Cut time at 5 years before death
dat.5yr <- dat.both[dat.both$years <= 5 & !is.na(dat.both$date_death), ]
dat.5yr$time <- dat.5yr$years

for (i in 1:2) {
  
  if (i==1) datt <- dat.5yr[ , c("time", grep("mean", colnames(dat.5yr), value = T))] else
    datt <- dat.5yr[ , c("time", grep("CV", colnames(dat.5yr), value = T))]
  colnames(datt) <- c("time", selected.biomarkersHD)
  
  # Calculate mean and 95%CI per time point for each variable
  time <- unique(dat.5yr$time)[order(unique(dat.5yr$time))]
  means <- list()
  cis <- list()
  for (var in 1:length(selected.biomarkersHD)) {
    means_raw <- tapply(datt[ , 1 + var], datt$time, mean)
    cis[[var]] <- tapply(datt[ , 1 + var], datt$time, 
                         function(x) qnorm(0.975) * sd(x) / sqrt(length(x)))
    
    # Center means so that values at 5 years start at 0
    means[[var]] <- means_raw - means_raw[which(names(means_raw)==5)]
  }
  
  # Find y-axis limits
  ylims <- range(rbind(mapply(function(x, y) x - y, means, cis), 
                       mapply(function(x, y) x + y, means, cis)))
  if (i==1) ylims <- c(ylims[1] - 0.1, ylims[2] + 0.1)
  
  plot(1:10, 1:10, type = "n", xlim = rev(range(dat.5yr$time)), 
       ylim = ylims, ylab = c("Levels (mean z-scores)", "Variability (CVs)")[i], 
       xlab = "Time before death (years)", frame.plot = F)
  mtext(c("A", "B")[i], side = 3, line = -0.5, at = 6.2, font = 2, cex = 1.3)
  
  for (var in 1:length(selected.biomarkersHD)) {
    points(x = time, y = means[[var]], pch = 16, col = myPal[var], cex = 0.9)
    segments(x0 = time, y0 = means[[var]] + cis[[var]], x1 = time, 
             y1 = means[[var]] - cis[[var]], col = myPal[var])
    segments(x0 = time - 0.02, y0 = means[[var]] + cis[[var]], x1 = time + 0.02, 
             y1 = means[[var]] + cis[[var]], col = myPal[var])
    segments(x0 = time - 0.02, y0 = means[[var]] - cis[[var]], x1 = time + 0.02, 
             y1 = means[[var]] - cis[[var]], col = myPal[var])
    for (j in 1:(length(time) - 1)) {
      segments(x0 = time[j], y0 = means[[var]][j], x1 = time[j + 1], 
               y1 = means[[var]][j + 1], lty = 3, col = myPal[var])
    }
  }
}

## Put legend in an empty plot
par(mar = c(0, 1, 0, 0))
plot(1:10, 1:10, type = "n", xlab = "", xaxt = "n", yaxt = "n", 
     ylab = "", frame.plot = F)
legend("topleft", legend = nms.4months, col = myPal, pch = 16, cex = 0.85, bty = "n")


### Panel C (PC trends)
par(mar = c(3, 3, 0.5, 0.5))
# Run PCA on all values to extract the loadings (which are going to be applied on means per year)
pca.all <- prcomp(dat[ , selected.biomarkersHD], center = T, scale. = T)

# Add "status" variable (for cox regressions)
dat.cox <- ddply(dat.means, .(id_no), status_fct)

# Calculate PC scores using loadings from all values
for (pc in 1:length(selected.biomarkersHD)) {
  dat.cox[ , paste("PC", pc, sep = "")] <- 
    as.matrix(dat.cox[ , paste(selected.biomarkersHD, "mean", sep = "_")]) %*% pca.all$rotation[ , pc]
}
for (pc in 1:length(selected.biomarkersHD)) {
  dat.cox <- pc_sign(dat.cox, paste("PC", pc, sep = ""), pca.all)
}
# Re-order columns (so that PCs are from 1 to 11)
dat.cox <- dat.cox[ , c(setdiff(colnames(dat.cox), 
                                paste("PC", 1:length(selected.biomarkersHD), sep = "")),
                        paste("PC", 1:length(selected.biomarkersHD), sep = ""))]


# Cut time at 5 years before death
dat.5yr <- dat.cox[dat.cox$years <= 5 & !is.na(dat.cox$date_death), ]
dat.5yr$time <- dat.5yr$years

# Calculate mean and 95%CI per time point for each CVPC
time <- unique(dat.5yr$time)[order(unique(dat.5yr$time))]
means <- list()
cis <- list()
for (var in 1:length(selected.biomarkersHD)) {
  means_raw <- tapply(dat.5yr[ , paste("PC", var, sep = "")], dat.5yr$time, mean)
  cis[[var]] <- tapply(dat.5yr[ , paste("PC", var, sep = "")], dat.5yr$time, 
                       function(x) qnorm(0.975) * sd(x) / sqrt(length(x)))
  
  # Center means so that values at 5 years start at 0
  means[[var]] <- means_raw - means_raw[which(names(means_raw)==5)]
}

# Find y-axis limits
ylims <- range(rbind(mapply(function(x, y) x - y, means, cis), 
                     mapply(function(x, y) x + y, means, cis)))
plot(1:10, 1:10, type = "n", xlim = rev(range(dat.5yr$time)), 
     ylim = ylims, ylab = "PCs (1-16)", 
     xlab = "Time before death (years)", frame.plot = F)
mtext("D", side = 3, line = -0.5, at = 6.2, font = 2, cex = 1.3)
for (var in 1:length(means)) {
  points(x = time, y = means[[var]], pch = 16, col = myPal[var], cex = 0.9)
  segments(x0 = time, y0 = means[[var]] + cis[[var]], x1 = time, 
           y1 = means[[var]] - cis[[var]], col = myPal[var])
  segments(x0 = time - 0.02, y0 = means[[var]] + cis[[var]], x1 = time + 0.02, 
           y1 = means[[var]] + cis[[var]], col = myPal[var])
  segments(x0 = time - 0.02, y0 = means[[var]] - cis[[var]], x1 = time + 0.02, 
           y1 = means[[var]] - cis[[var]], col = myPal[var])
  for (j in 1:(length(time) - 1)) {
    segments(x0 = time[j], y0 = means[[var]][j], x1 = time[j + 1], 
             y1 = means[[var]][j + 1], lty = 3, col = myPal[var])
  }
}

### Panel D (CVPC trends)
# Correct CVs for the number of visits included
dat.cvs[ , paste(selected.biomarkersHD, "CV", sep = "_")] <- 
  apply(dat.cvs[ , paste(selected.biomarkersHD, "CV", sep = "_")], 2, function(x) {
  residuals(nls(x ~ I(1 / sqrt(nb_vis) * a) + b, data = dat.cvs, start = list(a = 1, b = 1)))
})

# Add "status" variable (for cox regressions)
dat.cox.cv <- ddply(dat.cvs, .(id_no), status_fct)

# Calculate PCA on CVs
pca_cv <- prcomp(dat.cox.cv[ , paste(selected.biomarkersHD, "CV", sep = "_")], center = T, scale. = T)
dat.pc <- cbind(dat.cox.cv, pca_cv$x)
for (pc in 1:length(selected.biomarkersHD)) {
  dat.pc <- pc_sign(dat.pc, paste("PC", pc, sep = ""), pca_cv)
}

# Cut time at 5 years before death
dat.5yr <- dat.pc[dat.pc$years <= 5 & !is.na(dat.pc$date_death), ]
dat.5yr$time <- dat.5yr$years

# Calculate mean and 95%CI per time point for each CVPC
time <- unique(dat.5yr$time)[order(unique(dat.5yr$time))]
means <- list()
cis <- list()
for (var in 1:length(selected.biomarkersHD)) {
  means_raw <- tapply(dat.5yr[ , paste("PC", var, sep = "")], dat.5yr$time, mean)
  cis[[var]] <- tapply(dat.5yr[ , paste("PC", var, sep = "")], dat.5yr$time, 
                       function(x) qnorm(0.975) * sd(x) / sqrt(length(x)))
  
  # Center means so that values at 5 years start at 0
  means[[var]] <- means_raw - means_raw[which(names(means_raw)==5)]
}

# Find y-axis limits
ylims <- range(rbind(mapply(function(x, y) x - y, means, cis), 
                     mapply(function(x, y) x + y, means, cis)))
plot(1:10, 1:10, type = "n", xlim = rev(range(dat.5yr$time)), 
     ylim = ylims, ylab = "CVPCs (1-16)", 
     xlab = "Time before death (years)", frame.plot = F)
mtext("C", side = 3, line = -0.5, at = 6.2, font = 2, cex = 1.3)
for (var in 1:length(means)) {
  points(x = time, y = means[[var]], pch = 16, col = myPal[var], cex = 0.9)
  segments(x0 = time, y0 = means[[var]] + cis[[var]], x1 = time, 
           y1 = means[[var]] - cis[[var]], col = myPal[var])
  segments(x0 = time - 0.02, y0 = means[[var]] + cis[[var]], x1 = time + 0.02, 
           y1 = means[[var]] + cis[[var]], col = myPal[var])
  segments(x0 = time - 0.02, y0 = means[[var]] - cis[[var]], x1 = time + 0.02, 
           y1 = means[[var]] - cis[[var]], col = myPal[var])
  for (j in 1:(length(time) - 1)) {
    segments(x0 = time[j], y0 = means[[var]][j], x1 = time[j + 1], 
             y1 = means[[var]][j + 1], lty = 3, col = myPal[var])
  }
}

## Legend
par(mar = c(0, 1, 0, 0), xpd = T)
plot(1:10, 1:10, type = "n", xlab = "", xaxt = "n", yaxt = "n", 
     ylab = "", frame.plot = F)
legend("topleft", 
       legend = paste(paste("CVPC", 1:length(selected.biomarkersHD), sep = ""),
                      paste("PC", 1:length(selected.biomarkersHD), sep = ""),
                      sep = " / "), 
       col = myPal, pch = 16, cex = 0.85, bty = "n", xpd = T)

dev.off()
par(opar)


