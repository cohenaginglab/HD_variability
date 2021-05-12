library(survival)
library(splines)
library(mcp)
library(RColorBrewer)

rm(list = ls())

source("Functions_HD_variability_paper.R")

#### Load data ####
clean.HD$date_visit <- as.Date(clean.HD$date_visit)
clean.HD$date_death <- as.Date(clean.HD$date_death)

## variable list
vars.2wks <- c("wbc", "hb", "hct", "mch", "mchc", "mcv", "pltlt", "rdw", "k", 
               "sodium", "rbc")
nms.2weeks = c("WBC", "Hemoglobin", "Hematocrit", "MCH", "MCHC", "MCV", "Platelets", "RDW", 
               "Potassium", "Sodium", "RBC")
selected.biomarkersHD <- vars.2wks

# Transform variables not normally distributed
if ("pltlt" %in% selected.biomarkersHD) clean.HD$pltlt <- sqrt(clean.HD$pltlt)
clean.HD[ , which(colnames(clean.HD) %in% c("wbc", "rdw", "Glucose"))] <- 
  log(clean.HD[ , which(colnames(clean.HD) %in% c("wbc", "rdw", "Glucose"))])

# Add date of last contact
dat <- ddply(clean.HD, .(id_no), lastvis)

tiff("Figure_4.tiff", width = 180, height = 180, "mm", res = 300, compression = "lzw")
opar <- par(mfrow = c(2, 2), mar = c(3, 3, 0.5, 0.5), oma = c(0.2, 0.2, 0.2, 0.2), 
            mgp = c(1.5, 0.4, 0), tcl = -0.3, cex = 1.1)

# Set color palettes
myPal1 <- brewer.pal(3, "Set1")
myPal <- brewer.pal(length(selected.biomarkersHD), "Paired")

## Panel A - CVPC1 trend, calculated every 2, 3, and 4 months + change point analyses
plot(1:10, 1:10, type = "n", xlim = c(5, 0), 
     ylim = c(-0.9, 1.7), ylab = "CVPC1", 
     xlab = "Time before death (years)", frame.plot = F)
mtext("A", side = 3, line = -0.5, at = 6.2, font = 2, cex = 1.3)

cp.mods <- list()
time.windows <- c(2, 3, 4)
for (tp in 1:length(time.windows)) {
  
  # Calculate CVs every 2, 3, and 4 months
  dat.cvs <- ddply(dat, .(id_no), cv.4mth, selected.biomarkersHD = selected.biomarkersHD, 
                   cv = T, nb.month = time.windows[tp], min.vis.nb = 2, 
                   covars = c("id_no", "sex", "date_death", "date_birth", "date_last", 
                              "diabetes", "fu_length"))
  
  # Apply log transformation on CVs
  dat.cvs[ , selected.biomarkersHD] <- apply(dat.cvs[ , selected.biomarkersHD], 2, 
                                             log_zeros)
  
  # Some variables have no variation, leading to outliers
  dat.cvs[ , selected.biomarkersHD] <- apply(dat.cvs[ , selected.biomarkersHD], 2, 
                                             uniform_val)
  
  # Correct CVs for the number of visits included
  dat.cvs[ , selected.biomarkersHD] <- apply(dat.cvs[ , selected.biomarkersHD], 2, function(x) {
    residuals(nls(x ~ I(1 / sqrt(nb_vis) * a) + b, data = dat.cvs, start = list(a = 1, b = 1)))
  })
  
  # Add "status" variable (for cox regressions)
  dat.cox.cv <- ddply(dat.cvs, .(id_no), status_fct)
  
  # Calculate PCA on CVs
  pca_cv <- prcomp(dat.cox.cv[ , selected.biomarkersHD], center = T, scale. = T)
  dat.pc <- cbind(dat.cox.cv, pca_cv$x)
  for (pc in 1:length(selected.biomarkersHD)) {
    dat.pc <- pc_sign(dat.pc, paste("PC", pc, sep = ""), pca_cv)
  }
  
  # Add "CV" to column names
  colnames(dat.pc)[which(colnames(dat.pc) %in% selected.biomarkersHD)] <- 
    unlist(lapply(colnames(dat.pc)[which(colnames(dat.pc) %in% selected.biomarkersHD)], 
                  function (x) paste(x, "CV", sep = "_")))
  colnames(dat.pc)[grep("PC", colnames(dat.pc))] <- 
    unlist(lapply(colnames(dat.pc)[grep("PC", colnames(dat.pc))], 
                  function (x) paste("CV", x, sep = "")))
  
  # Cut time at 5 years before death
  dat.5yr <- dat.pc[dat.pc$years <= (5 * (12 / time.windows[tp])) & 
                      !is.na(dat.pc$date_death), ]
  dat.5yr$time <- dat.5yr$years / (12 / time.windows[tp])
  
  model = list(
    CVPC1 ~ 1 + time,          # intercept + slope
    1 + (1|id_no) ~ 0 + time   # joined slope, varying by id
  )
  dat.ordered <- dat.5yr[order(dat.5yr$time, decreasing = T), ]
  fit <- mcp(model, data = dat.ordered)
  cp.mods[[tp]] <- fit
  
  time <- unique(dat.5yr$time)[order(unique(dat.5yr$time))]
  means <- tapply(dat.5yr$CVPC1, dat.5yr$time, mean)
  cis <- tapply(dat.5yr$CVPC1, dat.5yr$time, 
                function(x) qnorm(0.975) * sd(x) / sqrt(length(x)))
  points(x = time, y = means, pch = 16, col = myPal1[tp])
  segments(x0 = time, y0 = means + cis, x1 = time, y1 = means - cis, 
           col = myPal1[tp])
  segments(x0 = time - 0.02, y0 = means + cis, x1 = time + 0.02, y1 = means + cis, 
           col = myPal1[tp])
  segments(x0 = time - 0.02, y0 = means - cis, x1 = time + 0.02, y1 = means - cis, 
           col = myPal1[tp])
  for (i in 1:(length(time) - 1)) {
    segments(x0 = time[i], y0 = means[i], x1 = time[i + 1], y1 = means[i + 1], 
             lty = 3, col = myPal1[tp])
  }
  abline(v = summary(fit)$mean[1], col = myPal1[tp], lty = 2)
  text(x = 2.7, y = c(1.7, 1.45, 1.2)[tp], 
       labels = paste("Change point =", 
                      round(as.numeric(summary(fit)$mean[1]) * 12, 1), "months"), 
       col = myPal1[tp], font = 2, cex = 0.9)
}
save(cp.mods, file = "Change point_CVPC1.RData")

## Legends 
par(mar = c(0.5, 0.5, 0.5, 0.5))
plot(1:10, 1:10, type = "n", xlab = "", xaxt = "n", yaxt = "n", 
     ylab = "", frame.plot = F)
# Legend for panel A
legend(x = 1, y = 10, legend = c("bimonthly", "trimonthly", "four-monthly"),   
       col = myPal1, pch = 16, cex = 0.9, bty = "n")
text(x = 2.3, y = 10.2, labels = "CVs calculated", font = 2, xpd = T, cex = 0.95)

# Legend for panels B and C
legend("bottomright", legend = nms.2weeks, col = myPal, pch = 16, cex = 0.95, 
       bty = "n", xpd = T)

##### Panels B and C (individual biomarker mean and CV trends) ####
par(mar = c(3, 3, 0.5, 0.5))
# Calculate mean for each biomarker per 3 months
dat0 <- dat
dat0[ , selected.biomarkersHD] <- scale(dat0[ , selected.biomarkersHD])
dat.means <- ddply(dat0, .(id_no), cv.4mth, 
                   selected.biomarkersHD = selected.biomarkersHD, 
                   cv = F, nb.month = 3, min.vis.nb = 2, 
                   covars = c("id_no", "sex", "date_death", "date_birth", "date_last", 
                              "diabetes", "fu_length"))

# Add "mean" to column names
colnames(dat.means)[which(colnames(dat.means) %in% selected.biomarkersHD)] <- 
  unlist(lapply(colnames(dat.means)[which(colnames(dat.means) %in% selected.biomarkersHD)], 
                function (x) paste(x, "mean", sep = "_")))


# Calculate CVs per 3 months
dat.cvs <- ddply(dat, .(id_no), cv.4mth, selected.biomarkersHD = selected.biomarkersHD, 
                 cv = T, nb.month = 3, min.vis.nb = 2, 
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
dat.5yr <- dat.both[dat.both$years <= (5 * (12/3)) & !is.na(dat.both$date_death), ]
dat.5yr$time <- dat.5yr$years / (12/3)

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
  
  # Set y-axis limits
  ylims <- range(rbind(mapply(function(x, y) x - y, means, cis), 
                       mapply(function(x, y) x + y, means, cis)))
  if (i==1) ylims <- c(ylims[1] - 0.1, ylims[2] + 0.1)
  
  if (i==1) {
    plot(1:10, 1:10, type = "n", xlim = rev(range(dat.5yr$time)), 
         ylim = ylims, ylab = "Levels (mean z-scores)", yaxp = c(-0.6, 0.9, 6),
         xlab = "Time before death (years)", frame.plot = F)
  } else {
    plot(1:10, 1:10, type = "n", xlim = rev(range(dat.5yr$time)), 
         ylim = ylims, ylab = "Variability (CVs)", 
         xlab = "Time before death (years)", frame.plot = F)
  }
  mtext(c("B", "C")[i], side = 3, line = -0.5, at = 6.2, font = 2, cex = 1.3)
  
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


dev.off()
par(opar)

