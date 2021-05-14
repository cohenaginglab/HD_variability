#######################################################################################################
#####                  This code is for the figure entitled "Trends before death of               ##### 
#####                              integrative multivariate indices"                              #####
#######################################################################################################
library(mcp)
library(RColorBrewer)

rm(list = ls())

source("Functions_HD_variability_paper.R")

#### Load data ####
clean.HD$date_visit <- as.Date(clean.HD$date_visit)
clean.HD$date_death <- as.Date(clean.HD$date_death)

## variable list
vars.2wks <- c("wbc", "hb", "hct", "mch", "mchc", "mcv", "pltlt", "rdw", "k", "sodium", "rbc")
nms.2weeks = c("WBC", "Hemoglobin", "Hematocrit", "MCH", "MCHC", "MCV", "Platelets", "RDW", 
               "Potassium", "Sodium", "RBC")
selected.biomarkersHD <- vars.2wks

# Transform variables not normally distributed
if ("pltlt" %in% selected.biomarkersHD) clean.HD$pltlt <- sqrt(clean.HD$pltlt)
clean.HD[ , which(colnames(clean.HD) %in% c("wbc", "rdw", "Glucose"))] <- 
  log(clean.HD[ , which(colnames(clean.HD) %in% c("wbc", "rdw", "Glucose"))])

# Add date of last contact
dat <- ddply(clean.HD, .(id_no), lastvis)


tiff("Mean_CV_trends.tiff", width = 180, height = 180, "mm", res = 300, compression = "lzw")
opar <- par(mfrow = c(2, 2), mar = c(3, 3, 0.5, 0.5), oma = c(0.2, 0.2, 0.2, 0.2), 
            mgp = c(1.5, 0.4, 0), tcl = -0.3, cex = 1.1)

# Set color palette
myPal2 <- brewer.pal(12, "Paired")[-11]

## Panel A - CVPC trends, calculated every 3 months, and change point analyses
# Calculate CVs every 3 months
dat.cvs <- ddply(dat, .(id_no), cv.4mth, selected.biomarkersHD = selected.biomarkersHD, cv = T, 
                 nb.month = 3, min.vis.nb = 2, 
                 covars = c("id_no", "sex", "date_death", "date_birth", "date_last", "diabetes", "fu_length"))

# Apply log transformation on CVs
dat.cvs[ , selected.biomarkersHD] <- apply(dat.cvs[ , selected.biomarkersHD], 2, log_zeros)

# Some variables have no variation, leading to outliers
dat.cvs[ , selected.biomarkersHD] <- apply(dat.cvs[ , selected.biomarkersHD], 2, uniform_val)

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
rm(dat.cox.cv, pca_cv)

# Add "CV" to column names
colnames(dat.pc)[which(colnames(dat.pc) %in% selected.biomarkersHD)] <- 
  unlist(lapply(colnames(dat.pc)[which(colnames(dat.pc) %in% selected.biomarkersHD)], 
                function (x) paste(x, "CV", sep = "_")))
colnames(dat.pc)[grep("PC", colnames(dat.pc))] <- 
  unlist(lapply(colnames(dat.pc)[grep("PC", colnames(dat.pc))], 
                function (x) paste("CV", x, sep = "")))


### Change point analysis 
# Cut time at 5 years before death
dat.5yr <- dat.pc[dat.pc$years <= (5 * (12/3)) & !is.na(dat.pc$date_death), ]
dat.5yr$time <- dat.5yr$years / (12/3)

rm(clean.HD, clean.HD2, dat, dat.cvs)

cp.mods <- list()
for (var in 1:length(selected.biomarkersHD)) {
  model = list(
    formula(paste(paste("CVPC", var, sep = ""), "~ 1 + time", sep = " ")),   
    1 + (1|id_no) ~ 0 + time   
  )
  dat.ordered <- dat.5yr[order(dat.5yr$time, decreasing = T), ]
  cp.mods[[var - 1]] <- mcp(model, data = dat.ordered)
  
}
save(cp.mods, file = "Change point_all_CVPCs.RData")

# Calculate mean and 95%CI per time point for each CVPC
time <- unique(dat.5yr$time)[order(unique(dat.5yr$time))]
means <- list()
cis <- list()
for (var in 1:length(selected.biomarkersHD)) {
  means_raw <- tapply(dat.5yr[ , paste("CVPC", var, sep = "")], dat.5yr$time, mean)
  cis[[var]] <- tapply(dat.5yr[ , paste("CVPC", var, sep = "")], dat.5yr$time, 
                       function(x) qnorm(0.975) * sd(x) / sqrt(length(x)))
  
  # Center means so that values at 5 years start at 0
  means[[var]] <- means_raw - means_raw[which(names(means_raw)==5)]
}

# Find y-axis limits
ylims <- range(rbind(mapply(function(x, y) x - y, means, cis), 
                     mapply(function(x, y) x + y, means, cis)))
plot(1:10, 1:10, type = "n", xlim = rev(range(dat.5yr$time)), 
     ylim = ylims, ylab = "CVPCs (1-11)", 
     xlab = "Time before death (years)", frame.plot = F)
mtext("A", side = 3, line = -0.5, at = 6.2, font = 2, cex = 1.3)
for (var in 1:length(means)) {
  points(x = time, y = means[[var]], pch = 16, col = myPal2[var], cex = 0.9)
  segments(x0 = time, y0 = means[[var]] + cis[[var]], x1 = time, 
           y1 = means[[var]] - cis[[var]], col = myPal2[var])
  segments(x0 = time - 0.02, y0 = means[[var]] + cis[[var]], x1 = time + 0.02, 
           y1 = means[[var]] + cis[[var]], col = myPal2[var])
  segments(x0 = time - 0.02, y0 = means[[var]] - cis[[var]], x1 = time + 0.02, 
           y1 = means[[var]] - cis[[var]], col = myPal2[var])
  for (j in 1:(length(time) - 1)) {
    segments(x0 = time[j], y0 = means[[var]][j], x1 = time[j + 1], 
             y1 = means[[var]][j + 1], lty = 3, col = myPal2[var])
  }
  
  # Add change point
  fit <- cp.mods[[var]]
  abline(v = summary(fit)$mean[1], col = myPal2[var], lty = 2)
  
}


### Panel B (PC trends)
# Run PCA on all values to extract the loadings (which are going to be applied on means per year)
pca.all <- prcomp(dat[ , selected.biomarkersHD], center = T, scale. = T)

# Calculate means every 3 months
dat.means <- ddply(dat, .(id_no), cv.4mth, selected.biomarkersHD = selected.biomarkersHD, 
                   cv = F, nb.month = 3, min.vis.nb = 2, 
                   covars = c("id_no", "sex", "date_death", "date_birth", "date_last", "diabetes", "fu_length"))

# Add "status" variable (for cox regressions)
dat.cox <- ddply(dat.means, .(id_no), status_fct)

# Calculate PC scores using loadings from all values
for (pc in 1:length(selected.biomarkersHD)) {
  dat.cox[ , paste("PC", pc, sep = "")] <- 
    as.matrix(dat.cox[ , selected.biomarkersHD]) %*% pca.all$rotation[ , pc]
}
for (pc in 1:length(selected.biomarkersHD)) {
  dat.cox <- pc_sign(dat.cox, paste("PC", pc, sep = ""), pca.all)
}
# Re-order columns (so that PCs are from 1 to 11)
dat.cox <- dat.cox[ , c(setdiff(colnames(dat.cox), 
                                paste("PC", 1:length(selected.biomarkersHD), sep = "")),
                        paste("PC", 1:length(selected.biomarkersHD), sep = ""))]


# Cut time at 5 years before death
dat.5yr <- dat.cox[dat.cox$years <= (5 * (12 / 3)) & 
                     !is.na(dat.cox$date_death), ]
dat.5yr$time <- dat.5yr$years / (12 / 3)

# Calculate mean and 95%CI per time point for each PC
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
     ylim = ylims, ylab = "PCs (1-11)", 
     xlab = "Time before death (years)", frame.plot = F)
mtext("B", side = 3, line = -0.5, at = 6.2, font = 2, cex = 1.3)
for (var in 1:length(means)) {
  points(x = time, y = means[[var]], pch = 16, col = myPal2[var], cex = 0.9)
  segments(x0 = time, y0 = means[[var]] + cis[[var]], x1 = time, 
           y1 = means[[var]] - cis[[var]], col = myPal2[var])
  segments(x0 = time - 0.02, y0 = means[[var]] + cis[[var]], x1 = time + 0.02, 
           y1 = means[[var]] + cis[[var]], col = myPal2[var])
  segments(x0 = time - 0.02, y0 = means[[var]] - cis[[var]], x1 = time + 0.02, 
           y1 = means[[var]] - cis[[var]], col = myPal2[var])
  for (j in 1:(length(time) - 1)) {
    segments(x0 = time[j], y0 = means[[var]][j], x1 = time[j + 1], 
             y1 = means[[var]][j + 1], lty = 3, col = myPal2[var])
  }
}

## Extract change point results into a table 
xx <- as.data.frame(matrix(NA, length(selected.biomarkersHD), 3))
colnames(xx) <- c("Change_point", "LCI", "UCI")
rownames(xx) <- paste("CVPC", 1:length(selected.biomarkersHD), sep = "")
for (pc in 1:length(selected.biomarkersHD)) {
  fit <- cp.mods[[pc]]
  xx$Change_point[pc] <- summary(fit)$mean[1]
  xx$LCI[pc] <- summary(fit)$lower[1]
  xx$UCI[pc] <- summary(fit)$upper[1]
}
plot(1:10, 1:10, type = "n", xlab = "", xaxt = "n", yaxt = "n", ylab = "", frame.plot = F)
text(x = 7.5, y = 9.7, labels = "Change point", font = 2, xpd = T)
text(x = rep(4, nrow(xx)), y = seq(8.5, -0.5, length.out = nrow(xx)), 
     labels = rownames(xx), adj = 1, cex = 0.9, xpd = T)
text(x = rep(7.5, nrow(xx)), y = seq(8.5, -0.5, length.out = nrow(xx)), 
     labels = paste(format(round(xx$Change_point, 2), nsmall = 2), " (",
                    format(round(xx$LCI, 2), nsmall = 2), " - ",  
                    format(round(xx$UCI, 2), nsmall = 2), ")", sep = ""), cex = 0.9, adj = 0.5, xpd = T)
segments(x0 = 1, y0 = -1.2, x1 = 10, y1 = -1.2, lwd = 2, xpd = T)
segments(x0 = 1, y0 = 10.3, x1 = 10, y1 = 10.3, lwd = 2, xpd = T)
segments(x0 = 1, y0 = 9.1, x1 = 10, y1 = 9.1)

## Legend
plot(1:10, 1:10, type = "n", xlab = "", xaxt = "n", yaxt = "n", ylab = "", frame.plot = F)
legend(x = 2, y = 9,
       legend = paste(paste("CVPC", 1:length(selected.biomarkersHD), sep = ""),
                      paste("PC", 1:length(selected.biomarkersHD), sep = ""), sep = " / "), 
       col = myPal2, pch = 16, cex = 0.9, bty = "n", xpd = T)



dev.off()
par(opar)

