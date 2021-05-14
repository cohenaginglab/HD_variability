library(survival)
library(splines)
library(RColorBrewer)

rm(list = ls())


source("Functions_HD_variability_paper.R")

#### Load data ####
clean.HD$date_visit <- as.Date(clean.HD$date_visit)
clean.HD$date_death <- as.Date(clean.HD$date_death)

## variable list
selected.biomarkersHD <- c("wbc", "hb", "hct", "mch", "mchc", "mcv", "pltlt", "rdw", 
                           "k", "sodium", "rbc")
nms.2weeks <- c("WBC", "Hemoglobin", "Hematocrit", "MCH", "MCHC", "MCV", "Platelets", 
                "RDW", "Potassium", "Sodium", "RBC")

# Transform variables not normally distributed
if ("pltlt" %in% selected.biomarkersHD) {
  clean.HD$pltlt <- sqrt(clean.HD$pltlt)
}
clean.HD[ , which(colnames(clean.HD) %in% c("wbc", "rdw", "Glucose"))] <- 
  log(clean.HD[ , which(colnames(clean.HD) %in% c("wbc", "rdw", "Glucose"))])

# Add date of last contact
dat <- ddply(clean.HD, .(id_no), lastvis)


## Cox models, sequentially adding further PCs/CVPCs
# Run PCA on all values to extract the loadings (which are going to be applied on means per year)
pca.all <- prcomp(dat[ , selected.biomarkersHD], center = T, scale. = T)

# Calculate mean for each biomarker per year 
dat0 <- dat
dat0[ , selected.biomarkersHD] <- scale(dat0[ , selected.biomarkersHD])
dat.means <- ddply(dat0, .(id_no), cv.year, vars = selected.biomarkersHD, 
                   cv = F, min.vis.nb = 2, covars = c("id_no", "sex", "date_death", 
                                                      "date_birth", "date_last", "diabetes", 
                                                      "fu_length"))
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


# Add "mean" to column names
colnames(dat.cox)[which(colnames(dat.cox) %in% selected.biomarkersHD)] <- 
  unlist(lapply(colnames(dat.cox)[which(colnames(dat.cox) %in% selected.biomarkersHD)], 
                function (x) paste(x, "mean", sep = "_")))


####### Calculate CVs per year and PCs on these CVs (CVPCs) #######
# Calculate CVs per year 
dat.cvs <- ddply(dat, .(id_no), cv.year, vars = selected.biomarkersHD, 
                 cv = T, min.vis.nb = 2, covars = c("id_no", "sex", "date_death", 
                                                    "date_birth", "date_last", 
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
# Re-order columns (so that CVPCs are from 1 to 11)
dat.pc <- dat.pc[ , c(setdiff(colnames(dat.pc), grep("CVPC", colnames(dat.pc), value = T)),
                      paste("CVPC", 1:length(selected.biomarkersHD), sep = ""))]


dat.both <- merge(dat.cox, 
                  dat.pc[ , -which(colnames(dat.pc) %in% c("sex","date_death", "date_birth",
                                                           "date_last", "diabetes", 
                                                           "fu_length", "age_visit", 
                                                           "nb_vis", "vis_sd", "status"))],
                  by = c("id_no", "years"))
# Cox models don't take zero values as time.
dat.both$years <- dat.both$years + 0.5

############# Cox models ###############
### Part 1 - Add other PCs one by one 
pcs <- c("PC1", "PC2", "PC6", "CVPC1", "CVPC3")
res_i <- list()
for (i in 1:length(pcs)) {         
  if (i < 4) vars <- paste("PC", 1:length(selected.biomarkersHD), sep = "") else 
    vars <- grep("CVPC", colnames(dat.both), value = T)
  
  res_k <- list()
  for (k in 1:2) {     # With and without covariates
    if (k==1) {
      covars <- ""
    } else {
      covars <- c("bs(age_visit, df = 5)", "sex", "diabetes", "fu_length")
    }
    
    res_j <- list()
    for (j in 1:length(selected.biomarkersHD)) {
      res <- as.data.frame(matrix(NA, 1, 12))
      colnames(res) <- c("Mean", "CV", "PCs", "controls", "Variable", "HR", "LCI", "UCI", 
                         "p_value", "HR95", "HR95_LCI", "HR95_UCI")
      res$Mean <- ifelse(nchar(pcs[i])==3, 1, 0)
      res$CV <- ifelse(nchar(pcs[i])==5, 1, 0)
      res$PCs <- c(1:length(selected.biomarkersHD))[j]
      res$controls <- ifelse(k==1, 0, 1)
      
      pc <- pcs[i]
      other.pcs <- setdiff(vars, pc)
      
      if (j==1) {
        fmla <- formula(paste("Surv(years, status) ~ ", 
                              ifelse(k==1, pc, paste(c(pc, covars), collapse = " + ")),
                              sep = ""))
      } else {
        fmla <- formula(paste("Surv(years, status) ~ ",
                              ifelse(k==1, paste(c(pc, other.pcs[1:(j - 1)]), 
                                                 collapse = " + "),
                                     paste(c(pc, other.pcs[1:(j - 1)], covars), 
                                           collapse = " + ")),
                              sep = ""))
      }
      
      # Run cox model
      cox.mod <- coxph(fmla, data = dat.both, weights = sqrt(nb_vis), cluster = id_no)
      res$Variable <- pc
      res[ , c("HR", "LCI", "UCI")] <- summary(cox.mod)$conf.int[1:nrow(res), c(1, 3, 4)]
      res[ , "p_value"] <- summary(cox.mod)$coefficients[1:nrow(res), 6]
      for (vv in 1:nrow(res)) {
        res[vv, c("HR95", "HR95_LCI", "HR95_UCI")] <- summary(cox.mod)$conf.int[vv, c(1, 3, 4)] ^ 
          (quantile(dat.both[ , rownames(summary(cox.mod)$coefficients)[vv]], 0.975) - 
             quantile(dat.both[ , rownames(summary(cox.mod)$coefficients)[vv]], 0.025))
        
      }
      res_j[[j]] <- res
    }
    res_k[[k]] <- bind_rows(res_j)
  }
  res_i[[i]] <- bind_rows(res_k)
}
res_part1 <- bind_rows(res_i)


### Part 2 - Rerun means with CVs, controlling for all CVPCs for each combination of PCs and vice-versa
vars <- list(paste("PC", 1:length(selected.biomarkersHD), sep = ""), 
             grep("CVPC", colnames(dat.both), value = T))
pcs <- c("PC1", "PC2", "PC6", "CVPC1", "CVPC3")
res_i <- list()
for (i in 1:length(pcs)) {         
  pc <- pcs[i]
  other.pcs <- setdiff(vars[[ifelse(i < 4, 1, 2)]], pc)
  
  res_k <- list()
  for (k in 1:2) {        # with and without covariates
    if (k==1) {
      covars <- ""
    } else {
      covars <- c("bs(age_visit, df = 5)", "sex", "diabetes", "fu_length")
    }
    
    res_j <- list()
    for (j in 1:length(selected.biomarkersHD)) {
      res <- as.data.frame(matrix(NA, 1, 12))
      colnames(res) <- c("Mean", "CV", "PCs", "controls", "Variable", "HR", "LCI", "UCI", 
                         "p_value", "HR95", "HR95_LCI", "HR95_UCI")
      res$Mean <- 1
      res$CV <- 1
      res$PCs <- c(1:length(selected.biomarkersHD))[j]
      res$controls <- ifelse(k==1, 0, 1)
      
      if (j==1) {
        fmla <- formula(paste("Surv(years, status) ~ ",
                              ifelse(k==1, paste(c(pc, vars[[ifelse(i < 4, 2, 1)]]), 
                                                 collapse = " + "),
                                     paste(c(pc, vars[[ifelse(i < 4, 2, 1)]], covars), 
                                           collapse = " + ")),
                              sep = ""))
      } else {
        fmla <- formula(paste("Surv(years, status) ~ ",
                              ifelse(k==1, paste(c(pc, other.pcs[1:(j - 1)], 
                                                   vars[[ifelse(i < 4, 2, 1)]]), 
                                                 collapse = " + "),
                                     paste(c(pc, other.pcs[1:(j - 1)], 
                                             vars[[ifelse(i < 4, 2, 1)]], covars), 
                                           collapse = " + ")),
                              sep = ""))
      }
      
      # Run cox model
      cox.mod <- coxph(fmla, data = dat.both, weights = sqrt(nb_vis), cluster = id_no)
      res$Variable <- rownames(summary(cox.mod)$coefficients)[1:nrow(res)]
      res[ , c("HR", "LCI", "UCI")] <- summary(cox.mod)$conf.int[1:nrow(res), c(1, 3, 4)]
      res[ , "p_value"] <- summary(cox.mod)$coefficients[1:nrow(res), 6]
      for (vv in 1:nrow(res)) {
        res[vv, c("HR95", "HR95_LCI", "HR95_UCI")] <- summary(cox.mod)$conf.int[vv, c(1, 3, 4)] ^ 
          (quantile(dat.both[ , rownames(summary(cox.mod)$coefficients)[vv]], 0.975) - 
             quantile(dat.both[ , rownames(summary(cox.mod)$coefficients)[vv]], 0.025))
        
      }
      res_j[[j]] <- res
    }
    res_k[[k]] <- bind_rows(res_j)
  }
  res_i[[i]] <- bind_rows(res_k)
}
res_part2 <- bind_rows(res_i)
res <- rbind(res_part1, res_part2)
write.csv(res, "ExtDatFigure3.csv", row.names = F)

#### Plot
tiff("Extended Data_Figure_3.tiff", width = 180, height = 270, "mm", res = 300, compression = "lzw")
opar <- par(mfrow = c(3, 2), mar = c(2.5, 2.5, 1.5, 0.5), oma = c(0.2, 0.2, 0.2, 0.2), 
            mgp = c(1.5, 0.4, 0), tcl = -0.3, cex = 0.9)

cols <- list(c("#4053d3", "#26317e"), c("#cd1719", "#880f10"))

ylims <- c(10, 8, 10, 14, 8)

pcs <- c("PC1", "PC2", "PC6", "CVPC1", "CVPC3")
for (pc in 1:length(pcs)) {
  
  if (pc < 4)  col.pc <- cols[[1]] else col.pc <- cols[[2]]
  
  # Select relevant results
  dat <- res[res$Variable==pcs[pc], ]
  
  plot(1:10, 1:10, type = "n", xlim = c(1, length(selected.biomarkersHD)), 
       xlab = "# of PCs included", ylim = c(min(dat$HR95_LCI), ylims[pc]),
       ylab = "HR95", frame.plot = F, cex.axis = 0.9)
  ## Add grid 
  abline(h = seq(2, 14, 2), col = "grey85", lwd = 2)
  abline(h = seq(3, 13, 2), col = "grey85")
  
  poss <- list(c(0, 0), c(0, 1), c(1, 0), c(1, 1))
  manuel_jit <- c(-0.15, 0.05, -0.05, 0.15)
  for (i in 1:4) {
    if (pc < 4) {
      dat.poss <- dat[dat$CV==poss[[i]][1] & dat$controls==poss[[i]][2], ]
    } else {
      dat.poss <- dat[dat$Mean==poss[[i]][1] & dat$controls==poss[[i]][2], ]
    }
    points(x = dat.poss$PCs + manuel_jit[i], y = dat.poss$HR95, pch = 16, 
           col = col.pc[ifelse(i < 3, 1, 2)])
    segments(x0 = dat.poss$PCs + manuel_jit[i], y0 = dat.poss$HR95_LCI, 
             x1 = dat.poss$PCs + manuel_jit[i], y1 = dat.poss$HR95_UCI, 
             col = col.pc[ifelse(i < 3, 1, 2)], lty = ifelse(i %in% c(1, 3), 1, 2))
    segments(x0 = dat.poss$PCs - 0.1 + manuel_jit[i], y0 = dat.poss$HR95_LCI, 
             x1 = dat.poss$PCs + 0.1 + manuel_jit[i], y1 = dat.poss$HR95_LCI, 
             col = col.pc[ifelse(i < 3, 1, 2)])
    segments(x0 = dat.poss$PCs - 0.1 + manuel_jit[i], y0 = dat.poss$HR95_UCI, 
             x1 = dat.poss$PCs + 0.1 + manuel_jit[i], y1 = dat.poss$HR95_UCI, 
             col = col.pc[ifelse(i < 3, 1, 2)])
    for (j in 1:(nrow(dat.poss) - 1)) {
      segments(x0 = dat.poss[j, "PCs"] + manuel_jit[i], y0 = dat.poss[j, "HR95"], 
               x1 = dat.poss[j + 1, "PCs"] + manuel_jit[i], y1 = dat.poss[j + 1, "HR95"], 
               col = col.pc[ifelse(i < 3, 1, 2)], lty = ifelse(i %in% c(1, 3), 1, 2))
    }
  }
  mtext(paste(c("A", "B", "C", "D", "E")[pc], pcs[pc], sep = " - "), side = 3, 
        line = 0.5, at = -0.5, adj = 0, font = 2, cex = 1.1)
  
}

# Legend
par(mar = c(2.5, 1, 2, 0.5), cex = 0.9)
plot(1:10, 1:10, type = "n", xlab = "", xaxt = "n", yaxt = "n", 
     ylab = "", frame.plot = F)
text(x = c(3, 3, 5, 5, 3), y = c(10, 9.4, 9.75, 9.05, 8.8), 
     labels = c("Primary", "variable", "Control", "variables", "(x-axis)"), 
     font = 2, cex = 0.95)
for (i in 1:length(unlist(cols))) {
  segments(x0 = 0, x1 = 1.5, y0 = c(8.2, 7.5, 6.8, 6.1)[i], 
           y1 = c(8.2, 7.5, 6.8, 6.1)[i], col = unlist(cols)[i], lwd = 2)
}
text(x = rep(3, 4), y = c(8.2, 7.5, 6.8, 6.1), 
     labels = c("PCs", "PCs", "CVPCs", "CVPCs"), cex = 0.9)
text(x = rep(5, 4), y = c(8.2, 7.5, 6.8, 6.1), 
     labels = c("none", "CVPCs (all)", "none", "PCs (all)"), cex = 0.9)

segments(x0 = 0, x1 = 1.5, y0 = 4.9, y1 = 4.9, lwd = 2)
segments(x0 = 0, x1 = 1.5, y0 = 4.2, y1 = 4.2, lty = 2, lwd = 2)
text(x = c(2.2, 2.2), y = c(4.9, 4.2), adj = 0, 
     labels = c("Without demographic variables", "With demographic variables"), 
     cex = 0.9)

dev.off()
par(opar)
