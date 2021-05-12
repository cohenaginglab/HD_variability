library(survival)
library(splines)
library(RColorBrewer)
library(pROC)

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


tiff("Figure_2.tiff", width = 180, height = 180, "mm", res = 300, compression = "lzw")
opar <- par(mfrow = c(2, 2), mar = c(3, 10, 3.5, 0.2), oma = c(0.2, 0.2, 0.2, 0.2), 
            mgp = c(1.5, 0.4, 0), tcl = -0.3)

### Panels A and B (HRs mean vs CVs and univariate vs multivariate) ###
cnms <- c("Model", "Variable", "Measure", "HR", "LCI", "UCI", "p", 
          "HR95", "HR95_LCI", "HR95_UCI", "zph_p")
res <- as.data.frame(matrix(NA, 1, length(cnms)))
colnames(res) <- cnms

res_ms <- list()
# Means vs CVs
for (ms in 1:2) {
  
  if (ms==1) {
    # Run PCA on all values to extract the loadings (which are going to be applied on means per year)
    pca.all <- prcomp(dat[ , selected.biomarkersHD], center = T, scale. = T)
    
    # Calculate mean for each biomarker per year 
    dat0 <- dat
    dat0[ , selected.biomarkersHD] <- scale(dat0[ , selected.biomarkersHD])
    dat.cvs <- ddply(dat0, .(id_no), cv.year, vars = selected.biomarkersHD, 
                     cv = F, min.vis.nb = 2, 
                     covars = c("id_no", "sex", "date_death", "date_birth", 
                                "date_last", "diabetes", "fu_length"))
    
    # Add "status" variable (for cox regressions)
    dat.cox <- ddply(dat.cvs, .(id_no), status_fct)
    
    # Calculate PC scores using loadings from all values
    for (pc in 1:length(selected.biomarkersHD)) {
      dat.cox[ , paste("PC", pc, sep = "")] <- 
        as.matrix(dat.cox[ , selected.biomarkersHD]) %*% pca.all$rotation[ , pc]
    }
    for (pc in 1:length(selected.biomarkersHD)) {
      dat.cox <- pc_sign(dat.cox, paste("PC", pc, sep = ""), pca.all)
    }
    
    dat.pc <- dat.cox
    # Reverse PC sign if necessary
    for (vv in 1:length(selected.biomarkersHD)) {
      dat.pc <- pc_sign(dat.pc, selected.biomarkersHD[vv])
    }
  }
  
  if (ms==2) {
    # Calculate CVs per year 
    dat.cvs <- ddply(dat, .(id_no), cv.year, vars = selected.biomarkersHD, 
                     cv = T, min.vis.nb = 2, 
                     covars = c("id_no", "sex", "date_death", "date_birth", 
                                "date_last", "diabetes", "fu_length"))
    
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
    dat.cox <- ddply(dat.cvs, .(id_no), status_fct)
    
    # Calculate PCA on CVs
    pca_cv <- prcomp(dat.cox[ , selected.biomarkersHD], center = T, scale. = T)
    dat.pc <- cbind(dat.cox, pca_cv$x)
    for (pc in 1:length(selected.biomarkersHD)) {
      dat.pc <- pc_sign(dat.pc, paste("PC", pc, sep = ""), pca_cv)
    }
  }
  
  res_um <- list()
  ## Univariate (means/CVs) vs multivariate (PCs/CVPCs)
  variables <- list("uni" = selected.biomarkersHD, 
                    "multi" = paste("PC", 1:length(selected.biomarkersHD), sep = ""))
  for (um in 1:2) {
    
    res1 <- list()
    for (vr in 1:length(selected.biomarkersHD)) {
      res2 <- res
      res2$Model <- names(variables)[um]
      res2$Variable <- variables[[um]][vr]
      res2$Measure <- c("mean", "CV")[ms]
      dat_var <- dat.pc[ , c("id_no", "years", "status", "age_visit", "sex", "diabetes", 
                             "fu_length", "nb_vis", variables[[um]][vr])]
      colnames(dat_var)[ncol(dat_var)] <- "predictor"
      mod <- coxph(Surv(years, status) ~ predictor + bs(age_visit, df = 5) + sex +  
                     diabetes + fu_length, data = dat_var, weights = sqrt(nb_vis),
                   cluster = id_no)
      res2[ , c("HR", "LCI", "UCI")] <- summary(mod)$conf.int[1, c(1, 3, 4)]
      res2[ , "p"] <- summary(mod)$coefficients[1, 6]
      res2[ , c("HR95", "HR95_LCI", "HR95_UCI")] <- summary(mod)$conf.int[1, c(1, 3, 4)] ^ 
        (quantile(na.omit(dat_var[ , "predictor"]), 0.975) - quantile(na.omit(dat_var[ , "predictor"]), 0.025))
      
      # Test the proportional hazards assumption of the cox model
      res2$zph_p <- cox.zph(mod)$table[1, 3]
      res1[[vr]] <- res2
    }
    res1 <- bind_rows(res1)
    
    res_um[[um]] <- res1
  }
  res_ms[[ms]] <- bind_rows(res_um)
}    
res_2wks <- bind_rows(res_ms)
write.csv(res_2wks, "Figure2_panelsA_B.csv", row.names = F)

# Set x-axis range 
x_range <- c(0, max(res_2wks$HR95_UCI))

# Set color palette
myPal1 <- brewer.pal(5, "Set1")[1:2]

for (um in 1:2) {
  x_txt <- c(-13.9, -11.5)[um]
  res_um <- res_2wks[res_2wks$Model==c("uni", "multi")[um], ]
  
  plot(1:10, 1:10, type = "n", xlim = x_range, ylim = c(0, 10.5), 
       ylab = "", xlab = "HR95", yaxt = "n", frame.plot = F)
  mtext(text = c("A - Biomarkers", 
                 "B - Integrative multivariate indices")[ifelse(um==1, 1, 2)], 
        side = 3, at = x_txt, line = 2.7, font = 2, adj = 0, cex = 1)
  abline(v = 1, lty = 3)
  
  text_labels <- list("uni" = nms.2weeks, 
                      "multi" = paste("PC", 1:(nrow(res_um) / 2), sep = ""))
  
  mtext(text = c("PH assumption", "P-value"), side = 3, at = -3.8,
        line = c(1.5, 0.7), cex = 0.7, adj = 0.5, font = 3)
  
  for (ms in 1:2) {
    res_ms <- res_um[res_um$Measure==c("mean", "CV")[ms], ]
    
    for (vr in 1:nrow(res_um)) {
      segments(x0 = res_ms[vr, "HR95_LCI"], x1 = res_ms[vr, "HR95_UCI"], 
               y0 = nrow(res_ms) - vr, y1 = nrow(res_ms) - vr, col = myPal1[ms])
      points(x = res_ms[vr, "HR95"], y = nrow(res_ms) - vr, pch = c(15, 17)[ms], 
             cex = 1.2, col = myPal1[ms])
      segments(x0 = res_ms[vr, "HR95_LCI"], x1 = res_ms[vr, "HR95_LCI"], 
               y0 = nrow(res_ms) - vr + 0.15, y1 = nrow(res_ms) - vr - 0.15, 
               col = myPal1[ms])
      segments(x0 = res_ms[vr, "HR95_UCI"], x1 = res_ms[vr, "HR95_UCI"], 
               y0 = nrow(res_ms) - vr + 0.15, y1 = nrow(res_ms) - vr - 0.15, 
               col = myPal1[ms])
    }
    
    # Add p values for PH assumption test
    mtext(text = c("Means", "CVs")[ms], side = 3, at = c(-5.6, -2)[ms],
          line = -0.25, cex = 0.6, adj = 0.5, col = myPal1[ms], font = 2)
    text(x = rep(c(-5.6, -2)[ms], nrow(res_ms)), y = (nrow(res_ms) - 1):0,
         labels = round_p_values(res_ms$zph_p), xpd = T, adj = 0.5, cex = 0.9, 
         col = myPal1[ms])
  }
  
  # Add variable names
  text(x = rep(x_txt, nrow(res_ms)), y = (nrow(res_ms) - 1):0,
       labels = text_labels[[um]], xpd = T, adj = 0)
  
  if (um==2) {
    legend(x = 4.8, y = 3, legend = c("Levels (means)", "Variability (CVs)"), 
           col = myPal1, pch = c(15, 17), cex = 0.9, pt.cex = 1.2, xpd = T)
  }
}


### Panel C (AUC) ###
####### Calculate means per year and PCs on these means #######
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

############# Cox and ROC models ###############
### Part 1 - Add other PCs one by one to PC1 and CVPC1
pcs <- c("PC1", "CVPC1")
res_i <- list()
for (i in 1:length(pcs)) {         
  if (i==1) vars <- paste("PC", 1:length(selected.biomarkersHD), sep = "")
  if (i==2) vars <- grep("CVPC", colnames(dat.both), value = T)
  
  res_k <- list()
  for (k in 1:2) {     # With and without covariates
    if (k==1) {
      covars <- ""
    } else {
      covars <- c("bs(age_visit, df = 5)", "sex", "diabetes", "fu_length")
    }
    
    res_j <- list()
    for (j in 1:length(selected.biomarkersHD)) {
      res <- as.data.frame(matrix(NA, 1, 6))
      colnames(res) <- c("Mean", "CV", "PCs", "controls", "Variable", "AUC")
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
      roc.mod <- roc(dat.both$status ~ predict(cox.mod))
      res$AUC <- as.numeric(auc(roc.mod))
      res$Variable <- pc
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
pcs <- c("PC1", "CVPC1")
res_i <- list()
for (i in 1:length(pcs)) {         
  pc <- pcs[i]
  other.pcs <- setdiff(vars[[ifelse(i==1, 1, 2)]], pc)
  
  res_k <- list()
  for (k in 1:2) {        # with and without covariates
    if (k==1) {
      covars <- ""
    } else {
      covars <- c("bs(age_visit, df = 5)", "sex", "diabetes", "fu_length")
    }
    
    res_j <- list()
    for (j in 1:length(selected.biomarkersHD)) {
      res <- as.data.frame(matrix(NA, 1, 6))
      colnames(res) <- c("Mean", "CV", "PCs", "controls", "Variable", "AUC")
      res$Mean <- 1
      res$CV <- 1
      res$PCs <- c(1:length(selected.biomarkersHD))[j]
      res$controls <- ifelse(k==1, 0, 1)
      
      if (j==1) {
        fmla <- formula(paste("Surv(years, status) ~ ",
                              ifelse(k==1, paste(c(pc, vars[[ifelse(i==1, 2, 1)]]), 
                                                 collapse = " + "),
                                     paste(c(pc, vars[[ifelse(i==1, 2, 1)]], covars), 
                                           collapse = " + ")),
                              sep = ""))
      } else {
        fmla <- formula(paste("Surv(years, status) ~ ",
                              ifelse(k==1, paste(c(pc, other.pcs[1:(j - 1)], 
                                                   vars[[ifelse(i==1, 2, 1)]]), 
                                                 collapse = " + "),
                                     paste(c(pc, other.pcs[1:(j - 1)], 
                                             vars[[ifelse(i==1, 2, 1)]], covars), 
                                           collapse = " + ")),
                              sep = ""))
      }
      
      # Run cox model
      cox.mod <- coxph(fmla, data = dat.both, weights = sqrt(nb_vis), cluster = id_no)
      roc.mod <- roc(dat.both$status ~ predict(cox.mod))
      res$AUC <- as.numeric(auc(roc.mod))
      res$Variable <- pc
      res_j[[j]] <- res
    }
    res_k[[k]] <- bind_rows(res_j)
  }
  res_i[[i]] <- bind_rows(res_k)
}
res_part2 <- bind_rows(res_i)


###### Run model with only control variables
res <- as.data.frame(matrix(NA, 1, 6))
colnames(res) <- c("Mean", "CV", "PCs", "controls", "Variable", "AUC")
res$Mean <- 0
res$CV <- 0
res$PCs <- 0
res$controls <- 1

cox.mod <- coxph(Surv(years, status) ~ bs(age_visit, df = 5) + sex + diabetes + fu_length, 
                 data = dat.both, weights = sqrt(nb_vis), cluster = id_no)
roc.mod <- roc(dat.both$status ~ predict(cox.mod))
res$AUC <- as.numeric(auc(roc.mod))


res_all <- rbind(res_part1, res_part2, res)
write.csv(res_all, "Figure2_panelC.csv", row.names = F)

## For PC1
dat.pc1 <- res_all[is.na(res_all$Variable) | res_all$Variable=="PC1", ]
# Add values for 0 PCs
dat.pc1 <- rbind(dat.pc1, 
                 res_all[res_all$Mean==0 & res_all$CV==1 & res_all$PCs==11 & 
                           res_all$Variable=="CVPC1", ])
dat.pc1[(nrow(dat.pc1) - 1):nrow(dat.pc1), "PCs"] <- 0
dat.pc1[nrow(dat.pc1) + 1, ] <- c(1, 0, 0, 0, "PC1", 0.5)
dat.pc1[ , -grep("Variable", colnames(dat.pc1))] <- 
  apply(dat.pc1[ , -grep("Variable", colnames(dat.pc1))], 2, function(x) as.numeric(x))

## For CVPC1
dat.cvpc1 <- res_all[is.na(res_all$Variable) | res_all$Variable=="CVPC1", ]
# Add values for 0 PCs
dat.cvpc1 <- rbind(dat.cvpc1, 
                   res_all[res_all$Mean==1 & res_all$CV==0 & res_all$PCs==11 & 
                             res_all$Variable=="PC1", ])
dat.cvpc1[(nrow(dat.cvpc1) - 1):nrow(dat.cvpc1), "PCs"] <- 0
dat.cvpc1[nrow(dat.cvpc1) + 1, ] <- c(0, 1, 0, 0, "CVPC1", 0.5)
dat.cvpc1[ , -grep("Variable", colnames(dat.cvpc1))] <- 
  apply(dat.cvpc1[ , -grep("Variable", colnames(dat.cvpc1))], 2, function(x) as.numeric(x))

# Plot
par(mar = c(2.5, 3, 2, 0.5), cex = 0.9)
cols <- list(c("#4053d3", "#26317e"), c("#cd1719", "#880f10"))

plot(1:10, 1:10, type = "n", xlim = c(0, length(selected.biomarkersHD)), 
     xlab = "# of PCs included", ylim = c(0.5, max(c(dat.pc1$AUC, dat.cvpc1$AUC)) + 0.02), 
     ylab = "AUC", frame.plot = F, cex.axis = 0.8, 
     xaxp = c(0, length(selected.biomarkersHD) - 1, 5))
## Add grid 
abline(h = c(0.5, 0.6, 0.7, 0.8), col = "grey85", lwd = 2)
abline(h = c(0.55, 0.65, 0.75, 0.85), col = "grey85")

poss <- list(c(0, 0), c(0, 1), c(1, 0), c(1, 1))
for (j in 1:2) {
  for (i in 1:4) {
    if (j==1) {
      dat.poss <- dat.pc1[dat.pc1$CV==poss[[i]][1] & dat.pc1$controls==poss[[i]][2], ]
    } else {
      dat.poss <- dat.cvpc1[dat.cvpc1$Mean==poss[[i]][1] & dat.cvpc1$controls==poss[[i]][2], ]
    }
    dat.poss <- dat.poss[order(dat.poss$PCs), ]  
    
    lines(x = dat.poss$PCs, y = dat.poss$AUC, col = cols[[j]][ifelse(i < 3, 1, 2)], 
          lty = ifelse(i %in% c(1, 3), 1, 2), lwd = 2) 
  }
}
mtext("C", side = 3, line = 1, at = -2.7, adj = 0, font = 2)

# Legend 
par(mar = c(0, 1, 4.5, 0.5), cex = 0.9)
plot(1:10, 1:10, type = "n", xlab = "", xaxt = "n", yaxt = "n", 
     ylab = "", frame.plot = F)
text(x = c(3, 3, 5, 5, 3), y = c(10, 9.4, 9.75, 9.05, 8.8), 
     labels = c("Primary", "variable", "Control", "variables", "(x-axis)"), 
     font = 2, cex = 0.8)
for (i in 1:length(unlist(cols))) {
  segments(x0 = 0, x1 = 1.5, y0 = c(8.2, 7.5, 6.8, 6.1)[i], 
           y1 = c(8.2, 7.5, 6.8, 6.1)[i], col = unlist(cols)[i], lwd = 2)
}
text(x = rep(3, 4), y = c(8.2, 7.5, 6.8, 6.1), 
     labels = c("PCs", "PCs", "CVPCs", "CVPCs"), cex = 0.75)
text(x = rep(5, 4), y = c(8.2, 7.5, 6.8, 6.1), 
     labels = c("none", "CVPCs (all)", "none", "PCs (all)"), cex = 0.75)

segments(x0 = 0, x1 = 1.5, y0 = 4.9, y1 = 4.9, lwd = 2)
segments(x0 = 0, x1 = 1.5, y0 = 4.2, y1 = 4.2, lty = 2, lwd = 2)
text(x = c(2.2, 2.2), y = c(4.9, 4.2), adj = 0, 
     labels = c("Without demographic variables", "With demographic variables"), 
     cex = 0.75)

dev.off()
par(opar)

