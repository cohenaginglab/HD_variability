#######################################################################################################
#####                 This code is for the figure entitled "Effect of the number of               ##### 
#####                          observations included in CV calculation"                           #####
#######################################################################################################

rm(list = ls())

source("Functions_HD_variability_paper.R")

#### Load data ####
clean.HD$date_visit <- as.Date(clean.HD$date_visit)
clean.HD$date_death <- as.Date(clean.HD$date_death)

## variable list
vars.2wks <- c("wbc", "hb", "hct", "mch", "mchc", "mcv", "pltlt", "rdw", "k", 
               "sodium", "rbc")

# Transform variables not normally distributed
if ("pltlt" %in% vars.2wks) {
  clean.HD$pltlt <- sqrt(clean.HD$pltlt)
}
clean.HD[ , which(colnames(clean.HD) %in% c("wbc", "rdw", "Glucose"))] <- 
  log(clean.HD[ , which(colnames(clean.HD) %in% c("wbc", "rdw", "Glucose"))])

# Add date of last contact
dat <- ddply(clean.HD, .(id_no), lastvis)


########### First part: plot CV vs number of observations included
# Calculate CVs per year 
dat.cvs <- ddply(dat, .(id_no), cv.year, vars = vars.2wks, 
                 cv = T, min.vis.nb = 2, covars = c("id_no", "sex", "date_death", 
                                                    "date_birth", "date_last", 
                                                    "diabetes", "fu_length"))
# Apply log transformation on CVs
dat.cvs[ , vars.2wks] <- apply(dat.cvs[ , vars.2wks], 2, log_zeros)

# Some variables have no variation, leading to outliers
dat.cvs[ , vars.2wks] <- apply(dat.cvs[ , vars.2wks], 2, uniform_val)

### Generate simulated data# Get means and SDs of biomarkers (2.5% to 97.5% percentiles)
biom_means <- apply(clean.HD[ , vars.2wks], 2, mean)
biom_SDs <- apply(clean.HD[ , vars.2wks], 2, sd)

dat.sim <- clean.HD[ , -which(colnames(clean.HD) %in% vars.2wks)]
for (v in 1:length(vars.2wks)) {
  set.seed(v)
  dat.sim[ , vars.2wks[v]] <- rnorm(nrow(dat.sim), biom_means[v], biom_SDs[v])
}

# Add date of last contact
dat2 <- ddply(dat.sim, .(id_no), lastvis)

# Calculate CVs per year 
dat.cvs2 <- ddply(dat2, .(id_no), cv.year, vars = vars.2wks, 
                 cv = T, min.vis.nb = 2, covars = c("id_no", "sex", "date_death", 
                                                    "date_birth", "date_last", 
                                                    "diabetes", "fu_length"))
# Apply log transformation on CVs
dat.cvs2[ , vars.2wks] <- apply(dat.cvs2[ , vars.2wks], 2, log_zeros)

# Some variables have no variation, leading to outliers
dat.cvs2[ , vars.2wks] <- apply(dat.cvs2[ , vars.2wks], 2, uniform_val)


########### Second part: replicate some analyses with at least 5 obs. in CVs
#### Association with mortality for all CVPCs
# Calculate CVs per year 
dat.cvs3 <- ddply(dat, .(id_no), cv.year, vars = vars.2wks, 
                  cv = T, min.vis.nb = 5, covars = c("id_no", "sex", "date_death", 
                                                     "date_birth", "date_last", 
                                                     "diabetes", "fu_length"))
# Apply log transformation on CVs
dat.cvs3[ , vars.2wks] <- apply(dat.cvs3[ , vars.2wks], 2, log_zeros)

# Some variables have no variation, leading to outliers
dat.cvs3[ , vars.2wks] <- apply(dat.cvs3[ , vars.2wks], 2, uniform_val)


# Add "status" variable (for cox regressions)
dat.cox.cv <- ddply(dat.cvs3, .(id_no), status_fct)

# Calculate PCA on CVs
pca_cv <- prcomp(dat.cox.cv[ , vars.2wks], center = T, scale. = T)
dat.pc <- cbind(dat.cox.cv, pca_cv$x)
for (pc in 1:length(vars.2wks)) {
  dat.pc <- pc_sign(dat.pc, paste("PC", pc, sep = ""), pca_cv)
}

# Add "CV" to column names
colnames(dat.pc)[which(colnames(dat.pc) %in% vars.2wks)] <- 
  unlist(lapply(colnames(dat.pc)[which(colnames(dat.pc) %in% vars.2wks)], 
                function (x) paste(x, "CV", sep = "_")))
colnames(dat.pc)[grep("PC", colnames(dat.pc))] <- 
  unlist(lapply(colnames(dat.pc)[grep("PC", colnames(dat.pc))], 
                function (x) paste("CV", x, sep = "")))
# Re-order columns (so that CVPCs are from 1 to 11)
dat.pc <- dat.pc[ , c(setdiff(colnames(dat.pc), grep("CVPC", colnames(dat.pc), value = T)),
                      paste("CVPC", 1:length(vars.2wks), sep = ""))]

### Cox models for all CVPCs
selected.biomarkersHD <- vars.2wks
res_list <- list()
for (vr in 1:length(selected.biomarkersHD)) {
  var <- paste("CVPC", 1:length(selected.biomarkersHD), sep = "")[vr]
  res <- as.data.frame(matrix(NA, 1, 8))
  colnames(res) <- c("Variable", "HR", "LCI", "UCI", "HR95", "HR95_LCI", 
                     "HR95_UCI", "zph_p")
  res$Variable <- var
  dat_var <- dat.pc[ , c("id_no", "years", "status", "age_visit", "sex", "diabetes", 
                         "fu_length", "nb_vis", var)]
  colnames(dat_var)[ncol(dat_var)] <- "predictor"
  mod <- coxph(Surv(years, status) ~ predictor + bs(age_visit, df = 5) + sex +  
                 diabetes + fu_length, data = dat_var, weights = sqrt(nb_vis),
               cluster = id_no)
  res[ , c("HR", "LCI", "UCI")] <- summary(mod)$conf.int[1, c(1, 3, 4)]
  res[ , c("HR95", "HR95_LCI", "HR95_UCI")] <- summary(mod)$conf.int[1, c(1, 3, 4)] ^ 
    (quantile(na.omit(dat_var[ , "predictor"]), 0.975) - 
       quantile(na.omit(dat_var[ , "predictor"]), 0.025))
  
  # Test the proportional hazards assumption of the cox model
  res$zph_p <- cox.zph(mod)$table[1, 3]
  res_list[[vr]] <- res
}
res_all_CVPCs <- do.call("rbind", res_list)


#### Association with mortality in different population subsets for CVPC1
res_list <- list()
for (ti in 1:8) {
  
  # select population subset
  if (ti==1) {
    dat4 <- dat 
  } else {
    if (ti==2) {
      dat4 <- dat[dat$diabetes==0, ]
    } else {
      if (ti==3) {
        dat4 <- dat[dat$diabetes==1, ]
      } else {
        if (ti==4) {
          dat4 <- dat[dat$sex=="M", ]
        } else {
          if (ti==5) {
            dat4 <- dat[dat$sex=="F", ]
          } else {
            if (ti==6) {
              dat4 <- dat[dat$age_visit < 60, ]
            } else {
              if (ti==7) {
                dat4 <- dat[dat$age_visit >= 60 & dat$age_visit < 75, ]
              } else {
                if (ti==8) {
                  dat4 <- dat[dat$age_visit >=75, ]
                }
              }
            }
          }
        }
      }
    }
  }
  
  # Calculate CVs per year 
  dat.cvs4 <- ddply(dat4, .(id_no), cv.year, vars = selected.biomarkersHD, 
                   cv = T, min.vis.nb = 5, covars = c("id_no", "sex", "date_death", 
                                                      "date_birth", "date_last", 
                                                      "diabetes", "fu_length"))
  # Apply log transformation on CVs
  dat.cvs4[ , selected.biomarkersHD] <- apply(dat.cvs4[ , selected.biomarkersHD], 2, 
                                             log_zeros)
  
  # Some variables have no variation, leading to outliers
  dat.cvs4[ , selected.biomarkersHD] <- apply(dat.cvs4[ , selected.biomarkersHD], 2, 
                                             uniform_val)
  
  # Add "status" variable (for cox regressions)
  dat.cox.cv <- ddply(dat.cvs4, .(id_no), status_fct)
  
  # Calculate PCA on CVs
  pca_cv <- prcomp(dat.cox.cv[ , selected.biomarkersHD], center = T, scale. = T)
  dat.pc <- cbind(dat.cox.cv, pca_cv$x)
  # Add "CV" to column names
  colnames(dat.pc)[which(colnames(dat.pc) %in% selected.biomarkersHD)] <- 
    unlist(lapply(colnames(dat.pc)[which(colnames(dat.pc) %in% selected.biomarkersHD)], 
                  function (x) paste(x, "CV", sep = "_")))
  colnames(dat.pc)[grep("PC", colnames(dat.pc))] <- 
    unlist(lapply(colnames(dat.pc)[grep("PC", colnames(dat.pc))], 
                  function (x) paste("CV", x, sep = "")))
  dat.pc <- pc_sign2(dat.pc, "CVPC1", pca_cv)
  
  # Re-order columns (so that CVPCs are from 1 to 11)
  dat.cox <- dat.pc[ , c(setdiff(colnames(dat.pc), grep("CVPC", colnames(dat.pc), value = T)),
                         paste("CVPC", 1:length(selected.biomarkersHD), sep = ""))]
  
  # Survival model
  # Prepare results dataframe
  res <- as.data.frame(matrix(NA, 1, 9))
  colnames(res) <- c("Subset", "HR", "LCI", "UCI", "p_value", 
                     "HR95", "HR95_LCI", "HR95_UCI", "zph_p")
  res$Subset <- c("All", "Not diabetic", "Diabetic", "Men", 
                  "Women", "< 60 yrs", "60-75 yrs", "75+ yrs")[ti]
  
  # Select relevant covariables for each population subset
  if (ti==1) covars <- c("bs(age_visit, df = 5)", "sex", "diabetes", "fu_length")
  if (ti %in% 2:3) covars <- c("bs(age_visit, df = 5)", "sex", "fu_length")
  if (ti %in% 4:5) covars <- c("bs(age_visit, df = 5)", "diabetes", "fu_length")
  if (ti > 5) covars <- c("bs(age_visit, df = 3)", "sex", "diabetes", "fu_length")
  fmla <- formula(paste("Surv(years, status) ~ ",
                        paste(c("CVPC1", covars), collapse = " + "),
                        sep = ""))
  
  # Run Cox model
  cox.mod <- coxph(fmla, data = dat.cox, weights = sqrt(nb_vis), cluster = id_no)
  
  res[ , c("HR", "LCI", "UCI")] <- summary(cox.mod)$conf.int[1:nrow(res), c(1, 3, 4)]
  res[ , "p_value"] <- summary(cox.mod)$coefficients[1, 6]
  res[ , c("HR95", "HR95_LCI", "HR95_UCI")] <- summary(cox.mod)$conf.int[1, c(1, 3, 4)] ^ 
    (quantile(dat.cox[ , "CVPC1"], 0.975) - quantile(dat.cox[ , "CVPC1"], 0.025))
  res$zph_p <- cox.zph(cox.mod)$table[1, 3]
  
  res_list[[ti]] <- res
}
res_all_CVPC1 <- do.call("rbind", res_list)


###### Plot
tiff("CV correct.tiff", width = 200, height = 160, "mm", 
     res = 300, compression = "lzw")
opar <- par(mfrow = c(2, 2), mar = c(2.5, 2.5, 1.5, 1.5), oma = c(0.2, 0.2, 0.2, 0.2), 
            mgp = c(1.5, 0.4, 0), tcl = -0.3) 

## Panels A and B: bias in CV, taking hemoglobin as an example
for (j in 1:2) {
  if (j==1) dat <- dat.cvs else dat <- dat.cvs2
  
  y <- dat[ , "hb"]
  nb <- dat$nb_vis
  plot(sqrt(nb), y, las = 1, type = "n", xlab = "Number of observations", 
       xaxt = "n", ylab = "CV (log)", ylim = c(-6, -1), frame.plot = F)
  abline(v = c(2, 4, 6, 8), col = "lightgrey")
  points(sqrt(nb), y, pch = 19, cex = 0.2)
  mtext(c("A", "B")[j], side = 3, line = 0.5, at = 0, adj = 0, font = 2)
  # Add x-axis manually
  segments(x0 = 2, x1 = 8, y0 = -6.2, y1 = -6.2, xpd = T)
  for (k in c(2, 4, 6, 8)) {
    segments(x0 = k, x1 = k, y0 = -6.2, y1 = -6.3, xpd = T)
  }
  text(x = c(2, 4, 6, 8), y = rep(-6.45, 4), labels = c(2, 4, 6, 8)^2, xpd = T)
  
  # Add regression line
  nb_sqrt_inv <- 1/sqrt(dat$nb_vis)
  regNLin <- lm(y ~ nb_sqrt_inv)
  pred <- predict.lm(regNLin, interval = "pred")
  nb_sort <- sort(nb[which(!is.na(nb))])
  nb_ord <- order(nb[which(!is.na(nb))])
  lines(sqrt(nb_sort), pred[nb_ord, 1], col = "indianred", lwd = 2)
}


## Panels C and D: using CVs calculated with at least 5 observations
par(mar = c(3, 10, 4, 1.5))

for (j in 1:2) {
  if (j==1) res_pc <- res_all_CVPCs else res_pc <- res_all_CVPC1
  
  # Set x-axis range 
  x_range <- c(0, max(res_pc$HR95_UCI))
  
  plot(1:10, 1:10, type = "n", xlim = x_range,  
       ylim = c(0, (nrow(res_pc) - 1)), ylab = "", xlab = "HR95", yaxt = "n", 
       frame.plot = F)
  abline(v = 1, lty = 3)
  mtext(c("C", "D")[j], side = 3, line = 1, at = c(-14.5, -27)[j], adj = 0, font = 2)
  
  for (vr in 1:nrow(res_pc)) {
    segments(x0 = res_pc[vr, "HR95_LCI"], x1 = res_pc[vr, "HR95_UCI"], 
             y0 = nrow(res_pc) - vr, y1 = nrow(res_pc) - vr)
    points(x = res_pc[vr, "HR95"], y = nrow(res_pc) - vr, pch = 15, cex = 1.2)
    segments(x0 = res_pc[vr, "HR95_LCI"], x1 = res_pc[vr, "HR95_LCI"], 
             y0 = nrow(res_pc) - vr + 0.15, y1 = nrow(res_pc) - vr - 0.15)
    segments(x0 = res_pc[vr, "HR95_UCI"], x1 = res_pc[vr, "HR95_UCI"], 
             y0 = nrow(res_pc) - vr + 0.15, y1 = nrow(res_pc) - vr - 0.15)
  }
  # Add variable names
  if (j==1) {
    text(x = rep(c(-12, -22)[j], nrow(res_pc)), y = (nrow(res_pc) - 1):0,
         labels = res_pc$Variable, xpd = T, adj = 0)
  } else {
    text(x = rep(c(-12, -22)[j], nrow(res_pc)), y = (nrow(res_pc) - 1):0,
         labels = res_pc$Subset, xpd = T, adj = 0)
  }
  
  # Add p values for PH assumption test
  mtext(text = c("PH", "assumption", "test"), side = 3, at = c(-2, -3.5)[j],
        line = c(2, 1.3, 0.5), cex = 0.6, adj = 0.5)
  text(x = rep(c(-2, -3.5)[j], nrow(res_pc)), y = (nrow(res_pc) - 1):0,
       labels = round_p_values(res_pc$zph_p), xpd = T, adj = 0.5, cex = 0.9)
  
}

dev.off()
par(opar)
