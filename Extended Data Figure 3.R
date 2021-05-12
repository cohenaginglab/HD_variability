library(survival)
library(splines)

rm(list=ls())

source("Functions_HD_variability_paper.R")

#### Load data ####
clean.HD$date_visit <- as.Date(clean.HD$date_visit)
clean.HD$date_death <- as.Date(clean.HD$date_death)

## variable lists
selected.biomarkersHD <-  c("wbc", "hb", "hct", "mch", "mchc", "mcv", "pltlt", "rdw", "k", "sodium", "rbc")

# Transform variables not normally distributed
if ("pltlt" %in% selected.biomarkersHD) {
  clean.HD$pltlt <- sqrt(clean.HD$pltlt)
}
clean.HD[ , which(colnames(clean.HD) %in% c("wbc", "rdw", "Glucose"))] <- 
  log(clean.HD[ , which(colnames(clean.HD) %in% c("wbc", "rdw", "Glucose"))])

# Add date of last contact
datt <- ddply(clean.HD, .(id_no), lastvis)

########### Run Cox models for different PCs/CVPCs across different population subsets
res_list <- list()
pcs <- c("PC1", "PC2", "PC6", "CVPC1", "CVPC3")
for (k in 1:length(pcs)) {
  res_pc <- list()
  for (ti in 1:8) {
    
    # select population subset
    if (ti==1) {
      dat <- datt 
    } else {
      if (ti==2) {
        dat <- datt[datt$diabetes==0, ]
      } else {
        if (ti==3) {
          dat <- datt[datt$diabetes==1, ]
        } else {
          if (ti==4) {
            dat <- datt[datt$sex=="M", ]
          } else {
            if (ti==5) {
              dat <- datt[datt$sex=="F", ]
            } else {
              if (ti==6) {
                dat <- datt[datt$age_visit < 60, ]
              } else {
                if (ti==7) {
                  dat <- datt[datt$age_visit >= 60 & datt$age_visit < 75, ]
                } else {
                  if (ti==8) {
                    dat <- datt[datt$age_visit >=75, ]
                  }
                }
              }
            }
          }
        }
      }
    }
    
    if (k < 4) {
      # Run PCA on raw biomarkers
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
      dat.cox <- pc_sign2(dat.cox, pcs[k], pca.all)
      # Re-order columns (so that PCs are from 1 to 11)
      dat.cox <- dat.cox[ , c(setdiff(colnames(dat.cox), 
                                      paste("PC", 1:length(selected.biomarkersHD), sep = "")),
                              paste("PC", 1:length(selected.biomarkersHD), sep = ""))]
    } else {
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
      # Add "CV" to column names
      colnames(dat.pc)[which(colnames(dat.pc) %in% selected.biomarkersHD)] <- 
        unlist(lapply(colnames(dat.pc)[which(colnames(dat.pc) %in% selected.biomarkersHD)], 
                      function (x) paste(x, "CV", sep = "_")))
      colnames(dat.pc)[grep("PC", colnames(dat.pc))] <- 
        unlist(lapply(colnames(dat.pc)[grep("PC", colnames(dat.pc))], 
                      function (x) paste("CV", x, sep = "")))
      dat.pc <- pc_sign2(dat.pc, pcs[k], pca_cv)
      
      # Re-order columns (so that CVPCs are from 1 to 11)
      dat.cox <- dat.pc[ , c(setdiff(colnames(dat.pc), grep("CVPC", colnames(dat.pc), value = T)),
                             paste("CVPC", 1:length(selected.biomarkersHD), sep = ""))]
    }
    
    # Survival model
    # Prepare results dataframe
    res <- as.data.frame(matrix(NA, 1, 10))
    colnames(res) <- c("Variable", "Subset", "HR", "LCI", "UCI", "p_value", 
                       "HR95", "HR95_LCI", "HR95_UCI", "zph_p")
    res$Variable <- pcs[k]
    res$Subset <- c("All", "Not diabetic", "Diabetic", "Men", 
                    "Women", "< 60 yrs", "60-75 yrs", "75+ yrs")[ti]
    
    # Select relevant covariables for each population subset
    if (ti==1) covars <- c("bs(age_visit, df = 5)", "sex", "diabetes", "fu_length")
    if (ti %in% 2:3) covars <- c("bs(age_visit, df = 5)", "sex", "fu_length")
    if (ti %in% 4:5) covars <- c("bs(age_visit, df = 5)", "diabetes", "fu_length")
    if (ti > 5) covars <- c("bs(age_visit, df = 3)", "sex", "diabetes", "fu_length")
    fmla <- formula(paste("Surv(years, status) ~ ",
                          paste(c(pcs[k], covars), collapse = " + "),
                          sep = ""))
    
    # Run Cox and ROC models
    cox.mod <- coxph(fmla, data = dat.cox, weights = sqrt(nb_vis), cluster = id_no)
    
    res[ , c("HR", "LCI", "UCI")] <- summary(cox.mod)$conf.int[1:nrow(res), c(1, 3, 4)]
    res[ , "p_value"] <- summary(cox.mod)$coefficients[1, 6]
    res[ , c("HR95", "HR95_LCI", "HR95_UCI")] <- summary(cox.mod)$conf.int[1, c(1, 3, 4)] ^ 
      (quantile(dat.cox[ , pcs[k]], 0.975) - quantile(dat.cox[ , pcs[k]], 0.025))
    res$zph_p <- cox.zph(cox.mod)$table[1, 3]
    
    res_pc[[ti]] <- res
  }
  res_list[[k]] <- do.call("rbind", res_pc)
  
}
res_all <- do.call("rbind", res_list)
write.csv(res_all, "ExtDatFigure3.csv", row.names = F)

#### Plot results #####
tiff("Extended Data_Figure_3.tiff", width = 180, height = 240, "mm", res = 300, compression = "lzw")
opar <- par(mfrow = c(3, 2), mar = c(3, 11, 5, 1.5), oma = c(0.2, 0.2, 0.2, 0.2), 
            mgp = c(1.5, 0.4, 0), tcl = -0.3) 

# Set x-axis range 
x_range <- c(0, max(res_all$HR95_UCI))

pcs <- c("PC1", "PC2", "PC6", "CVPC1", "CVPC3")
for (pc in 1:length(pcs)) {
  res_pc <- res_all[res_all$Variable==pcs[pc], ]
  
  plot(1:10, 1:10, type = "n", xlim = x_range,  
       ylim = c(0, (nrow(res_pc) - 1)), ylab = "", xlab = "HR95", yaxt = "n", 
       frame.plot = F)
  abline(v = 1, lty = 3)
  mtext(paste(c("A", "B", "C", "D", "E")[pc], pcs[pc], sep = " - "), side = 3, 
              line = 3, at = -21, adj = 0, font = 2, xpd = T)
  
  for (vr in 1:nrow(res_pc)) {
    segments(x0 = res_pc[vr, "HR95_LCI"], x1 = res_pc[vr, "HR95_UCI"], 
             y0 = nrow(res_pc) - vr, y1 = nrow(res_pc) - vr)
    points(x = res_pc[vr, "HR95"], y = nrow(res_pc) - vr, pch = 15, cex = 1.2)
    segments(x0 = res_pc[vr, "HR95_LCI"], x1 = res_pc[vr, "HR95_LCI"], 
             y0 = nrow(res_pc) - vr + 0.15, y1 = nrow(res_pc) - vr - 0.15)
    segments(x0 = res_pc[vr, "HR95_UCI"], x1 = res_pc[vr, "HR95_UCI"], 
             y0 = nrow(res_pc) - vr + 0.15, y1 = nrow(res_pc) - vr - 0.15)
  }
  # Add subset names
  text(x = rep(-18, nrow(res_pc)), y = (nrow(res_pc) - 1):0,
       labels = res_pc$Subset, xpd = T, adj = 0, cex = 1.2)
  
  # Add p values for PH assumption test
  mtext(text = c("PH", "assumption", "test"), side = 3, at = -3.5,
        line = c(2, 1.3, 0.5), cex = 0.7, adj = 0.5)
  text(x = rep(-3.5, nrow(res_pc)), y = (nrow(res_pc) - 1):0,
       labels = round_p_values(res_pc$zph_p), xpd = T, adj = 0.5, cex = 1.1)
  
}

dev.off()
par(opar)



