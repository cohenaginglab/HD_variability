#######################################################################################################
#####        This code is for the figure entitled "Relative biomarker contributions to PC1,       ##### 
#####                     PC2, PC7, and CVPC3 for the 4-month variable list"                      #####
#######################################################################################################
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(grid)
library(tidyverse)

rm(list=ls())

source("Functions_HD_variability_paper.R")

#### Load data ####
clean.HD$date_visit <- as.Date(clean.HD$date_visit)
clean.HD$date_death <- as.Date(clean.HD$date_death)

## variable lists
selected.biomarkersHD <- c("wbc", "hb", "hct", "mch", "mchc", "mcv", "pltlt", "rdw", "ca", "rbc", "creat", 
                           "Glucose", "k", "Phosphates.sériques", "sodium", "alb")
vars.nms <- c("WBC", "Hemoglobin", "Hematocrit", "MCH", "MCHC", "MCV", "Platelets", "RDW", "Calcium", 
              "RBC", "Creatinine", "Glucose", "Potassium", "Phosphate","Sodium", "Albumin")

# Select complete observations
clean.HD2 <- clean.HD[complete.cases(clean.HD[ , selected.biomarkersHD]), ]

# Transform variables not normally distributed
if ("pltlt" %in% selected.biomarkersHD) {
  clean.HD2$pltlt <- sqrt(clean.HD2$pltlt)
}
clean.HD2[ , which(colnames(clean.HD2) %in% c("wbc", "rdw", "Glucose"))] <- 
  log(clean.HD2[ , which(colnames(clean.HD2) %in% c("wbc", "rdw", "Glucose"))])

# Add date of last contact
datt <- ddply(clean.HD2, .(id_no), lastvis)

# Select time intervals (in months)
time.ints <- list(c(0, 6), c(12, 24), c(36, 60), c(84, 120))         
res_list <- list()
for (ti in 1:(length(time.ints) + 8)) {
  
  # select observations according to selected time frame
  if (ti==1) {
    dat <- datt 
  } else {
    if (ti==(length(time.ints) + 2)) {
      dat <- datt[datt$diabetes==0, ]
    } else {
      if (ti==(length(time.ints) + 3)) {
        dat <- datt[datt$diabetes==1, ]
      } else {
        if (ti==(length(time.ints) + 4)) {
          dat <- datt[datt$sex=="M", ]
        } else {
          if (ti==(length(time.ints) + 5)) {
            dat <- datt[datt$sex=="F", ]
          } else {
            if (ti==(length(time.ints) + 6)) {
              dat <- datt[datt$age_visit < 60, ]
            } else {
              if (ti==(length(time.ints) + 7)) {
                dat <- datt[datt$age_visit >= 60 & datt$age_visit < 75, ]
              } else {
                if (ti==(length(time.ints) + 8)) {
                  dat <- datt[datt$age_visit >=75, ]
                } else {
                  dat <- ddply(datt[!is.na(datt$date_death), ], .(id_no), pca_interval, time.interval = time.ints[[ti-1]])
                }
              }
            }
          }
        }
      }
    }
  }
  
  # Run PCA on raw biomarkers
  pca.all <- prcomp(dat[ , selected.biomarkersHD], center = T, scale. = T)
  
  # Calculate CVs per year 
  dat <- ddply(dat, .(id_no), cv.year, vars = selected.biomarkersHD, cv = T,
               min.vis.nb = 2, covars = c("id_no", "sex", "date_death", "date_birth", 
                                          "date_last", "diabetes", "fu_length"))
  # Apply log transformation on CVs
  dat[ , selected.biomarkersHD] <- apply(dat[ , selected.biomarkersHD], 2, log_zeros)
  # Remove outliers (due to only two equal values included in CV)
  for (j in 1:length(selected.biomarkersHD)) {
    dat[ , selected.biomarkersHD[j]] <- 
      ifelse(dat[ , selected.biomarkersHD[j]] %in% boxplot(dat[ , selected.biomarkersHD[j]], range = 6, plot = F)$out==T,
             NA, dat[ , selected.biomarkersHD[j]])
  }
  dat <- dat[complete.cases(dat[ , selected.biomarkersHD]), ]
  
  # Correct CVs for the number of visits included
  dat[ , selected.biomarkersHD] <- apply(dat[ , selected.biomarkersHD], 2, function(x) {
    residuals(nls(x ~ I(1 / sqrt(nb_vis) * a) + b, data = dat, start = list(a = 1, b = 1)))
  })
  
  # perform PCA
  pca_cv <- prcomp(dat[ , selected.biomarkersHD], center = T, scale. = T)
  
  res_list[[ti]] <- list("means" = pca.all, "cvs" = pca_cv)
}

cols <- colorRampPalette(brewer.pal(12, "Paired"))
myPal <- rev(cols(length(selected.biomarkersHD)))
names(myPal) <- vars.nms

##### Plot
tiff("Biom contrib_4mth.tiff", 220, 180, "mm", res = 600, compression = "lzw")
opar <- par(mfrow = c(2, 3), mar = c(5, 2.5, 1.5, 0.5), oma = c(0.2, 0.2, 0.2, 0.2), 
            mgp = c(1.5, 0.4, 0), tcl = -0.3, cex = 0.9)

pcs <- c(1, 2, 7, 3)
for (pc in 1:length(pcs)) {
  # Get variable contributions
  contribs <- list()
  for (ti in 1:length(res_list)) {
    xx <- as.data.frame(matrix(NA, length(selected.biomarkersHD), 3))
    colnames(xx) <- c("Variable", "Contribution", "Time")
    xx[ , "Variable"] <- vars.nms
    if (pc < 4) {
      xx[ , "Contribution"] <- abs(res_list[[ti]]$means$rotation[ , pcs[pc]]) / 
        sum(abs(res_list[[ti]]$means$rotation[ , pcs[pc]]))
    } else {
      xx[ , "Contribution"] <- abs(res_list[[ti]]$cvs$rotation[ , pcs[pc]]) / 
        sum(abs(res_list[[ti]]$cvs$rotation[ , pcs[pc]]))
    }
    xx[ , "Time"] <- ti
    contribs[[ti]] <- xx
  }
  pc.contrib <- bind_rows(contribs)
  pc.contrib$Time <- factor(pc.contrib$Time, levels = 1:(length(time.ints) + 8), 
                            labels = c("All", "0-6 months", "1-2 years", "3-5 years", 
                                       "7-10 years", "Not diabetic", "Diabetic", "Men", 
                                       "Women", "< 60 yrs", "60-75 yrs", "75+ yrs"))
  
  # Order variables based on PCA "All" loadings (abs)
  pca.all <- pc.contrib[pc.contrib$Time=="All", ]
  order.all <- pca.all[order(pca.all$Contribution, decreasing = T), "Variable"]
  
  plot(seq(0.7, length(levels(pca.all$Time)), length.out = length(levels(pc.contrib$Time))), 
       seq(0, 1, length.out = length(levels(pc.contrib$Time))), type = "n", frame.plot = F,
       xlab = "", ylab = "Variable contribution", xaxt = "n")
  for (ti in 1:length(levels(pca.all$Time))) {
    
    for (vv in 1:length(selected.biomarkersHD)) {
      dat_time <- pc.contrib[pc.contrib$Time==levels(pc.contrib$Time)[ti], ]
      dat.ordered <- dat_time[match(order.all, dat_time$Variable), ]
      rect(xleft = (ti - 0.35), xright = (ti + 0.35), 
           ybottom = ifelse(vv==1, 0, sum(dat.ordered$Contribution[1:(vv - 1)])), 
           ytop = sum(dat.ordered$Contribution[1:vv]),
           col = myPal[which(names(myPal)==dat.ordered$Variable[vv])], 
           border = myPal[which(names(myPal)==dat.ordered$Variable[vv])])
    }
  }
  text(x = seq(1, length(levels(pca.all$Time)), length.out = length(levels(pc.contrib$Time))),
       y = rep(-0.02, length(levels(pca.all$Time))),
       labels = c("All", "0-0.5", "1-2", "3-5", "7-10", "No", "Yes", "M", "W", "<60", 
                  "60-75", "75+"), srt = 90, xpd = T, cex = 0.9, adj = 1)
  segments(2, -0.205, 5, -0.205, lwd = 2, xpd = T)
  segments(6, -0.205, 7, -0.205, lwd = 2, xpd = T)
  segments(8, -0.205, 9, -0.205, lwd = 2, xpd = T)
  segments(10, -0.205, 12, -0.205, lwd = 2, xpd = T)
  mtext(c("Years", "Diabetes", "Sex", "Age group"), side = 1, at = c(3.5, 6.5, 8.5, 11),
        line = 2.2, cex = 0.7, font = 2)
  mtext(c("before death", "(years)"), side = 1, at = c(3.5, 11), line = 3, cex = 0.7, 
        font = 2)
  mtext(paste(c("A", "B", "C", "D")[pc], " - ", ifelse(pc < 4, "PC", "CVPC"),
              pcs[pc], sep = ""), side = 3, at = -1.5, adj = 0, line = 0.5, 
        font = 2)
  
  if (pc==2) {
    par(mar = c(0.5, 0.5, 1.5, 0.5))
    plot(1:10, 1:10, type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "",
         frame.plot = F)
    legend("topleft", legend = vars.nms, col = myPal[which(names(myPal)==vars.nms)], 
           pch = 15, pt.cex = 2, bty = "n", cex = 0.9)
    par(mar = c(4, 2.5, 1.5, 0.5))
  }
}
dev.off()
par(opar)
