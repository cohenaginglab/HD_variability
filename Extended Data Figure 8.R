library(ggplot2)
library(ggpubr)

rm(list = ls())

source("Functions_HD_variability_paper.R")

#### Load biomarker data ####
clean.HD$date_visit <- as.Date(clean.HD$date_visit)
clean.HD$date_death <- as.Date(clean.HD$date_death)
#### Load hospitalisation dates ####
dat.HO$Date_début_visit_HD <- as.Date(dat.HO$Date_début_visit_HD)
dat.HO$Date_fin_visit_HD <- as.Date(dat.HO$Date_fin_visit_HD)


## variable list
selected.biomarkersHD <- c("wbc", "hb", "hct", "mch", "mchc", "mcv", "pltlt", "rdw", "k", "sodium", "rbc")
vars.nms = c("WBC", "Hemoglobin", "Hematocrit", "MCH", "MCHC", "MCV", "Platelets", "RDW", 
             "Potassium", "Sodium", "RBC")

# Transform variables not normally distributed
if ("pltlt" %in% selected.biomarkersHD) clean.HD$pltlt <- sqrt(clean.HD$pltlt)
clean.HD[ , which(colnames(clean.HD) %in% c("wbc", "rdw", "Glucose"))] <- 
  log(clean.HD[ , which(colnames(clean.HD) %in% c("wbc", "rdw", "Glucose"))])

# Add date of last contact
dat <- ddply(clean.HD, .(id_no), lastvis)

# Calculate CVs per year
dat.cvs <- ddply(dat, .(id_no), cv.year, vars = selected.biomarkersHD, cv = T, min.vis.nb = 2, 
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
dat.cox <- ddply(dat.cvs, .(id_no), status_fct)

# Calculate PCA on CVs
pca_cv <- prcomp(dat.cox[ , selected.biomarkersHD], center = T, scale. = T)
dat.pc <- cbind(dat.cox, pca_cv$x)
dat.pc <- pc_sign(dat.pc, "PC1", pca_cv)
colnames(dat.pc)[grep("PC", colnames(dat.pc))] <- 
  unlist(lapply(colnames(dat.pc)[grep("PC", colnames(dat.pc))], 
                function (x) paste("CV", x, sep = "")))

# Select random patients with at least 6 CVPC1 values
dat.pc$date_birth <- as.Date(dat.pc$date_birth)
# Cut time at 10 years before death/censoring
dat.10yr <- dat.pc[dat.pc$years <= 10, ]

datt <- list("dead" = dat.10yr[!is.na(dat.10yr$date_death), ], 
             "alive" = dat.10yr[is.na(dat.10yr$date_death), ])
set.seed(123)
for (i in 1:2) {
  dat <- ddply(datt[[i]], .(id_no), function(x) {
    x$vis_nb <- nrow(x)
    if (nrow(x) > 6)
    return(x) } )
  
  dat.db.diag <- split(dat, dat$diabetes)
  dat.db.diag2 <- lapply(dat.db.diag, function(x) {
    ids <- sample(x[!duplicated(x$id_no), "id_no"], 8)
    x <- x[x$id_no %in% ids, ]
  })
  dat1 <- bind_rows(dat.db.diag2)
  # Attribute dates to CVPC1 values using "years" and "date_last"
  dat1$date_CVPC1 <- as.Date(ymd(dat1$date_last) - months(dat1$years * 12 + 6), origin = "1970-01-01")
  dat1$yrs_CVPC1 <- as.numeric(dat1$date_last - dat1$date_CVPC1) / 365.25
  
  # get HO dates for these individuals
  hos <- list()
  for (id in 1:length(unique(dat1$id_no))) {
    dat.id <- dat1[dat1$id_no==unique(dat1$id_no)[id], ]
    ID <- as.character(unique(clean.HD[clean.HD$id_no==dat.id$id_no[1], "ID_patient"]))
    ho.id <- dat.HO[dat.HO$ID_patien==ID, ]
    ho.id$yrs_HO_st <- as.numeric(dat.id$date_last[1] - ho.id$Date_début_visit_HD) / 365.25
    ho.id$yrs_HO_end <- as.numeric(dat.id$date_last[1] - ho.id$Date_fin_visit_HD) / 365.25
    ho.id <- ho.id[which(ho.id$yrs_HO_st < max(dat.id$yrs_CVPC1)), ]
    ho.id$id_no <- dat.id$id_no[1]
    hos[[id]] <- ho.id
  }
  hos <- bind_rows(hos)
  
  max_ho <- max(tapply(hos$Date_début_visit_HD, hos$id_no, function(x) length(unique(x))))
  dat1[ , paste("HO_st", 1:max_ho, sep = "_")] <- NA
  dat1[ , paste("HO_end", 1:max_ho, sep = "_")] <- NA
  
  # Merge HOs info with data
  for (id in 1:length(unique(dat1$id_no))) {
    ho.id <- hos[hos$id_no==unique(dat1$id_no)[id], ]
    for (k in 1:nrow(ho.id)) {
      dat1[dat1$id_no==unique(dat1$id_no)[id], paste("HO_st", k, sep = "_")] <- ho.id[k, "yrs_HO_st"]
      dat1[dat1$id_no==unique(dat1$id_no)[id], paste("HO_end", k, sep = "_")] <- ho.id[k, "yrs_HO_end"]
    }
  }
  datt[[i]] <- dat1
}

#### Plot
cols <- c("#be2409", "#00baff")

# Set axis limits
datt2 <- lapply(datt, function(x) {
  x <- ddply(x, .(id_no), function(y) {
    y$yrs_follow <- ifelse(!is.na(y$date_death), 
                           as.numeric(unique(y$date_death) - unique(y$date_birth)) / 365.25 - min(y$age_visit),
                           as.numeric(unique(y$date_last) - unique(y$date_birth)) / 365.25 - min(y$age_visit))
    return(y)
  })
  return(x)
})
xlims <- c(0, 10)
ylims <- c(min(unlist(lapply(datt, function(x) min(x$CVPC1)))), max(unlist(lapply(datt, function(x) max(x$CVPC1)))))


plots <- list()
for (dd in 1:2) {
  for (id in 1:16) {
    dat.id <- datt[[dd]][datt[[dd]]$id_no==unique(datt[[dd]]$id_no)[id], ]
    dat.id <- dat.id[order(dat.id$age_visit), ]
    
    plots[[ifelse(dd==1, id, id + 16)]] <- ggplot(dat.id, aes(x = yrs_CVPC1, y = CVPC1)) +  
      # The number of "geom_rect" elements has to match the maximum number of HO event (depends on random sample) 
      geom_rect(aes(xmin = HO_st_1, xmax = HO_end_1, ymin = ylims[1], ymax = ylims[2]), fill = "#96ba75") +
      geom_rect(aes(xmin = HO_st_2, xmax = HO_end_2, ymin = ylims[1], ymax = ylims[2]), fill = "#96ba75") +
      geom_rect(aes(xmin = HO_st_3, xmax = HO_end_3, ymin = ylims[1], ymax = ylims[2]), fill = "#96ba75") +
      geom_rect(aes(xmin = HO_st_4, xmax = HO_end_4, ymin = ylims[1], ymax = ylims[2]), fill = "#96ba75") +
      geom_rect(aes(xmin = HO_st_5, xmax = HO_end_5, ymin = ylims[1], ymax = ylims[2]), fill = "#96ba75") +
      geom_rect(aes(xmin = HO_st_6, xmax = HO_end_6, ymin = ylims[1], ymax = ylims[2]), fill = "#96ba75") +
      geom_rect(aes(xmin = HO_st_7, xmax = HO_end_7, ymin = ylims[1], ymax = ylims[2]), fill = "#96ba75") +
      geom_rect(aes(xmin = HO_st_8, xmax = HO_end_8, ymin = ylims[1], ymax = ylims[2]), fill = "#96ba75") +
      geom_rect(aes(xmin = HO_st_9, xmax = HO_end_9, ymin = ylims[1], ymax = ylims[2]), fill = "#96ba75") +
      geom_rect(aes(xmin = HO_st_10, xmax = HO_end_10, ymin = ylims[1], ymax = ylims[2]), fill = "#96ba75") +
      geom_rect(aes(xmin = HO_st_11, xmax = HO_end_11, ymin = ylims[1], ymax = ylims[2]), fill = "#96ba75") +
      geom_rect(aes(xmin = HO_st_12, xmax = HO_end_12, ymin = ylims[1], ymax = ylims[2]), fill = "#96ba75") +
      geom_rect(aes(xmin = HO_st_13, xmax = HO_end_13, ymin = ylims[1], ymax = ylims[2]), fill = "#96ba75") +
      geom_rect(aes(xmin = HO_st_14, xmax = HO_end_14, ymin = ylims[1], ymax = ylims[2]), fill = "#96ba75") +
      #geom_rect(aes(xmin = HO_st_15, xmax = HO_end_15, ymin = ylims[1], ymax = ylims[2]), fill = "#96ba75") +
      #geom_rect(aes(xmin = HO_st_16, xmax = HO_end_16, ymin = ylims[1], ymax = ylims[2]), fill = "#96ba75") +
      #geom_rect(aes(xmin = HO_st_17, xmax = HO_end_17, ymin = ylims[1], ymax = ylims[2]), fill = "#96ba75") +
      #geom_rect(aes(xmin = HO_st_18, xmax = HO_end_18, ymin = ylims[1], ymax = ylims[2]), fill = "#96ba75") +
      geom_point(color = cols[dd], size = 0.7) + 
      geom_smooth(color = cols[dd], size = 0.5, method = "auto", se = F, 
                  linetype = unique(dat.id$diabetes) + 1) +
      scale_x_reverse(limits = rev(xlims)) + 
      scale_y_continuous(limits = ylims) + 
      theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks.length = unit(0.07, "cm"),
        axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5), axis.ticks = element_line(size = 0.5), 
        legend.position = "none", plot.margin = unit(c(0.12, 0.12, 0.12, 0.12), "lines"), 
        panel.background = element_rect(fill = "white", colour = "gray"))
  }
}

tiff("ExtDatFigure_8.tiff", 200, 120, "mm", res = 600, compression = "lzw", bg = "transparent")
ggarrange(plotlist = plots, ncol = 8, nrow = 4)
dev.off()

