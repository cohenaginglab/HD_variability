library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(grid)
library(ggcorrplot)
library(Hmisc)


rm(list = ls())

source("Functions_HD_variability_paper.R")

### Load data
clean.HD$date_visit <- as.Date(clean.HD$date_visit)
clean.HD$date_death <- as.Date(clean.HD$date_death)

## variables
selected.biomarkersHD <- c("wbc", "hb", "hct", "mch", "mchc", "mcv", "pltlt", "rdw", 
                           "k", "sodium", "rbc")
vars.nms <- c("WBC", "Hemoglobin", "Hematocrit", "MCH", "MCHC", "MCV", "Platelets", 
              "RDW", "Potassium", "Sodium", "RBC")

# Transform variables not normally distributed
if ("pltlt" %in% selected.biomarkersHD) {
  clean.HD$pltlt <- sqrt(clean.HD$pltlt)
}
clean.HD[ , which(colnames(clean.HD) %in% c("wbc", "rdw", "Glucose"))] <- 
  log(clean.HD[ , which(colnames(clean.HD) %in% c("wbc", "rdw", "Glucose"))])

# Add date of last contact
datt <- ddply(clean.HD, .(id_no), lastvis)


##### Panels A and B: PCA at different time intervals before death
# Select time intervals
time.ints <- list(c(0, 6), c(12, 24), c(36, 60), c(84, 120))         

# Perform PCA (on raw biomarkers and CVs per year)
pcas <- list("raw" = list(), "CVs" = list())
for (i in 1:length(pcas)) {
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
    
    if (i==2) {
      # Calculate CVs per year 
      dat <- ddply(dat, .(id_no), cv.year, vars = selected.biomarkersHD, cv = T)
      # Apply log transformation on CVs
      dat[ , selected.biomarkersHD] <- apply(dat[ , selected.biomarkersHD], 2, log_zeros)
      # Some variables have no variation, leading to outliers
      dat[ , selected.biomarkersHD] <- apply(dat[ , selected.biomarkersHD], 2, 
                                                 uniform_val)
      
      # Correct CVs for the number of visits included
      dat[ , selected.biomarkersHD] <- apply(dat[ , selected.biomarkersHD], 2, function(x) {
        residuals(nls(x ~ I(1 / sqrt(nb_vis) * a) + b, data = dat, start = list(a = 1, b = 1)))
      })
    }
    
    print(nrow(dat))
    print(length(unique(dat$id_no)))
    
    # perform PCA
    pcas[[i]][[ti]] <- prcomp(dat[ , selected.biomarkersHD], center = T, scale. = T)
    print(paste(c("raw", "CVs")[i], " - ", c("All", "0-6 months", "1-2 years", "3-5 years", "7-10 years", "Not diabetic", 
                                             "Diabetic", "Men", "Women", "< 60 yrs", "60-75 yrs", "75+ yrs")[ti]), sep = "")
    print(summary(pcas[[i]][[ti]]))
  }
}

######## Prepare panels #########
#### 1. Panel "A" - Plot of explained variance by PC 
# Combine proportion of variance by PC
variances <- list()
for (i in 1:length(pcas)) {
  for (ti in 1:length(pcas[[1]])) {
    xx <- as.data.frame(matrix(NA, length(selected.biomarkersHD), 4))
    colnames(xx) <- c("PC", "Variance", "Variable.list", "Time")
    xx[ , "PC"] <- 1:length(selected.biomarkersHD)
    xx[ , "Variance"] <- summary(pcas[[i]][[ti]])$importance[2, ]
    xx[ , "Variable.list"] <- names(pcas)[i]
    xx[ , "Time"] <- ti
    variances[[ifelse(i==1, ti, length(pcas[[1]]) + ti)]] <- xx
  }
}
pc.var <- bind_rows(variances)
pc.var$Time <- factor(pc.var$Time, levels = 1:(length(time.ints) + 8), 
                      labels = c("All", "0-6 months", "1-2 years", "3-5 years", "7-10 years", "Not diabetic", "Diabetic", 
                                 "Men", "Women", "< 60 yrs", "60-75 yrs", "75+ yrs"))
plot_a <- ggplot(pc.var, aes(x = PC, y = Variance, color = Time, shape = Variable.list)) + 
  geom_point(aes(x = PC, y = Variance), size = 2) +
  geom_line(aes(x = PC, y = Variance)) +
  scale_color_brewer(palette = "Paired") +
  labs(y = "Proportion of variance", shape = "PCA", color = "Subsets") + 
  scale_x_continuous(breaks = seq(2, length(selected.biomarkersHD), 2)) +
  scale_y_continuous(limits = c(-0.01, 0.55)) +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 10, face = "bold"), 
        axis.title.y = element_text(size = 10, face = "bold"), 
        axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), 
        legend.position = c(0.74, 0.6), legend.direction = "vertical", 
        legend.text = element_text(size = 8), legend.box = "horizontal",
        legend.title =  element_text(size = 9, face = "bold"), 
        legend.key.height =  unit(3.5, "mm"), legend.key.width =  unit(4, "mm"), 
        legend.margin = margin(1, 1, 0.5, 1), legend.key = element_rect(fill = NA, color = NA), 
        plot.margin = unit(c(0.25, 0.5, 0.5, 0.25), "lines"), 
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        axis.line = element_line(colour = "black")) +
  guides(color = guide_legend(order = 0), shape = guide_legend(order = 1))

#### 2. Panel "b" (barplot)
# Get variable contributions
contribs <- list()
for (ti in 1:length(pcas$CVs)) {
  xx <- as.data.frame(matrix(NA, length(selected.biomarkersHD), 3))
  colnames(xx) <- c("Variable", "Contribution", "Time")
  xx[ , "Variable"] <- vars.nms
  xx[ , "Contribution"] <- abs(pcas$CVs[[ti]]$rotation[ , 1]) / sum(abs(pcas$CVs[[ti]]$rotation[ , 1]))
  xx[ , "Time"] <- ti
  contribs[[ti]] <- xx
}
pc.contrib <- bind_rows(contribs)
pc.contrib$Time <- factor(pc.contrib$Time, levels = 1:(length(time.ints) + 8), 
                          labels = c("All", "0-6 months", "1-2 years", "3-5 years", "7-10 years", "Not diabetic", "Diabetic", 
                                     "Men", "Women", "< 60 yrs", "60-75 yrs", "75+ yrs"))

# Order variables based on PCA "All" loadings (abs)
pca.all <- pc.contrib[pc.contrib$Time=="All", ]
order.all <- pca.all[order(pca.all$Contribution, decreasing = T), "Variable"]
pc.contrib$Variable <- factor(pc.contrib$Variable, levels = rev(order.all))

myPal <- rev(brewer.pal(length(selected.biomarkersHD), "Paired"))

plot_b <- ggplot(pc.contrib, aes(x = Time, y = Contribution)) +
  geom_col(aes(fill = Variable), width = 0.7) + 
  scale_fill_manual(values = myPal, na.value = "white") +
  scale_x_discrete(labels = c("All", "0-0.5", "1-2", "3-5", "7-10", "No", "Yes", "M", "W", "<60", "60-75", "75+")) +
  labs(x = "", y = "Variable contribution") +
  theme(axis.title.y = element_text(size = 10, face = "bold"), axis.text.x = element_text(size = 7), 
        axis.text.y = element_text(size = 7), axis.title.x = element_text(size = 10, face = "bold", margin = margin(b = 5)),
        legend.text = element_text(size = 8), legend.title =  element_text(size = 9, face = "bold"), 
        legend.position = "right", legend.margin = margin(0, 0, 0, 0.05), legend.key.size =  unit(4, "mm"),
        plot.margin = unit(c(0.7, 0.2, 1, 0.1), "lines"), panel.background = element_rect(fill = "white"),
        panel.border = element_blank(), axis.line = element_line(size = 0.25, linetype = 1)) +
  annotation_custom(grob = segmentsGrob(x0 = unit(0, "npc"), x1 = unit(1, "npc"), y0 = unit(0, "npc"), 
                                        y1 = unit(1, "npc"), gp = gpar(lwd = 1.5)),
                    xmin = 2, xmax = 5, ymin = -0.16, ymax = -0.16) +
  annotation_custom(grob = segmentsGrob(x0 = unit(0, "npc"), x1 = unit(1, "npc"), y0 = unit(0, "npc"), 
                                        y1 = unit(1, "npc"), gp = gpar(lwd = 1.5)),
                    xmin = 6, xmax = 7, ymin = -0.16, ymax = -0.16) +
  annotation_custom(grob = segmentsGrob(x0 = unit(0, "npc"), x1 = unit(1, "npc"), y0 = unit(0, "npc"), 
                                        y1 = unit(1, "npc"), gp = gpar(lwd = 1.5)),
                    xmin = 8, xmax = 9, ymin = -0.16, ymax = -0.16) +
  annotation_custom(grob = segmentsGrob(x0 = unit(0, "npc"), x1 = unit(1, "npc"), y0 = unit(0, "npc"), 
                                        y1 = unit(1, "npc"), gp = gpar(lwd = 1.5)),
                    xmin = 10, xmax = 12, ymin = -0.16, ymax = -0.16) +
  annotation_custom(grob = textGrob(label = "Years", hjust = 0.5, gp = gpar(cex = 0.7, font = 2)),
                    xmin = 3.5, xmax = 3.5, ymin = -0.21, ymax = -0.21) +
  annotation_custom(grob = textGrob(label = "before death", hjust = 0.5, gp = gpar(cex = 0.7, font = 2)),
                    xmin = 3.5, xmax = 3.5, ymin = -0.28, ymax = -0.28) +
  annotation_custom(grob = textGrob(label = "Diabetes", hjust = 0.5, gp = gpar(cex = 0.7, font = 2)),
                    xmin = 6.5, xmax = 6.5, ymin = -0.21, ymax = -0.21) +
  annotation_custom(grob = textGrob(label = "Sex", hjust = 0.5, gp = gpar(cex = 0.7, font = 2)),
                    xmin = 8.5, xmax = 8.5, ymin = -0.21, ymax = -0.21) +
  annotation_custom(grob = textGrob(label = "Age group", hjust = 0.5, gp = gpar(cex = 0.7, font = 2)), 
                    xmin = 11, xmax = 11, ymin = -0.21, ymax = -0.21) +
  annotation_custom(grob = textGrob(label = "(years)", hjust = 0.5, gp = gpar(cex = 0.7, font = 2)), 
                    xmin = 11, xmax = 11, ymin = -0.28, ymax = -0.28)

gt <- ggplot_gtable(ggplot_build(plot_b))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
plot_b2 <- as_ggplot(gt)


######### 3. Panel "c" (corrplot) 
## Calculate CVs, means, and PCs
# Calculate mean biomarker values per year
dat.means <- ddply(datt, .(id_no), .fun = cv.year, 
                   vars = selected.biomarkersHD, cv = F, min.vis.nb = 2) 

## Calculate CVs per year 
dat.cvs <- ddply(datt, .(id_no), cv.year, vars = selected.biomarkersHD, 
                 cv = T, min.vis.nb = 2) 

# Apply log transformation on CVs
dat.cvs[ , selected.biomarkersHD] <- apply(dat.cvs[ , selected.biomarkersHD], 2, log_zeros)

# Some variables have no variation, leading to outliers
dat.cvs[ , selected.biomarkersHD] <- apply(dat.cvs[ , selected.biomarkersHD], 2, 
                                           uniform_val)
# Correct CVs for the number of visits included
dat.cvs[ , selected.biomarkersHD] <- apply(dat.cvs[ , selected.biomarkersHD], 2, function(x) {
  residuals(nls(x ~ I(1 / sqrt(nb_vis) * a) + b, data = dat.cvs, start = list(a = 1, b = 1)))
})

# Calculate PCA on CVs
pca_cv <- prcomp(dat.cvs[ , selected.biomarkersHD], center = T, scale. = T)
pca.cv <- cbind(dat.cvs, pca_cv$x) 
# Reverse CVPC1 sign 
pca.cv$PC1 <- -1 * as.matrix(scale(dat.cvs[ , selected.biomarkersHD])) %*% pca_cv$rotation[ , 1]
# add "CV" to names of PCs and biomarkers
colnames(pca.cv)[which(colnames(pca.cv) %in% selected.biomarkersHD)] <- 
  unlist(lapply(colnames(pca.cv)[which(colnames(pca.cv) %in% selected.biomarkersHD)], function (x) paste(x, "CV", sep = "_")))
colnames(pca.cv)[grep("PC", colnames(pca.cv))] <- unlist(lapply(colnames(pca.cv)[grep("PC", colnames(pca.cv))], 
                                                                function (x) paste("CV", x, sep = "")))

dat.means <- merge(dat.means, pca.cv[ , c("id_no", "years")], by = c("id_no", "years"))
pca.cv <- merge(pca.cv, dat.means[ , c("id_no", "years")], by = c("id_no", "years")) 
dat.means2 <- merge(dat.means, 
                    pca.cv[ , c("id_no", "years", "CVPC1", "CVPC2", "CVPC3")], 
                    by = c("id_no", "years"))

## Calculate correlations
# Create empty correlation matrix
corrs <- matrix(NA, length(selected.biomarkersHD) + 3, length(selected.biomarkersHD) + 3)
rownames(corrs) <- c(vars.nms, "CVPC1", "CVPC2", "CVPC3")
colnames(corrs) <- c(vars.nms, "CVPC1", "CVPC2", "CVPC3")
corps <- corrs

# Correlations with raw biomarkers
# get correlation coefficients
corr.raw <- rcorr(as.matrix(dat.means2[ , c(selected.biomarkersHD, "CVPC1", "CVPC2", "CVPC3")]), 
                  type="pearson")$r
corrs[lower.tri(corrs)] <- corr.raw[lower.tri(corr.raw)]

# get p values
corp.raw <- rcorr(as.matrix(dat.means2[ , c(selected.biomarkersHD, "CVPC1", "CVPC2", "CVPC3")]), 
                  type="pearson")$P
corps[lower.tri(corps)] <- corp.raw[lower.tri(corp.raw)]

# Correlations with CVs per year
# get correlation coefficients
corr.cv <- rcorr(as.matrix(pca.cv[ , c(grep("_CV", colnames(pca.cv), value = T), 
                                       "CVPC1", "CVPC2", "CVPC3")]), type="pearson")$r
corrs[upper.tri(corrs)] <- corr.cv[upper.tri(corr.cv)]

# get p values
corp.cv <- rcorr(as.matrix(pca.cv[ , c(grep("_CV", colnames(pca.cv), value = T), 
                                       "CVPC1", "CVPC2", "CVPC3")]), type="pearson")$P
corps[upper.tri(corps)] <- corp.cv[upper.tri(corp.cv)]

plot_c <- ggcorrplot_custom(corrs, corr_size = c(1, 5.5), method = "square", outline.col = "white", show.diag = F, 
                            colors = c("#E46726", "white", "#6D9EC1"), insig = "pch", tl.cex = 8, pch.cex = 2, p.mat = corps) +
  theme(legend.text = element_text(size = 7), legend.title =  element_text(size = 9, face = "bold"), 
        legend.position = "right", legend.margin = margin(0.1, 0.1, 0.1, 0.1), legend.key.size =  unit(4, "mm"),
        plot.margin = unit(c(0, 0, 0, 0), "lines"), panel.grid = element_blank())

#### 4.  Panel "d" (histogram of CVs)
# Plot histogram of CVs
cor.cvs <- rcorr(as.matrix(pca.cv[ , grep("_CV", colnames(pca.cv), value = T)]), type="pearson")$r
cor.cvs <- cor.cvs[upper.tri(cor.cvs)]
plot_d <- qplot(cor.cvs, geom = "histogram") +
  labs(x = "Correlation coefficient", y = "Count") + 
  theme_minimal() +
  theme(axis.title.x = element_text(size = 10, face = "bold"), 
        axis.title.y = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7), 
        plot.margin = unit(c(1, 1.5, 0.5, 1), "cm"), legend.text.align = 1, 
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), 
        axis.line = element_line(colour = "black"))
mean(cor.cvs)             
sd(cor.cvs)               
range(cor.cvs)           


####### Save Figure ########
tiff("Figure_3.tiff", 210, 150, "mm", res = 600, compression = "lzw", bg = "transparent")
ggarrange(plotlist = list(plot_a, plot_b2, plot_c, plot_d), ncol = 2, nrow = 2,  
          labels = c("A", "B", "C", "D"),  font.label = list(size = 13))    
dev.off()

