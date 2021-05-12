require(plyr)
require(dplyr)
require(lubridate)
require(tsibble)
require(DescTools)
require(splines)
require(survival)

# Function to calculate log with zeros
log_zeros <- function(x) {
  sign.val <- sign(x)
  min.val <- min(abs(x[x!=0]))
  log.x <- log(ifelse(x==0, min.val/100, abs(x)))
  return(log.x)
}

# Calculate date of last contact
lastvis <- function (dat.id) {
  dat.id$date_last <- rep(max(c(max(dat.id$date_visit), max(dat.id$date_death)), na.rm = T), nrow(dat.id))
  dat.id$fu_length <- as.numeric(dat.id$date_last[1] - min(dat.id$date_visit)) / 365.25
  return (dat.id)
}

# Function to calculate CVs by year from last contact
cv.year <- function(dat.id, vars = selected.biomarkersHD, cv = T, min.vis.nb = 7,
                    covars = c("id_no", "sex", "date_death", "date_birth", "date_last", "diabetes", "fu_length")) {
  
  dat.id$years <- floor((as.numeric(dat.id$date_last - dat.id$date_visit) / 365))
  
  # remove years with insufficient visit number, otherwise CVPC1 is biased
  counts <- table(dat.id$years)
  dat.id2 <- dat.id[dat.id$years %in% names(counts[which(counts >= min.vis.nb)]), ]
  
  if (nrow(dat.id2) > 0) {
    # calculate CVs
    ID_by_yr <- split(dat.id2, dat.id2$years)
    ID_cv <- as.data.frame(matrix(NA, length(ID_by_yr), length(vars) + length(covars) + 4))
    colnames(ID_cv) <- c(covars, "age_visit", "years", "nb_vis", "vis_sd", vars)
    ID_cv[ , covars] <- rdply(length(ID_by_yr), ID_by_yr[[1]][1, covars], .id = NULL)
    ID_cv[ , "age_visit"] <- unlist(lapply(ID_by_yr, function(x) min(x$age_visit)))
    ID_cv[ , "years"] <- as.numeric(names(ID_by_yr))
    ID_cv[ , "nb_vis"] <- unlist(lapply(ID_by_yr, nrow))
    ID_cv[ , "vis_sd"] <- unlist(lapply(ID_by_yr, function (z) {
      z <- z[order(z$date_visit), ]
      if (nrow(z)<3) {
        as.numeric(z$date_visit[2] - z$date_visit[1]) / 7
      } else {
        sd(unlist(Map(function(x, y) as.numeric(y - x) / 7, z$date_visit[-length(z$date_visit)], z$date_visit[-1])))
      }
      
    }))
    
    # Calculate CVs (or means)
    if (cv==T) {
      cvs <- lapply(ID_by_yr, function(x) { 
        apply(x[ , vars], 2, function(z) {
          cv <- as.numeric(as.character(sd(z))) / as.numeric(as.character(mean(z)))
          return(ifelse(cv=="NaN", 0, cv))
        }) })
    } else {
      cvs <- lapply(ID_by_yr, function(x) { 
        apply(x[ , vars], 2, function(z) as.numeric(as.character(mean(z)))) })
    }
    ID_cv[ , vars] <- t(bind_cols(cvs))
    ID_by_yr2 <- as.data.frame(ID_cv)
    
    
    return(ID_by_yr2)
  }
}


# Function to calculate CVs every 4 months (or other time window)
cv.4mth <- function(dat.id, selected.biomarkersHD, cv = T, nb.month = 4, min.vis.nb = 7,
                    covars = c("id_no", "sex", "date_death", "date_birth", "date_last", "diabetes", "fu_length")) {
  
  dat.id$years <- floor(as.numeric(dat.id$date_last - dat.id$date_visit) / (365 / (12/nb.month)))
  
  # remove years with insufficient visit number, otherwise CVPC1 is biased
  counts <- table(dat.id$years)
  dat.id2 <- dat.id[dat.id$years %in% names(counts[which(counts >= min.vis.nb)]), ]
  
  if (nrow(dat.id2) > 0) {
    # calculate CVs
    ID_by_yr <- split(dat.id2, dat.id2$years)
    ID_cv <- as.data.frame(matrix(NA, length(ID_by_yr), length(selected.biomarkersHD) + length(covars) + 4))
    colnames(ID_cv) <- c(covars, "age_visit", "years", "nb_vis", "vis_sd", selected.biomarkersHD)
    ID_cv[ , covars] <- rdply(length(ID_by_yr), ID_by_yr[[1]][1, covars], .id = NULL)
    ID_cv[ , "age_visit"] <- unlist(lapply(ID_by_yr, function(x) min(x$age_visit)))
    ID_cv[ , "years"] <- as.numeric(names(ID_by_yr))
    ID_cv[ , "nb_vis"] <- unlist(lapply(ID_by_yr, nrow))
    ID_cv[ , "vis_sd"] <- unlist(lapply(ID_by_yr, function (z) {
      z <- z[order(z$date_visit), ]
      if (nrow(z)<3) {
        as.numeric(z$date_visit[2] - z$date_visit[1]) / 7
      } else {
        sd(unlist(Map(function(x, y) as.numeric(y - x) / 7, z$date_visit[-length(z$date_visit)], z$date_visit[-1])))
      }
      
    }))
    
    # Calculate CVs (or means)
    if (cv==T) {
      cvs <- lapply(ID_by_yr, function(x) { 
        apply(x[ , selected.biomarkersHD], 2, function(z) {
          cv <- as.numeric(as.character(sd(z))) / as.numeric(as.character(mean(z)))
          return(ifelse(cv=="NaN", 0, cv))
        }) })
    } else {
      cvs <- lapply(ID_by_yr, function(x) { 
        apply(x[ , selected.biomarkersHD], 2, function(z) as.numeric(as.character(mean(z)))) })
    }
    ID_cv[ , selected.biomarkersHD] <- t(bind_cols(cvs))
    ID_by_yr2 <- as.data.frame(ID_cv)
    
    
    return(ID_by_yr2)
  }
}



# Deal with outliers due to uniform values (mainly in sodium)
uniform_val <- function(var) {
  outlier <- var[which(var %in% boxplot(var, range = 6, plot = F)$out==T)]
  if (length(outlier) > 0) {
    var[var %in% outlier] <- min(var[!var %in% outlier]) * 1.1
  }
  return(var)
}


# Change sign of PC
pc_sign <- function(data, pc.name = "PC1", pca) {
  datt <- data
  colnames(datt)[which(colnames(datt)==pc.name)] <- "predictor"
  mod <- coxph(Surv(years, status) ~ predictor + bs(age_visit, df = 5) + sex +  
                 diabetes + fu_length, data = datt, weights = sqrt(nb_vis),
               cluster = id_no)
  hr <- as.numeric(summary(mod)$conf.int[1, 1])
  if (hr < 1) {
    PC <- -1 * datt$predictor
    datt.2 <- cbind(datt[ , -which(colnames(datt)=="predictor")], PC)
    colnames(datt.2)[which(colnames(datt.2)=="PC")] <- pc.name
    return(datt.2)
  } else return(data)
}


# Change sign of PC (in specific subsets of population)
pc_sign2 <- function(data, pc.name = "PC1", pca) {
  datt <- data
  colnames(datt)[which(colnames(datt)==pc.name)] <- "predictor"
  mod <- coxph(Surv(years, status) ~ predictor + bs(age_visit, df = 3)+ fu_length, 
               data = datt)
  hr <- as.numeric(summary(mod)$conf.int[1, 1])
  if (hr < 1) {
    PC <- -1 * datt$predictor
    datt.2 <- cbind(datt[ , -which(colnames(datt)=="predictor")], PC)
    colnames(datt.2)[which(colnames(datt.2)=="PC")] <- pc.name
    return(datt.2)
  } else return(data)
}

round_p_values <- function(p_value) {
  xx <- round(p_value, 2)
  yy <- unlist(lapply(xx, function(x) {
    if (x==0) {
      x <- round(x, 3)
      if (x==0) x <- "<0.001"
    }
    return(x)
  }))
  zz <- unlist(lapply(yy, function(x) ifelse(nchar(x) < 4, paste(x, "0", sep = ""), x)))
  return(yy)
}



# Function to select a time frame before death
pca_interval <- function (dat.id, time.interval) {
  death <- ymd(unique(dat.id$date_death))
  
  if (!is.na(death)) {
    min.time <- death %m-% months(time.interval[1])
    max.time <- death %m-% months(time.interval[2])
    dat.int <- dat.id[dat.id$date_visit<=min.time & dat.id$date_visit>=max.time, ]
    return(dat.int)
  }
}


# Function to select a time frame before death with defined time 
find_interval_v2 <- function (dat.id, time, increase = 1, selected.biomarkersHD) {
  death <- ymd(unique(dat.id$date_death))
  
  if (!is.na(death)) {
    tis <- list()
    for (ti in 1:(length(time)-1)) {
      time.interval <- time[ti]
      
      # Define T1 and T2 (using ± 1 week)
      time1 <- range(ymd(death %m-% months(time.interval)) - weeks(1), ymd(death %m-% months(time.interval)) + weeks(1))  
      time2 <- range(ymd(death %m-% months(time.interval + increase)) - weeks(1), 
                     ymd(death %m-% months(time.interval + increase)) + weeks(1))
      
      # Select visits corresponding to T1 and T2
      dat.T1 <- dat.id[dat.id$date_visit >= time1[1] & dat.id$date_visit <= time1[2], ]
      colnames(dat.T1)[which(colnames(dat.T1) %in% selected.biomarkersHD)] <- 
        unlist(lapply(colnames(dat.T1)[which(colnames(dat.T1) %in% selected.biomarkersHD)], 
                      function (x) paste(x, "T1", sep = "_")))
      if (nrow(dat.T1) > 1) {
        dat.T1 <- dat.T1[which.min(ymd(death %m-% months(time.interval)) - dat.T1$date_visit), ]
      }
      
      dat.T2 <- dat.id[dat.id$date_visit >= time2[1] & dat.id$date_visit <= time2[2], ]
      colnames(dat.T2)[which(colnames(dat.T2) %in% selected.biomarkersHD)] <- 
        unlist(lapply(colnames(dat.T2)[which(colnames(dat.T2) %in% selected.biomarkersHD)], 
                      function (x) paste(x, "T2", sep = "_")))
      if (nrow(dat.T2) > 1) {
        dat.T2 <- dat.T2[which.min(abs(ymd(death %m-% months(time.interval + increase)) - dat.T2$date_visit)), ]
      }
      
      # Re-organize into one dataset
      if (nrow(dat.T1) > 0 & nrow(dat.T2) > 0) {
        tis[[ti]] <- cbind(dat.T1[ , colnames(dat.T1)[grep("T1", colnames(dat.T1))]], 
                           dat.T2[ , colnames(dat.T2)[grep("T2", colnames(dat.T2))]])
      }
    }
    dat.id2 <- bind_rows(tis)
    return(dat.id2)
  }
}


find_interval_v3 <- function (dat.id, time, increase = 1, selected.biomarkersHD) {
  
  # Determine elapsed time between last HD visit and death 
  death <- ymd(unique(dat.id$date_death))
  lastvis <- max(dat.id$date_visit)
  mth.diff <- round(as.numeric(death - lastvis) / 365 * 12)
  #week.diff <- as.numeric(death - lastvis)/7
  
  if (!is.na(lastvis)) {
    tis <- list()
    for (ti in (mth.diff + 1):(length(time)-1)) {
      time.interval <- time[ti]
      
      # Define T1 and T2 (using ± 1 week)
      time1 <- range(ymd(lastvis %m-% months(time.interval)) - weeks(1), ymd(lastvis %m-% months(time.interval)) + weeks(1))  
      time2 <- range(ymd(lastvis %m-% months(time.interval + increase)) - weeks(1), 
                     ymd(lastvis %m-% months(time.interval + increase)) + weeks(1))
      
      # Select visits corresponding to T1 and T2
      dat.T1 <- dat.id[dat.id$date_visit >= time1[1] & dat.id$date_visit <= time1[2], ]
      colnames(dat.T1)[which(colnames(dat.T1) %in% selected.biomarkersHD)] <- 
        unlist(lapply(colnames(dat.T1)[which(colnames(dat.T1) %in% selected.biomarkersHD)], 
                      function (x) paste(x, "T1", sep = "_")))
      if (nrow(dat.T1) > 1) {
        dat.T1 <- dat.T1[which.min(ymd(death %m-% months(time.interval)) - dat.T1$date_visit), ]
      }
      
      dat.T2 <- dat.id[dat.id$date_visit >= time2[1] & dat.id$date_visit <= time2[2], ]
      colnames(dat.T2)[which(colnames(dat.T2) %in% selected.biomarkersHD)] <- 
        unlist(lapply(colnames(dat.T2)[which(colnames(dat.T2) %in% selected.biomarkersHD)], 
                      function (x) paste(x, "T2", sep = "_")))
      if (nrow(dat.T2) > 1) {
        dat.T2 <- dat.T2[which.min(abs(ymd(death %m-% months(time.interval + increase)) - dat.T2$date_visit)), ]
      }
      
      # Re-organize into one dataset
      if (nrow(dat.T1) > 0 & nrow(dat.T2) > 0) {
        tis[[ti]] <- cbind(dat.T1[ , colnames(dat.T1)[grep("T1", colnames(dat.T1))]], 
                           dat.T2[ , colnames(dat.T2)[grep("T2", colnames(dat.T2))]])
      }
    }
    dat.id2 <- bind_rows(tis)
    return(dat.id2)
  }
}




# Add "status" variable
status_fct <- function(dat.id) {
  dat.id$status <- ifelse(!is.na(dat.id$date_death) & dat.id$years==min(dat.id$years), 1, 0)
  return(dat.id)
}



# Exclude dead individuals whom we lost during follow-up
no.fu <- function(dat.id, time = 3) {
  if (!is.na(dat.id$date_death[1])) {
    max.vis <- max(dat.id$date_visit)
    death <- dat.id$date_death[1]
    if (as.numeric(death - max.vis) <= (time * 365/12)) {
      dat.id$fu.until.death <- 1
    } else {
      dat.id$fu.until.death <- 0
    }
  } else {
    dat.id$fu.until.death <- 0
  }
  return(dat.id)
}


############## For ggcorrplot
ggcorrplot_custom <- function(corr, corr_size = c(4, 10),
                              method = c("square", "circle"),
                              type = c("full", "lower", "upper"),
                              ggtheme = ggplot2::theme_minimal,
                              title = "",
                              show.legend = TRUE,
                              legend.title = "Corr",
                              show.diag = FALSE,
                              colors = c("blue", "white", "red"),
                              outline.color = "gray",
                              hc.order = FALSE,
                              hc.method = "complete",
                              lab = FALSE,
                              lab_col = "black",
                              lab_size = 4,
                              p.mat = NULL,
                              sig.level = 0.05,
                              insig = c("pch", "blank"),
                              pch = 4,
                              pch.col = "black",
                              pch.cex = 5,
                              tl.cex = 12,
                              tl.col = "black",
                              tl.srt = 45,
                              digits = 2) {
  type <- match.arg(type)
  method <- match.arg(method)
  insig <- match.arg(insig)
  
  if(inherits(corr, "cor_mat")){
    # cor_mat object from rstatix
    cor.mat <- corr
    corr <- .tibble_to_matrix(cor.mat)
    p.mat <- .tibble_to_matrix(attr(cor.mat, "pvalue"))
  }
  
  if (!is.matrix(corr) & !is.data.frame(corr)) {
    stop("Need a matrix or data frame!")
  }
  corr <- as.matrix(corr)
  
  corr <- base::round(x = corr, digits = digits)
  
  if (hc.order) {
    ord <- .hc_cormat_order(corr)
    corr <- corr[ord, ord]
    if (!is.null(p.mat)) {
      p.mat <- p.mat[ord, ord]
      p.mat <- base::round(x = p.mat, digits = digits)
    }
  }
  
  # Get lower or upper triangle
  if (type == "lower") {
    corr <- .get_lower_tri(corr, show.diag)
    p.mat <- .get_lower_tri(p.mat, show.diag)
  }
  else if (type == "upper") {
    corr <- .get_upper_tri(corr, show.diag)
    p.mat <- .get_upper_tri(p.mat, show.diag)
  }
  
  # Melt corr and pmat
  corr <- reshape2::melt(corr, na.rm = TRUE)
  colnames(corr) <- c("Var1", "Var2", "value")
  corr$pvalue <- rep(NA, nrow(corr))
  corr$signif <- rep(NA, nrow(corr))
  
  if (!is.null(p.mat)) {
    p.mat <- reshape2::melt(p.mat, na.rm = TRUE)
    corr$coef <- corr$value
    corr$pvalue <- p.mat$value
    corr$signif <- as.numeric(p.mat$value <= sig.level)
    p.mat <- subset(p.mat, p.mat$value > sig.level)
    if (insig == "blank") {
      corr$value <- corr$value * corr$signif
    }
  }
  
  
  corr$abs_corr <- abs(corr$value) * 10
  
  # heatmap
  p <-
    ggplot2::ggplot(
      data = corr,
      mapping = ggplot2::aes_string(x = "Var1", y = "Var2", fill = "value")
    )
  
  # modification based on method
  if (method == "square") {
    p <- p +
      ggplot2::geom_tile(color = outline.color)
  } else if (method == "circle") {
    p <- p +
      ggplot2::geom_point(
        color = outline.color,
        shape = 21,
        ggplot2::aes_string(size = "abs_corr")
      ) +
      ggplot2::scale_size(range = corr_size) +
      ggplot2::guides(size = FALSE)
  } 
  
  # adding colors
  p <-
    p + ggplot2::scale_fill_gradient2(
      low = colors[1],
      high = colors[3],
      mid = colors[2],
      midpoint = 0,
      limit = c(-1, 1),
      space = "Lab",
      name = legend.title
    )
  
  # depending on the class of the object, add the specified theme
  if (class(ggtheme)[[1]] == "function") {
    p <- p + ggtheme()
  } else if (class(ggtheme)[[1]] == "theme") {
    p <- p + ggtheme
  }
  
  
  p <- p +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = tl.srt,
        vjust = 1,
        size = tl.cex,
        hjust = 1
      ),
      axis.text.y = ggplot2::element_text(size = tl.cex)
    ) +
    ggplot2::coord_fixed()
  
  label <- round(x = corr[, "value"], digits = digits)
  if(!is.null(p.mat) & insig == "blank"){
    ns <- corr$pvalue > sig.level
    if(sum(ns) > 0) label[ns] <- " "
  }
  
  # matrix cell labels
  if (lab) {
    p <- p +
      ggplot2::geom_text(
        mapping = ggplot2::aes_string(x = "Var1", y = "Var2"),
        label = label,
        color = lab_col,
        size = lab_size
      )
  }
  
  # matrix cell glyphs
  if (!is.null(p.mat) & insig == "pch") {
    p <- p + ggplot2::geom_point(
      data = p.mat,
      mapping = ggplot2::aes_string(x = "Var1", y = "Var2"),
      shape = pch,
      size = pch.cex,
      color = pch.col
    )
  }
  
  # add titles
  if (title != "") {
    p <- p +
      ggplot2::ggtitle(title)
  }
  
  # removing legend
  if (!show.legend) {
    p <- p +
      ggplot2::theme(legend.position = "none")
  }
  
  # removing panel
  p <- p +
    .no_panel()
  p
}



#' Compute the matrix of correlation p-values
#'
#' @param x numeric matrix or data frame
#' @param ... other arguments to be passed to the function cor.test.
#' @rdname ggcorrplot
#' @export

cor_pmat <- function(x, ...) {
  
  # initializing values
  mat <- as.matrix(x)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  
  # creating the p-value matrix
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- stats::cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  
  # name rows and columns of the p-value matrix
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  
  # return the final matrix
  p.mat
}



#+++++++++++++++++++++++
# Helper Functions
#+++++++++++++++++++++++

# Get lower triangle of the correlation matrix
.get_lower_tri <- function(cormat, show.diag = FALSE) {
  if (is.null(cormat)) {
    return(cormat)
  }
  cormat[upper.tri(cormat)] <- NA
  if (!show.diag) {
    diag(cormat) <- NA
  }
  return(cormat)
}

# Get upper triangle of the correlation matrix
.get_upper_tri <- function(cormat, show.diag = FALSE) {
  if (is.null(cormat)) {
    return(cormat)
  }
  cormat[lower.tri(cormat)] <- NA
  if (!show.diag) {
    diag(cormat) <- NA
  }
  return(cormat)
}

# hc.order correlation matrix
.hc_cormat_order <- function(cormat, hc.method = "complete") {
  dd <- stats::as.dist((1 - cormat) / 2)
  hc <- stats::hclust(dd, method = hc.method)
  hc$order
}

.no_panel <- function() {
  ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank()
  )
}


# Convert a tbl to matrix
.tibble_to_matrix <- function(x){
  x <-  as.data.frame(x)
  rownames(x) <- x[, 1]
  x <- x[, -1]
  as.matrix(x)
}
