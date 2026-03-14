
library(ape)
library(Biostrings)
library(dplyr)
library(ggplot2)
library(HTSSIP)
library(magrittr)
library(phyloseq)
library(randomcoloR)
library(regex)
library(scales)
library(stringr)
library(tidyr)
library(vegan)


setwd("C:/Users/User/Desktop/Svalbard/Svalbard_Data/Data_2023/Svalbard23_Data_qSIP/AWI2")

Ordered_fig_EAF_13C <- read.table("Ordered_OTU_greater_than_5_13C.csv", header = TRUE, sep=";",dec = ",", fill = TRUE, check.names = FALSE)
Ordered_fig_EAF_CN <- read.table("Ordered_OTU_greater_than_5_CN.csv", header = TRUE, sep=";",dec = ",", fill = TRUE, check.names = FALSE)
Ordered_fig_EAF_15N <- read.table("Ordered_OTU_greater_than_5_15N.csv", header = TRUE, sep=";",dec = ",", fill = TRUE, check.names = FALSE)

unique_phylum <- unique(Ordered_fig_EAF_CN$Phylum)
unique_class <- unique(Ordered_fig_EAF_CN$Class)

class_names <- c("Pezizomycetes", "Leotiomycetes", "unidentified", "Dothideomycetes",
                 "Basidiobolomycetes", "Microbotryomycetes", "Saccharomycetes",
                 "Eurotiomycetes", "Rhizophydiomycetes", "Agaricomycetes",
                 "Orbiliomycetes", "Sordariomycetes", "Malasseziomycetes")

colors_13 <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
               "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
               "#aec7e8", "#ffbb78", "#98df8a")

master_taxa_colors <- setNames(colors_13, class_names)


### qSIP plot

### 13C + 15N

Ordered_fig_EAF_CN <- read.table("Ordered_OTU_greater_than_5_CN.csv", header = TRUE, sep=";",dec = ",", fill = TRUE, check.names = FALSE)

unique_phylum <- unique(Ordered_fig_EAF_CN$Phylum)
unique_class <- unique(Ordered_fig_EAF_CN$Class)

qSIP.plot.maker <- function(DATA, global_y_min, global_y_max) {
  DATA$Class <- factor(DATA$Class)
  phyla <- levels(DATA$Class)
  taxa.id <- data.frame(taxon = DATA$OTU, code = DATA$Class)
  ranks <- order(DATA$A, decreasing = FALSE)
  tax.order <- data.frame(axis.loc = seq(1, dim(DATA)[1]), ranks = ranks)
  
  plot(y = 1:dim(DATA)[1], x = DATA$A[tax.order$ranks], type = "n", bty = "l",
       xlab = "", ylab = "", xaxt = "n", yaxt = "n", xlim = c(global_y_min, global_y_max), main = "")
  counter <- 1
  for (p in seq_along(phyla)) {
    curr.comp.ranked <- DATA[tax.order$ranks, ]
    mids <- curr.comp.ranked$EAF_13C_15N[as.character(curr.comp.ranked$OTU) %in% as.character(taxa.id$taxon[taxa.id$code == phyla[p]])]
    lowers <- curr.comp.ranked$EAF_CI_low_13C_15N[as.character(curr.comp.ranked$OTU) %in% as.character(taxa.id$taxon[taxa.id$code == phyla[p]])]
    uppers <- curr.comp.ranked$EAF_CI_high_13C_15N[as.character(curr.comp.ranked$OTU) %in% as.character(taxa.id$taxon[taxa.id$code == phyla[p]])]
    tax.nums <- counter:(counter + length(mids) - 1)
    counter <- counter + length(mids)
    arrows(x0 = lowers, y0 = tax.nums, x1 = uppers, y1 = tax.nums, length = 0, angle = 90, code = 3, col = as.character(phyla.cols[p]), lwd = 1)
    points(x = mids, y = tax.nums, pch = 21, cex = 0.8, col = NA, bg = as.character(phyla.cols[p]))
    par(xpd = NA)
    text(x = global_y_max, y = mean(tax.nums), labels = phyla[p], adj = c(0, 0.5), pos = 2.2, offset = 2.0, cex = 1.1, col = as.character(phyla.cols[p]))
    par(xpd = FALSE)
  }
  abline(v = 0, col = "black", lty = 2)
  par(mgp = c(3, 0, 0))
  axis(side = 2, at = 1:dim(DATA)[1], labels = NA, tck = -0.015, las = 1, cex.axis = 0.6)
  mtext("Ranked OTUs", side = 2, line = 0.4, cex = 0.75)
  par(mgp = c(3, 0.45, 0))
  axis(side = 1, tck = -0.015, las = 1, cex.axis = 0.6)
  mtext(expression(paste("Excess atom fraction "^13, "C", sep = "")), side = 1, line = 1.7, at = NA, cex = 0.75)
}

svg(file="Fungi_AWI2_CN.svg", width=10, height=10)
qSIP.plot.maker(DATA=Ordered_fig_EAF_CN, -0.6, 1)
dev.off()


### 15N

Ordered_fig_EAF_15N <- read.table("Ordered_OTU_greater_than_5_15N.csv", header = TRUE, sep=";",dec = ",", fill = TRUE, check.names = FALSE)

unique_phylum <- unique(Ordered_fig_EAF_15N$Phylum)
unique_class <- unique(Ordered_fig_EAF_15N$Class)

qSIP.line.plot.maker <- function(DATA, global_y_min, global_y_max) {
  DATA$Class <- factor(DATA$Class)
  phyla <- levels(DATA$Class)
  phyla.cols <- master_phyla_colors[phyla]
  taxa.id <- data.frame(taxon = DATA$OTU, code = DATA$Class)
  ranks <- order(DATA$A, decreasing = FALSE)
  tax.order <- data.frame(axis.loc = seq(1, dim(DATA)[1]), ranks = ranks)
  
  plot(y = 1:dim(DATA)[1], x = DATA$A[tax.order$ranks], type = "n", bty = "l",
       xlab = "", ylab = "", xaxt = "n", yaxt = "n", xlim = c(global_y_min, global_y_max), main = "")
  counter <- 1
  for (p in seq_along(phyla)) {
    curr.comp.ranked <- DATA[tax.order$ranks, ]
    mids <- curr.comp.ranked$EAF_15N[as.character(curr.comp.ranked$OTU) %in% as.character(taxa.id$taxon[taxa.id$code == phyla[p]])]
    tax.nums <- counter:(counter + length(mids) - 1)
    if (length(mids) > 1) {
      lines(x = mids, y = tax.nums, col = as.character(phyla.cols[p]), lwd = 2)
    }
    points(x = mids, y = tax.nums, pch = 16, cex = 0.2, col = as.character(phyla.cols[p]))
    par(xpd = NA)
    text(x = global_y_max, y = mean(tax.nums), labels = phyla[p], adj = c(0, 0.5), pos = 2.2, offset = 2.0, cex = 1.1, col = as.character(phyla.cols[p]))
    par(xpd = FALSE)
    counter <- counter + length(mids)
  }
  abline(v = 0, col = "black", lty = 2)
  par(mgp = c(3, 0, 0))
  axis(side = 2, at = 1:dim(DATA)[1], labels = NA, tck = -0.015, las = 1, cex.axis = 0.6)
  mtext("Ranked OTUs", side = 2, line = 0.4, cex = 0.75)
  par(mgp = c(3, 0.45, 0))
  axis(side = 1, tck = -0.015, las = 1, cex.axis = 0.6)
  mtext(expression(paste("Excess atom fraction "^15, "N", sep = "")), side = 1, line = 1.7, at = NA, cex = 0.75)
}

svg(file="Fungi_AWI2_15N.svg", width=10, height=10)
qSIP.line.plot.maker(DATA=Ordered_fig_EAF_15N, -0.6, 1)
dev.off()
