library(dplyr)
library(magrittr)
library(ggplot2)
library(HTSSIP)
library(phyloseq)
library(ggplot2)
library(scales)
library("ape")
library(vegan)
#library(regex)
library(Biostrings)
library(stringr)
library(tidyverse)

setwd("/home/worsi/Svalbard_ITS_2024/qSIP_HTSSIP/qSIP calculations/K8-40m-hot-qSIP/")
OTU_table <- read.csv("OTU_table_fjord.csv", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

# ---- SELECTION OF FRACTIONS  ----

OTU_table$Unlabeled_11a <- (apply(OTU_table[, grepl("^11afraction[0-9]+$", colnames(OTU_table))],1,sum))
OTU_table$Unlabeled_11b <- (apply(OTU_table[, grepl("^11bfraction[0-9]+$", colnames(OTU_table))],1,sum))
OTU_table$Unlabeled_11c <- (apply(OTU_table[, grepl("^11cfraction[0-9]+$", colnames(OTU_table))],1,sum))

OTU_table$Labeled_12a <- (apply(OTU_table[, grepl("^12afraction[0-9]+$", colnames(OTU_table))],1,sum))
OTU_table$Labeled_12b <- (apply(OTU_table[, grepl("^12bfraction[0-9]+$", colnames(OTU_table))],1,sum))
OTU_table$Labeled_12c <- (apply(OTU_table[, grepl("^12cfraction[0-9]+$", colnames(OTU_table))],1,sum))

# ---- REMOVE MISSING SAMPLE  ---- !! SOMETIMES YOU HAVE TO START THIS COMMAND AGAIN and RUN THE REST TO GET IT TO WORK NO IDEA WHY
OTU_table <- OTU_table[, !(colnames(OTU_table) == "")]

OTU_table <- OTU_table %>%
  filter(
    rowSums(select(., Unlabeled_11a, Unlabeled_11b, Unlabeled_11c) >= 4) >= 2,
    rowSums(select(., Labeled_12a, Labeled_12b, Labeled_12c) >= 4) >= 2
  )
  

  



# ----  Dataset creation ----
OTU_table_samples <- OTU_table %>%
  select(matches("^#OTU ID$|^1[12][a-c]fraction[0-9]+$"))
colnames(OTU_table_samples)[1] <- "OTU"
rownames(OTU_table_samples) <- OTU_table_samples[,1]
OTU_table_samples[,1] <- NULL
OTU_table_samples <- as.matrix.data.frame(OTU_table_samples)
OTU_table_phyloseq_dataset = otu_table(OTU_table_samples, taxa_are_rows = TRUE)
write.csv("Dataset.csv", x=OTU_table_samples)

##Glucose taxa information
OTU_table_taxa <- select(OTU_table, matches("OTU|Domain|Phylum|Class|Subclass|Order|Family|Genus|Species"))
colnames(OTU_table_taxa)[1] <- "OTU"
rownames(OTU_table_taxa) <- OTU_table_taxa[,1]
OTU_table_taxa[,1] <- NULL

## Change the taxa information 
t1 <- sapply(OTU_table_taxa, grepl, pattern="^uncultured|^Ambiguous|^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)$")
t2 <- apply(t1, 1, any)
t3 <- apply(t1, 1, which.max)
for(i in seq_len(nrow(OTU_table_taxa))) {
  if(t2[i]) {OTU_table_taxa[i, t3[i]:ncol(OTU_table_taxa)]  <- paste("uncultured", OTU_table_taxa[i, t3[i]-1])}
}

## Create Phyloseq Taxa table object 
OTU_table_taxa <- as.matrix.data.frame(OTU_table_taxa)
TAX_OTU = tax_table(OTU_table_taxa)


### Glucose qPCR and sample data###
con <- file("qPCR_data_fjord.csv", encoding = "UTF-7")
qPCRdata_OTU <- data.frame(read.csv(con, header = TRUE, sep = ",",check.names = FALSE,  row.names = 1, as.is=TRUE))
con1 <- file("qSIP_Sampledata_fjord.csv", encoding = "UTF-7")
sampledata_OTU <- data.frame(read.csv(con1, header = TRUE, sep = ",",check.names = FALSE,  row.names = 1, as.is=TRUE))


sampledata_OTU = sample_data(data.frame(sampledata_OTU))
physeq_analysis_OTU = phyloseq(OTU_table_phyloseq_dataset, TAX_OTU, sampledata_OTU)
random_tree_OTU = rtree(ntaxa(physeq_analysis_OTU), rooted = TRUE, tip.label = taxa_names(physeq_analysis_OTU))
physeq_analysis_OTU = merge_phyloseq(physeq_analysis_OTU, random_tree_OTU)

### Tundra incubations  ###
params = get_treatment_params(physeq_analysis_OTU, c('Substrate', 'Day'))
params
params = dplyr::filter(params, Substrate!='12CCon')
params
ex = "(Substrate=='12CCon' & Day=='${Day}') | (Substrate=='${Substrate}' & Day == '${Day}')"
physeq_analysis_OTU_l = phyloseq_subset(physeq_analysis_OTU, params, ex)

doParallel::registerDoParallel(2)
physeq_analysis_OTU_l_df = SIP_betaDiv_ord(physeq_analysis_OTU_l, parallel=TRUE)
phyloseq_ord_plot(physeq_analysis_OTU_l_df)

physeq_analysis_OTU_t = OTU_qPCR_trans(physeq_analysis_OTU, qPCRdata_OTU)
phyloseq::otu_table(physeq_analysis_OTU_t)

df_OTU <- HTSSIP:::qSIP_atom_excess_format(physeq_analysis_OTU_t, control_expr='Substrate=="12CCon"', treatment_rep='Microcosm_replicate')


atomX_OTU = qSIP_atom_excess(physeq_analysis_OTU_t,
                            control_expr='Substrate=="12CCon"',
                            treatment_rep='Microcosm_replicate', isotope = "13C", df_OTU_W = NULL)

atomX_OTU_W <- as.data.frame(atomX_OTU$W)
write.csv("atomX_OTU_W_greater_than_10.csv", x= atomX_OTU_W)

df_atomX_boot_OTU = qSIP_bootstrap(atomX_OTU, isotope = "13C", n_boot=100)
CI_threshold = 0
df_atomX_boot_OTU = df_atomX_boot_OTU %>%
  mutate(Incorporator = A_CI_low > CI_threshold,
         OTU = reorder(OTU, -A))

n_incorp = df_atomX_boot_OTU %>%
  filter(Incorporator == TRUE) %>%
  nrow
cat('Number of incorporators:', n_incorp, '\n')

# ---- Writing the qSIP outputs to the spreadsheets. ----

Ordered_OTU <- as.data.frame(df_atomX_boot_OTU)
as.data.frame(OTU_table_taxa)
Ordered_OTU <- merge(Ordered_OTU, as.data.frame(OTU_table_taxa), by.x="OTU", by.y="row.names", all.x=TRUE)
rownames(OTU_table_taxa)
unique(Ordered_OTU$Phylum)
write.csv("Ordered_OTU_greater_than_5.csv", x=Ordered_OTU)
Ordered_OTU <- read.csv("Ordered_OTU_greater_than_5.csv",  sep = ",", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)



