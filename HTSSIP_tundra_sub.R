
library(ape)
library(doParallel)
library(dplyr)
library(foreach)
library(HTSSIP)
library(phyloseq)
library(reshape2)
library(stringr)
library(tidyr)
library(vegan)


setwd("C:/Users/User/Desktop/Svalbard/Svalbard_Data/Data_2023/Svalbard23_Data_qSIP/Tundra")

OTU_table <- read.csv("OTU_table_tundra.csv", sep=";", header = TRUE, check.names = FALSE)

OTU_table$Unlabeled_C <- (apply(OTU_table[, grepl("^Cf[0-9]+$", colnames(OTU_table))],1,sum))
OTU_table$Unlabeled_E <- (apply(OTU_table[, grepl("^Ef[0-9]+$", colnames(OTU_table))],1,sum))
OTU_table$Unlabeled_G <- (apply(OTU_table[, grepl("^Gf[0-9]+$", colnames(OTU_table))],1,sum))

OTU_table$Labeled_D <- (apply(OTU_table[, grepl("^Df[0-9]+$", colnames(OTU_table))],1,sum))
OTU_table$Labeled_B <- (apply(OTU_table[, grepl("^Ff[0-9]+$", colnames(OTU_table))],1,sum))
OTU_table$Labeled_H <- (apply(OTU_table[, grepl("^Hf[0-9]+$", colnames(OTU_table))],1,sum))

OTU_table <- OTU_table %>%
  filter(
    rowSums(select(., Unlabeled_C, Unlabeled_E, Unlabeled_G) >= 5) >= 2,
    rowSums(select(., Labeled_D, Labeled_F, Labeled_H) >= 5) >= 2)

OTU_table_samples <- select(OTU_table, matches("OTU|^[C-H]f[0-9]+$"))
colnames(OTU_table_samples)[1] <- "OTU"
rownames(OTU_table_samples) <- OTU_table_samples[,1]
OTU_table_samples[,1] <- NULL
OTU_table_samples <- as.matrix.data.frame(OTU_table_samples)
OTU_table_phyloseq_dataset = otu_table(OTU_table_samples, taxa_are_rows = TRUE)
write.csv("Dataset.csv", x=OTU_table_samples)


### Tundra taxa information

OTU_table_taxa <- select(OTU_table, matches("OTU|Kingdom|Phylum|Class|Subclass|Order|Suborder|Family|Genus|Species"))
OTU_table_taxa <- OTU_table_taxa %>%
  select(-Subclass, -Suborder)
colnames(OTU_table_taxa)[1] <- "OTU"
rownames(OTU_table_taxa) <- OTU_table_taxa[,1]
OTU_table_taxa[,1] <- NULL

t1 <- sapply(OTU_table_taxa, grepl, pattern="^unidentified|^uncultured|^Ambiguous|^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)$|^$")
t2 <- apply(t1, 1, any)
t3 <- apply(t1, 1, which.max)
for(i in seq_len(nrow(OTU_table_taxa))) {
  if(t2[i]) {OTU_table_taxa[i, t3[i]:ncol(OTU_table_taxa)]  <- paste("unidentified", OTU_table_taxa[i, t3[i]-1])}
}


### Create Phyloseq Taxa table object

OTU_table_taxa <- as.matrix.data.frame(OTU_table_taxa)
TAX_OTU = tax_table(OTU_table_taxa)


### Tundra qPCR and sample data ##

con <- file("qPCR_summary_tundra.csv", encoding = "UTF-8")
qPCRdata_OTU <- data.frame(read.csv(con, header = TRUE, sep = ";", check.names = FALSE,  row.names = 1, as.is=TRUE))
con1 <- file("qSIP_Sampledata_tundra.csv", encoding = "UTF-8")
sampledata_OTU <- data.frame(read.csv(con1, header = TRUE, sep = ";", check.names = FALSE,  row.names = 1, as.is=TRUE))

sampledata_OTU = sample_data(data.frame(sampledata_OTU))
physeq_analysis_OTU = phyloseq(OTU_table_phyloseq_dataset, TAX_OTU, sampledata_OTU)
random_tree_OTU = rtree(ntaxa(physeq_analysis_OTU), rooted = TRUE, tip.label = taxa_names(physeq_analysis_OTU))
physeq_analysis_OTU = merge_phyloseq(physeq_analysis_OTU, random_tree_OTU)


### Tundra incubations

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


### EAF calculations

### Inclusion of 15N in HTSSIP package

### 1) Function to unlock & replace functions inside the HTSSIP namespace
replace_in_ns <- function(fun_name, fun) {
  ns <- asNamespace("HTSSIP")
  if (bindingIsLocked(fun_name, ns)) unlockBinding(fun_name, ns)
  assign(fun_name, fun, envir = ns)
  lockBinding(fun_name, ns)
  invisible(TRUE)
}

### 2) calc_atom_excess: adding 15N natural abundance
new_calc_atom_excess <- function(Mlab, Mlight, Mheavymax, isotope = "13C") {
  isotope <- toupper(isotope)
  if (isotope == "13C") {
    x <- 0.01111233
  } else if (isotope == "18O") {
    x <- 0.002000429
  } else if (isotope == "15N") {
    x <- 0.003663004   # natural abundance of 15N
  } else {
    stop("isotope not recognized")
  }
  (Mlab - Mlight) / (Mheavymax - Mlight) * (1 - x)
}

### 3) calc_Mheavymax: adding Mheavymax for 15N
new_calc_Mheavymax <- function(Mlight, isotope = "13C", Gi = Gi) {
  isotope <- toupper(isotope)
  if (isotope == "13C") {
    Mhm <- -0.4987282 * Gi + 9.974564 + Mlight
  } else if (isotope == "18O") {
    Mhm <- 12.07747 + Mlight
  }  else if (isotope == "15N") {
    Mhm <- (0.5024851 * Gi) + 3.517396 + Mlight
  } else {
    stop("isotope not recognized")
  }
  
  Mhm
}

### 4) Replace HTSSIP package
replace_in_ns("calc_atom_excess", new_calc_atom_excess)
replace_in_ns("calc_Mheavymax",  new_calc_Mheavymax)


### 13C

atomX_OTU_13C = qSIP_atom_excess(physeq_analysis_OTU_t,
                                 control_expr='Substrate=="12CCon"',
                                 treatment_rep='Microcosm_replicate', isotope = "13C", df_OTU_W = NULL)
atomX_OTU_A_13C <- as.data.frame(atomX_OTU_13C$A)
atomX_OTU_W_13C <- as.data.frame(atomX_OTU_13C$W)
write.csv("atomX_OTU_W_13C.csv", x= atomX_OTU_W_13C)

df_atomX_boot_OTU_13C = qSIP_bootstrap(atomX_OTU_13C, isotope = "13C", n_boot=1000)
CI_threshold = 0
df_atomX_boot_OTU_13C = df_atomX_boot_OTU_13C %>%
  mutate(Incorporator = A_CI_low > CI_threshold,
         OTU = reorder(OTU, -A))

n_incorp = df_atomX_boot_OTU_13C %>%
  filter(Incorporator == TRUE) %>%
  nrow
cat('Number of incorporators:', n_incorp, '\n')

write.csv("13C_EAF_OTU.csv", x= df_atomX_boot_OTU_13C)


### 15N

atomX_OTU_15N = qSIP_atom_excess(physeq_analysis_OTU_t,
                                 control_expr='Substrate=="12CCon"',
                                 treatment_rep='Microcosm_replicate', isotope = "15N", df_OTU_W = NULL)
atomX_OTU_A_15N <- as.data.frame(atomX_OTU_15N$A)
atomX_OTU_W_15N <- as.data.frame(atomX_OTU_15N$W)
write.csv("atomX_OTU_W_15N.csv", x= atomX_OTU_W_15N)

df_atomX_boot_OTU_15N = qSIP_bootstrap(atomX_OTU_15N, isotope = "15N", n_boot=1000)
CI_threshold = 0
df_atomX_boot_OTU_15N = df_atomX_boot_OTU_15N %>%
  mutate(Incorporator = A_CI_low > CI_threshold,
         OTU = reorder(OTU, -A))

n_incorp = df_atomX_boot_OTU_15N %>%
  filter(Incorporator == TRUE) %>%
  nrow
cat('Number of incorporators:', n_incorp, '\n')

write.csv("15N_EAF_OTU.csv", x= df_atomX_boot_OTU_15N)


### Writing qSIP outputs

Ordered_OTU_13C <- as.data.frame(df_atomX_boot_OTU_13C)
as.data.frame(OTU_table_taxa)
Ordered_OTU_13C <- merge(Ordered_OTU_13C, as.data.frame(OTU_table_taxa), by.x="OTU", by.y="row.names", all.x=TRUE)
rownames(OTU_table_taxa)
names(Ordered_OTU_13C)
unique(Ordered_OTU_13C$Phylum)
unique(Ordered_OTU_13C$Class)
unique(Ordered_OTU_13C$Order)
unique(Ordered_OTU_13C$Family)

write.csv("Ordered_OTU_greater_than_5_13C.csv", x=Ordered_OTU_13C)
Ordered_OTU_13C <- read.csv("Ordered_OTU_greater_than_5_13C.csv",  sep = ",", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

Ordered_OTU_15N <- as.data.frame(df_atomX_boot_OTU_15N)
as.data.frame(OTU_table_taxa)
Ordered_OTU_15N <- merge(Ordered_OTU_15N, as.data.frame(OTU_table_taxa), by.x="OTU", by.y="row.names", all.x=TRUE)
rownames(OTU_table_taxa)
names(Ordered_OTU_15N)
unique(Ordered_OTU_15N$Phylum)
unique(Ordered_OTU_15N$Class)
unique(Ordered_OTU_15N$Order)
unique(Ordered_OTU_15N$Family)

write.csv("Ordered_OTU_greater_than_5_15N.csv", x=Ordered_OTU_15N)
Ordered_OTU_15N <- read.csv("Ordered_OTU_greater_than_5_15N.csv",  sep = ",", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)


### Notes

### From here calculations were done manually in Excel
### We sumed up the incorporation due to 13C and 15N proportional to their average relative content in an amino acid molecule (C:N = 3.78:1)
### Figures were produced with a separate code
