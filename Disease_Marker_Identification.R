# Load libraries
library(dplyr)
library(tidyverse)
library(GEOquery)
library(ggplot2)
library(tibble)

# Read the data
dat <- read.csv(file = "GSE183947_fpkm.csv")
dim(dat)

# Get metadata / rename / cut
gse <- getGEO(GEO = 'GSE183947', GSEMatrix = TRUE)
metadata <- pData(phenoData(gse[[1]]))
head(metadata)
metadata.subset <- metadata %>%
  select(1, 10, 11, 17)
metadata.modified <- metadata.subset %>%
  rename(tissue = characteristics_ch1) %>%
  rename(metastasis = characteristics_ch1.1) %>%
  mutate(tissue = gsub("tissue: ", "", tissue)) %>%
  mutate(metastasis = gsub("metastasis: ", "", metastasis))

head(dat)

# Reshape data
dat.long <- dat %>%
  rename(gene = X) %>%
  gather(key = 'samples', value = 'FPKM', -gene)

# Join dataframes: dat.long + metadata.modified
dat.long <- dat.long %>%
  left_join(metadata.modified, by = c("samples" = "description"))

write.csv(dat.long, "DATA-Long.csv")

# Eliminate the entries with medium expressions and separate into cancer and normal
dat.Normal <- subset(dat.long, tissue != "normal breast tissue")
dat.Tumor <- subset(dat.long, tissue != "breast tumor")

# Mean of same FPKM from Normal and Tumor
Mean_Normal <- setNames(aggregate(dat.Normal$FPKM, list(dat.Normal$gene), mean),
                        c("Gene", "FPKM"))
Mean_Tumor <- setNames(aggregate(dat.Tumor$FPKM, list(dat.Tumor$gene), mean),
                       c("Gene", "FPKM"))

# Separate depending on Expression
Genes_TUX_FPKM <- subset(Mean_Tumor, FPKM <= 20)
Genes_TOX_FPKM <- subset(Mean_Tumor, FPKM >= 100)
Genes_NUX_FPKM <- subset(Mean_Normal, FPKM <= 20)
Genes_NOX_FPKM <- subset(Mean_Normal, FPKM >= 100)

# Eliminate column of FPKM
Genes_TUX <- select(Genes_TUX_FPKM, -FPKM)
Genes_TOX <- select(Genes_TOX_FPKM, -FPKM)
Genes_NUX <- select(Genes_NUX_FPKM, -FPKM)
Genes_NOX <- select(Genes_NOX_FPKM, -FPKM)

# Export the columns that are repeated
NUTO <- inner_join(Genes_NUX, Genes_TOX, by = "Gene")
NOTU <- inner_join(Genes_NOX, Genes_TUX, by = "Gene")

NUTO <- left_join(NUTO, Genes_NUX_FPKM, by = "Gene") %>%
  left_join(Genes_TOX_FPKM, by = "Gene") %>%
  rename(Normal = FPKM.x, Tumor = FPKM.y)
NOTU <- left_join(NOTU, Genes_NOX_FPKM, by = "Gene") %>%
  left_join(Genes_TUX_FPKM, by = "Gene") %>%
  rename(Normal = FPKM.x, Tumor = FPKM.y)

NUTO$Overexpression <- NUTO$Tumor / NUTO$Normal
NOTU$Underexpression <- NOTU$Normal / NOTU$Tumor

write.csv(NUTO, "NUTO.csv")
write.csv(NOTU, "NOTU.csv")

ggplot(NUTO, aes(x = Gene, y = Overexpression)) +
  geom_col()

ggplot(NOTU, aes(x = Gene, y = Underexpression)) +
  geom_col()