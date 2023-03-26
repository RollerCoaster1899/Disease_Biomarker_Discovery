setwd("~/5.0 Proyectos/9.- Análisis de RNA-Seq/data")

# load libraries
library(dplyr)
library(tidyverse)
library(GEOquery)
library(ggplot2)
library(tibble)

# Read the data 
dat <- read.csv(file = "../data/DATA-1.csv")
dim(dat)

#get metadata / rename / cut 
gse <- getGEO(GEO = 'GSE183947', GSEMatrix = TRUE)

gse

  metadata <- pData(phenoData(gse[[1]]))
head(metadata)

metadata.subset <- select(metadata, c(1,10,11,17))

metadata.modified <- metadata %>% 
  select(1,10,11,17) %>%
  rename(tissue = characteristics_ch1) %>%
  rename(metastasis = characteristics_ch1.1) %>%
  mutate(tissue = gsub("tissue: ", "", tissue)) %>%
  mutate(metastasis = gsub("metastasis: ", "", metastasis))

head(dat)  

# reshaping data
dat.long <- dat %>%
  rename (gene = X) %>%
  gather(key = 'samples', value = 'FPKM', -gene)
  
# join dataframes = dat.long + metadata.modified

dat.long <- dat.long %>%
  left_join(.,metadata.modified, by = c("samples" = "description"))

write.csv(dat.long, "../data/DATA-Long.csv")

# Eliminate the entries with medium expressions and separated in cancer and normal
dat.Normal <- subset(dat.long, tissue != "normal breast tissue")
dat.Tumor <- subset(dat.long,tissue != "breast tumor") 

# Mean of same FPKM from Normal and Tumor 
Mean_Normal <- setNames(aggregate(dat.Normal$FPKM, list(dat.Normal$gene), mean),
                        c("Gene", "FPKM"))
Mean_Tumor <- setNames(aggregate(dat.Tumor$FPKM, list(dat.Tumor$gene), mean),
                        c("Gene", "FPKM"))

#Separate depending on Expression
Genes_TUX_FPKM <- subset(Mean_Tumor, FPKM <= 20)
Genes_TOX_FPKM <- subset(Mean_Tumor, FPKM >= 100)
Genes_NUX_FPKM <- subset(Mean_Normal, FPKM <=20)
Genes_NOX_FPKM <- subset(Mean_Normal, FPKM >= 100)

#Eliminate column of FPKM
Genes_TUX <- Genes_TUX_FPKM
Genes_TOX <- Genes_TOX_FPKM
Genes_NUX <- Genes_NUX_FPKM
Genes_NOX <- Genes_NOX_FPKM

Genes_TUX$FPKM <- NULL
Genes_TOX$FPKM <- NULL
Genes_NUX$FPKM <- NULL
Genes_NOX$FPKM <- NULL

# Exportar las columnas que se repiten
NUTO <- inner_join(Genes_NUX, Genes_TOX)
NOTU <- inner_join(Genes_NOX, Genes_TUX)

NUTO <- setNames(list(NUTO, Genes_NUX_FPKM, Genes_TOX_FPKM) %>% reduce(left_join, by = "Gene"), 
         c("Gene", "Normal", "Tumor"))
NOTU <- setNames(list(NOTU, Genes_NOX_FPKM, Genes_TUX_FPKM) %>% reduce(left_join, by = "Gene"), 
         c("Gene", "Normal", "Tumor"))

NUTO <- transform(NUTO, Overexpression = NUTO$Tumor/NUTO$Normal)
NOTU <- transform(NOTU, Underexpression = NOTU$Normal/NOTU$Tumor)

write.csv(NUTO, "../data/NUTO.csv")
write.csv(NOTU, "../data/NOTU.csv")

NUTO %>%
  ggplot(., aes(x = Gene, y = Overexpression)) +
  geom_col()

NOTU %>%
  ggplot(., aes(x = Gene, y = Underexpression)) +
  geom_col()

