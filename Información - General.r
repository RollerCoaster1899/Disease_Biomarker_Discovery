# script to manipulate gene expression data 
# setwd("~/Folder with the File")
# file.exists("../Folder with the File/file")
# dim(variable) te dice el tamaño del database
# Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 1000) Para aumentar el size del buffer
#if (!require("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")
#install.packages('R.utils')
# Nombre <- variable
# head(variable) para ver lo de arriba de las variables
# para seleccionar varias columnas (numero) de una tabla select(variable, c(1,10,11,17))
# para seleccionar varias columnas (texto) de una tabla select(variable, c("1,10,11,17"))
# rename(NombreNuevo = NombreViejo) %>%
#  mutate(columna = gsub("Texto Antes: ", "Texto Despues", columna)) cambiar en la columna
# gather(key = 'samples', value = 'FPKM', -gene)
# ggplot(data, aes(x = variable, y = variable1))

library(data.table)
library(pheatmap)
library(ggsignif)
library(ggpubr)
library(data.table) 
library(ggpubr)
library(reshape2)
library(RColorBrewer)
library(corrplot)
library(ggpubr)
library(reshape2)
library(edgeR)
library(sva)
library(stringr)

dat.M10Y <- subset(dat.long, FPKM <= 10 & tissue != "normal breast tissue")
dat.M10N <- subset(dat.long, FPKM <= 10 & tissue != "breast tumor") 
dat.M70Y <- subset(dat.long, FPKM >= 70 & tissue != "normal breast tissue")
dat.M70N <- subset(dat.long, FPKM >= 70 & tissue != "breast tumor") 



NUTO <- NUTO %>%
  left_join(.,Genes_NUX_FPKM, Genes_TOX_FPKM, by = 'Gene')

NOTU <- NOTU %>%
  left_join(.,Genes_NOX_FPKM, Genes_TUX_FPKM, by = 'Gene')

NUTO %>%
  ggplot(., aes(x = Normal, y = Tumor)) + 
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE)



# explore data / graphing barplot
dat.long %>%
  filter(gene == 'BRCA1') %>%
  ggplot(., aes(x = samples, y = FPKM, fill = tissue)) +
  geom_col()

# 2. density plot
dat.long %>%
  filter(gene == 'BRCA1') %>%
  ggplot(., aes(x = FPKM, fill = tissue)) +
  geom_density(alpha = 0.3)

# 3. boxplot
dat.long %>%
  filter(gene == 'BRCA1') %>%
  ggplot(., aes(x = metastasis, y = FPKM)) + 
  geom_boxplot()

# 4. Scatterplot
NUTO %>%
  spread(key = Gene, value = FPKM) %>%
  ggplot(., aes(x = Tumor, y = BRCA2, color = tissue)) + 
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE)

# 5. heatmap
genes.of.interest <- c('BRCA1', 'BRCA2', 'TP53', 'ALK', 'MYCN')

p <- dat.long %>%
  filter(gene %in% genes.of.interest) %>%
  ggplot(., aes(x = samples, y = gene, fill = FPKM)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_gradient(low = 'white', high = 'red') +
  geom_tile()

ggsave(p, filename = 'heatmap_save2.png', width = 10,height = 8)

NUTO <- transform(NUTO, OVEREX = NUTO$Tumor/NUTO$Normal)
NOTU <- transform(NOTU, OVEREX = NOTU$Tumor/NOTU$Normal)






