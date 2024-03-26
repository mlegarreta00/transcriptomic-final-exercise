# Instalación y activación de librerías:
BiocManager::install("apeglm")
install.packages("VennDiagram")
library("DESeq2")
library("tidyverse")
library("VennDiagram")

# Versiones de los paquetes utilizados 

R.version.string
packageVersion("DESeq2")
packageVersion("tidyverse")
packageVersion("VennDiagram")

# Creación del archivo .rnk

dds <- readRDS("input/dds3.rds")
res <- results(dds, alpha = 0.05, contrast = c("group", "DPN24h", "Control24h"))
summary(res)
## Reducción de los datos
res.ape <- lfcShrink(dds = dds, coef = "group_DPN24h_vs_Control24h", type = "apeglm",
                         res = res)
summary(res.ape)
## Visualización datos reducidos frente a los no reducidos
par(mfrow = c(1,2))
plotMA(res, ylim = c(-3, 3))
plotMA(res.ape, ylim = c(-3, 3))
## Guardar archivo .rnk
rnk <- data.frame(Feature = rownames(res.ape), LFC = res.ape$log2FoldChange)
write.table(rnk, file = "input/DPN_ranked.rnk", sep = "\t", quote = FALSE,
            col.names = FALSE, row.names = FALSE)

# Análisis GSEA en consola BASH

# Resultados

DPN_perturbed <- read.table(file = "resultados/GSEA/my_analysis.GseaPreranked.1711477851437/gsea_report_for_na_pos_1711477851437.tsv",
                            header = T, sep = "\t")
DPN_perturbed_refined <- DPN_perturbed[c(1,7:8)]
DPN_perturbed_refined

DPN_unperturbed <- read.table(file = "resultados/GSEA/my_analysis.GseaPreranked.1711477851437/gsea_report_for_na_neg_1711477851437.tsv",
                              header = T, sep = "\t")
DPN_unperturbed_refined <- DPN_unperturbed[c(1,7:8)]
DPN_unperturbed_refined

