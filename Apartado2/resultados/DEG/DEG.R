# Instalación y activación de librerías:

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("vsn", force = TRUE)
library("DESeq2")
library("tidyverse")
library("pheatmap")
library("RColorBrewer")
library("vsn")

# Versiones de los paquetes utilizados 

R.version.string
packageVersion("DESeq2")
packageVersion("tidyverse")
packageVersion("pheatmap")
packageVersion("RColorBrewer")
packageVersion("vsn")

# Preprocesado de los datos

## Cargamos los archivos
raw_data <- read.csv(file = "input/rawcounts.tsv", sep = "\t", row.names = 1)
colnames(raw_data)
experiment_data <- read.csv(file = "input/metadata.tsv", sep = "\t")
rownames(experiment_data) <- colnames(raw_data)
experiment_data <- mutate(.data = experiment_data,
                          X = NULL,
                          patient = as.factor(patient),
                          agent = as.factor(agent),
                          time = as.factor(time))
all(colnames(raw_data) %in% rownames(experiment_data))
all(colnames(raw_data) == rownames(experiment_data))
## Creación nueva variable
experiment_data$group <- as.factor(paste0(experiment_data$agent, experiment_data$time))
levels(experiment_data$group)

# Creación de DESeqDataSet

dds <- DESeqDataSetFromMatrix(countData = raw_data,
                              colData = experiment_data,
                              design = ~ patient + group)
dds
## Eliminación de genes con cuenta menor que 10 
keep <- rowSums(counts(dds)) >= 10
dds2 <- dds[keep, ]
dim(dds)
dim(dds2)
dds3 <- DESeq(dds2, test = "Wald")
saveRDS(object = dds3, file = "input/dds3.rds")
# Análisis exploratorio

## VST
vsd <- vst(dds2, blind = TRUE)
## Visualización gráfica de VST
### Normalización de los datos de expresión
normal_data <- normTransform(dds2)
### Comparación efectos de transformación de la varianza
### frente a las muestras normalizadas
normal_graph <- meanSdPlot(assay(normal_data))
vsd_graph <- meanSdPlot(assay(vsd))
### Análisis de componentes principales (PCA)
plotPCA(vsd, intgroup = "patient")
plotPCA(vsd, intgroup = "group")
plotPCA(vsd, intgroup = "agent")
plotPCA(vsd, intgroup = "time")
### Matriz de distancias
#### Cálculo de las distancias entre muestras
sampleDists <- dist(t(assay(vsd)))
#### Creación de matriz de distancias
sampleDistMatrix <- as.matrix( sampleDists )
#### Asignación de nombres a las filas de la matriz de distancias
rownames(sampleDistMatrix) <- paste( vsd$patient, vsd$group, sep = " - " )
colnames(sampleDistMatrix) <- NULL
#### Asignación de colores para el heatmap
colors <- colorRampPalette( rev(brewer.pal(9, "Greens")) )(255)
#### Generación del heatmap
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

# Análisis de expresión diferencial
## Eliminación del outlier
patient_outlier <- which(experiment_data$patient == 4)
agent_outlier <- which(experiment_data$agent == "Control")
time_outlier <- which(experiment_data$time == "24h")
patient_agent <- intersect(patient_outlier, agent_outlier)
patient_time <- intersect(patient_outlier, time_outlier)
outlier <- intersect(patient_agent, patient_time)
experiment_data_clean <- experiment_data[-outlier,]
raw_data_clean <- raw_data[ ,-outlier]
## Creación nuevo DESeqDataSet
dds_outlier <- DESeqDataSetFromMatrix(countData = raw_data_clean,
                                      colData = experiment_data_clean,
                                      design = ~ patient + group)
dds_outlier
keep <- rowSums(counts(dds_outlier)) >= 10
dds2_outlier <- dds_outlier[keep, ]
dim(dds_outlier)
dim(dds2_outlier)
## VST
vsd_outlier <- vst(dds2_outlier, blind = TRUE)
## Visualización gráfica de VST
### Normalización de los datos de expresión
normal_data_outlier <- normTransform(dds2_outlier)
### Comparación efectos de transformación de la varianza
### frente a las muestras normalizadas
normal_graph_outlier <- meanSdPlot(assay(normal_data_outlier))
vsd_graph_outlier <- meanSdPlot(assay(vsd_outlier))
### Análisis de componentes principales (PCA)
plotPCA(vsd_outlier, intgroup = "patient")
plotPCA(vsd_outlier, intgroup = "group")
plotPCA(vsd_outlier, intgroup = "agent")
plotPCA(vsd_outlier, intgroup = "time")
### Matriz de distancias
#### Cálculo de las distancias entre muestras
sampleDists_outlier <- dist(t(assay(vsd_outlier)))
#### Creación de matriz de distancias
sampleDistMatrix_outlier <- as.matrix( sampleDists_outlier )
#### Asignación de nombres a las filas de la matriz de distancias
rownames(sampleDistMatrix_outlier) <- paste( vsd_outlier$patient, vsd_outlier$group, sep = " - " )
colnames(sampleDistMatrix_outlier) <- NULL
#### Asignación de colores para el heatmap
colors <- colorRampPalette( c("#E91E63", "#EC407A", "#F06292", "#F48FB1", "#FCE4EC"))(255)
#### Generación del heatmap
pheatmap(sampleDistMatrix_outlier,
         clustering_distance_rows = sampleDists_outlier,
         clustering_distance_cols = sampleDists_outlier,
         col = colors)
### Generación nuevo DESeqDataSet
dds3_outlier <- DESeq(dds2_outlier, test = "Wald")
### Generación de gráficos
plotDispEsts(dds3_outlier)
plotMA(dds3_outlier)
### Tratamiento DPN
res_DPN_outlier <- results(object = dds3_outlier,
                           contrast = c("group", "Control24h", "DPN24h"),
                           alpha = 0.05,
                           pAdjustMethod = "BH"
                           )
summary(res_DPN_outlier)
res_DPN_lfc1_outlier <- results(object = dds3_outlier,
                           contrast = c("group", "Control24h", "DPN24h"),
                           alpha = 0.05,
                           lfcThreshold = 1,
                           pAdjustMethod = "BH"
)
summary(res_DPN_lfc1_outlier)
### Tratamiento con OHT
res_OHT_outlier <- results(object = dds3_outlier,
                           contrast = c("group", "Control24h", "OHT24h"),
                           alpha = 0.05,
                           pAdjustMethod = "BH"
)
summary(res_OHT_outlier)
res_OHT_lfc1_outlier <- results(object = dds3_outlier,
                                contrast = c("group", "Control24h", "OHT24h"),
                                alpha = 0.05,
                                lfcThreshold = 1,
                                pAdjustMethod = "BH"
)
summary(res_OHT_lfc1_outlier)


# Visualización de resultados
mat_DPN_outlier <- assay(vsd_outlier)[head(order(res_DPN_outlier$padj), 30), ]
mat_DPN <- mat_DPN_outlier - rowMeans(mat_DPN_outlier)
anno <- as.data.frame(colData(vsd_outlier)[ , c("patient", "agent", "time")])
ann_colors <- list(
    patient = c("1" ="blue", "2"= "purple", "3" = "green", "4" = "yellow"),
    agent = c(Control ="turquoise", DPN ="pink", OHT = "orange"),
    time = c("24h" = "aquamarine", "48h" ="red")
)

mat_OHT_outlier <- assay(vsd_outlier)[head(order(res_OHT_outlier$padj), 30), ]
mat_OHT <- mat_OHT_outlier - rowMeans(mat_OHT_outlier)
anno <- as.data.frame(colData(vsd_outlier)[ , c("patient", "agent", "time")])
ann_colors <- list(
    patient = c("1" ="blue", "2"= "purple", "3" = "green", "4" = "yellow"),
    agent = c(Control ="turquoise", DPN ="pink", OHT = "orange"),
    time = c("24h" = "aquamarine", "48h" ="red")
)

par(mfrow = c(1, 2))
pheatmap(mat = mat_DPN_outlier, annotation_col = anno, show_colnames = F, 
         annotation_colors = ann_colors, main = "DPN_outlier")
pheatmap(mat = mat_OHT_outlier, annotation_col = anno, show_colnames = F, 
         annotation_colors = ann_colors, main = "OHT_outlier")


