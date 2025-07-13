# Este script realiza análisis de RNA-seq en lenguaje de R, incluyendo
# la importación de metadatos y de la matriz de recuentos.

# Establecer directorio de trabajo 
setwd(dir= "/home/laltadilll/Documentos/TFM/")

#### DE LA CARGA DE METADATOS Y MATRIZ DE RECUENTOS AL OBJETO DGEList ####

### 1. IMPORTAR METADATOS

## 1.1. Carga de librerías
library(readr) # para lectura de archivos .tsv
library(limma) # necesario para cargar paquete edgeR y análisis GO
library(edgeR) # para análisis transcriptómico
library(tidyr) # necesario para cargar paquete dplyr
library(dplyr) # para edición de columnas en el dataframe
library(statmod) # para modelización estadística
library(ggplot2) # para representación de gráficos (Volcano Plot)
library(ggrepel) # para evitar solapamientos de texto en gráficos
library(patchwork) # para representar varios gráficos ggplot2 juntos
library(BioVenn) # diagramas de Venn
library(pheatmap) # para representación del heatmap
library(org.Hs.eg.db) # anotación génica de Homo sapiens

## 1.2. Cargar metadatos
metadata <- read_tsv("Metadata/Metadata_with_sex_age.txt")

## 1.3. Editar columna "Sex_Age" metadatos
metadata <- metadata %>% 
  separate(Sex_Age,
           into = c("Age","Sex"),
           sep  = 1) %>% # separa justo después del primer carácter 
  mutate(Age = recode(Age, Y="Young", O="Old"),
         Sex = recode(Sex, F="Female", M="Male"),
         Timepoint = recode(Timepoint, PreTraining="Pre", PostTraining="Post")) %>%
  filter(Age != "Old")

### 2. IMPORTAR Y PREPRAR MATRIZ DE RECUENTOS. CONVERSIÓN A OBJETO DGEList

## 2.1. Cargar matriz de recuentos
counts1 = read.delim(file = "Results/GSE97084_GeneCount_raw.tsv.gz", 
                     stringsAsFactors = FALSE)

## 2.2. Guardar anotaciones (GeneID + Chr) para añadirlas en objeto DGEList
annotation_org <- AnnotationDbi::select(org.Hs.eg.db, keys = counts1$GeneID, 
                                       keytype = "ENSEMBL", 
                                       columns = c("ENTREZID", "SYMBOL", "GENENAME"))

annotation_org_uniq <- annotation_org[!duplicated(annotation_org$ENSEMBL), ]
all(counts1$GeneID == annotation_org_uniq$ENSEMBL)
counts1$EntrezID <- annotation_org_uniq$ENTREZID

## 2.3. Editar matriz de recuentos y metadatos previo a conversión

# Simplificar el nombre del ID de cada participante
colnames(counts1) <- gsub("^s_([0-9]+[A-Za-z])\\.180_GeneCount$", "\\1", 
                         colnames(counts1))

# Establecer como encabezado de las filas la columna "GeneID"
nrow(counts1[is.na(counts1$EntrezID),]) # 24778 IDs vacíos (NA)
nrow(counts1[duplicated(counts1$EntrezID),]) # 60 IDs duplicados
counts1_clean <- counts1 %>% filter(!is.na(counts1$EntrezID) &
                                    !duplicated(counts1$EntrezID))
rownames(counts1_clean) <- counts1_clean$EntrezID

## 2.2. Guardar anotaciones (Chr) para añadirlas en objeto DGEList
annotation_chr <- counts1_clean[, 1]
annotation_org_clean <- annotation_org_uniq[!duplicated(annotation_org_uniq$ENTREZID) &
                                                 !is.na(annotation_org_uniq$ENTREZID),]
all(rownames(counts1_clean) == annotation_org_clean$ENTREZID)
annotation_org_filtered <- annotation_org_clean %>% dplyr::select(-ENSEMBL)
annotation <- cbind(annotation_chr,annotation_org_filtered)
names(annotation)[1] <- "Chr"

# Eliminar las columnas "Chr", "GeneName", "Start", "Stop" y "CodingLength"
counts1_clean <- counts1_clean[,-c(1:6,ncol(counts1_clean))]

# Cambiar a mayúsculas (18a) los IDs para comparar con los de metadata
colnames(counts1_clean) <- toupper(colnames(counts1_clean))

# Eliminar en los conteos los IDs descartados previamente en metadata (27B, Old)
counts1_filtered <- counts1_clean[, colnames(counts1_clean) %in% metadata$ID]

# Eliminar en metadata los IDs no incluidos en los conteos
metadata1 <- metadata[metadata$ID %in% colnames(counts1_filtered),]

## 2.4. Categorizar grupos (según tipo de entrenamiento, sexo_edad y timepoint)
group1 <- paste(metadata1$Sex,metadata1$Training, metadata1$Timepoint, sep = ".")
group1 <- factor(group1)
table(group1)

#### FILTRADO Y NORMALIZACIÓN OBJETO DGEList ####

### 1. CREAR OBJETO DGEList  
dge1 <- DGEList(counts1_filtered)
dge1$samples$group <- group1 # añadir los grupos previamente definidos 
dge1$genes <- annotation # añadir las anotaciones guardadas

### 2. FILTRAR GENES DE BAJA EXPRESIÓN 

## 2.1. Filtrar los genes con un mínimo de n lecturas por muestra y total
genes_filtered_1 <- filterByExpr(dge1) 
summary(genes_filtered_1) #perdemos 15.650 genes de 32934 (un 48%, valor esperable)

# Realizar el subsetting de genes filtrados en el objeto DGEList.
# Restablecer el tamaño de librería con los datos filtrados. 
dge1 <- dge1[genes_filtered_1,keep.lib.sizes=FALSE] 

### 3. NORMALIZAR CONTEO DE GENES

# Normalizar de la composición de la librería
dge1 <- calcNormFactors(dge1) 
head(dge1$samples)

#### MULTIDIMENSIONAL SCALING (MDS) PLOT ####

# Generar 3 archivos diferentes que contengan:
# MDS PLOT GENERAL (2 plots: con y sin genes chrY)
# MDS PLOT POR TIPO DE ENTRENAMIENTO (6 plots: con y sin genes chrY)
# MDS PLOT POR TIMEPOINT (4 plots: con y sin genes chrY)

### 1. SUBSETTING GENES (SIN GENES CHRY) ####

## 1.1. Creación, filtrado y normalización objeto DGEList (sin genes chrY)
# Subset de la matriz de recuentos filtrando las filas sin genes chrY
annotation_no_chrY <- annotation %>% filter(Chr != "chrY")
counts1_no_chrY <- counts1_filtered[rownames(counts1_filtered) %in% annotation_no_chrY$ENTREZID,]

# Crear objeto DGEList
dge1_no_chrY <- DGEList(counts1_no_chrY)
dge1_no_chrY$samples$group <- group1 # añadir los grupos previamente definidos 
dge1_no_chrY$genes <- annotation_no_chrY # añadir las anotaciones guardadas

# Filtrar genes de baja expresión 
genes_filtered_1 <- filterByExpr(dge1_no_chrY) 
dge1_no_chrY <- dge1_no_chrY[genes_filtered_1,keep.lib.sizes=FALSE] 

# Normalizar conteo de genes
dge1_no_chrY <- calcNormFactors(dge1_no_chrY)

### 2. SUBSETTING PARTICIPANTES (POR ENTRENAMIENTO Y POR TIMEPOINT, ####
                                 #CON O SIN GENES CHRY)

# Por entrenamiento (H =HIIT, R = RESISTANCE, C = COMBINED)
# Por timepoint (Pre = PreTraining, Post = PostTraining)
# Con/sin genes chrY (no_chrY = sin genes chrY)

# Subset del objeto group
group1_H <- grep("HIIT", dge1$samples$group, value = TRUE)
group1_R <- grep("Resistance", dge1$samples$group, value = TRUE)
group1_C <- grep("Combined", dge1$samples$group, value = TRUE)
group1_Pre <- grep("Pre", dge1$samples$group, value = TRUE)
group1_Post <- grep("Post", dge1$samples$group, value = TRUE)

group1_H <- factor(group1_H)
group1_R <- factor(group1_R)
group1_C <- factor(group1_C)
group1_Pre <- factor(group1_Pre)
group1_Post <- factor(group1_Post)

# Subset del objeto DGEList (por entrenamiento y por timepoint)
dge1_H <- dge1[, dge1$samples$group %in% group1_H]
dge1_R <- dge1[, dge1$samples$group %in% group1_R]
dge1_C <- dge1[, dge1$samples$group %in% group1_C]
dge1_H_no_chrY <- dge1_no_chrY[, dge1_no_chrY$samples$group %in% group1_H]
dge1_R_no_chrY <- dge1_no_chrY[, dge1_no_chrY$samples$group %in% group1_R]
dge1_C_no_chrY <- dge1_no_chrY[, dge1_no_chrY$samples$group %in% group1_C]

dge1_Pre <- dge1[, dge1$samples$group %in% group1_Pre]
dge1_Post <- dge1[, dge1$samples$group %in% group1_Post]
dge1_Pre_no_chrY <- dge1_no_chrY[, dge1_no_chrY$samples$group %in% group1_Pre]
dge1_Post_no_chrY <- dge1_no_chrY[, dge1_no_chrY$samples$group %in% group1_Post]

### 3. DEFINICIÓN DE FORMAS, COLORES Y ETIQUETAS PARA CADA PLOT MDS ####
# HIIT: Círculos. Mujeres: "#7F3C8D", Hombres: "#00BFA0"
# Resistance: Triángulos. Mujeres: "#E03C31", Hombres: "#F0E442"
# Combined: Cuadrados. Mujeres: "#56B4E9", Hombres: "#CC79A7"
# PreTraining: forma vacía
# PostTraining: forma rellena

## 3.1. Formas
pch_all <- rep(c(0,15,1,16,2,17),2)
pch_H <- rep(c(16,1),2) 
pch_R <- rep(c(17,2),2) 
pch_C <- rep(c(15,0),2) 
pch_Pre <- rep(c(0,1,2),2)
pch_Post <- rep(c(15,16,17),2)

## 3.2. Colores
colors_all <- rep(c("#56B4E9","#7F3C8D","#E03C31",
                    "#CC79A7","#00BFA0","#F0E442"),
                  each=2) # paleta okabe_ito y Tol (compatible con daltonismo)
colors_H <- rep(c("#7F3C8D","#00BFA0"),each=2) 
colors_R <- rep(c("#E03C31","#F0E442"),each=2) 
colors_C <- rep(c("#56B4E9","#CC79A7"),each=2) 
colors_Pre <- c("#56B4E9","#7F3C8D","#E03C31",
                "#CC79A7","#00BFA0","#F0E442")
colors_Post <- colors_Pre

## 3.3. Etiquetas (en plots con muchos datos, seleccionar las más relevantes)
labels_all <- rownames(dge1$samples) %in% c("2C","2B","18B","14B","23B","24B")
labels_all_no_chrY <- rownames(dge1$samples) %in% c("2C","2B","18B","14B","23B","24B")
labels_Pre <- rownames(dge1_Pre$samples) %in% c("2B", "18B","14B")
labels_Post <- rownames(dge1_Post$samples) %in% c("23B", "24B")
labels_Pre_no_chrY <- rownames(dge1_Pre_no_chrY$samples) %in% c("2B","18B", "14B")
labels_Post_no_chrY <- rownames(dge1_Post_no_chrY$samples) %in% c("2C", "23B", "24B")

### 4. GUARDAR MDS PLOT en objeto mds: ####
# 1. MDS PLOT GENERAL (con y sin genes chrY)
# 2. MDS POR ENTRENAMIENTO (con y sin genes chrY)
# 3. MDS PLOT POR TIMEPOINT (con y sin genes chrY)

## 4.1. MDS PLOT GENERAL (con y sin genes chrY) 
mds1 <- plotMDS(dge1, plot=FALSE) 
mds1_no_chrY <- plotMDS(dge1_no_chrY, plot=FALSE)

## 4.2. MDS POR ENTRENAMIENTO (con y sin genes chrY) 
mds1_H <- plotMDS(dge1_H, plot=FALSE) 
mds1_R <- plotMDS(dge1_R, plot=FALSE) 
mds1_C <- plotMDS(dge1_C, plot=FALSE)
mds1_H_no_chrY <- plotMDS(dge1_H_no_chrY, plot=FALSE) 
mds1_R_no_chrY <- plotMDS(dge1_R_no_chrY, plot=FALSE) 
mds1_C_no_chrY <- plotMDS(dge1_C_no_chrY, plot=FALSE)

## 4.3. MDS PLOT POR TIMEPOINT (con y sin genes chrY)
mds1_Pre <- plotMDS(dge1_Pre, plot=FALSE) 
mds1_Post <- plotMDS(dge1_Post, plot=FALSE) 
mds1_Pre_no_chrY <- plotMDS(dge1_Pre_no_chrY, plot=FALSE) 
mds1_Post_no_chrY <- plotMDS(dge1_Post_no_chrY, plot=FALSE) 

### 5. GUARDAR MDS PLOT EN TRES ARCHIVOS ####

## 5.1. MDS PLOT GENERAL (2 plots: con y sin genes chrY) 
pdf("Plots/mds_All_wn_chrY.pdf", width = 7, height = 6)
par(mfrow = c(2, 1),mar=c(5, 5, 1, 1))

# MDS Plot General (con genes chrY)
plot(mds1, 
     col=colors_all[group1], 
     pch=pch_all[group1], 
     xlim=c(-1,2.85),
     ylim=c(-1,2.1),
     xlab="", 
     ylab="Leading logFC dim 2",
     xaxt="n",
     cex.axis=0.8)
text(mds1$x[labels_all],
     mds1$y[labels_all], 
     labels=rownames(dge1$samples)[labels_all],
     pos=4,      # 1=abajo,2=izquierda,3=arriba,4=derecha
     cex=0.7,
     font = 2,
     offset=0.4)

# MDS Plot general (sin genes chrY)
plot(mds1_no_chrY, 
     col=colors_all[group1], 
     pch=pch_all[group1], 
     xlim=c(-1,2.85),
     ylim=c(-1,2.1),
     xlab="Leading logFC dim 1", 
     ylab="Leading logFC dim 2",
     cex.axis=0.8)
text(mds1_no_chrY$x[labels_all_no_chrY],
     mds1_no_chrY$y[labels_all_no_chrY], 
     labels=rownames(dge1$samples)[labels_all_no_chrY],
     pos=4,      
     cex=0.7,  
     font=2,
     offset=0.4)
legend("topright", 
       legend=levels(group1), 
       pch=pch_all, 
       col=colors_all, 
       inset=c(0.0,0.01), 
       cex=0.55, 
       text.font=2,
       ncol=2)

dev.off()

## 5.2. MDS PLOT POR TIPO DE ENTRENAMIENTO (6 plots: con y sin genes chrY)
pdf("Plots/mds_HRC_wn_chrY.pdf", width = 7, height = 4)
par(mfrow = c(2, 3),mar=c(2, 2, 0.5, 0.5))

# HIIT (con genes chrY)
plot(mds1_H, 
     col=colors_H[group1_H], 
     pch=pch_H[group1_H], 
     xlim=c(-1.4, 3),
     ylim=c(-1.5, 1.5),
     xlab="", 
     ylab="Leading logFC dim 2",
     xaxt="n",
     cex.axis=0.8)
text(mds1_H, 
     labels=rownames(dge1_H$samples),
     pos=4,    
     cex=0.6,
     font=2, # negrita
     offset=0.4)
legend("top", 
       legend=levels(group1_H), 
       pch=pch_H, 
       col=colors_H, 
       cex=0.6,
       text.font=2,
       inset=(c(0.01,0.01)),
       ncol=2)

# Resistance (con genes chrY)
plot(mds1_R, 
     col=colors_R[group1_R], 
     pch=pch_R[group1_R], 
     xlim=c(-1.4,3),
     ylim=c(-1.5, 1.5),
     xlab="",
     ylab="",
     xaxt="n",
     yaxt="n")
text(mds1_R, 
     labels=rownames(dge1_R$samples),
     pos=4,      
     cex=0.6, 
     font=2,
     offset=0.4)
legend("top", 
       legend=levels(group1_R), 
       pch=pch_R, 
       col=colors_R, 
       cex=0.6, 
       text.font=2,
       inset=c(0.01,0.01),
       ncol=2)

# Combined (con genes chrY)
plot(mds1_C, 
     col=colors_C[group1_C], 
     pch=pch_C[group1_C], 
     xlim=c(-1.4,3),
     ylim=c(-1.5, 1.5),
     xlab="", 
     ylab="",
     xaxt="n",
     yaxt="n")
text(mds1_C, 
     labels=rownames(dge1_C$samples),
     pos=4,      
     cex=0.6, 
     font = 2,
     offset=0.4)
legend("right", 
       legend=levels(group1_C), 
       pch=pch_C, 
       col=colors_C, 
       cex=0.6,
       text.font=2,
       inset=c(0.01,0.01),
       ncol=1)

# HIIT (sin genes chrY)
plot(mds1_H_no_chrY, 
     col=colors_H[group1_H], 
     pch=pch_H[group1_H], 
     xlim=c(-1.4, 3),
     ylim=c(-1.5, 1.5),
     xlab="Leading logFC dim 1", 
     ylab="",
     cex.axis=0.8)
text(mds1_H_no_chrY, 
     labels=rownames(dge1_H_no_chrY$samples),
     pos=4,      
     cex=0.6,
     font=2,
     offset=0.4)

# Resistance (sin genes chrY)
plot(mds1_R_no_chrY, 
     col=colors_R[group1_R], 
     pch=pch_R[group1_R], 
     xlim=c(-1.4,3),
     ylim=c(-1.5, 1.5),
     xlab="Leading logFC dim 1", 
     ylab="",
     yaxt="n",
     cex.axis=0.8)
text(mds1_R_no_chrY, 
     labels=rownames(dge1_R_no_chrY$samples),
     pos=4,    
     cex=0.6, 
     font=2,
     offset=0.4)

# Combined (sin genes chrY)
plot(mds1_C_no_chrY, 
     col=colors_C[group1_C], 
     pch=pch_C[group1_C],
     xlim=c(-1.4,3),
     ylim=c(-1.5,1.5),
     xlab="Leading logFC dim 1", 
     ylab="",
     yaxt="n",
     cex.axis=0.8)
text(mds1_C_no_chrY, 
     labels=rownames(dge1_C_no_chrY$samples),
     pos=4,      
     cex=0.6, 
     font=2,
     offset=0.4)

dev.off()

## 5.3. MDS PLOT POR TIMEPOINT (4 plots: con y sin genes chrY)
pdf("Plots/mds_PrePost_wn_chrY.pdf", width = 7, height = 4)
par(mfrow = c(2, 2),mar=c(2, 2, 0.5, 0.5))

# PreTraining (con genes chrY)
plot(mds1_Pre, 
     col=colors_Pre[group1_Pre], 
     pch=pch_Pre[group1_Pre], 
     xlim=c(-2.8,2.8),
     ylim=c(-1.2,2),
     xlab="", 
     ylab="Leading logFC dim 1",
     xaxt="n",
     cex.axis=0.8)
text(mds1_Pre$x[labels_Pre],
     mds1_Pre$y[labels_Pre], 
     labels=rownames(dge1_Pre$samples)[labels_Pre],
     pos=4,      
     cex=0.7, 
     font=2,
     offset=0.4)

# PostTraining (con genes chrY)
plot(mds1_Post, 
     col=colors_Post[group1_Post], 
     pch=pch_Post[group1_Post], 
     xlim=c(-2.8,2.8),
     ylim=c(-1.2,2),
     xlab="", 
     ylab="",
     xaxt="n",
     yaxt="n")
text(mds1_Post$x[labels_Post], 
     mds1_Post$y[labels_Post],
     labels=rownames(dge1_Post$samples)[labels_Post],
     pos=4,     
     cex=0.7, 
     font=2,
     offset=0.4)

# Pretraining (sin genes chrY)
plot(mds1_Pre_no_chrY, 
     col=colors_Pre[group1_Pre], 
     pch=pch_Pre[group1_Pre], 
     xlim=c(-2.8,2.8),
     ylim=c(-1.2,2),
     xlab="Leading logFC dim 1", 
     ylab="Leading logFC dim 2",
     cex.axis=0.8)
text(mds1_Pre_no_chrY$x[labels_Pre_no_chrY],
     mds1_Pre_no_chrY$y[labels_Pre_no_chrY], 
     labels=rownames(dge1_Pre_no_chrY$samples)[labels_Pre_no_chrY],
     pos=4,      
     cex=0.7, 
     font=2,
     offset=0.4)
legend("topright", 
       legend=levels(group1_Pre), 
       pch=pch_Pre, 
       col=colors_Pre, 
       cex=0.6, 
       text.font=2,
       inset=c(0.005,0.01),
       ncol=1)

# Posttraining (sin genes chrY)
plot(mds1_Post_no_chrY, 
     col=colors_Post[group1_Post], 
     pch=pch_Post[group1_Post], 
     xlim=c(-2.8,2.8),
     ylim=c(-1.2,2),
     xlab="Leading logFC dim 1", 
     ylab="",
     yaxt="n",
     cex.axis=0.8)
text(mds1_Post_no_chrY$x[labels_Post_no_chrY], 
     mds1_Post_no_chrY$y[labels_Post_no_chrY],
     labels=rownames(dge1_Post_no_chrY$samples)[labels_Post_no_chrY],
     pos=4,      
     cex=0.7, 
     font=2,
     offset=0.4)
legend("topright", 
       legend=levels(group1_Post), 
       pch=pch_Post, 
       col=colors_Post, 
       cex=0.6, 
       text.font=2, 
       inset=c(0.005,0.01),
       ncol=1)

dev.off()        

#### ESTIMACIÓN DE DISPERSIONES Y AJUSTE DE MODELO LINEAR GENERALIZADO (GLM) ####

# 1. Estimation of dispersions (common, trended, and tagwise)
# 2. Generalized Linear Model (GLM) fitting using quasi-likelihood methods
# 3. Visualization of dispersions (BCV and quasi-likelihood dispersion plots)

### 1. CREACIÓN DE LA MATRIZ DE DISEÑO

## 1.1. Establecer la relación de agrupación de las muestras, sin intercepto
design1 <- model.matrix(~ 0 + group1) 

## 1.2. Asociar nombres de la columnas a los grupos
colnames(design1) <- levels(group1) 

### 2. ESTIMACIÓN DISPERSIONES (common, trended y tagwise)

## 2.1. Calcular los tres tipos de dispersión
dge1 <- estimateDisp(dge1, design1, robust = TRUE) 

## 2.2. Repesentar gráficamente la variabilidad estimada
plotBCV(dge1)

### 3. AJUSTE MODELO GLM CON MÉTODO CUASI-VEROSIMILITUD 

## 3.1. Ajustar los datos para minimizar el impacto de genes altamente variables
fit1 <- glmQLFit(dge1, design = design1, robust = TRUE)

## 3.2. Representar gráficamente el ajuste
plotQLDisp(fit1)

#### PRUEBA DE SIGNIFICANCIA (F-TEST) ####

### 1. DEFINICIÓN DEL CONTRASTE: Female vs. Male HIIT PostTraining

## 1.1. Definir grupos a comparar
F_vs_M.H.Post <- makeContrasts(Female.HIIT.Post-Male.HIIT.Post, 
                              levels = design1)
F_vs_M.R.Post <- makeContrasts(Female.Resistance.Post-Male.Resistance.Post, 
                              levels = design1)
F_vs_M.C.Post <- makeContrasts(Female.Combined.Post-Male.Combined.Post, 
                              levels = design1)
F_vs_M.H.Pre <- makeContrasts(Female.HIIT.Pre-Male.HIIT.Pre, 
                              levels = design1)
F_vs_M.R.Pre <- makeContrasts(Female.Resistance.Pre-Male.Resistance.Pre, 
                              levels = design1)
F_vs_M.C.Pre <- makeContrasts(Female.Combined.Pre-Male.Combined.Pre, 
                              levels = design1)

### 2. ANÁLISIS DE EXPRESIÓN DIFERENCIAL

## 2.1. Aplicar el F-test del método de cuasi-verosimilitud (QLF)
res1_H.Post <- glmQLFTest(fit1, contrast = F_vs_M.H.Post)
res1_R.Post <- glmQLFTest(fit1, contrast = F_vs_M.R.Post)
res1_C.Post <- glmQLFTest(fit1, contrast = F_vs_M.C.Post)
res1_H.Pre <- glmQLFTest(fit1, contrast = F_vs_M.H.Pre)
res1_R.Pre <- glmQLFTest(fit1, contrast = F_vs_M.R.Pre)
res1_C.Pre <- glmQLFTest(fit1, contrast = F_vs_M.C.Pre)

## 2.2. Corregir el testeo múltiple usando el método Benjamini-Hochberg (BH) 
res1_corrected_H.Post <- topTags(res1_H.Post, n = Inf) 
res1_corrected_R.Post <- topTags(res1_R.Post, n = Inf)
res1_corrected_C.Post <- topTags(res1_C.Post, n = Inf)
res1_corrected_H.Pre <- topTags(res1_H.Pre, n = Inf) 
res1_corrected_R.Pre <- topTags(res1_R.Pre, n = Inf)
res1_corrected_C.Pre <- topTags(res1_C.Pre, n = Inf)

## 2.3. Filtrar genes por significancia estadística (FDR) y log Fold Change (FC)
is.de1_H.Post <- decideTests(res1_H.Post, adjust.method = "BH", p.value = 0.05, lfc = 0.5)
is.de1_R.Post <- decideTests(res1_R.Post, adjust.method = "BH", p.value = 0.05, lfc = 0.5)
is.de1_C.Post <- decideTests(res1_C.Post, adjust.method = "BH", p.value = 0.05, lfc = 0.5)
is.de1_H.Pre <- decideTests(res1_H.Pre, adjust.method = "BH", p.value = 0.05, lfc = 0.5)
is.de1_R.Pre <- decideTests(res1_R.Pre, adjust.method = "BH", p.value = 0.05, lfc = 0.5)
is.de1_C.Pre <- decideTests(res1_C.Pre, adjust.method = "BH", p.value = 0.05, lfc = 0.5)

summary(is.de1_H.Post) # Down: 86, NotSig: 17150, Up: 48 
summary(is.de1_R.Post) # Down: 18, NotSig: 17253, Up: 13
summary(is.de1_C.Post) # Down: 15, NotSig: 17262, Up: 7
summary(is.de1_H.Pre) # Down: 16, NotSig: 17259, Up: 9
summary(is.de1_R.Pre) # Down: 52, NotSig: 17188, Up: 44
summary(is.de1_C.Pre) # Down: 17, NotSig: 17263, Up: 4

#### VISUALIZACIÓN DE RESULTADOS DEL ANÁLISIS DIFERENCIAL ####

### 1. VOLCANO PLOT 
# Gráfico que representa la relación entre la magnitud de cambio de 
# expresión génica (logFC) y su significancia estadística (-log10(FDR))

## 1.1. Preparación de datos para representación (up, down, notsign)

# Recoger resumen estadístico en nuevo objeto data
data1_H.Post <- res1_corrected_H.Post$table
data1_R.Post <- res1_corrected_R.Post$table
data1_C.Post <- res1_corrected_C.Post$table
data1_H.Pre <- res1_corrected_H.Pre$table
data1_R.Pre <- res1_corrected_R.Pre$table
data1_C.Pre <- res1_corrected_C.Pre$table

# Crear nueva columna para diferenciar genes infra y sobreexpresados y del chrY
data1_H.Post$DE = "No significant"    
data1_H.Post$DE[data1_H.Post$logFC > 0.5 & data1_H.Post$FDR < 0.05] = "Upregulated"     
data1_H.Post$DE[data1_H.Post$logFC < -0.5 & data1_H.Post$FDR < 0.05] = "Downregulated"
data1_H.Post$DE[abs(data1_H.Post$logFC) > 0.5 & data1_H.Post$FDR < 0.05
                & data1_H.Post$Chr == "chrY"] = "chrY genes"

data1_R.Post$DE = "No significant"    
data1_R.Post$DE[data1_R.Post$logFC > 0.5 & data1_R.Post$FDR < 0.05] = "Upregulated"     
data1_R.Post$DE[data1_R.Post$logFC < -0.5 & data1_R.Post$FDR < 0.05] = "Downregulated"
data1_R.Post$DE[abs(data1_R.Post$logFC) > 0.5 & data1_R.Post$FDR < 0.05
                & data1_R.Post$Chr == "chrY"] = "chrY genes"

data1_C.Post$DE = "No significant"    
data1_C.Post$DE[data1_C.Post$logFC > 0.5 & data1_C.Post$FDR < 0.05] = "Upregulated"     
data1_C.Post$DE[data1_C.Post$logFC < -0.5 & data1_C.Post$FDR < 0.05] = "Downregulated"
data1_C.Post$DE[abs(data1_C.Post$logFC) > 0.5 & data1_C.Post$FDR < 0.05
                & data1_C.Post$Chr == "chrY"] = "chrY genes"

data1_H.Pre$DE = "No significant"    
data1_H.Pre$DE[data1_H.Pre$logFC > 0.5 & data1_H.Pre$FDR < 0.05] = "Upregulated"     
data1_H.Pre$DE[data1_H.Pre$logFC < -0.5 & data1_H.Pre$FDR < 0.05] = "Downregulated"
data1_H.Pre$DE[abs(data1_H.Pre$logFC) > 0.5 & data1_H.Pre$FDR < 0.05
               & data1_H.Pre$Chr == "chrY"] = "chrY genes"

data1_R.Pre$DE = "No significant"    
data1_R.Pre$DE[data1_R.Pre$logFC > 0.5 & data1_R.Pre$FDR < 0.05] = "Upregulated"     
data1_R.Pre$DE[data1_R.Pre$logFC < -0.5 & data1_R.Pre$FDR < 0.05] = "Downregulated"
data1_R.Pre$DE[abs(data1_R.Pre$logFC) > 0.5 & data1_R.Pre$FDR < 0.05
               & data1_R.Pre$Chr == "chrY"] = "chrY genes"

data1_C.Pre$DE = "No significant"    
data1_C.Pre$DE[data1_C.Pre$logFC > 0.5 & data1_C.Pre$FDR < 0.05] = "Upregulated"     
data1_C.Pre$DE[data1_C.Pre$logFC < -0.5 & data1_C.Pre$FDR < 0.05] = "Downregulated"
data1_C.Pre$DE[abs(data1_C.Pre$logFC) > 0.5 & data1_C.Pre$FDR < 0.05
               & data1_C.Pre$Chr == "chrY"] = "chrY genes"

table(data1_H.Post$DE) # chrY genes: 14, Down: 72, NotSig: 17150, Up: 48 
table(data1_R.Post$DE) # chrY genes: 12, Down: 6, NotSig: 17253, Up: 13 
table(data1_C.Post$DE) # chrY genes: 12, Down: 3, NotSig: 17262, Up: 7
table(data1_H.Pre$DE) # chrY genes: 12, Down: 4, NotSig: 17259, Up: 9 
table(data1_R.Pre$DE) # chrY genes: 14, Down: 38, NotSig: 17188, Up: 44
table(data1_C.Pre$DE) # chrY genes: 13, Down: 4, NotSig: 17263, Up: 4

# Subset de los DEGs más significativos para etiquetarlos en el plot
top_genes1_H.Post <- data1_H.Post %>% 
  filter(DE != "chrY genes" & DE != "No significant" &
         abs(logFC) > 0.5 & FDR < 0.01) %>%
  dplyr::select(SYMBOL)

top_genes1_R.Post <- data1_R.Post %>% 
  filter(DE != "chrY genes" & DE != "No significant" &
         abs(logFC) > 0.5 & FDR < 0.01) %>%
  dplyr::select(SYMBOL)

top_genes1_C.Post <- data1_C.Post %>% 
  filter(DE != "chrY genes" & DE != "No significant" &
         abs(logFC) > 0.5 & FDR < 0.01) %>%
  dplyr::select(SYMBOL)

top_genes1_H.Pre <- data1_H.Pre %>% 
  filter(DE != "chrY genes" & DE != "No significant" &
           abs(logFC) > 0.5 & FDR < 0.01) %>%
  dplyr::select(SYMBOL)

top_genes1_R.Pre <- data1_R.Pre %>% 
  filter(DE != "chrY genes" & DE != "No significant" &
           abs(logFC) > 0.5 & FDR < 0.01) %>%
  dplyr::select(SYMBOL)

top_genes1_C.Pre <- data1_C.Pre %>% 
  filter(DE != "chrY genes" & DE != "No significant" &
           abs(logFC) > 0.5 & FDR < 0.01) %>%
  dplyr::select(SYMBOL)

## 1.2. Creación y representación del Volcano plot

# Crear gráfico de Volcano Plot (parte inferior)
p1_H.Post_bottom <- ggplot(data1_H.Post,aes(x=logFC, y=-log10(FDR), col=DE)) + 
  geom_point(size=2) + 
  theme_classic() +
  geom_text_repel(data=data1_H.Post[data1_H.Post$SYMBOL %in% top_genes1_H.Post$SYMBOL,], 
                  aes(label=top_genes1_H.Post$SYMBOL), 
                  vjust=0.5, 
                  hjust=0.5, 
                  size=3.5, 
                  color="black",
                  fontface="bold") +
  scale_y_continuous(limits = c(0, 5)) +
  labs(x = "Mujeres vs Hombres Post HIIT (logFC)") +
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold",size = 16),
        axis.title.y=element_text(face="bold"))

p1_R.Post_bottom <- ggplot(data1_R.Post,aes(x=logFC, y=-log10(FDR), col=DE)) + 
  geom_point(size=2) + 
  theme_classic() +
  geom_text_repel(data=data1_R.Post[data1_R.Post$SYMBOL %in% top_genes1_R.Post$SYMBOL,], 
                  aes(label=top_genes1_R.Post$SYMBOL), 
                  vjust=-0.5, 
                  hjust=0.5, 
                  size=3.5, 
                  color="black",
                  fontface="bold") +
  scale_y_continuous(limits = c(0, 5)) +
  labs(x = "Mujeres vs Hombres Post Fuerza (logFC)") +
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold",size = 16),
        axis.title.y=element_text(face="bold"))

p1_C.Post_bottom <- ggplot(data1_C.Post,aes(x=logFC, y=-log10(FDR), col=DE)) + 
  geom_point(size=2) + 
  theme_classic() +
  geom_text_repel(data=data1_C.Post[data1_C.Post$SYMBOL %in% top_genes1_C.Post$SYMBOL,], 
                  aes(label=top_genes1_C.Post$SYMBOL), 
                  vjust=-0.5, 
                  hjust=0.5, 
                  size=3.5, 
                  color="black",
                  fontface="bold") +
  scale_y_continuous(limits=c(0, 5)) +
  labs(x="Mujeres vs Hombres Post Combinado (logFC)") +
  guides(color=guide_legend(title = NULL)) +
  theme(legend.position.inside=c(0.85, 0.65),
        legend.text=element_text(face="bold",size = 12),
        axis.title.x=element_text(face="bold",size = 16),
        axis.title.y=element_text(face="bold"))

p1_H.Pre_bottom <- ggplot(data1_H.Pre,aes(x=logFC, y=-log10(FDR), col=DE)) + 
  geom_point(size=2) + 
  theme_classic() +
  geom_text_repel(data=data1_H.Pre[data1_H.Pre$SYMBOL %in% top_genes1_H.Pre$SYMBOL,], 
                  aes(label=top_genes1_H.Pre$SYMBOL), 
                  vjust=0.5, 
                  hjust=0.5, 
                  size=3.5, 
                  color="black",
                  fontface="bold") +
  scale_y_continuous(limits = c(0, 5)) +
  labs(x = "Mujeres vs Hombres Pre HIIT (logFC)") +
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold",size = 16),
        axis.title.y=element_text(face="bold"))

p1_R.Pre_bottom <- ggplot(data1_R.Pre,aes(x=logFC, y=-log10(FDR), col=DE)) + 
  geom_point(size=2) + 
  theme_classic() +
  geom_text_repel(data=data1_R.Pre[data1_R.Pre$SYMBOL %in% top_genes1_R.Pre$SYMBOL,], 
                  aes(label=top_genes1_R.Pre$SYMBOL), 
                  vjust=-0.5, 
                  hjust=0.5, 
                  size=3.5, 
                  color="black",
                  fontface="bold") +
  scale_y_continuous(limits = c(0, 5)) +
  labs(x = "Mujeres vs Hombres Pre Fuerza (logFC)") +
  theme(legend.position = "none",
        axis.title.x=element_text(face="bold",size = 16),
        axis.title.y=element_text(face="bold"))

p1_C.Pre_bottom <- ggplot(data1_C.Pre,aes(x=logFC, y=-log10(FDR), col=DE)) + 
  geom_point(size=2) + 
  theme_classic() +
  geom_text_repel(data=data1_C.Pre[data1_C.Pre$SYMBOL %in% top_genes1_C.Pre$SYMBOL,], 
                  aes(label=top_genes1_C.Pre$SYMBOL), 
                  vjust=-0.5, 
                  hjust=0.5, 
                  size=3.5, 
                  color="black",
                  fontface="bold") +
  scale_y_continuous(limits=c(0, 5)) +
  labs(x="Mujeres vs Hombres Pre Combinado (logFC)") +
  guides(color=guide_legend(title = NULL)) +
  theme(legend.position.inside=c(0.85, 0.65),
        legend.text=element_text(face="bold",size = 12),
        axis.title.x=element_text(face="bold",size = 16),
        axis.title.y=element_text(face="bold"))

# Crear gráfico Volcano Plot (parte superior)
p1_H.Post_top <- ggplot(data1_H.Post,aes(x=logFC, y=-log10(FDR), col=DE)) + 
  geom_point(size=2) + 
  theme_classic() +
  geom_text_repel(data=data1_H.Post[data1_H.Post$SYMBOL %in% top_genes1_H.Post$SYMBOL,], 
                  aes(label=top_genes1_H.Post$SYMBOL), 
                  vjust=-0.5, 
                  hjust=0.5, 
                  size=3.5, 
                  color="black",
                  fontface="bold") +
  scale_y_continuous(limits = c(5, 20)) +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank())

p1_R.Post_top <- ggplot(data1_R.Post,aes(x=logFC, y=-log10(FDR), col=DE)) + 
  geom_point(size=2) + 
  theme_classic() +
  geom_text_repel(data=data1_R.Post[data1_R.Post$SYMBOL %in% top_genes1_R.Post$SYMBOL,], 
                  aes(label=top_genes1_R.Post$SYMBOL), 
                  vjust=-0.5, 
                  hjust=0.5, 
                  size=3.5, 
                  color="black",
                  fontface="bold") +
  scale_y_continuous(limits = c(5, 20)) +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank())

p1_C.Post_top <- ggplot(data1_C.Post,aes(x=logFC, y=-log10(FDR), col=DE)) + 
  geom_point(size=2) + 
  theme_classic() +
  geom_text_repel(data=data1_C.Post[data1_C.Post$SYMBOL %in% top_genes1_C.Post$SYMBOL,], 
                  aes(label=top_genes1_C.Post$SYMBOL), 
                  vjust=-0.5, 
                  hjust=0.5, 
                  size=3.5, 
                  color="black",
                  fontface="bold") +
  scale_y_continuous(limits = c(5, 20)) +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank()) 

p1_H.Pre_top <- ggplot(data1_H.Pre,aes(x=logFC, y=-log10(FDR), col=DE)) + 
  geom_point(size=2) + 
  theme_classic() +
  geom_text_repel(data=data1_H.Pre[data1_H.Pre$SYMBOL %in% top_genes1_H.Pre$SYMBOL,], 
                  aes(label=top_genes1_H.Pre$SYMBOL), 
                  vjust=-0.5, 
                  hjust=0.5, 
                  size=3.5, 
                  color="black",
                  fontface="bold") +
  scale_y_continuous(limits = c(5, 20)) +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank())

p1_R.Pre_top <- ggplot(data1_R.Pre,aes(x=logFC, y=-log10(FDR), col=DE)) + 
  geom_point(size=2) + 
  theme_classic() +
  geom_text_repel(data=data1_R.Pre[data1_R.Pre$SYMBOL %in% top_genes1_R.Pre$SYMBOL,], 
                  aes(label=top_genes1_R.Pre$SYMBOL), 
                  vjust=-0.5, 
                  hjust=0.5, 
                  size=3.5, 
                  color="black",
                  fontface="bold") +
  scale_y_continuous(limits = c(5, 20)) +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank())

p1_C.Pre_top <- ggplot(data1_C.Pre,aes(x=logFC, y=-log10(FDR), col=DE)) + 
  geom_point(size=2) + 
  theme_classic() +
  geom_text_repel(data=data1_C.Pre[data1_C.Pre$SYMBOL %in% top_genes1_C.Pre$SYMBOL,], 
                  aes(label=top_genes1_C.Pre$SYMBOL), 
                  vjust=-0.5, 
                  hjust=0.5, 
                  size=3.5, 
                  color="black",
                  fontface="bold") +
  scale_y_continuous(limits = c(5, 20)) +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank()) 

# Guardar cada Volcano Plot de forma individual (parte inferior + superior)
p1_H.Post <- p1_H.Post_top/p1_H.Post_bottom + plot_layout(heights = c(0.2, 0.8))
p1_R.Post <- p1_R.Post_top/p1_R.Post_bottom + plot_layout(heights = c(0.2, 0.8))
p1_C.Post <- p1_C.Post_top/p1_C.Post_bottom + plot_layout(heights = c(0.2, 0.8))
p1_H.Pre <- p1_H.Pre_top/p1_H.Pre_bottom + plot_layout(heights = c(0.2, 0.8))
p1_R.Pre <- p1_R.Pre_top/p1_R.Pre_bottom + plot_layout(heights = c(0.2, 0.8))
p1_C.Pre <- p1_C.Pre_top/p1_C.Pre_bottom + plot_layout(heights = c(0.2, 0.8))

# Fusionar los Volcano Plot guardados en uno solo (como tres columnas)
p1_HRC.Post <- p1_H.Post | p1_R.Post | p1_C.Post
p1_HRC.Pre <- p1_H.Pre | p1_R.Pre | p1_C.Pre

# Crear un archivo .pdf con el plot final
ggsave("Plots/p1_HRC.Post.pdf", plot = p1_HRC.Post, width = 18, height = 8)
ggsave("Plots/p1_HRC.Pre.pdf", plot = p1_HRC.Pre, width = 18, height = 8)

#### DIAGRAMAS DE VENN ####

### 1. Extraer DEGs (por sexo, training y timepoint)
## 1.1. Post-entrenamiento
degs1_H.Post.F <- data1_H.Post[data1_H.Post$DE == "Upregulated",] %>% pull(SYMBOL)
degs1_H.Post.M <- data1_H.Post[data1_H.Post$DE == "Downregulated" | data1_H.Post$DE == "chrY genes",] %>% pull(SYMBOL)
degs1_R.Post.F <- data1_R.Post[data1_R.Post$DE == "Upregulated",] %>% pull(SYMBOL)
degs1_R.Post.M <- data1_R.Post[data1_R.Post$DE == "Downregulated" | data1_R.Post$DE == "chrY genes",] %>% pull(SYMBOL)
degs1_C.Post.F <- data1_C.Post[data1_C.Post$DE == "Upregulated",] %>% pull(SYMBOL)
degs1_C.Post.M <- data1_C.Post[data1_C.Post$DE == "Downregulated" | data1_C.Post$DE == "chrY genes",] %>% pull(SYMBOL)

## 1.2. Pre-entrenamiento
degs1_H.Pre.F <- data1_H.Pre[data1_H.Pre$DE == "Upregulated",] %>% pull(SYMBOL)
degs1_H.Pre.M <- data1_H.Pre[data1_H.Pre$DE == "Downregulated" | data1_H.Pre$DE == "chrY genes",] %>% pull(SYMBOL) 
degs1_R.Pre.F <- data1_R.Pre[data1_R.Pre$DE == "Upregulated",] %>% pull(SYMBOL)
degs1_R.Pre.M <- data1_R.Pre[data1_R.Pre$DE == "Downregulated" | data1_R.Pre$DE == "chrY genes",] %>% pull(SYMBOL) 
degs1_C.Pre.F <- data1_C.Pre[data1_C.Pre$DE == "Upregulated",] %>% pull(SYMBOL)
degs1_C.Pre.M <- data1_C.Pre[data1_C.Pre$DE == "Downregulated" | data1_C.Pre$DE == "chrY genes",] %>% pull(SYMBOL) 

### 2. Extraer DEGs únicos para post-entrenamiento
degs1_uniq_H.Post.F <- setdiff(degs1_H.Post.F, degs1_H.Pre.F)
degs1_uniq_R.Post.F <- setdiff(degs1_R.Post.F, degs1_R.Pre.F)
degs1_uniq_C.Post.F <- setdiff(degs1_C.Post.F, degs1_C.Pre.F)
degs1_uniq_H.Post.M <- setdiff(degs1_H.Post.M, degs1_H.Pre.M)
degs1_uniq_R.Post.M <- setdiff(degs1_R.Post.M, degs1_R.Pre.M)
degs1_uniq_C.Post.M <- setdiff(degs1_C.Post.M, degs1_C.Pre.M)

### 3. Representar Diagrama de Venn 
## 3.1. PosTraining (DEGs únicos)
draw.venn(degs1_uniq_H.Post.F,
          degs1_uniq_R.Post.F,
          degs1_uniq_C.Post.F,
          title = "Mujeres",
          subtitle = NULL,
          xtitle = NULL,
          ytitle = NULL,
          ztitle = NULL,
          x_c = "#7F3C8D", 
          y_c = "#E03C31", 
          z_c = "#56B4E9",
          output = "pdf",
          filename = "Plots/biovenn1_uniq_Post_F")

draw.venn(degs1_uniq_H.Post.M,
          degs1_uniq_R.Post.M,
          degs1_uniq_C.Post.M,
          title = "Hombres",
          subtitle = NULL,
          xtitle = NULL,
          ytitle = NULL,
          ztitle = NULL,
          x_c = "#00BFA0", 
          y_c = "#F0E442", 
          z_c = "#CC79A7",
          output = "pdf",
          filename = "Plots/biovenn1_uniq_Post_M")

## 3.2. PosTraining (todos los DEGs)
draw.venn(degs1_H.Post.F,
          degs1_R.Post.F,
          degs1_C.Post.F,
          title = "Mujeres",
          subtitle = NULL,
          xtitle = NULL,
          ytitle = NULL,
          ztitle = NULL,
          x_c = "#7F3C8D", 
          y_c = "#E03C31", 
          z_c = "#56B4E9",
          output = "pdf",
          filename = "Plots/biovenn1_Post_F")

draw.venn(degs1_H.Post.M,
          degs1_R.Post.M,
          degs1_C.Post.M,
          title = "Hombres",
          subtitle = NULL,
          xtitle = NULL,
          ytitle = NULL,
          ztitle = NULL,
          x_c = "#00BFA0", 
          y_c = "#F0E442", 
          z_c = "#CC79A7",
          output = "pdf",
          filename = "Plots/biovenn1_Post_M")

## 3.3. PreTraining (todos los DEGs)
draw.venn(degs1_H.Pre.F,
          degs1_R.Pre.F,
          degs1_C.Pre.F,
          title = "Mujeres",
          subtitle = NULL,
          xtitle = NULL,
          ytitle = NULL,
          ztitle = NULL,
          x_c = "#7F3C8D", 
          y_c = "#E03C31", 
          z_c = "#56B4E9",
          output = "pdf",
          filename = "Plots/biovenn1_Pre_F")

draw.venn(degs1_H.Pre.M,
          degs1_R.Pre.M,
          degs1_C.Pre.M,
          title = "Hombres",
          subtitle = NULL,
          xtitle = NULL,
          ytitle = NULL,
          ztitle = NULL,
          x_c = "#00BFA0", 
          y_c = "#F0E442", 
          z_c = "#CC79A7",
          output = "pdf",
          filename = "Plots/biovenn1_Pre_M")

#### HEATMAP ####

# 1. NORMALIZACIÓN logCPM DE LOS DATOS DGEList (PostTraining)
# Normalizar conteos de genes de la matriz de conteos PostTraining
logCPM1 = cpm(dge1_Post, log = TRUE) 

# Sustituir nombres filas (entrez id por símbolo)
rownames(logCPM1) = dge1_Post$genes$SYMBOL 

# Sustituir nombres columnas (id participante por grupo numerado)
colnames(logCPM1) <- paste0(dge1_Post$samples$group, 
                            "_", 
                            ave(seq_along(dge1_Post$samples$group), 
                                dge1_Post$samples$group, 
                                FUN = seq_along))

# 2. FILTRADO DE DATOS POR DEGs
degs1_uniq_Post = unique(c(degs1_uniq_H.Post.F,degs1_uniq_H.Post.M,
                           degs1_uniq_R.Post.F,degs1_uniq_R.Post.M,
                           degs1_uniq_C.Post.F,degs1_uniq_C.Post.M)) #121 genes

heatmap1_Post = logCPM1[degs1_uniq_Post,]

# 3. VISUALIZACIÓN DEL HEATMAP

pheatmap(heatmap1_Post, 
         scale = "row",
         cluster_rows = T,
         cluster_cols = T,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D2",
         cutree_rows = 2,
         cutree_cols = 3,
         fontsize_row = 6,
         show_rownames = FALSE,
         display_numbers = F) #escalamos las muestras

#### ANÁLISIS DE ONTOLOGÍA GÉNICA (GO) ####

### 1. Realizar análisis de enriquecimiento de Ontología Génica (GO) 

## 1.1. Análisis GO a partir del ENTREZIDs 
go1_uniq_H.Post.F <- goana(res1_H.Post$genes$ENTREZID[res1_H.Post$genes$SYMBOL %in% degs1_uniq_H.Post.F], species="Hs")
go1_uniq_H.Post.M <- goana(res1_H.Post$genes$ENTREZID[res1_H.Post$genes$SYMBOL %in% degs1_uniq_H.Post.M], species="Hs")
go1_uniq_R.Post.F <- goana(res1_R.Post$genes$ENTREZID[res1_R.Post$genes$SYMBOL %in% degs1_uniq_R.Post.F], species="Hs")
go1_uniq_R.Post.M <- goana(res1_R.Post$genes$ENTREZID[res1_R.Post$genes$SYMBOL %in% degs1_uniq_R.Post.M], species="Hs")
go1_uniq_C.Post.F <- goana(res1_C.Post$genes$ENTREZID[res1_C.Post$genes$SYMBOL %in% degs1_uniq_C.Post.F], species="Hs")
go1_uniq_C.Post.M <- goana(res1_C.Post$genes$ENTREZID[res1_C.Post$genes$SYMBOL %in% degs1_uniq_C.Post.M], species="Hs")

## 1.2. Extraer los términos GO más relevantes en data frames
topgo1_df_uniq_H.Post.F = topGO(go1_uniq_H.Post.F, number = 10, ontology = "BP") %>%
                       as.data.frame()
topgo1_df_uniq_H.Post.M = topGO(go1_uniq_H.Post.M, number = 10, ontology = "BP") %>%
                       as.data.frame()
topgo1_df_uniq_R.Post.F = topGO(go1_uniq_R.Post.F, number = 10, ontology = "BP") %>%
                       as.data.frame()
topgo1_df_uniq_R.Post.M = topGO(go1_uniq_R.Post.M, number = 10, ontology = "BP") %>%
                       as.data.frame()
topgo1_df_uniq_C.Post.F = topGO(go1_uniq_C.Post.F, number = 10, ontology = "BP") %>%
                       as.data.frame()
topgo1_df_uniq_C.Post.M = topGO(go1_uniq_C.Post.M, number = 10, ontology = "BP") %>%
                       as.data.frame()

# 1.3. Guardar cada data frame en formato .csv
write.csv(topgo1_df_uniq_H.Post.F, "GO/GO1_H_Post_F.csv", row.names = FALSE)
write.csv(topgo1_df_uniq_H.Post.M, "GO/GO1_H_Post_M.csv", row.names = FALSE)
write.csv(topgo1_df_uniq_R.Post.F, "GO/GO1_R_Post_F.csv", row.names = FALSE)
write.csv(topgo1_df_uniq_R.Post.M, "GO/GO1_R_Post_M.csv", row.names = FALSE)
write.csv(topgo1_df_uniq_C.Post.F, "GO/GO1_C_Post_F.csv", row.names = FALSE)
write.csv(topgo1_df_uniq_C.Post.M, "GO/GO1_C_Post_M.csv", row.names = FALSE)

