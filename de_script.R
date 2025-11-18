
# ------------------------------------------------------------------------------
# 0) Preparación del entorno
#    - BiocManager instala paquetes de Bioconductor (como DESeq2)
#    - library(DESeq2) carga el paquete para usar sus funciones
# ------------------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")     # se instala solo si no está disponible

BiocManager::install("DESeq2")        # instala DESeq2 (puede tardar)
library(DESeq2)                       # carga DESeq2 para este script


# ------------------------------------------------------------------------------
# 1) Descarga de datos del experimento COVID-19 (E-ENAD-46) - VERSIÓN CORREGIDA
# ------------------------------------------------------------------------------

# Set working directory
setwd("E:/Bioinfosofiu/DEanalysis")

# ✅ CORREGIDO: URLs correctas y nombres de archivo consistentes
counts_original <- read.delim(
  "https://www.ebi.ac.uk/gxa/experiments-content/E-ENAD-46/resources/DifferentialSecondaryDataFiles.RnaSeq/raw-counts"
)

# ✅ CORREGIDO: Metadata desde el diseño experimental, NO analytics
metadata_original <- read.delim(
  "https://www.ebi.ac.uk/gxa/experiments-content/E-ENAD-46/resources/ExperimentDesignFile.RnaSeq/experiment-design"
)

# ✅ CORREGIDO: Guardar con nombres CORRECTOS (E-ENAD-46,)
write.table(counts_original, "data/raw_counts_E-ENAD-46.tsv", 
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(metadata_original, "metadata/experiment_design_E-ENAD-46.tsv", 
            sep = "\t", row.names = FALSE, quote = FALSE)
#-------------------------------------------------------------------------------
# 1.1) Tomar solo muestras de pulmón
#-------------------------------------------------------------------------------

# 1. Filtrar metadata para solo pulmón
metadata_pulmon <- metadata_original[metadata_original$Factor.Value.organism.part. == "lung", ]
cat("Muestras de pulmón:", nrow(metadata_pulmon), "\n")
print(table(metadata_pulmon$Factor.Value.disease.))

write.table(metadata_pulmon, "metadata/experiment_design_E-ENAD-46_pulmon.tsv", 
            sep = "\t", row.names = FALSE, quote = FALSE)

# 2. Filtrar counts para solo muestras de pulmón
muestras_pulmon <-  metadata_pulmon$Run
counts <- counts_original[, colnames(counts_original) %in% muestras_pulmon]
write.table(counts, "data/raw_counts_E-ENAD-46_pulmon.tsv", 
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("Dimensión de counts pulmonar:", dim(counts), "\n")

cat("=== VERIFICACIÓN DE DESCARGA ===\n")
cat("Dimensión de conteos:", dim(counts), "\n")
cat("Dimensión de metadatos:", dim(metadata_pulmon), "\n")

# ------------------------------------------------------------------------------
# 2) Verificación rápida
# Acomodar los datos al formato que DESeq2 espera
#    - Filas de 'counts' = genes (con sus IDs en rownames)
#    - Columnas de 'counts' = muestras (los nombres deben coincidir con metadata)
#    - Rownames de 'metadata' = IDs de muestra


# 2.1) PRIMERO: Asignar los IDs como rownames a counts
rownames(counts) <- counts_original$Gene.ID  # 

# 2.2) SEGUNDO: Extraer la información de genes
genes_info <- counts_original[, c("Gene.ID", "Gene.Name")] 
cat("genes_info - dimensión:", dim(genes_info), "\n")
print(head(genes_info))


# 2.3) Eliminar de 'counts' las columnas que no son conteos (las 2 primeras)
#      (Esto deja solo números de conteo por muestra.)
#No hacer porque no están en mis datos


# 2.4) Poner los IDs de muestra como rownames en 'metadata'
#      'Run' suele ser el identificador único de cada muestra en EBI/GXA.
rownames(metadata_pulmon) <- metadata_pulmon$Run
head(metadata_pulmon)

#Encontrar columna clave del estudio
columna_enfermedad <- grep("Factor.Value.disease", colnames(metadata_pulmon), value = TRUE)
cat("Columna de enfermedad encontrada:", columna_enfermedad, "\n")

# 2.5) Quedarnos solo con la columna del factor experimental de interés
#      Aquí nos interesa el genotipo. 'drop=FALSE' asegura que siga siendo data.frame.
metadata_pulmon <- metadata_pulmon[, c("Factor.Value.disease."), drop = FALSE]
metadata_pulmon

# 2.6) Renombrar la columna a algo simple (evita errores y hace el código legible)
colnames(metadata_pulmon) <- "condicion"
metadata_pulmon

# 2.7) Limpiar etiquetas para que no haya espacios raros:
metadata_pulmon$condicion[metadata_pulmon$condicion == "COVID-19"] <- "covid"
metadata_pulmon$condicion[metadata_pulmon$condicion == "normal"] <- "control"


# 2.8) Declarar 'condicion' como factor y fijar el orden (referencia primero)
metadata_pulmon$condicion <- factor(metadata_pulmon$condicion, levels = c("control", "covid"))

cat("\n🎯 VERIFICACIÓN FINAL DEL CONTRASTE:\n")
cat("Niveles del factor:", levels(metadata_pulmon$condicion), "\n")
cat("Distribución de muestras:\n")
print(table(metadata_pulmon$condicion))

# ------------------------------------------------------------------------------
# 3) 3) Chequeo rápido: EXPRESIÓN DE ACE2 (receptor del SARS-CoV-2)
#    - Buscamos el ID del gen usando su nombre (Gene.Name == 'ACE2')
#    - Extraemos los conteos de ese gen en todas las muestras
#    - Creamos una tablita con metadatos + conteos para graficar
# ------------------------------------------------------------------------------
gene_name <- "ACE2"
gene_id <- genes_info$Gene.ID[genes_info$Gene.Name == gene_name]
cat("ID de ACE2 encontrado:", gene_id, "\n")

gene_counts <- counts[gene_id, ]                           # vector de conteos por muestra
gene_counts


# cbind: "pega" columnas. as.numeric asegura que el vector sea numérico simple.
gene_data <- cbind(metadata_pulmon, counts = as.numeric(gene_counts))
gene_data

#Box plot  por grupo (wildtype vs knockout)
library(ggplot2)

gene_data <- cbind(metadata_pulmon, counts = as.numeric(gene_counts))

ggplot(gene_data, aes(x = condicion, y = counts, fill = condicion)) +
  geom_boxplot() +
  labs(title = "Conteos crudos de ACE2 por condición",
       x = "Condición",
       y = "Conteos (sin normalizar)")

# ------------------------------------------------------------------------------
# 4) Crear el objeto DESeq y correr el análisis
#    - design = ~ genotipo indica que queremos probar diferencias por genotipo
#    - Filtramos genes muy poco expresados (ruido) para evitar falsos positivos
# ------------------------------------------------------------------------------

#Aquí estamos buscando el patrón objetivo que queremos que analice
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData   = metadata_pulmon,  
                              design    = ~ condicion)

# Filtrar genes con baja suma de conteos (aquí umbral > 10 es un ejemplo simple)
dds <- dds[rowSums(counts(dds)) > 10, ]

# Ejecuta la estimación de tamaños, dispersión y el modelo (Wald por defecto)
dds <- DESeq(dds)

# ------------------------------------------------------------------------------
# 5) Extraer resultados del contraste de interés
# ------------------------------------------------------------------------------
res <- results(dds, contrast = c("condicion", "covid", "control"), alpha = 1e-5)
cat("📊 Resultados obtenidos para", nrow(results), "genes\n")
head(res)

#Ahora tomamos todos los genes significativos con valor p por comparar
sig_genes_p <- res[!is.na(res$padj) & res$pvalue < 0.05 & abs(res$log2FoldChange) > 0.5, ]
cat("GENES SIGNIFICATIVOS (pvalue < 0.05 & |log2FC| > 0.5):", nrow(sig_genes_p), "\n")
head (sig_genes_p)

#Ahora tomamos todos los genes significativos con padj que son los que queremos
sig_genes <- res[!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) > 0.5, ]
cat("GENES SIGNIFICATIVOS (pvalue < 0.05 & |log2FC| > 0.5):", nrow(sig_genes), "\n")
head (sig_genes)


# VERIFICACIÓN genes con val p: Esto incluye UP y DOWN
if(nrow(sig_genes_p) > 0) {
  genes_up <- sum(sig_genes_p$log2FoldChange > 0)
  genes_down <- sum(sig_genes_p$log2FoldChange < 0)
  
  cat("📈 Genes UP regulados:", genes_up, "\n")
  cat("📉 Genes DOWN regulados:", genes_down, "\n")
  cat("📊 Total genes significativos:", genes_up + genes_down, "\n")
}


# VERIFICACIÓN genes con p adj: Esto incluye UP y DOWN
if(nrow(sig_genes) > 0) {
  genes_up <- sum(sig_genes$log2FoldChange > 0)
  genes_down <- sum(sig_genes$log2FoldChange < 0)
  
  cat("📈 Genes UP regulados:", genes_up, "\n")
  cat("📉 Genes DOWN regulados:", genes_down, "\n")
  cat("📊 Total genes significativos:", genes_up + genes_down, "\n")
}

# ------------------------------------------------------------------------------
# 7) Spot checks: comparar algunos genes con la interfaz de GXA
#    - Convertimos 'res' a data.frame y unimos los nombres de genes
#    - Buscamos genes de interés y revisamos sus resultados
#    - De aquí determinamos en cuantos genes se tuvo un efecto y cuales son significativos
# ------------------------------------------------------------------------------
# 1. genes_info debe contener SOLO los genes que están en res
genes_info_filtrado <- genes_info[genes_info$Gene.ID %in% rownames(res), ]
cat("genes_info filtrado - dimensión:", dim(genes_info_filtrado), "\n")
cat("genes_info original - dimensión:", dim(genes_info), "\n")

# 2. Ahora hacer el merge con la versión filtrada
res_df <- as.data.frame(res)
res_df$Gene.ID <- rownames(res_df)

res_df <- merge(res_df, genes_info_filtrado, by = "Gene.ID")
cat("res_df después del merge CORREGIDO - dimensión:", dim(res_df), "\n")

print(head(res_df))

write.table(res_df, "results/ DE_results.csv", 
            sep = "\t", row.names = FALSE, quote = FALSE)

#Ahora tomamos todos los genes significativos pero con nombres
total_genes_verificados <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 0.5, ]
head(total_genes_verificados)

#Ahora vamos a crear la tabla de resultados
# Separar en UP y DOWN regulados
genes_up <- total_genes_verificados[total_genes_verificados$log2FoldChange > 0, ]
genes_down <- total_genes_verificados[total_genes_verificados$log2FoldChange < 0, ]

cat("Genes UP regulados en COVID:", nrow(genes_up), "\n")
cat("Genes DOWN regulados en COVID:", nrow(genes_down), "\n")

# 4. Mostrar los TOP 10 en cada categoría
cat("\n🔝 TOP 10 GENES MÁS UP REGULADOS:\n")
top_up <- genes_up[order(genes_up$log2FoldChange, decreasing = TRUE), 
                   c("Gene.Name","Gene.ID", "log2FoldChange","pvalue", "padj")]
print(head(top_up, 10))

cat("\n🔻 TOP 10 GENES MÁS DOWN REGULADOS:\n")
top_down <- genes_down[order(genes_down$log2FoldChange, decreasing = TRUE), 
                       c("Gene.Name","Gene.ID","log2FoldChange","pvalue", "padj")]
print(head(top_down, 10))


# Crear tabla combinada solo con genes UP y DOWN
#Union de los genes up y down más significativos
genes_verificados <- rbind(
  cbind(head(top_up, 3)[, c("Gene.Name","Gene.ID","log2FoldChange", "pvalue", "padj")], Regulacion = "UP"),
  cbind(head(top_down,3)[, c("Gene.Name","Gene.ID", "log2FoldChange", "pvalue", "padj")], Regulacion = "DOWN")
)

write.table(genes_verificados, "results/genes_verificados_significancia.csv", 
            sep = "\t", row.names = FALSE, quote = FALSE)

# Mostrar resultados
cat("TABLA COMBINADA - GENES UP Y DOWN:\n")
print(genes_verificados)
cat("Total en tabla combinada:", nrow(genes_verificados), "\n")


# Genes COVID-19 relevantes para verificar
genes_a_verificar <- c("ACE2", "TMPRSS2","IFNG")
genes_a_verificar_1 <- res_df[res_df$Gene.Name %in% genes_a_verificar, ]
print(genes_a_verificar_1)

write.table(genes_a_verificar_1, "results/genes_verificados_deinteres.csv", 
            sep = "\t", row.names = FALSE, quote = FALSE)

# Interpretación rápida de los genes top
cat("\nINTERPRETACION:\n")
for(i in 1:nrow(genes_verificados)) {
  gen <- genes_verificados$Gene.Name[i]
  lfc <- genes_verificados$log2FoldChange[i]
  padj <- genes_verificados$padj[i]
  
  if(!is.na(padj) & padj < 0.05) {
    direccion <- ifelse(lfc > 0, "UP en COVID", "DOWN en COVID")
    cat("-", gen, ":", direccion, "(log2FC =", round(lfc, 2), ")\n")
  } else {
    cat("-", gen, ": No significativo\n")
  }
}

# Interpretación rápida de los genes de interes a verificar

for(i in 1:nrow(genes_a_verificar_1)) {
  gen <- genes_a_verificar_1$Gene.Name[i]
  lfc <- genes_a_verificar_1$log2FoldChange[i]
  padj <- genes_a_verificar_1$padj[i]
  
  if(!is.na(padj) & padj < 0.05) {
    direccion <- ifelse(lfc > 0, "UP en COVID", "DOWN en COVID")
    cat("-", gen, ":", direccion, "(log2FC =", round(lfc, 2), ")\n")
  } else {
    cat("-", gen, ": No significativo\n")
  }
}
# ------------------------------------------------------------------------------
# 8) Visualizaciones rápidas
#    - MA plot: muestra relación entre abundancia media y cambio de expresión
#    - Volcano plot: log2FC vs significancia (p-value)
# ------------------------------------------------------------------------------
plotMA(res)  # puntos rojos suelen ser genes significativos

#Limpieza previa a hacer el volcano
res_df$unique_id <- ifelse(res_df$Gene.Name == "" | duplicated(res_df$Gene.Name),
                           paste0("gene_", 1:nrow(res_df)),
                           res_df$Gene.Name)


#Volcano plot de todos los resultados
BiocManager::install("EnhancedVolcano")

library(EnhancedVolcano)
volcano_padj <- EnhancedVolcano(res_df,
                lab = res_df$Gene.Name,  # Mostrar nombres originales
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Volcano Plot: COVID-19 vs Control')

volcano_padj
