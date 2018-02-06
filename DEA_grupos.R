# ANALISIS DIFERENCIAL DE LA EXPRESIÓN EN CÁNCER DE PROSTATA

# Cargar paquetes necesarios. Como "attach" paquetes se cargarán otros automáticamente 
library("TCGAbiolinks")
library("SummarizedExperiment")

# Resumen sobre los datos de mi proyecto disponibles en TCGA
# resumen <- getDataCategorySummary("TCGA-PRAD", legacy = FALSE)

# Buscar en TCGA datos de RNA-seq
prost_query <- GDCquery(project = "TCGA-PRAD", 
                        data.category = "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification",
                        workflow.type = "HTSeq - FPKM")

#View(getResults(prost_query))

# Obtener bcr (identificador) de las muestras para generar grupos (1.-Muestras tumor primario; 2.- Muestas de tejido soĺido normal)
prost_muestras <- getResults(prost_query, cols = c("cases"))
prost_TP <- TCGAquery_SampleTypes(barcode = prost_muestras, typesample = "TP")  # TP -> Tumor Primary --------> 498 samples
prost_NT <- TCGAquery_SampleTypes(barcode = prost_muestras, typesample = "NT")  # NT -> Solid Tissue Normal --> 52 samples


# Buscar los datos específicos de mis grupos y descargarlos.
datos_query <- GDCquery(project = "TCGA-PRAD",
                        data.category = "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification",
                        workflow.type = "HTSeq - FPKM",
                        barcode = c(prost_TP, prost_NT))

#GDCdownload(query = datos_query, directory = "/home/cone/Escritorio/Proyect_Final/data/")

# Crear objeto SummarizedExperiment y guardar una copia 
data_se <- GDCprepare(query = datos_query,
                      save = TRUE,
                      directory = "/home/cone/Escritorio/Proyect_Final/data/",
                      save.filename = "se_prostata.rda")

# A partir del objeto "data_se", podemos extraer diferentes datos:
# 1.- Matriz de expresión 
# express_matrix <- assay(data_se)
# 2.- Información clínica
# clin_inf <- colData(data_se)
# 3.- Información sobre genes
# gen_informacion <- rowRanges(data_se)


# PREPROCESADO DE LOS DATOS: Permite encontrar muestras con baja correlacion que pueden ser identificadas como outliers
# 1.-AAIC --> Plot con la intensidad array-array. Matriz cuadrada simétrica con las correlaciones de Pearson entre muestras
# 2.-Boxplot -->  de las correlaciones muestra a muestra
prost_pre <- TCGAanalyze_Preprocessing(object = data_se[,c(1:20)],
                                       cor.cut = 0.6,  # Threshold to filter samples according their spearman
                                       # correlation in samples by samples. Default cor.cut = 0.
                                       filename = "AAIC_boxplot_correlation.png",
                                       width = 1000,
                                       height = 1000,
                                       datatype = names(assays(data_se))[1])


# NORMALIZACIÓN --> con el paquete "EDASeq"
# Whithinlane normalizacion para ajustar los efectos del contenido en GC en los read counts: loess robust local regression,
# global-scaling, y full-quantile normalization.
# Between-lan normalización para ajustar las diferencias de distribuciones entre lineas: global-scaling y full-quantile normalization.
prost_norm<- TCGAanalyze_Normalization(tabDF = data_se,
                                       geneInfo = geneInfoHT, # Information matrix of 20531 genes about geneLength and gcContent
                                       method = "gcContent")  

dim(prost_norm)
boxplot(prost_pre, outline = FALSE)
boxplot(prost_norm[,1:20], outline = FALSE)


# FILTRADO --> con el paquete "geneFilter". Devuelve todos los mRNA con una media (a través de todas las muestras), superior
# al umbral definido por la media en el cuartil a traves de todas las muestras. 
prost_filt <- TCGAanalyze_Filtering(tabDF = prost_norm,
                                    method = "quantile", 
                                    qnt.cut = 0.80); dim(prost_filt)  # qnt.cut = umbral para el filtrado

# ANÁLISIS DIFERENCIAL DE LA EXPRESIÓN --> paquete edgeR. Identifica genes expresados diferencialmente entre dos grupos
prost_GED <- TCGAanalyze_DEA(mat1 = prost_filt[,prost_NT],
                             mat2 = prost_filt[,prost_TP],
                             Cond1type = "Normal",
                             Cond2type = "Tumor",
                             fdr.cut = 0.01,      # umbral para filtrar DEGs según p-value corrected
                             logFC.cut = 1,       # umbral para filtrar DEGs según su logFC
                             method = "glmLRT"); dim(prost_GED)


# AÑADIR INFORMACIÓN RELACIONADA CON DEGs (e.g. media en dos condiciones). TCGAanalyze_LevelTab 
prost_GED_inf <- TCGAanalyze_LevelTab(prost_GED,"Tumor","Normal", prost_filt[,prost_TP],prost_filt[,prost_NT]); dim(prost_GED_inf)

# ANÁLISIS DE ENRIQUECIMIENTO --> recuperar o generar un perfil funcional de un grupo de genes. Dado un grupo de genes
#                                 regulados bajo ciertas condiciones, este análisis identificará clases de genes que son 
#                                 sobrerepresentados utilizando anotaciones. Gene Ontology y Pathways

Genelist <-  prost_GED_inf$external_gene_name # grupo de genes para realizar el EA
an_enr <- TCGAanalyze_EAcomplete(TFname="DEA genes Normal Vs Tumor",Genelist)
TCGAvisualize_EAbarplot(tf = rownames(an_enr$ResBP),
                        GOBPTab = an_enr$ResBP,
                        GOCCTab = an_enr$ResCC,
                        GOMFTab = an_enr$ResMF,
                        PathTab = an_enr$ResPat,
                        nRGTab = Genelist,
                        nBar = 15)

##############################################################################
######## ANÁLISIS DIFENCIAL DE LA EXPRESIÓN (según el gleason score)  ########
##############################################################################
# 1.- GLEASON < 7 vs GLEASON = 7 (gleason bajo vs gleason medio)
# 2.- GLEASON < 7 vs GLEASON > 7 (gleason bajo vs gleason alto)
# 3.- GLEASON = 7  vs GLEASON > 7 (gleason medio vs gleason alto)
##############################################################################
##############################################################################

# 1.- GLEASON < 7 vs GLEASON = 7 (gleason bajo vs gleason medio)
# Obtenemos los barcodes para las muestras que presenten un gleason < 7
gl6 <- subset(data_se, data_se$barcode, data_se$subtype_Reviewed_Gleason_sum == 6)
gl6_lista <- as.factor(gl6$barcode)
gl6_group <- TCGAquery_SampleTypes(barcode = gl6_lista, typesample = "TP")

# Obtenemos barcodes para muestras que presenten un gleason = 7
gl7 <- subset(data_se, data_se$barcode, data_se$subtype_Reviewed_Gleason_sum == 7)
gl7_lista <- as.factor(gl7$barcode)
gl7_group <- TCGAquery_SampleTypes(barcode = gl7_lista, typesample = "TP")

query_gl6_7 <- GDCquery(project = "TCGA-PRAD",
                        data.category = "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification",
                        workflow.type = "HTSeq - FPKM",
                        barcode = c(gl6_group, gl7_group))
#GDCdownload(query = query_gl6_7,  directory = "/home/cone/Escritorio/Proyect_Final/data/")

se_gl6_7 <- GDCprepare(query = query_gl6_7,
                       save = TRUE,
                       directory = "/home/cone/Escritorio/Proyect_Final/data/",
                       save.filename = "se_gl6_7.rda")

pre_gl6_7 <- TCGAanalyze_Preprocessing(object = se_gl6_7,
                                       cor.cut = 0.6,  
                                       filename = "AAIC_gl6_7.png",
                                       width = 2500,
                                       height = 3000,
                                       datatype = names(assays(se_gl6_7))[1])

gl6_7_norm<- TCGAanalyze_Normalization(tabDF = se_gl6_7,
                                       geneInfo = geneInfoHT, # Information matrix of 20531 genes about
                                       # geneLength and gcContent
                                       method = "gcContent")  # "gcContent" or "geneLength"

boxplot(pre_gl6_7, outline = FALSE)
boxplot(gl6_7_norm, outline = FALSE)

gl6_7_filter <- TCGAanalyze_Filtering(tabDF = gl6_7_norm,
                                      method = "quantile",
                                      qnt.cut = 0.25); dim(gl6_7_filter)

gl6_7DEA <- TCGAanalyze_DEA(mat1 = gl6_7_filter[,gl6_group],
                            mat2 = gl6_7_filter[,gl7_group],
                            Cond1type = "Gleason score < 6",
                            Cond2type = "Gleason score = 7",
                            fdr.cut = 0.01,
                            logFC.cut = 1,
                            method = "glmLRT"); dim(gl6_7DEA)


##############################################################################
##############################################################################

# 2.- GLEASON < 7 vs GLEASON > 7 (gleason bajo vs gleason alto)
# Obtenemos barcodes para muestras que presenten un gleason > 7
as.factor(data_se$subtype_Reviewed_Gleason_sum)
nums = c(8,9,10)
gl_mas7 <- subset(data_se, data_se$barcode, data_se$subtype_Reviewed_Gleason_sum >= 8)
gl_mas7_lista <-as.factor(gl_mas7$barcode)
#gl6_group <- TCGAquery_SampleTypes(barcode = gl6_lista, typesample = "TP")
gl_mas7_group <- TCGAquery_SampleTypes(barcode = gl_mas7_lista, typesample = "TP")
#gl_mas7_group

query_gl_bajo_alto<- GDCquery(project = "TCGA-PRAD",
                              data.category = "Transcriptome Profiling",
                              data.type = "Gene Expression Quantification",
                              workflow.type = "HTSeq - FPKM",
                              barcode = c(gl6_group, gl_mas7_group))
#GDCdownload(query = query_gl_bajo_alto, directory = "/home/cone/Escritorio/Proyect_Final/data/")
gl_bajo_alto <- GDCprepare(query = query_gl_bajo_alto,
                           save = TRUE,
                           directory = Directorio_sub,
                           save.filename = "se_gl_bajo_alto.rda")

gl_bajo_altoAAIC <- TCGAanalyze_Preprocessing(object = gl_bajo_alto,
                                              cor.cut = 0.6,  
                                              filename = "AAICgl_bajo_alto.png",
                                              width = 2500,
                                              height = 3000,
                                              datatype = names(assays(gl_bajo_alto))[1])

gl_bajo_alto_norm<- TCGAanalyze_Normalization(tabDF = gl_bajo_alto,
                                              geneInfo = geneInfoHT, 
                                              method = "gcContent")  

boxplot(gl_bajo_altoAAIC, outline = FALSE)
boxplot(gl_bajo_alto_norm, outline = FALSE)

gl_bajo_alto_filt <- TCGAanalyze_Filtering(tabDF = gl_bajo_alto_norm,
                                           method = "quantile", 
                                           qnt.cut = 0.25)     

prost_GED <- TCGAanalyze_DEA(mat1 = gl_bajo_alto_filt[,gl6_group],
                             mat2 = gl_bajo_alto_filt[,gl_mas7_group],
                             Cond1type = "Gleason score < 7",
                             Cond2type = "Gleason score > 7",
                             fdr.cut = 0.01,
                             logFC.cut = 1,
                             method = "glmLRT")



##############################################################################
##############################################################################

# 3.- GLEASON = 7  vs GLEASON > 7 (gleason medio vs gleason alto)

query_gl7_mas7 <- GDCquery(project = "TCGA-PRAD",
                           data.category = "Transcriptome Profiling",
                           data.type = "Gene Expression Quantification",
                           workflow.type = "HTSeq - FPKM",
                           barcode = c(gl7_group, gl_mas7_group))
#GDCdownload(query = query_gl7_mas7, directory = "/home/cone/Escritorio/Proyect_Final/data/")

gl7_mas7 <- GDCprepare(query = query_gl7_mas7,
                       save = TRUE,
                       directory = Directorio_sub,
                       save.filename = "se_gl7_mas7.rda")

gl7_mas7AAICC <- TCGAanalyze_Preprocessing(object = gl7_mas7,
                                           cor.cut = 0.6, 
                                           filename = "AAIC7_mas7.png",
                                           width = 2500,
                                           height = 3000,
                                           datatype = names(assays(gl7_mas7))[1])

gl7_mas7_norm<- TCGAanalyze_Normalization(tabDF = gl7_mas7,
                                          geneInfo = geneInfoHT, 
                                          method = "gcContent")  

gl7_mas7_filt <- TCGAanalyze_Filtering(tabDF = gl_7_mas7_norm,
                                       method = "quantile", 
                                       qnt.cut = 0.25)      

gl7_mas7GED <- TCGAanalyze_DEA(mat1 = gl7_mas7_filt[,gl7_group],
                               mat2 = gl7_mas7_filt[,gl_mas7_group],
                               Cond1type = "Gleason score = 7",
                               Cond2type = "Gleason score > 7",
                               fdr.cut = 0.01,
                               logFC.cut = 1,
                               method = "glmLRT")


load("~/Escritorio/Proyect_Final/se_prostata.rda")
se_prost <- data
tumor <- se_prost[,se_prost$shortLetterCode == "TP"]


