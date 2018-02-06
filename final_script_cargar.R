# ANALISIS DIFERENCIAL DE LA EXPRESIÓN EN CÁNCER DE PROSTATA

# Cargar paquetes necesarios. 
library("TCGAbiolinks")
library("SummarizedExperiment")

load("~/Escritorio/Proyect_Final/se_prostata.rda")
data_se <- data
load("~/Escritorio/Proyect_Final/model.rda")
model_se <- data

# Resumen sobre los datos de mi proyecto disponibles en TCGA
# resumen <- getDataCategorySummary("TCGA-PRAD", legacy = FALSE)

# Buscar en TCGA datos de RNA-seq
#prost_query <- GDCquery(project = "TCGA-PRAD", 
#                        data.category = "Transcriptome Profiling",
#                        data.type = "Gene Expression Quantification",
#                        workflow.type = "HTSeq - FPKM")
#GDCdownload(query = prost_query, directory = "/home/cone/Escritorio/Proyect_Final/data/")
#View(getResults(prost_query))

# Obtener bcr (identificador) de las muestras para generar grupos (1.-Muestras tumor primario; 2.- Muestas de tejido soĺido normal)
#prost_muestras <- getResults(prost_query, cols = c("cases"))
#prost_TP <- TCGAquery_SampleTypes(barcode = prost_muestras, typesample = "TP")  # TP -> Tumor Primary --------> 498 samples
#prost_NT <- TCGAquery_SampleTypes(barcode = prost_muestras, typesample = "NT")  # NT -> Solid Tissue Normal --> 52 samples


# Buscar los datos específicos de mis grupos y descargarlos.
#datos_query <- GDCquery(project = "TCGA-PRAD",
#                        data.category = "Transcriptome Profiling",
#                        data.type = "Gene Expression Quantification",
#                        workflow.type = "HTSeq - FPKM",
#                        barcode = c(prost_TP, prost_NT))

# Crear objeto SummarizedExperiment y guardar una copia 
#data_se <- GDCprepare(query = datos_query,
#                      save = TRUE,
#                      directory = "/home/cone/Escritorio/Proyect_Final/data/",
#                      save.filename = "se_prostata.rda")

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
prost_pre <- TCGAanalyze_Preprocessing(object = data_se[,c(1:5,515:520)],
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
                                       method = "gcContent");  dim(prost_norm)

boxplot(prost_pre, outline = FALSE)
boxplot(prost_norm[,1:20], outline = FALSE)


# FILTRADO --> con el paquete "geneFilter". Devuelve todos los mRNA con una media (a través de todas las muestras), superior
# al umbral definido por la media en el cuartil a traves de todas las muestras. 
prost_filt <- TCGAanalyze_Filtering(tabDF = prost_norm,
                                    method = "quantile", 
                                    qnt.cut = 0.66); dim(prost_filt)  # qnt.cut = umbral para el filtrado

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
                        nBar = 10)

#################################################################################################
#################################################################################################
# Necesitamos los barcodes, según el grupo, para generar nuestro objeto "SummarizedExperiment" 
# y poder hacer un DEA (Análisis Diferencial de la Expresión)

# Barcodes para las muestras que presenten un gleason < 7 (Gleason Bajo)
gl6 <- subset(data_se, data_se$barcode, data_se$subtype_Reviewed_Gleason_sum < 7)
gl6_ident <- TCGAquery_SampleTypes(barcode = gl6$barcode, typesample = "TP")
#gl6$subtype_Reviewed_Gleason_sum

# Barcodes para muestras que presenten un gleason = 7 (Gleason Medio)
gl7 <- subset(data_se, data_se$barcode, data_se$subtype_Reviewed_Gleason_sum == 7)
gl7_ident <- TCGAquery_SampleTypes(barcode = gl7$barcode, typesample = "TP")
#gl7$subtype_Reviewed_Gleason_sum

# Barcodes para muestras que presentan gleason > 7 (Gleason Alto)
gl_mas7 <- subset(data_se, data_se$barcode, data_se$subtype_Reviewed_Gleason_sum > 7)
gl_mas7_ident <- TCGAquery_SampleTypes(barcode = gl_mas7$barcode, typesample = "TP")
#gl_mas7$subtype_Reviewed_Gleason_sum

# Barcodes para muestras que presental gleason <= 7 (Gleason Bajo-Medio)
gl_bm <- subset(data_se, data_se$barcode, data_se$subtype_Reviewed_Gleason_sum <= 7)
gl_bm_ident <- TCGAquery_SampleTypes(barcode = gl_bm$barcode, typesample = "TP")
#gl_bm$subtype_Reviewed_Gleason_sum

##################################################################################################
##################################################################################################

# MODELOS PREDICTIVOS
# 1.- Realizar DEA entre grupos gleason medio-bajo VS gleason alto

#model_query <- GDCquery(project = "TCGA-PRAD",
#                        data.category = "Transcriptome Profiling",
#                        data.type = "Gene Expression Quantification",
#                        workflow.type = "HTSeq - FPKM",
#                        barcode = c(gl_mas7_ident, gl_bm_ident))

#model_se <- GDCprepare(query = model_query,
#                       save = TRUE,
#                       save.filename = "model.rda",
#                       directory = "/home/cone/Escritorio/Proyect_Final/data/")

model_norm<- TCGAanalyze_Normalization(tabDF = model_se,
                                       geneInfo = geneInfoHT, 
                                       method = "gcContent")  ; dim(model_norm)

model_filt <- TCGAanalyze_Filtering(tabDF = model_norm,
                                    method = "quantile", 
                                    qnt.cut = 0.75); dim(model_filt)

model_GED <- TCGAanalyze_DEA(mat1 = model_filt[,gl_mas7_ident],
                             mat2 = model_filt[,gl_bm_ident],
                             Cond1type = "Gleason score <= 7",
                             Cond2type = "Gleason score > 7",
                             fdr.cut = 0.01,
                             logFC.cut = 1,
                             method = "glmLRT") ;dim(model_GED)

model_GED_inf <- TCGAanalyze_LevelTab(model_GED,"gl <=7","gl>7", model_filt[,gl_mas7_ident],model_filt[,gl_bm_ident])

# Generar tabla con variables predictoras (GEDs) y variable respuesta (gleason medio-bajo ó gleason alto)
genes_DEA <- model_GED_inf$ensembl_gene_id
name_genes <- model_GED_inf$external_gene_name
tabla_glmnet <- model_filt[genes_DEA,]
rownames(tabla_glmnet) <- name_genes
tabla_glmnet <- t(tabla_glmnet)
gleason <- model_se$subtype_Reviewed_Gleason_sum
tabla_glmnet2 <- cbind(gleason,tabla_glmnet)
tabla_final <- as.data.frame(tabla_glmnet2)
tabla_final$G <- NA
tabla_final$G [tabla_final$gleason <= 7] <- "NO"
tabla_final$G [tabla_final$gleason > 7] <- "SI"
as.factor(tabla_final$G)
#View(tabla_final)
glmnet_datos <- tabla_final[,-1]
dim(glmnet_datos)
#View(glmnet_datos)

# MODELO DE REGRESION PENALIZADA 
xtabs(~ G + gleason, tabla_final) # Crea tabla (opcionalx matriz de dispersión), a partir de factores
                                  # contenidos en el data.frame

set.seed(7)                       
ss <- sample(1:2, size=nrow(glmnet_datos), replace = T, prob = c(0.66, 0.33))
train <-glmnet_datos[ss==1,]
test <- glmnet_datos[ss==2,]

# Model matrix de los predictores en el set de training
X<-model.matrix(~., train)
Y<-X[,35]
X<-X[,-c(35,1)]

# Model matrix de los predictores en el set de validacion
Xpred <- model.matrix(~., test)
Ypred<-Xpred[,35]
Xpred<-Xpred[,-c(1, 35)]


library(glmnet)
# Validación Cruzada
# Permite evaluar los modelos de regresión penalizados con validación cruzada para optimizar el parámetro lambda.
# alpha --> alpha = 1 (Lasso), alpha = 0 (Ridge), alpha = intermedio (Elastic Net).

# ALPHA = 0.1
cv01<-cv.glmnet(X, Y, family="binomial", alpha=0.1) 
cv01
plot(cv01)
lambda <- cv01$lambda.1se  # factor de penalización
glmnet01 <- glmnet(X, Y, family="binomial", alpha=0.1)
predict(glmnet01, s=lambda, type="coef")

# cv1$lambda # valores de lambda utilizados para los ajustes
# cv1$cvm # Media para el error de la validación cruzada
# cv1$cvsd # Estimación del error estándar de cvm
# cv01$lambda.min # Valor de lambda para obtener el mín cvm
# cv01$lambda.1se # Mayor valor de lambda cuyo error está en 1 y el error estándar del mín


predict(glmnet01, s=lambda, type="coef")
predict(glmnet01, s=median(lambdas), type="coef")


#Área bajo la curva aparente (sobre los datos de training)
library(pROC)
pred_test01 <- predict(glmnet01, s=lambda, type="response", newx=X)
roc_tes01 <- roc(Y ~ pred_test01[,1])
roc_tes01
ci(roc_tes01)
plot(roc_tes01)
#Área bajo la curva validada (sobre los datos de validación)
pred_val01<-predict(glmnet01, s=lambda, type="response", newx=Xpred)
roc_val01 <- roc(Ypred ~ pred_val01[,1])
roc_val01
ci(roc_val01)
plot(roc_val01)




# MODELOS PREDICTIVOS MACHINE LEARNING
glmnet_datos
dim(glmnet_datos)

nw_tabla <- glmnet_datos
library(caret)
set.seed(123)
semillas <- vector(mode = "list", length = 11)
for(i in 1:11) semillas[[i]] <- sample.int(1000,10)
semillas[[11]] <- sample.int(1000, 1)


fitcontrol <- trainControl(method = "cv",number=10,
                           returnResamp = "final",
                           summaryFunction = twoClassSummary,
                           classProbs=TRUE,
                           seeds = semillas,
                           verboseIter=TRUE)

set.seed(342)
glmnettFIT <- train(G~.,data=nw_tabla,
                    method="glmnet",
                    tuneLength=10,
                    trControl=fitcontrol,
                    metric = "ROC")

summary(glmnettFIT)
glmnet.prob <- predict(glmnettFIT, newdata=nw_tabla,type = "prob")
glmnet.roc <- roc(nw_tabla$G, glmnet.prob$SI, dataGrid = TRUE, gridLength=100)
glmnet.roc
plot(glmnet.roc)


set.seed(342)
C50Fit <- train(G~.,data=nw_tabla,
                method="C5.0",
                preProc= "range",
                trControl=fitcontrol,
                metric = "ROC")
summary(C50Fit)
C50.prob <- predict(C50Fit, newdata = nw_tabla, type = "prob")
C50.prob
C50.roc <- roc(nw_tabla$G, C50.prob$SI, dataGrid = TRUE, gridLength=100)
C50.roc
plot(C50.roc)

set.seed(342)
rpartFit <- train(G~.,data=nw_tabla,
                  method="rpart",
                  tuneLength=10,
                  trControl=fitcontrol,
                  metric = "ROC")
summary(rpartFit)
rpartfit.prob <- predict(rpartFit, newdata = nw_tabla, type = "prob")
rpartfit.roc <- roc(nw_tabla$G, rpartFit.prob$SI, dataGrid = TRUE, gridLength=100)
rpartFit.roc
plot(rpartFit.roc)
