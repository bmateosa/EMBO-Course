############ Public transcriptomics S. Vermeire
library(DESeq2)
library(biomaRt)
library(clusterProfiler)
library(limma)
library(ggplot2)
library(ggrepel)
library(stringr)


                  #####################################
                  ############## CD DATA ##############
                  #####################################

####### 1. PREPARO LOS DATOS ####### 

######### Quito los datos que no me interesan

matrix<-read.delim("Proyecto_Monogenic_IBD/Public_transcriptomics/GSE75214_series_matrix.txt/GSE75214_series_matrix.txt")
matrix_2<-matrix[-(1:33),]
colnames(matrix_2)
matrix_3<-matrix_2[,118:ncol(matrix_2)]
matrix_3<-matrix_3[-1,]

muestras<-data.frame(tipo=colnames(matrix_3), nombre=as.character(unlist(matrix_3[1,])))
rownames(muestras)<-muestras$nombre
colnames(matrix_3)<-unlist(matrix_3[1,])

genes<-matrix_2[,1]
genes<-genes[-1]
rownames(matrix_3)<-genes

######### Preparo COLData con los tipos de muestras que tengo (CD + control)

muestras$cond<-rep("CD_active",nrow(muestras))
muestras[52:67,3] <- "CD_inactive"
muestras[68:78,3]<- "control"

########## Nombres de los genes

mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL") ### para definer el mart (en este caso ensemble human genome build 37 genes)
mart <- useDataset("hsapiens_gene_ensembl", mart)

listAttributes(mart)  ### lista de informacion que puedes pedir
listFilters(mart) ### lista de filtros por los que filtrar ,suelen ser iguales que los atributos

values <- rownames(matrix_3)
affy.to.hgnc <- getBM(mart=mart, attributes = c("affy_hugene_1_0_st_v1","hgnc_symbol"), filter = "affy_hugene_1_0_st_v1", values = values, uniqueRows = TRUE) ### correr biomart, permite sacar varios attributes pero solo con un fltro, en tu caso son probes de afimetrix, en values especificas el nombre de esas probes (22222_at)
affy.to.ensembl <- getBM(mart=mart, attributes = c("affy_hugene_1_0_st_v1","ensembl_gene_id"), filter = "affy_hugene_1_0_st_v1", values = values, uniqueRows = TRUE) ### correr biomart, permite sacar varios attributes pero solo con un fltro, en tu caso son probes de afimetrix, en values especificas el nombre de esas probes (22222_at)

####### 2. LIMMA #######
design<-model.matrix(~ muestras$cond -1)
rownames(design)<-rownames(muestras)
colnames(design)<-c("CD_active", "CD_inactive", "control")
head(design)
colSums(design)

contrast.matrix<-makeContrasts(CD_active-control, CD_active-CD_inactive, levels=design)
matrix_3<-apply(matrix_3, 2, as.numeric)
rownames(matrix_3)<-genes
matrix_3[is.na(matrix_3)]<-0

a<-voom(matrix_3, design, plot=TRUE)
fit<-lmFit(matrix_3, design)
fit2<-contrasts.fit(fit, contrast.matrix)
fit2.2<-eBayes(fit2)
topTable(fit2.2)
topTable(fit2.2, coef="CD_active - control")
results_active.control<-topTable(fit2.2, coef="CD_active - control", number=Inf, sort.by = "logFC")
results_active.inactive<-topTable(fit2.2, coef="CD_active - CD_inactive", number=Inf, sort.by = "logFC")


write.table(results_active.control, "Proyecto_Monogenic_IBD/Resultados_LIMMA/Results_CDactive_control.txt", quote=F, sep = "\t", row.names = T, col.names = T)
write.table(results_active.inactive, "Proyecto_Monogenic_IBD/Resultados_LIMMA/Results_CDactive_CDinactive.txt", quote=F, sep = "\t", row.names = T, col.names = T)




####### 3. ANÁLISIS #######

########## CD activo vs control

active_control<-read.table("Proyecto_Monogenic_IBD/Resultados_LIMMA/Results_CDactive_control.txt",  header = T)

pvalor<-0.05/nrow(active_control)
active_control$label<-as.numeric(row.names(active_control))
active_control<-active_control[order(active_control$label),]
active_control<-active_control[-nrow(active_control),]

    ###### Añado gene symbol (label)

        affy.to.hgnc <- getBM(mart=mart, attributes = c("affy_hugene_1_0_st_v1","hgnc_symbol"), filter = "affy_hugene_1_0_st_v1", values = values, uniqueRows = TRUE) ### correr biomart, permite sacar varios attributes pero solo con un fltro, en tu caso son probes de afimetrix, en values especificas el nombre de esas probes (22222_at)
        affy.to.hgnc <- affy.to.hgnc[order(affy.to.hgnc$affy_hugene_1_0_st_v1),]
        active_control$gen<-affy.to.hgnc[match(active_control$label, affy.to.hgnc$affy_hugene_1_0_st_v1),2]

    ###### Añado CD genes (type)

        id_CD<-read.table("Proyecto_Monogenic_IBD/Resultados_DESeq_monogenic_IBD_genes/Listas_con_genes_importantes/Genes_CD_symbol.txt", header=T)
        active_control$Type<-rep("other",nrow(active_control))
        active_control$Type[match(id_CD$x,active_control$gen)]<-"CD"
        
    ###### Añado CD genes (para colorear)
        
        cd<-str_subset(active_control$Type, "CD")
        cd_g<-active_control[active_control$Type %in% cd,]

    ###### Añado transparencia a los puntos
        
        active_control$alpha<-rep(1,nrow(active_control))
        active_control$alpha[active_control$Type=="other"]<-0.2
        
        
    ###### volcano ggplot
       pdf("Proyecto_Monogenic_IBD/Plots/LIMMA_CD_genes/Volcano_CD_genes_LIMMA.pdf") 
        active_control_plot<-ggplot(active_control, aes(x = logFC, y = -log10(P.Value), label=gen, color=Type)) +
          geom_point(size = 2.5, alpha=active_control$alpha, )+
          geom_point(data=cd_g)+
          geom_text_repel(data=cd_g, aes(label=gen))+
          geom_hline(yintercept = -log10(pvalor), linetype = "dashed") +
          geom_vline(xintercept = c(-1,1), linetype = "dashed")+
          xlab("log2 fold change") +
          ylab("-log10 p value")+
          theme_classic()+
          scale_color_manual(values=c("red","grey"))+
          ggtitle("CD genes in public data \n")+
          theme(plot.title=element_text(size=20,face= "bold", hjust=0.5))
        plot(active_control_plot)
        dev.off()
        
########## CD activo vs CD inactivo
        
        active_inactive<-read.table("Proyecto_Monogenic_IBD/Resultados_LIMMA/Results_CDactive_CDinactive.txt", header= T)
        
        pvalor<-0.05/nrow(active_inactive)
        active_inactive$label<-as.numeric(row.names(active_inactive))
        active_inactive<-active_inactive[order(active_inactive$label),]
        active_inactive<-active_inactive[-nrow(active_inactive),]       
        
        ###### Añado gene symbol (label)
        
        affy.to.hgnc <- getBM(mart=mart, attributes = c("affy_hugene_1_0_st_v1","hgnc_symbol"), filter = "affy_hugene_1_0_st_v1", values = values, uniqueRows = TRUE) ### correr biomart, permite sacar varios attributes pero solo con un fltro, en tu caso son probes de afimetrix, en values especificas el nombre de esas probes (22222_at)
        affy.to.hgnc <- affy.to.hgnc[order(affy.to.hgnc$affy_hugene_1_0_st_v1),]
        active_inactive$gen<-affy.to.hgnc[match(active_inactive$label, affy.to.hgnc$affy_hugene_1_0_st_v1),2]

        ###### Añado CD genes (type)
        
        id_CD<-read.table("Proyecto_Monogenic_IBD/Resultados_DESeq_monogenic_IBD_genes/Listas_con_genes_importantes/Genes_CD_symbol.txt", header=T)
        active_inactive$Type<-rep("other",nrow(active_inactive))
        active_inactive$Type[match(id_CD$x,active_inactive$gen)]<-"CD"
        
        ###### Añado CD genes (para colorear)
        
        cd<-str_subset(active_inactive$Type, "CD")
        cd_g<-active_inactive[active_inactive$Type %in% cd,]
        
        ###### Añado transparencia a los puntos
        
        active_inactive$alpha<-rep(1,nrow(active_inactive))
        active_inactive$alpha[active_inactive$Type=="other"]<-0.2
        
        ###### volcano ggplot
        pdf("Proyecto_Monogenic_IBD/Plots/LIMMA_CD_genes/Volcano_CD_genes_LIMMA_inactive.pdf") 
        active_inactive_plot<-ggplot(active_inactive, aes(x = logFC, y = -log10(P.Value), label=gen, color=Type)) +
          geom_point(size = 2.5, alpha=active_inactive$alpha, )+
          geom_point(data=cd_g)+
          geom_text_repel(data=cd_g, aes(label=gen))+
          geom_hline(yintercept = -log10(pvalor), linetype = "dashed") +
          geom_vline(xintercept = c(-1,1), linetype = "dashed")+
          xlab("log2 fold change") +
          ylab("-log10 p value")+
          theme_classic()+
          scale_color_manual(values=c("red","grey"))+
          ggtitle("CD genes in public data \n")+
          theme(plot.title=element_text(size=20,face= "bold", hjust=0.5))
        plot(active_inactive_plot)
        dev.off()
        
        
        
                  #####################################
                  ############## UC DATA ##############
                  #####################################
        
####### 1. PREPARO LOS DATOS ####### 
        
        ######### Quito los datos que no me interesan
        
        matrix<-read.delim("Proyecto_Monogenic_IBD/Public_transcriptomics/GSE75214_series_matrix.txt/GSE75214_series_matrix.txt")
        matrix_2<-matrix[-(1:33),]
        colnames(matrix_2)
        matrix_3<-matrix_2[,10:117]
        matrix_3<-matrix_3[-1,]
        
        muestras<-data.frame(tipo=colnames(matrix_3), nombre=as.character(unlist(matrix_3[1,])))
        rownames(muestras)<-muestras$nombre
        colnames(matrix_3)<-unlist(matrix_3[1,])
        
        genes<-matrix_2[,1]
        genes<-genes[-1]
        rownames(matrix_3)<-genes
        
        ######### Preparo COLData con los tipos de muestras que tengo (CD + control)
        
        muestras$cond<-rep("UC_active",nrow(muestras))
        muestras[86:nrow(muestras),3] <- "UC_inactive"
        muestras[1:11,3]<- "control"
        
        ########## Nombres de los genes
        
        mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL") ### para definer el mart (en este caso ensemble human genome build 37 genes)
        mart <- useDataset("hsapiens_gene_ensembl", mart)
        
        listAttributes(mart)  ### lista de informacion que puedes pedir
        listFilters(mart) ### lista de filtros por los que filtrar ,suelen ser iguales que los atributos
        
        values <- rownames(matrix_3)
        affy.to.hgnc <- getBM(mart=mart, attributes = c("affy_hugene_1_0_st_v1","hgnc_symbol"), filter = "affy_hugene_1_0_st_v1", values = values, uniqueRows = TRUE) ### correr biomart, permite sacar varios attributes pero solo con un fltro, en tu caso son probes de afimetrix, en values especificas el nombre de esas probes (22222_at)

####### 2. LIMMA #######
        
        design<-model.matrix(~ muestras$cond -1)
        rownames(design)<-rownames(muestras)
        colnames(design)<-c("control","UC_active", "UC_inactive")
        head(design)
        colSums(design)
        
        contrast.matrix<-makeContrasts(UC_active-control, UC_active-UC_inactive, levels=design)
        matrix_3<-apply(matrix_3, 2, as.numeric)
        rownames(matrix_3)<-genes
        matrix_3[is.na(matrix_3)]<-0
        
        a<-voom(matrix_3, design, plot=TRUE)
        fit<-lmFit(matrix_3, design)
        fit2<-contrasts.fit(fit, contrast.matrix)
        fit2.2<-eBayes(fit2)
        topTable(fit2.2)
        topTable(fit2.2, coef="UC_active - control")
        
        results_active.control<-topTable(fit2.2, coef="UC_active - control", number=Inf, sort.by = "logFC")
        results_active.inactive<-topTable(fit2.2, coef="UC_active - UC_inactive", number=Inf, sort.by = "logFC")
        
        
        write.table(results_active.control, "Proyecto_Monogenic_IBD/Resultados_LIMMA/Results_UCactive_control.txt", quote=F, sep = "\t", row.names = T, col.names = T)
        write.table(results_active.inactive, "Proyecto_Monogenic_IBD/Resultados_LIMMA/Results_UCactive_UCinactive.txt", quote=F, sep = "\t", row.names = T, col.names = T)
        
        affy.to.ensembl <- getBM(mart=mart, attributes = c("affy_hugene_1_0_st_v1","ensembl_gene_id"), filter = "affy_hugene_1_0_st_v1", values = values, uniqueRows = TRUE) ### correr biomart, permite sacar varios attributes pero solo con un fltro, en tu caso son probes de afimetrix, en values especificas el nombre de esas probes (22222_at)
        
####### 3. ANÁLISIS #######
        
  ########## UC activo vs control   ########## 
        
        active_control<-read.table("Proyecto_Monogenic_IBD/Resultados_LIMMA/Results_UCactive_control.txt",  header = T)
        
        pvalor<-0.05/nrow(active_control)
        active_control$label<-as.numeric(row.names(active_control))
        active_control<-active_control[order(active_control$label),]
        active_control<-active_control[-nrow(active_control),]
        
        ###### Añado gene symbol (label)
        
        affy.to.hgnc <- getBM(mart=mart, attributes = c("affy_hugene_1_0_st_v1","hgnc_symbol"), filter = "affy_hugene_1_0_st_v1", values = values, uniqueRows = TRUE) ### correr biomart, permite sacar varios attributes pero solo con un fltro, en tu caso son probes de afimetrix, en values especificas el nombre de esas probes (22222_at)
        affy.to.hgnc <- affy.to.hgnc[order(affy.to.hgnc$affy_hugene_1_0_st_v1),]
        active_control$gen<-affy.to.hgnc[match(active_control$label, affy.to.hgnc$affy_hugene_1_0_st_v1),2]
        
        ###### Añado UC genes (type)
        
        id_UC<-read.table("Proyecto_Monogenic_IBD/Resultados_DESeq_monogenic_IBD_genes/Listas_con_genes_importantes/Genes_UC_symbol.txt", header=T)
        active_control$Type<-rep("other",nrow(active_control))
        active_control$Type[match(id_UC$x,active_control$gen)]<-"UC"
        
        ##### Busco cuántos DEGs de UC genes tengo
        UC_active_control<-active_control[active_control$Type=="UC",]
        DEGs_UC_active_control<-UC_active_control[UC_active_control$logFC>1 & UC_active_control$P.Value<pvalor | UC_active_control$logFC<(-1) & UC_active_control$P.Value<pvalor,]
        
        ###### Añado UC genes (para colorear)
        
        uc<-str_subset(active_control$Type, "UC")
        uc_g<-active_control[active_control$Type %in% uc,]
        
        ###### Añado transparencia a los puntos
        
        active_control$alpha<-rep(1,nrow(active_control))
        active_control$alpha[active_control$Type=="other"]<-0.2
        
        
        ###### volcano ggplot
        pdf("Proyecto_Monogenic_IBD/Plots/LIMMA_UC_genes/Volcano_UC_genes_LIMMA_Severine.pdf") 
        active_control_plot<-ggplot(active_control, aes(x = logFC, y = -log10(P.Value), label=gen, color=Type)) +
          geom_point(size = 2.5, alpha=active_control$alpha, )+
          geom_point(data=uc_g)+
          geom_text_repel(data=uc_g, aes(label=gen))+
          geom_hline(yintercept = -log10(pvalor), linetype = "dashed") +
          geom_vline(xintercept = c(-1,1), linetype = "dashed")+
          xlab("log2 fold change") +
          ylab("-log10 p value")+
          theme_classic()+
          scale_color_manual(values=c("grey", "red"))+
          ggtitle("UC genes in public data \n")+
          theme(plot.title=element_text(size=20,face= "bold", hjust=0.5))
        plot(active_control_plot)
        dev.off()     
        
  ########## UC activo vs UC inactivo   ########## 
        
        active_inactive<-read.table("Proyecto_Monogenic_IBD/Resultados_LIMMA/Results_UCactive_UCinactive.txt", header= T)
        
        pvalor<-0.05/nrow(active_inactive)
        active_inactive$label<-as.numeric(row.names(active_inactive))
        active_inactive<-active_inactive[order(active_inactive$label),]
        active_inactive<-active_inactive[-nrow(active_inactive),]       
        
        ###### Añado gene symbol (label)
        values<-rownames(active_inactive)
        affy.to.hgnc <- getBM(mart=mart, attributes = c("affy_hugene_1_0_st_v1","hgnc_symbol"), filter = "affy_hugene_1_0_st_v1", values = values, uniqueRows = TRUE) ### correr biomart, permite sacar varios attributes pero solo con un fltro, en tu caso son probes de afimetrix, en values especificas el nombre de esas probes (22222_at)
        affy.to.hgnc <- affy.to.hgnc[order(affy.to.hgnc$affy_hugene_1_0_st_v1),]
        active_inactive$gen<-affy.to.hgnc[match(active_inactive$label, affy.to.hgnc$affy_hugene_1_0_st_v1),2]
        
        ###### Añado UC genes (type)
        
        id_UC<-read.table("Proyecto_Monogenic_IBD/Resultados_DESeq_monogenic_IBD_genes/Listas_con_genes_importantes/Genes_UC_symbol.txt", header=T)
        active_inactive$Type<-rep("other",nrow(active_inactive))
        active_inactive$Type[match(id_UC$x,active_inactive$gen)]<-"UC"
        
        ##### Busco cuántos DEGs de UC genes tengo
        UC_active_inactive<-active_inactive[active_inactive$Type=="UC",]
        DEGs_UC_active_inactive<-UC_active_inactive[UC_active_inactive$logFC>1 & UC_active_inactive$P.Value<pvalor | UC_active_inactive$logFC<(-1) & UC_active_inactive$P.Value<pvalor,]
        
        ###### Añado UC genes (para colorear)
        
        uc<-str_subset(active_inactive$Type, "UC")
        uc_g<-active_inactive[active_inactive$Type %in% uc,]
        
        ###### Añado transparencia a los puntos
        
        active_inactive$alpha<-rep(1,nrow(active_inactive))
        active_inactive$alpha[active_inactive$Type=="other"]<-0.2
        
        ###### volcano ggplot
        pdf("Proyecto_Monogenic_IBD/Plots/LIMMA_UC_genes/Volcano_UC_genes_LIMMA_inactive.pdf") 
        active_inactive_plot<-ggplot(active_inactive, aes(x = logFC, y = -log10(P.Value), label=gen, color=Type)) +
          geom_point(size = 2.5, alpha=active_inactive$alpha, )+
          geom_point(data=uc_g)+
          geom_text_repel(data=uc_g, aes(label=gen))+
          geom_hline(yintercept = -log10(pvalor), linetype = "dashed") +
          geom_vline(xintercept = c(-1,1), linetype = "dashed")+
          xlab("log2 fold change") +
          ylab("-log10 p value")+
          theme_classic()+
          scale_color_manual(values=c("grey","red"))+
          ggtitle("UC genes in public data - active vs inactive \n")+
          theme(plot.title=element_text(size=20,face= "bold", hjust=0.5))
        plot(active_inactive_plot)
        dev.off()
        
        