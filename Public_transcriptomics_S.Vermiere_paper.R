############ Public transcriptomics S. Vermeire

################## Data --> https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75214

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

####### 1. PROCESSING DATA ####### 

      ######### Deleting data I don't need
      
      matrix<-read.delim("Proyecto_Monogenic_IBD/Public_transcriptomics/GSE75214_series_matrix.txt/GSE75214_series_matrix.txt")
      matrix_2<-matrix[-(1:33),]
      matrix_3<-matrix_2[,118:ncol(matrix_2)]
      matrix_3<-matrix_3[-1,]
      
      muestras<-data.frame(tipo=colnames(matrix_3), nombre=as.character(unlist(matrix_3[1,])))
      rownames(muestras)<-muestras$nombre
      colnames(matrix_3)<-unlist(matrix_3[1,])
      
      genes<-matrix_2[,1]
      genes<-genes[-1]
      rownames(matrix_3)<-genes

      ######### Preparing COLData with the type of data I have (CD + control)
      
      muestras$cond<-rep("CD_active",nrow(muestras))
      muestras[52:67,3] <- "CD_inactive"
      muestras[68:78,3]<- "control"

      ########## Gene names
      
      mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL") ### to define mart (for this case ensemble human genome build 37 genes)
      mart <- useDataset("hsapiens_gene_ensembl", mart)
      
      listAttributes(mart)  ### list of info you could ask for
      listFilters(mart) ### list of filters
      
      values <- rownames(matrix_3)
      affy.to.hgnc <- getBM(mart=mart, attributes = c("affy_hugene_1_0_st_v1","hgnc_symbol"), filter = "affy_hugene_1_0_st_v1", values = values, uniqueRows = TRUE) 
      affy.to.ensembl <- getBM(mart=mart, attributes = c("affy_hugene_1_0_st_v1","ensembl_gene_id"), filter = "affy_hugene_1_0_st_v1", values = values, uniqueRows = TRUE) 

      
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


####### 3. ANALYSIS PLOT #######

  ########## CD active vs control########## 
    
    active_control<-read.table("Proyecto_Monogenic_IBD/Resultados_LIMMA/Results_CDactive_control.txt",  header = T)
    
    pvalor<-0.05/nrow(active_control)
    active_control$label<-as.numeric(row.names(active_control))
    active_control<-active_control[order(active_control$label),]
    active_control<-active_control[-nrow(active_control),]

    ###### Add gene symbol (label)

    affy.to.hgnc <- getBM(mart=mart, attributes = c("affy_hugene_1_0_st_v1","hgnc_symbol"), filter = "affy_hugene_1_0_st_v1", values = values, uniqueRows = TRUE) ### correr biomart, permite sacar varios attributes pero solo con un fltro, en tu caso son probes de afimetrix, en values especificas el nombre de esas probes (22222_at)
    affy.to.hgnc <- affy.to.hgnc[order(affy.to.hgnc$affy_hugene_1_0_st_v1),]
    active_control$gen<-affy.to.hgnc[match(active_control$label, affy.to.hgnc$affy_hugene_1_0_st_v1),2]

    ###### Add CD genes (type)

    id_CD<-read.table("Proyecto_Monogenic_IBD/Resultados_DESeq_monogenic_IBD_genes/Listas_con_genes_importantes/Genes_CD_symbol.txt", header=T)
    active_control$Type<-rep("other",nrow(active_control))
    active_control$Type[match(id_CD$x,active_control$gen)]<-"CD"
        
    ###### Add CD genes (for color)
        
    cd<-str_subset(active_control$Type, "CD")
    cd_g<-active_control[active_control$Type %in% cd,]

    ###### Add transparency to points
        
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
        

        