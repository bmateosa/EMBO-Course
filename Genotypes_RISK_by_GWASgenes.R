################### HIGH-LOW - GENOTYPES FROM RISK study ##################

####### We selected top genes from IBD GWAS, and we want to see if there are transcriptomics differences between
####### homocygous patients with high risk and homocygous patients with no risk.
####### We use two methods to select SNPs: Top CS2G and the closest SNP to the gene.
###### Genotypes and expression was analyzed from the RISK cohort

library(stringr)
library(DESeq2)
library(tidyverse)
library(data.table)


###### 1. RUN SNPS BY TWO DIFFERENT METHODS ###### 

      closest_SNPs<-read.table("Proyecto_Monogenic_IBD/PRS_RISK/SNPs_to_GWASgenes/Closest_snp_toRISK_genes.txt", header=T)
      closest_SNPs<-closest_SNPs[-c(7),]
      
      TopCS2G_SNPs<-read.table("Proyecto_Monogenic_IBD/PRS_RISK/SNPs_to_GWASgenes/TopCS2G_snp_toRISK_genes.txt", header=T)
      
      genotipos<-read.table("RISK/genotipos/02f_risk_ichip_Ileum_365ind_101K_MAF5_GENO1perc.bed")
      SNPs_risk<- read.table("RISK/genotipos/02f_risk_ichip_Ileum_365ind_101K_MAF5_GENO1perc.bim")
      individuals_risk<-read.table("RISK/genotipos/02f_risk_ichip_Ileum_365ind_101K_MAF5_GENO1perc.fam")
      
      closest_SNPs_risk<-SNPs_risk[which(SNPs_risk$V2 %in% closest_SNPs$SNP),]
      TopCS2G_SNPs_risk<-SNPs_risk[which(SNPs_risk$V2 %in% TopCS2G_SNPs$snp1),]
      
      closest_SNPs_risk$GENE<-closest_SNPs[match(closest_SNPs_risk$V2,closest_SNPs$SNP),6]
      TopCS2G_SNPs_risk$GENE<- TopCS2G_SNPs[match(TopCS2G_SNPs_risk$V2,TopCS2G_SNPs$snp1),2]
      
      lista_closest_SNPs<-as.vector(closest_SNPs_risk$V2)
      write.table(lista_closest_SNPs, "Proyecto_Monogenic_IBD/PRS_RISK/SNPs_to_GWASgenes/Closest_snp_toRISK_genes_list.txt", col.names = F, row.names = F, quote = F, sep = "\t")
      
      lista_TopCS2G_SNPs<-as.vector(TopCS2G_SNPs_risk$V2)
      write.table(lista_TopCS2G_SNPs, "Proyecto_Monogenic_IBD/PRS_RISK/SNPs_to_GWASgenes/TopCS2G_snp_toRISK_genes_list.txt", col.names = F, row.names = F, quote = F, sep = "\t")

###### 2.RUN PLINK ###### 

      pathtoplink<-"PLINK/"
      system(paste0(pathtoplink, "plink --bfile RISK/genotipos/02f_risk_ichip_Ileum_365ind_101K_MAF5_GENO1perc",
      " --extract Proyecto_Monogenic_IBD/PRS_RISK/SNPs_to_GWASgenes/Closest_snp_toRISK_genes_list.txt",
      " --recode",
      " --out Proyecto_Monogenic_IBD/PRS_RISK/SNPs_to_GWASgenes/Risk_genotypes_interestingSNP_GWASgenes"))
      
      
      pathtoplink<-"PLINK/"
      system(paste0(pathtoplink, "plink --bfile RISK/genotipos/02f_risk_ichip_Ileum_365ind_101K_MAF5_GENO1perc",
                    " --extract Proyecto_Monogenic_IBD/PRS_RISK/SNPs_to_GWASgenes/TopCS2G_snp_toRISK_genes_list.txt",
                    " --recode",
                    " --out Proyecto_Monogenic_IBD/PRS_RISK/SNPs_to_GWASgenes/Risk_genotypes_interestingSNP_GWASgenes_TopCS2G"))


###### 3. READING GENERATED DATA AND SELECTION OF CD INDIVIDUALS (EXCLUDING CONTROLS) ###### 

      risk_closest<-read.table("Proyecto_Monogenic_IBD/PRS_RISK/SNPs_to_GWASgenes/Risk_genotypes_interestingSNP_GWASgenes_closest.ped")
      risk_TopCS2G<-read.table("Proyecto_Monogenic_IBD/PRS_RISK/SNPs_to_GWASgenes/Risk_genotypes_interestingSNP_GWASgenes_TopCS2G.ped")
      
      fenotipos<-read.table("Proyecto_Monogenic_IBD/01_NewAnalysis_MiOkProb_Phen.txt", header=T)
      
      fenotipos_id<-fenotipos$label
      fenotipos_id<-sub("-","",fenotipos_id, fixed = T)
      fenotipos_id<-paste0("RS",fenotipos_id)
      risk_closest$V2<-str_sub(risk_closest$V2,1,7)
      risk_TopCS2G$V2<-str_sub(risk_TopCS2G$V2,1,7)
      
      mis_genotipos_closest<-risk_closest[which(risk_closest$V2 %in% fenotipos_id),]
      mis_genotipos_TopCS2G<-risk_TopCS2G[which(risk_TopCS2G$V2 %in% fenotipos_id),]
      
      controles<-fenotipos[fenotipos$TO_TAKE_MIOK.1=="CTRL",1]
      controles_id<-sub("-","", controles, fixed=T)
      controles_id<-paste0("RS", controles_id)
      
      mis_genotipos_closest<-mis_genotipos_closest[which(!mis_genotipos_closest$V2 %in% controles_id),]
      mis_genotipos_TopCS2G<-mis_genotipos_TopCS2G[which(!mis_genotipos_TopCS2G$V2 %in% controles_id),]
      
      
      write.table(mis_genotipos_closest,"Proyecto_Monogenic_IBD/PRS_RISK/SNPs_to_GWASgenes/Genotipos_CDpatients_RISK_misSNPs_closest.txt", row.names = F, quote = F, sep="\t")
      write.table(mis_genotipos_TopCS2G,"Proyecto_Monogenic_IBD/PRS_RISK/SNPs_to_GWASgenes/Genotipos_CDpatients_RISK_misSNPs_TopCS2G.txt", row.names = F, quote = F, sep="\t")

###### 4. ADDING SNPS NAMES ######

      lista_closest_SNPs<-read.table("Proyecto_Monogenic_IBD/PRS_RISK/SNPs_to_GWASgenes/Closest_snp_toRISK_genes_list.txt")
      lista_TopCS2G_SNPs<-read.table("Proyecto_Monogenic_IBD/PRS_RISK/SNPs_to_GWASgenes/TopCS2G_snp_toRISK_genes_list.txt")
      
      lista_closest_SNPs<-as.vector(rep(lista_closest_SNPs$V1,each=2))
      lista_TopCS2G_SNPs<-as.vector(rep(lista_TopCS2G_SNPs$V1,each=2))
      
      mis_genotipos_closest<-read.table("Proyecto_Monogenic_IBD/PRS_RISK/SNPs_to_GWASgenes/Genotipos_CDpatients_RISK_misSNPs_closest.txt", header=T)
      mis_genotipos_TopCS2G<-read.table("Proyecto_Monogenic_IBD/PRS_RISK/SNPs_to_GWASgenes/Genotipos_CDpatients_RISK_misSNPs_TopCS2G.txt", header=T)


###### 5. LOOKING FOR THE SNPS IN GWAS TO GET THE RISK ALLELE ######

      GWAS<-fread("IBDGWAS/EUR.CD.gwas_info03_filtered.assoc", header=T,stringsAsFactors = F)
      GWAS_closest<-GWAS[which(GWAS$SNP %in% lista_closest_SNPs),]
      GWAS_TopCS2G<-GWAS[which(GWAS$SNP %in% lista_TopCS2G_SNPs),]
      
      GWAS_closest$beta<-log(GWAS_closest$OR)
      GWAS_TopCS2G$beta<-log(GWAS_TopCS2G$OR)
      
      colnames(GWAS_closest)[4]<-"alelo_efecto"
      colnames(GWAS_TopCS2G)[4]<-"alelo_efecto"
      
      GWAS_closest$riesgo[GWAS_closest$beta<0]<-"protector"
      GWAS_closest$riesgo[GWAS_closest$beta>0]<-"riesgo"
      
      GWAS_TopCS2G$riesgo[GWAS_TopCS2G$beta<0]<-"protector"
      GWAS_TopCS2G$riesgo[GWAS_TopCS2G$beta>0]<-"riesgo"

###### 6. DEFINE RISK INDIVIDUALS ######

################ 1. Changing colnames from mis_genotipos to later select columns by SNP

      ############## CLOSEST
        colnames_closest_SNPs<-vector("character", length(lista_closest_SNPs))
        for(i in seq_along(lista_closest_SNPs)){
          if(i %% 2==0){
            colnames_closest_SNPs[i]<-paste0(lista_closest_SNPs[i],"_A2")
          } else{
            colnames_closest_SNPs[i]<- paste0(lista_closest_SNPs[i],"_A1")
          }
        }
        
        colnames(mis_genotipos_closest)[7:ncol(mis_genotipos_closest)]<-colnames_closest_SNPs

      ############## TopCS2G
        lista_TopCS2G_SNPs<-read.table("Proyecto_Monogenic_IBD/PRS_RISK/SNPs_to_GWASgenes/TopCS2G_snp_toRISK_genes_list.txt")
        lista_TopCS2G_SNPs<-as.vector(rep(lista_TopCS2G_SNPs$V1,each=2))
        mis_genotipos_TopCS2G<-read.table("Proyecto_Monogenic_IBD/PRS_RISK/SNPs_to_GWASgenes/Genotipos_CDpatients_RISK_misSNPs_TopCS2G.txt", header=T)
        
        colnames_TopCS2G_SNPs<-vector("character", length(lista_TopCS2G_SNPs))
        for(i in seq_along(lista_TopCS2G_SNPs)){
          if(i %% 2==0){
            colnames_TopCS2G_SNPs[i]<-paste0(lista_TopCS2G_SNPs[i],"_A2")
          } else{
            colnames_TopCS2G_SNPs[i]<- paste0(lista_TopCS2G_SNPs[i],"_A1")
          }
        }
        
        colnames(mis_genotipos_TopCS2G)[7:ncol(mis_genotipos_TopCS2G)]<-colnames_TopCS2G_SNPs

################ 2. Loop to define who is homozygous with risk, heterocygous or homozygous with no risk.
        
    ############## Closest_SNPs
        
      for(i in 1:length(lista_closest_SNPs)){
        mi_snp<-mis_genotipos_closest[,grep(lista_closest_SNPs[i], colnames(mis_genotipos_closest))]
        row.names(mi_snp)<-mis_genotipos_closest$V2
        GWAS_mi_snp<-GWAS_closest[GWAS_closest$SNP==lista_closest_SNPs[i],]
        if (GWAS_mi_snp$beta>0)
          {alelo.riesgo<-GWAS_mi_snp$alelo_efecto
          cuantos<-rowSums(mi_snp==alelo.riesgo)
          mi_snp$cuantos<-cuantos
        }else{
          alelo.riesgo<-GWAS_mi_snp$A2
          cuantos<-rowSums(mi_snp==alelo.riesgo)
          mi_snp$cuantos<-cuantos
        }
        write.table(mi_snp, paste0("Proyecto_Monogenic_IBD/PRS_RISK/SNPs_to_GWASgenes/Individuos_RISK_homocigotos_paraSNPs/",lista_closest_SNPs[i],"_closest_SNPs.txt"), col.names = T, row.names = T, quote = F, sep = "\t")
      }

    ############## TopCS2G

      for(i in 1:length(lista_TopCS2G_SNPs)){
        mi_snp<-mis_genotipos_TopCS2G[,grep(lista_TopCS2G_SNPs[i], colnames(mis_genotipos_TopCS2G))]
        row.names(mi_snp)<-mis_genotipos_TopCS2G$V2
        GWAS_mi_snp<-GWAS_TopCS2G[GWAS_TopCS2G$SNP==lista_TopCS2G_SNPs[i],]
        if (GWAS_mi_snp$beta>0)
        {alelo.riesgo<-GWAS_mi_snp$alelo_efecto
        cuantos<-rowSums(mi_snp==alelo.riesgo)
        mi_snp$cuantos<-cuantos
        }else{
          alelo.riesgo<-GWAS_mi_snp$A2
          cuantos<-rowSums(mi_snp==alelo.riesgo)
          mi_snp$cuantos<-cuantos
        }
        write.table(mi_snp, paste0("Proyecto_Monogenic_IBD/PRS_RISK/SNPs_to_GWASgenes/Individuos_RISK_homocigotos_paraSNPs/",lista_TopCS2G_SNPs[i],"_TopCS2G_SNPs.txt"), col.names = T, row.names = T, quote = F, sep = "\t")
      }

###### 7. LET's DO DESeq2 BETWEEN INDIVIDUALS WITH RISK (HOMOCYGOUS) AND INDIVIDUALS WITH NO RISK (HOMOCYGOUS) ###### 

  ##############CLOSEST
        
      counts<-read.table("Proyecto_Monogenic_IBD/01_NewAnalysis_MiOkProb_Norm.txt",header=T)
      closest_SNPs<-read.table("Proyecto_Monogenic_IBD/PRS_RISK/SNPs_to_GWASgenes/Closest_snp_toRISK_genes_list.txt")

      ind_closest<-read.table("Proyecto_Monogenic_IBD/PRS_RISK/SNPs_to_GWASgenes/Genotipos_CDpatients_RISK_misSNPs_closest.txt", header=T)
      individuals<-ind_closest$V2
      
      colnames(counts)[3:ncol(counts)]<-sub("norm.","RS",colnames(counts)[3:ncol(counts)], fixed = T)
      colnames(counts)[3:ncol(counts)]<-sub(".","",colnames(counts)[3:ncol(counts)], fixed = T)
      
      counts_data<-counts[,which(colnames(counts) %in% individuals)]
      counts_data<-counts_data[,order(colnames(counts_data), decreasing = F)]
      
      all(colnames(counts_data_misind)==rownames(ColData_SNPs))
      
      DEGs_genes<-data.frame()
      for(i in 1:length(closest_SNPs$V1)){
        SNP<-closest_SNPs$V1[i]
        colData_closest<-read.table(paste0("Proyecto_Monogenic_IBD/PRS_RISK/SNPs_to_GWASgenes/Individuos_RISK_homocigotos_paraSNPs/",SNP,"_closest_SNPs.txt"), header=T)
        colData_closest<-colData_closest[order(rownames(colData_closest), decreasing = F),]
        riesgo<-colData_closest[which(colData_closest$cuantos=="2"),]
        noriesgo<-colData_closest[which(colData_closest$cuantos=="0"),]
        ColData_SNPs<-rbind(riesgo,noriesgo)
        ColData_SNPs<-ColData_SNPs[order(rownames(ColData_SNPs),decreasing=F),]
        counts_data_misind<-counts_data[,colnames(counts_data)%in%rownames(ColData_SNPs)]
        row.names(counts_data_misind)<-counts$Id
        counts_data_misind[is.na(counts_data_misind)]<-0
        ColData_SNPs$cuantos<-as.factor(ColData_SNPs$cuantos)
        dds_SNPs<-DESeqDataSetFromMatrix(countData = counts_data_misind,
                                               colData = ColData_SNPs,
                                               design = ~cuantos)
        dds_SNPs$cuantos <- relevel(dds_SNPs$cuantos, ref="0")
        
        dds_SNPs_res<-DESeq(dds_SNPs)
        res0.05<- results(dds_SNPs_res, alpha = 0.05 )
        res.0.05_naomit<-na.omit(res0.05)
        res.0.05_naomit$SNP<-rep(SNP, nrow(res.0.05_naomit))
        res.0.05_naomit<-as.data.frame(res.0.05_naomit)
        DEGs_genes<-rbind(DEGs_genes,res.0.05_naomit)
      }
      
      write.table(DEGs_genes, "Proyecto_Monogenic_IBD/PRS_RISK/SNPs_to_GWASgenes/Resultados_DESeq/Resultados_SNPs_closest_homocigotosRiesgo_homocigotosNORiesgo.txt", row.names = T, col.names = T, quote=F, sep="\t")
      bonferroni<-0.05/nrow(DEGs_genes)
      DEGs_signif<-DEGs_genes[DEGs_genes$log2FoldChange>1 & DEGs_genes$padj<bonferroni | DEGs_genes$log2FoldChange<(-1) & DEGs_genes$padj<bonferroni,]
      DEGS_delosSPS<-data.frame(table(DEGs_signif$SNP))
      
      
      Closest_snp<-read.table("Proyecto_Monogenic_IBD/PRS_RISK/SNPs_to_GWASgenes/Closest_snp_toRISK_genes.txt")
      genes_Adri<-read.table("Datos_Adri/genes_RISKyMAGMA_para_SNPs.txt", header=T)


  ##############TopCS2G

      counts<-read.table("Proyecto_Monogenic_IBD/01_NewAnalysis_MiOkProb_Norm.txt",header=T)
      TopCS2G_SNPs<-read.table("Proyecto_Monogenic_IBD/PRS_RISK/SNPs_to_GWASgenes/TopCS2G_snp_toRISK_genes_list.txt")
      
      ind_TopCS2G<-read.table("Proyecto_Monogenic_IBD/PRS_RISK/SNPs_to_GWASgenes/Genotipos_CDpatients_RISK_misSNPs_TopCS2G.txt", header=T)
      individuals<-ind_TopCS2G$V2
      
      colnames(counts)[3:ncol(counts)]<-sub("norm.","RS",colnames(counts)[3:ncol(counts)], fixed = T)
      colnames(counts)[3:ncol(counts)]<-sub(".","",colnames(counts)[3:ncol(counts)], fixed = T)
      
      counts_data<-counts[,which(colnames(counts) %in% individuals)]
      counts_data<-counts_data[,order(colnames(counts_data), decreasing = F)]
      
      all(colnames(counts_data_misind)==rownames(ColData_SNPs))
      
      DEGs_genes<-data.frame()
      for(i in 1:length(TopCS2G_SNPs$V1)){
        SNP<-TopCS2G_SNPs$V1[i]
        colData_TopCS2G<-read.table(paste0("Proyecto_Monogenic_IBD/PRS_RISK/SNPs_to_GWASgenes/Individuos_RISK_homocigotos_paraSNPs/",SNP,"_TopCS2G_SNPs.txt"), header=T)
        colData_TopCS2G<-colData_TopCS2G[order(rownames(colData_TopCS2G), decreasing = F),]
        riesgo<-colData_TopCS2G[which(colData_TopCS2G$cuantos=="2"),]
        noriesgo<-colData_TopCS2G[which(colData_TopCS2G$cuantos=="0"),]
        ColData_SNPs<-rbind(riesgo,noriesgo)
        ColData_SNPs<-ColData_SNPs[order(rownames(ColData_SNPs),decreasing=F),]
        counts_data_misind<-counts_data[,colnames(counts_data)%in%rownames(ColData_SNPs)]
        row.names(counts_data_misind)<-counts$Id
        counts_data_misind[is.na(counts_data_misind)]<-0
        ColData_SNPs$cuantos<-as.factor(ColData_SNPs$cuantos)
        dds_SNPs<-DESeqDataSetFromMatrix(countData = counts_data_misind,
                                         colData = ColData_SNPs,
                                         design = ~cuantos)
        dds_SNPs$cuantos <- relevel(dds_SNPs$cuantos, ref="0")
        
        dds_SNPs_res<-DESeq(dds_SNPs)
        res0.05<- results(dds_SNPs_res, alpha = 0.05 )
        res.0.05_naomit<-na.omit(res0.05)
        res.0.05_naomit$SNP<-rep(SNP, nrow(res.0.05_naomit))
        res.0.05_naomit<-as.data.frame(res.0.05_naomit)
        DEGs_genes<-rbind(DEGs_genes,res.0.05_naomit)
      }
      
      write.table(DEGs_genes, "Proyecto_Monogenic_IBD/PRS_RISK/SNPs_to_GWASgenes/Resultados_DESeq/Resultados_SNPs_TopCS2G_homocigotosRiesgo_homocigotosNORiesgo.txt", row.names = T, col.names = T, quote=F, sep="\t")
      bonferroni<-0.05/nrow(DEGs_genes)
      DEGs_signif<-DEGs_genes[DEGs_genes$log2FoldChange>1 & DEGs_genes$padj<bonferroni| DEGs_genes$log2FoldChange<(-1) & DEGs_genes$padj<bonferroni,]
      DEGS_delosSPS<-data.frame(table(DEGs_signif$SNP))
      
