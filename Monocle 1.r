### Fige.g.f
counts<-pbmc@assays$SCT@counts  # raw data counts 更适合于monocle 的分析
pd <-  pbmc@meta.data
pd$leiden<-gsub("b'","",pd$leiden)
pd$leiden<-gsub("'","",pd$leiden)

selcol<-which(pd$leiden %in% c("pMN > Motor neuron_7","pMN_2","OPC_11"))
counts<-counts[,selcol]
pd<-pd[selcol,]

fData <- data.frame(gene_short_name = row.names(counts), row.names = row.names(counts))
colnames(pd)

pd <- new("AnnotatedDataFrame", data = pd) #
fd <- new("AnnotatedDataFrame", data =fData)

Obj <- newCellDataSet(
as(counts, "sparseMatrix"),
phenoData = pd,
featureData = fd,lowerDetectionLimit = 0.1,
expressionFamily = negbinomial.size() #gaussianff()
)

Obj <- estimateSizeFactors(Obj)
Obj <- estimateDispersions(Obj)


###Unsupervised cell clustering 
disp_table <- dispersionTable(Obj)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
Obj <- setOrderingFilter(Obj, unsup_clustering_genes$gene_id)
plot_ordering_genes(Obj)
Obj <- detectGenes(Obj, min_expr = 0.1)
expressed_genes <-  row.names(subset(fData(Obj),num_cells_expressed >= 10))
Obj <- Obj[expressed_genes,]

Obj<-reduceDimension(Obj, max_components=2, norm_method = 'log', reduction_method = 'tSNE', num_dim =3,  verbose = F) #reduction_method
Obj <- clusterCells(Obj,num_clusters=6) # 
plot_cell_clusters(Obj, color_by = 'as.factor(Cluster)')
plot_cell_clusters(Obj, color_by = 'as.factor(leiden)')
num<-table(pData(Obj)$leiden,pData(Obj)$Cluster)

#We also provide the decision plot for users to check the Ρ, Δ for each cell and decide the threshold for defining the cell clusters.
plot_rho_delta(Obj, rho_threshold = 2, delta_threshold = 4 )
Obj <- clusterCells(Obj,
                 rho_threshold = 28,
                 delta_threshold = 5,
                 skip_rho_sigma = T,
                 verbose = F)
#save(pheno,sraw,snorm,file="refine.cluster.Rdata")
plot_cell_clusters(Obj, color_by = 'as.factor(Cluster)')
plot_cell_clusters(Obj, color_by = 'as.factor(leiden)')

CDE <- differentialGeneTest(Obj[expressed_genes], fullModelFormulaStr = '~leiden', cores = detectCores()/ 2)
CDE<-arrange(CDE,pval)

TF<-read_tsv("/home/dell/database/up_stream/dre_TF.txt")

sCDE<-CDE[1:500,]
row.names(sCDE)[row.names(sCDE) %in% TF$Symbol]
seltf<- row.names(sCDE)[row.names(sCDE) %in% TF$Symbol]
sel<-sCDE[row.names(sCDE) %in% TF$Symbol,]
sel$TF<-row.names(seltfde)

ordering_genes <- row.names(CDE)[order(CDE$qval)][1:1000]
Obj1 <- setOrderingFilter(Obj, ordering_genes = c(ordering_genes))
Obj1 <- reduceDimension(Obj1, method = 'ICA') # ICA or DDRTREE # method 很重要。
Obj1 <- orderCells(Obj1)
#num<-table(pData(Obj)$leiden,pData(Obj)$State)


GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$leiden)[,"pMN_2"]
    return(as.numeric(names(T0_counts)[which
          (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
Obj1<- orderCells(Obj1, root_state =3)


sel<-sel[1:45,]

source("plot_genes_branched_heatmap.R")

pdf(str_c("plots/","figG.pdf"),family="Arial",width=7,height=7)
plot_genes_branched_heatmap(Obj1[sel$gene_short_name,],
                                          branch_point = 1,
                                          num_clusters = 4,
                                          cores = 1,
                                          use_gene_short_name = T,
                                          show_rownames = T,fontsize=12,branch_labels = c("OPC_11", "pMN>Motor neuron_7"))
dev.off()


blast_genes<-c("myt1b", "myt1a", "mnx1", "neurod1", "isl1", "olig1", "sox10", "nkx2.2a")

pdf(str_c("plots/","figE.1.pdf"),family="Arial",width=7,height=7)
plot_genes_branched_pseudotime(Obj1[blast_genes,],
                       branch_point = 1,
                       color_by = "leiden",panel_order=blast_genes,
                       ncol = 3)+mytheme+theme(legend.position="none",legend.title= element_blank(),strip.background=element_blank(),panel.border=element_rect(fill = NA,colour='black',size=0.8),strip.text=element_text(size=18,colour="black"),title=element_text(face="bold",size=12))+xlab("Pseudotime")

dev.off()


pdf(str_c("plots/","figE.2.pdf"),family="Arial",width=8,height=7)
plot_genes_branched_pseudotime(Obj1[blast_genes,],
                       branch_point = 1,
                       color_by = "leiden",
                       ncol = 3)+mytheme+theme(legend.position="right",legend.title= element_blank(),strip.background=element_blank(),panel.border=element_rect(fill = NA,colour='black',size=0.8),strip.text=element_text(size=18,colour="black"),title=element_text(face="bold",size=12))+xlab("Pseudotime")

dev.off()

blast_genes <- row.names(subset(fData(Obj1),gene_short_name %in% c("myt1a","myt1b")))

pdf(str_c("plots/","figF.pdf"),family="Arial",width=4,height=4)

plot_genes_violin(Obj1[blast_genes,],
                  grouping = "leiden",
                 min_expr = 1.0,color_by = "leiden")+mytheme+theme(legend.position="none",legend.title= element_blank(),strip.background=element_blank(),panel.border=element_rect(fill = NA,colour='black',size=0.8),strip.text=element_text(size=12,colour="black"))+xlab("")
dev.off()