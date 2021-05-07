library(monocle)
library(Seurat)
library(tidyverse)
library('extrafont')

data <- Read10X(data.dir = str_c("2.2.filtered_feature_bc_matrix/"))
flabel<-"1c"

pbmc <- CreateSeuratObject(counts = data, project = flabel, min.cells = 20 ,min.genes = 0)
print(tail(pbmc@meta.data))
mito.genes <- grep(pattern = "^mt-", x = rownames(x = pbmc@assays$RNA), value = TRUE)
mito.genes
pbmc <- PercentageFeatureSet(pbmc, pattern = "^mt-", col.name = "percent.mt")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < ifelse(max(pbmc$nFeature_RNA)*0.9>=4000,4000,max(pbmc$nFeature_RNA)*0.9) & percent.mt < 5)
# Normalize
pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs=50)
pbmc <- RunUMAP(pbmc, dims = 1:30,verbose = FALSE)#,umap.method = "umap-learn",metric = "correlation")#min.dist  controls how tightly the embedding is allowed compress points together
pheno<-read_csv("pheno.csv")
pheno$index<-gsub("b'","",pheno$index)
pheno$index<-gsub("'","",pheno$index)
meta<-pbmc@meta.data
meta2<-meta
meta2$si<-row.names(meta2)
meta2<-merge(meta2,pheno,by.x="si",by.y="index",all.x=TRUE)
meta2$leiden[is.na(meta2$leiden)]<-"unk"
meta2<-meta2[match(row.names(meta),meta2$si),]
row.names(meta2)<-row.names(meta)
pbmc@meta.data<-meta2

write_tsv(meta2,"sample.info.tsv")
p2<- DimPlot(pbmc, reduction = "umap",label=TRUE,group.by = "leiden")
ggsave(p2,device="pdf",filename=str_c("plots/",flabel,".UMAP.plot2.pdf"),width=6,height=6)

counts<-pbmc@assays$SCT@counts  
pd <-  pbmc@meta.data
pd$leiden<-gsub("b'","",pd$leiden)
pd$leiden<-gsub("'","",pd$leiden)

selcol<-which(pd$leiden %in% c("pMN_4","pMN > Motor neuron_7","pMN_2","OPC_11"))
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

## trajectory 

CDE <- differentialGeneTest(Obj[expressed_genes], fullModelFormulaStr = '~leiden', cores = detectCores()/ 2)
CDE<-arrange(CDE,pval)

TF<-read_tsv("dre_TF.txt")

sCDE<-CDE[1:500,]
row.names(sCDE)[row.names(sCDE) %in% TF$Symbol]
seltf<- row.names(sCDE)[row.names(sCDE) %in% TF$Symbol]
seltfde<-sCDE[row.names(sCDE) %in% TF$Symbol,]
seltfde$TF<-row.names(seltfde)
write_tsv(seltfde,"sel.tf.de.info.tsv")

ordering_genes <- row.names(CDE)[order(CDE$qval)][1:1000]
Obj1 <- setOrderingFilter(Obj, ordering_genes = c(ordering_genes))
Obj1 <- reduceDimension(Obj1, method = 'ICA') #
Obj1 <- orderCells(Obj1)
#num<-table(pData(Obj)$leiden,pData(Obj)$State)


GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$leiden)[,"pMN_4"]
    return(as.numeric(names(T0_counts)[which
          (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
Obj1.1<- orderCells(Obj1, root_state = GM_state(Obj1))
pdf("figcd.pdf",width=6,height=4,family="Arial")  
plot_cell_trajectory(Obj1, color_by = 'leiden',show_tree=FALSE)
plot_cell_trajectory(Obj1.1, color_by = 'Pseudotime',show_tree=FALSE)
dev.off()


