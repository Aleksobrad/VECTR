
###add samhita metadata and CD4/CD8 metadata
####three turnover analyses (early, late, very late), three treatment analyses (10to15, 11to17, 11to14)
vp_metadat_cd4cd8=read.csv("/Users/aleksandar/Dropbox/Turnover_Treatment_Paper/objects_metadata/viper_object_metadata.csv")
vp_metadat=read.csv("/Users/aleksandar/Dropbox/Turnover_Treatment_Paper/objects_metadata/comb.vp.metadata.annotated.for.ao.csv")
rownames(vp_metadat)=vp_metadat[,1]
rownames(vp_metadat_cd4cd8)=vp_metadat_cd4cd8[,1]

vp.meta.seurat=readRDS("~/Documents/Documents/MED SCHOOL/summer 2019/arnold_han_collaboration/vp.meta.seurat.alldata.mc38.batchcorrected.v2.rds")
vp.meta.seurat$TO1_time=vp_metadat$TO1_time
vp.meta.seurat$TO2_type=vp_metadat$TO2_type
vp.meta.seurat$TX1_condition=vp_metadat$TX1_condition
vp.meta.seurat$TX2_time=vp_metadat$TX2_time
vp.meta.seurat$ao_10to15=vp_metadat$ao_10to15
vp.meta.seurat$ao_11to14=vp_metadat$ao_11to14
vp.meta.seurat$ao_11to17=vp_metadat$ao_11to17
vp.meta.seurat$ao_IMP_TX=vp_metadat$ao_IMP_TX
vp.meta.seurat$ao_IMP_TO=vp_metadat$ao_IMP_TO
vp.meta.seurat$CD4_CD8_assignment=vp_metadat_cd4cd8$CD4_CD8_assignment

vp.meta.seurat.CD4=vp.meta.seurat[,which(vp.meta.seurat$CD4_CD8_assignment == "CD4")]
vp.meta.seurat.CD8=vp.meta.seurat[,which(vp.meta.seurat$CD4_CD8_assignment == "CD8")]

vp.meta.seurat.CD4 <- RunPCA(vp.meta.seurat.CD4,features=rownames(vp.meta.seurat.CD4))
vp.meta.seurat.CD4 <- RunUMAP(vp.meta.seurat.CD4, dims = 1:50, verbose = FALSE,umap.method="umap-learn",metric="correlation")
vp.meta.seurat.CD4 <- FindNeighbors(vp.meta.seurat.CD4, dims = 1:50, verbose = FALSE)
vp.meta.seurat.CD4 <- FindClusters(vp.meta.seurat.CD4, resolution=seq(0.01,1,by=0.01), verbose = FALSE,algorithm=1) 
clust=vp.meta.seurat.CD4@meta.data[,which(grepl("RNA_snn_res.",colnames(vp.meta.seurat.CD4@meta.data)))]
mat=vp.meta.seurat.CD4@assays$RNA@scale.data
out=sil_subsample(mat,clust)
means=out[[1]]
sd=out[[2]]
x=seq(0.01,1,by=0.01)
errbar(x,means,means+sd,means-sd,ylab="mean silhouette score",xlab="resolution parameter")
lines(x,means)
x=x[20:length(x)]
means=means[20:length(means)]
best=tail(x[which(means==max(means))],n=1)
legend("topright",paste("Best",best,sep = " = "))
vp.meta.seurat.CD4$seurat_clusters=vp.meta.seurat.CD4@meta.data[,which(colnames(vp.meta.seurat.CD4@meta.data)==paste("RNA_snn_res.",best,sep=""))]
Idents(vp.meta.seurat.CD4) <- "seurat_clusters"
plot(DimPlot(vp.meta.seurat.CD4, reduction = "umap",label = TRUE) + NoLegend())
saveRDS(vp.meta.seurat.CD4, file = "~/Documents/Documents/MED SCHOOL/summer 2019/arnold_han_collaboration/vp.meta.seurat.alldata.mc38.batchcorrected.CD4.v3.rds")
#LDA
model <- lda(t(as.matrix(vp.meta.seurat.CD4@assays$RNA@counts)),vp.meta.seurat.CD4$seurat_clusters)
lda_cd4=predict(model)$x
plot(lda_cd4[,1:2],pch=16,col=vp.meta.seurat.CD4$seurat_clusters)
ggplot(data.frame(UMAP_1=lda_cd4[,1],UMAP_2=lda_cd4[,2],cluster=vp.meta.seurat.CD4$seurat_clusters),aes(x=UMAP_1,y=UMAP_2,color=cluster))+geom_point(size=0.01)+geom_density_2d(color="black")+theme_bw()+theme(panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))+xlim(-7,15)+ylim(-7,7)
colnames(lda_cd4)=paste("LD",colnames(lda_cd4),sep="_")
vp.meta.seurat.CD4@reductions[["lda"]]=CreateDimReducObject(embeddings = lda_cd4, key = "LD_", assay = DefaultAssay(vp.meta.seurat.CD4))
DimPlot(vp.meta.seurat.CD4,reduction = "lda",label = T,repel=T)
saveRDS(vp.meta.seurat.CD4, file = "~/Documents/Documents/MED SCHOOL/summer 2019/arnold_han_collaboration/vp.meta.seurat.alldata.mc38.batchcorrected.CD4.v3.rds")


vp.meta.seurat.CD8 <- RunPCA(vp.meta.seurat.CD8,features=rownames(vp.meta.seurat.CD8))
vp.meta.seurat.CD8 <- RunUMAP(vp.meta.seurat.CD8, dims = 1:50, verbose = FALSE,umap.method="umap-learn",metric="correlation")
vp.meta.seurat.CD8 <- FindNeighbors(vp.meta.seurat.CD8, dims = 1:50, verbose = FALSE)
vp.meta.seurat.CD8 <- FindClusters(vp.meta.seurat.CD8, resolution=seq(0.01,1,by=0.01), verbose = FALSE,algorithm=1) 
clust=vp.meta.seurat.CD8@meta.data[,which(grepl("RNA_snn_res.",colnames(vp.meta.seurat.CD8@meta.data)))]
mat=vp.meta.seurat.CD8@assays$RNA@scale.data
out=sil_subsample(mat,clust)
means=out[[1]]
sd=out[[2]]
x=seq(0.01,1,by=0.01)
errbar(x,means,means+sd,means-sd,ylab="mean silhouette score",xlab="resolution parameter")
lines(x,means)
x=x[20:length(x)]
means=means[20:length(means)]
best=tail(x[which(means==max(means))],n=1)
legend("topright",paste("Best",best,sep = " = "))
vp.meta.seurat.CD8$seurat_clusters=vp.meta.seurat.CD8@meta.data[,which(colnames(vp.meta.seurat.CD8@meta.data)==paste("RNA_snn_res.",best,sep=""))]
Idents(vp.meta.seurat.CD8) <- "seurat_clusters"
plot(DimPlot(vp.meta.seurat.CD8, reduction = "umap",label = TRUE) + NoLegend())
saveRDS(vp.meta.seurat.CD8, file = "~/Documents/Documents/MED SCHOOL/summer 2019/arnold_han_collaboration/vp.meta.seurat.alldata.mc38.batchcorrected.CD8.v3.rds")
#LDA
model <- lda(t(as.matrix(vp.meta.seurat.CD8@assays$RNA@counts)),vp.meta.seurat.CD8$seurat_clusters)
lda_cd8=predict(model)$x
plot(lda_cd8[,1:2],pch=16,col=vp.meta.seurat.CD8$seurat_clusters)
ggplot(data.frame(UMAP_1=lda_cd8[,1],UMAP_2=lda_cd8[,2],cluster=vp.meta.seurat.CD8$seurat_clusters),aes(x=UMAP_1,y=UMAP_2,color=cluster))+geom_point(size=0.01)+geom_density_2d(color="black")+theme_bw()+theme(panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))
colnames(lda_cd8)=paste("LD",colnames(lda_cd8),sep="_")
vp.meta.seurat.CD8@reductions[["lda"]]=CreateDimReducObject(embeddings = lda_cd8, key = "LD_", assay = DefaultAssay(vp.meta.seurat.CD8))
DimPlot(vp.meta.seurat.CD8,reduction = "lda",label = T,repel=T)
saveRDS(vp.meta.seurat.CD8, file = "~/Documents/Documents/MED SCHOOL/summer 2019/arnold_han_collaboration/vp.meta.seurat.alldata.mc38.batchcorrected.CD8.v3.rds")



vp.meta.seurat.CD4=readRDS("~/Documents/Documents/MED SCHOOL/summer 2019/arnold_han_collaboration/vp.meta.seurat.alldata.mc38.batchcorrected.CD4.v3.rds")
vp.meta.seurat.CD8=readRDS("~/Documents/Documents/MED SCHOOL/summer 2019/arnold_han_collaboration/vp.meta.seurat.alldata.mc38.batchcorrected.CD8.v3.rds")
vp.meta.seurat.CD4.v4=vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$seurat_clusters %in% c(0,1,2,3,4,6,8))]
vp.meta.seurat.CD8.v4=merge(vp.meta.seurat.CD8,vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$seurat_clusters %in% c(5,7))])
vp.meta.seurat.CD8.v4@assays$RNA@scale.data=as.matrix(vp.meta.seurat.CD8.v4@assays$RNA@counts)
vp.meta.seurat.CD8.v4 <- RunPCA(vp.meta.seurat.CD8.v4,features=rownames(vp.meta.seurat.CD8.v4))
vp.meta.seurat.CD8.v4 <- RunUMAP(vp.meta.seurat.CD8.v4, dims = 1:50, verbose = FALSE,umap.method="umap-learn",metric="correlation")
vp.meta.seurat.CD8.v4 <- FindNeighbors(vp.meta.seurat.CD8.v4, dims = 1:50, verbose = FALSE)
vp.meta.seurat.CD8.v4 <- FindClusters(vp.meta.seurat.CD8.v4, resolution=seq(0.01,1,by=0.01), verbose = FALSE,algorithm=1) 
clust=vp.meta.seurat.CD8.v4@meta.data[,which(grepl("RNA_snn_res.",colnames(vp.meta.seurat.CD8.v4@meta.data)))]
mat=vp.meta.seurat.CD8.v4@assays$RNA@scale.data
out=sil_subsample(mat,clust)
means=out[[1]]
sd=out[[2]]
x=seq(0.01,1,by=0.01)
errbar(x,means,means+sd,means-sd,ylab="mean silhouette score",xlab="resolution parameter")
lines(x,means)
x=x[20:length(x)]
means=means[20:length(means)]
best=tail(x[which(means==max(means))],n=1)
legend("topright",paste("Best",best,sep = " = "))
vp.meta.seurat.CD8.v4$seurat_clusters=vp.meta.seurat.CD8.v4@meta.data[,which(colnames(vp.meta.seurat.CD8.v4@meta.data)==paste("RNA_snn_res.",best,sep=""))]
Idents(vp.meta.seurat.CD8.v4) <- "seurat_clusters"
plot(DimPlot(vp.meta.seurat.CD8.v4, reduction = "umap",label = TRUE) + NoLegend())
model <- lda(t(as.matrix(vp.meta.seurat.CD8.v4@assays$RNA@counts)),vp.meta.seurat.CD8.v4$seurat_clusters)
lda_cd8=predict(model)$x
plot(lda_cd8[,1:2],pch=16,col=vp.meta.seurat.CD8.v4$seurat_clusters)
ggplot(data.frame(UMAP_1=lda_cd8[,1],UMAP_2=lda_cd8[,2],cluster=vp.meta.seurat.CD8.v4$seurat_clusters),aes(x=UMAP_1,y=UMAP_2,color=cluster))+geom_point(size=0.01)+geom_density_2d(color="black")+theme_bw()+theme(panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))
colnames(lda_cd8)=paste("LD",colnames(lda_cd8),sep="_")
vp.meta.seurat.CD8.v4@reductions[["lda"]]=CreateDimReducObject(embeddings = lda_cd8, key = "LD_", assay = DefaultAssay(vp.meta.seurat.CD8.v4))
DimPlot(vp.meta.seurat.CD8.v4,reduction = "lda",label = T,repel=T)
saveRDS(vp.meta.seurat.CD4.v4,"~/Documents/Documents/MED SCHOOL/summer 2019/arnold_han_collaboration/vp.meta.seurat.alldata.mc38.batchcorrected.CD4.v4.rds")
saveRDS(vp.meta.seurat.CD8.v4,"~/Documents/Documents/MED SCHOOL/summer 2019/arnold_han_collaboration/vp.meta.seurat.alldata.mc38.batchcorrected.CD8.v4.rds")


##subsampling boxplots-- frequency and tcr-sharing...



###generate all the plots!
arnoldhan_data.integrated=readRDS("~/Documents/Documents/MED SCHOOL/summer 2019/arnold_han_collaboration/combined_batch_corrected_arnoldhan_data.integrated.mc38.rds")
arnoldhan_data.integrated=arnoldhan_data.integrated[,which(!(arnoldhan_data.integrated$sample %in% c("low_count","multiplet","noisy_hash")))]
vp.meta.seurat.CD8=readRDS("~/Documents/Documents/MED SCHOOL/summer 2019/arnold_han_collaboration/vp.meta.seurat.alldata.mc38.batchcorrected.CD8.v4.rds")
vp.meta.seurat.CD4=readRDS("~/Documents/Documents/MED SCHOOL/summer 2019/arnold_han_collaboration/vp.meta.seurat.alldata.mc38.batchcorrected.CD4.v4.rds")
vp.meta.seurat.CD4=vp.meta.seurat.CD4[,which(colnames(vp.meta.seurat.CD4) %in% colnames(arnoldhan_data.integrated))]
vp.meta.seurat.CD8=vp.meta.seurat.CD8[,which(colnames(vp.meta.seurat.CD8) %in% colnames(arnoldhan_data.integrated))]
saveRDS(arnoldhan_data.integrated,"~/Documents/Documents/MED SCHOOL/summer 2019/arnold_han_collaboration/combined_batch_corrected_arnoldhan_data.integrated.mc38.rds")
saveRDS(vp.meta.seurat.CD4,"~/Documents/Documents/MED SCHOOL/summer 2019/arnold_han_collaboration/vp.meta.seurat.alldata.mc38.batchcorrected.CD4.v4.rds")
saveRDS(vp.meta.seurat.CD8,"~/Documents/Documents/MED SCHOOL/summer 2019/arnold_han_collaboration/vp.meta.seurat.alldata.mc38.batchcorrected.CD8.v4.rds")
vp.meta.seurat.CD4$model="mc38"
vp.meta.seurat.CD8$model="mc38"
vp.meta.seurat.CD4$cdr3b_nt=arnoldhan_data.integrated$cdr3b_nt[colnames(vp.meta.seurat.CD4)]
vp.meta.seurat.CD8$cdr3b_nt=arnoldhan_data.integrated$cdr3b_nt[colnames(vp.meta.seurat.CD8)]
vp.meta.seurat.CD4$cdr3b=arnoldhan_data.integrated$cdr3b[colnames(vp.meta.seurat.CD4)]
vp.meta.seurat.CD8$cdr3b=arnoldhan_data.integrated$cdr3b[colnames(vp.meta.seurat.CD8)]
vp.meta.seurat.CD4$timepoint=arnoldhan_data.integrated$timepoint[colnames(vp.meta.seurat.CD4)]
vp.meta.seurat.CD8$timepoint=arnoldhan_data.integrated$timepoint[colnames(vp.meta.seurat.CD8)]
vp.meta.seurat.CD4$experiment=arnoldhan_data.integrated$experiment[colnames(vp.meta.seurat.CD4)]
vp.meta.seurat.CD8$experiment=arnoldhan_data.integrated$experiment[colnames(vp.meta.seurat.CD8)]
vp.meta.seurat.CD4$condition=arnoldhan_data.integrated$treatment[colnames(vp.meta.seurat.CD4)]
vp.meta.seurat.CD8$condition=arnoldhan_data.integrated$treatment[colnames(vp.meta.seurat.CD8)]
vp.meta.seurat.CD4$sample=arnoldhan_data.integrated$sample[colnames(vp.meta.seurat.CD4)]
vp.meta.seurat.CD8$sample=arnoldhan_data.integrated$sample[colnames(vp.meta.seurat.CD8)]
vp.meta.seurat.CD4$revised_timepoint=arnoldhan_data.integrated$revised_timepoint[colnames(vp.meta.seurat.CD4)]
vp.meta.seurat.CD8$revised_timepoint=arnoldhan_data.integrated$revised_timepoint[colnames(vp.meta.seurat.CD8)]
Idents(vp.meta.seurat.CD4)=vp.meta.seurat.CD4$seurat_clusters
Idents(vp.meta.seurat.CD8)=vp.meta.seurat.CD8$seurat_clusters
#early: day10-11 PRE, day13-14 POST, day16-17 POST, day19-20 POST, day 22-23 POST
#late: day16-17 PRE, day19-20 POST, day 22-23 POST, day 25-26 POST, day 28-29 POST
vp.meta.seurat.CD4$revised_timepoint[which(vp.meta.seurat.CD4$experiment=="eto" & vp.meta.seurat.CD4$condition=="D0")]="DAY10-11"
vp.meta.seurat.CD4$revised_timepoint[which(vp.meta.seurat.CD4$experiment=="eto" & vp.meta.seurat.CD4$condition=="D3")]="DAY13-14"
vp.meta.seurat.CD4$revised_timepoint[which(vp.meta.seurat.CD4$experiment=="eto" & vp.meta.seurat.CD4$condition=="D6")]="DAY16-17"
vp.meta.seurat.CD4$revised_timepoint[which(vp.meta.seurat.CD4$experiment=="eto" & vp.meta.seurat.CD4$condition=="D9")]="DAY19-20"
vp.meta.seurat.CD4$revised_timepoint[which(vp.meta.seurat.CD4$experiment=="eto" & vp.meta.seurat.CD4$condition=="D12")]="DAY22-23"
vp.meta.seurat.CD4$revised_timepoint[which(vp.meta.seurat.CD4$experiment=="lto" & vp.meta.seurat.CD4$condition=="D0")]="DAY16-17"
vp.meta.seurat.CD4$revised_timepoint[which(vp.meta.seurat.CD4$experiment=="lto" & vp.meta.seurat.CD4$condition=="D3")]="DAY19-20"
vp.meta.seurat.CD4$revised_timepoint[which(vp.meta.seurat.CD4$experiment=="lto" & vp.meta.seurat.CD4$condition=="D6")]="DAY22-23"
vp.meta.seurat.CD4$revised_timepoint[which(vp.meta.seurat.CD4$experiment=="lto" & vp.meta.seurat.CD4$condition=="D9")]="DAY25-26"
vp.meta.seurat.CD4$revised_timepoint[which(vp.meta.seurat.CD4$experiment=="lto" & vp.meta.seurat.CD4$condition=="D12")]="DAY28-29"
vp.meta.seurat.CD8$revised_timepoint[which(vp.meta.seurat.CD8$experiment=="eto" & vp.meta.seurat.CD8$condition=="D0")]="DAY10-11"
vp.meta.seurat.CD8$revised_timepoint[which(vp.meta.seurat.CD8$experiment=="eto" & vp.meta.seurat.CD8$condition=="D3")]="DAY13-14"
vp.meta.seurat.CD8$revised_timepoint[which(vp.meta.seurat.CD8$experiment=="eto" & vp.meta.seurat.CD8$condition=="D6")]="DAY16-17"
vp.meta.seurat.CD8$revised_timepoint[which(vp.meta.seurat.CD8$experiment=="eto" & vp.meta.seurat.CD8$condition=="D9")]="DAY19-20"
vp.meta.seurat.CD8$revised_timepoint[which(vp.meta.seurat.CD8$experiment=="eto" & vp.meta.seurat.CD8$condition=="D12")]="DAY22-23"
vp.meta.seurat.CD8$revised_timepoint[which(vp.meta.seurat.CD8$experiment=="lto" & vp.meta.seurat.CD8$condition=="D0")]="DAY16-17"
vp.meta.seurat.CD8$revised_timepoint[which(vp.meta.seurat.CD8$experiment=="lto" & vp.meta.seurat.CD8$condition=="D3")]="DAY19-20"
vp.meta.seurat.CD8$revised_timepoint[which(vp.meta.seurat.CD8$experiment=="lto" & vp.meta.seurat.CD8$condition=="D6")]="DAY22-23"
vp.meta.seurat.CD8$revised_timepoint[which(vp.meta.seurat.CD8$experiment=="lto" & vp.meta.seurat.CD8$condition=="D9")]="DAY25-26"
vp.meta.seurat.CD8$revised_timepoint[which(vp.meta.seurat.CD8$experiment=="lto" & vp.meta.seurat.CD8$condition=="D12")]="DAY28-29"

table(vp.meta.seurat.CD8$TO1_time)
table(vp.meta.seurat.CD8$TO2_type)
table(vp.meta.seurat.CD8$TX2_time)
table(vp.meta.seurat.CD8$TX1_condition)

##density plots
DimPlot(vp.meta.seurat.CD4,reduction = "lda",label = T)
DimPlot(vp.meta.seurat.CD8,reduction = "lda",label=T)
ggplot(data.frame(LD_1=vp.meta.seurat.CD8@reductions$lda@cell.embeddings[,1],LD_2=vp.meta.seurat.CD8@reductions$lda@cell.embeddings[,2],cluster=vp.meta.seurat.CD8$seurat_clusters),aes(x=LD_1,y=LD_2,color=cluster))+geom_point(size=0.01)+geom_density_2d(color="black",binwidth=0.0005)+theme_bw()+theme(panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))
ggplot(data.frame(LD_1=vp.meta.seurat.CD4@reductions$lda@cell.embeddings[,1],LD_2=vp.meta.seurat.CD4@reductions$lda@cell.embeddings[,2],cluster=vp.meta.seurat.CD4$seurat_clusters),aes(x=LD_1,y=LD_2,color=cluster))+geom_point(size=0.01)+geom_density_2d(color="black")+theme_bw()+theme(panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))
markers.cd8.viper <- FindAllMarkers(vp.meta.seurat.CD8, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1,test.use = "t")
top10.CD8.viper <- markers.cd8.viper %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
geneHeatmap_plot(vp.meta.seurat.CD8@assays$RNA@data,vp.meta.seurat.CD8$seurat_clusters,top10.CD8.viper$gene,genes_by_cluster=T,scaled=T,n_top_genes_per_cluster = 5)
markers.cd4.viper <- FindAllMarkers(vp.meta.seurat.CD4, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1,test.use = "t")
top10.CD4.viper <- markers.cd4.viper %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
geneHeatmap_plot(vp.meta.seurat.CD4@assays$RNA@data,vp.meta.seurat.CD4$seurat_clusters,top10.CD4.viper$gene,genes_by_cluster=T,scaled=T,n_top_genes_per_cluster = 5)


exp.meta.seurat.CD4=arnoldhan_data.integrated[,colnames(vp.meta.seurat.CD4)]
exp.meta.seurat.CD4$seurat_clusters=vp.meta.seurat.CD4$seurat_clusters
Idents(exp.meta.seurat.CD4)="seurat_clusters"
markers.cd4.exp <- FindAllMarkers(exp.meta.seurat.CD4, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1,test.use = "MAST")
top10.CD4.viper <- markers.cd4.exp %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
geneHeatmap_plot(exp.meta.seurat.CD4@assays$SCT@counts,exp.meta.seurat.CD4$seurat_clusters,top10.CD4.viper$gene,genes_by_cluster=T,scaled=F,n_top_genes_per_cluster = 5)
exp.meta.seurat.CD8=arnoldhan_data.integrated[,colnames(vp.meta.seurat.CD8)]
exp.meta.seurat.CD8$seurat_clusters=vp.meta.seurat.CD8$seurat_clusters
Idents(exp.meta.seurat.CD8)="seurat_clusters"
markers.cd8.exp <- FindAllMarkers(exp.meta.seurat.CD8, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1,test.use = "MAST")
top10.CD8.viper <- markers.cd8.exp %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
geneHeatmap_plot(exp.meta.seurat.CD8@assays$SCT@counts,exp.meta.seurat.CD8$seurat_clusters,top10.CD8.viper$gene,genes_by_cluster=T,scaled=F,n_top_genes_per_cluster = 5)
saveRDS(exp.meta.seurat.CD4,"~/Dropbox/Turnover_Treatment_Paper/data_and_code_for_samhita/gene_expression_cd4.rds")
saveRDS(exp.meta.seurat.CD8,"~/Dropbox/Turnover_Treatment_Paper/data_and_code_for_samhita/gene_expression_cd8.rds")
saveRDS(markers.cd4.exp,"~/Dropbox/Turnover_Treatment_Paper/data_and_code_for_samhita/differential_gene_expression_table_cd4.rds")
saveRDS(markers.cd8.exp,"~/Dropbox/Turnover_Treatment_Paper/data_and_code_for_samhita/differential_gene_expression_table_cd8.rds")
saveRDS(markers.cd4.viper,"~/Dropbox/Turnover_Treatment_Paper/data_and_code_for_samhita/differential_viper_table_cd4.rds")
saveRDS(markers.cd8.viper,"~/Dropbox/Turnover_Treatment_Paper/data_and_code_for_samhita/differential_viper_table_cd8.rds")

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
s.genes=human2mouse(s.genes, db = homologene::homologeneData)[,2]
g2m.genes=human2mouse(g2m.genes, db = homologene::homologeneData)[,2]
exp.meta.seurat.CD4 <- CellCycleScoring(exp.meta.seurat.CD4, s.features = s.genes, g2m.features = g2m.genes)
exp.meta.seurat.CD8 <- CellCycleScoring(exp.meta.seurat.CD8, s.features = s.genes, g2m.features = g2m.genes)
VlnPlot(exp.meta.seurat.CD4,"G2M.Score",pt.size=0)
VlnPlot(exp.meta.seurat.CD8,"G2M.Score",pt.size=0)

arnoldhan_data.integrated<-CellCycleScoring(arnoldhan_data.integrated, s.features = s.genes, g2m.features = g2m.genes)
vp.meta.seurat.CD4$G2M.Score=arnoldhan_data.integrated$G2M.Score[colnames(vp.meta.seurat.CD4)]
vp.meta.seurat.CD8$G2M.Score=arnoldhan_data.integrated$G2M.Score[colnames(vp.meta.seurat.CD8)]

VlnPlot(vp.meta.seurat.CD4[,vp.meta.seurat.CD4$TO1_time %in% c("DAY13-14","DAY16-17","DAY19-20","DAY22-23","DAY25-26")],"G2M.Score",pt.size=0,group.by="TO1_time")+NoLegend()
VlnPlot(vp.meta.seurat.CD8[,vp.meta.seurat.CD8$TO1_time %in% c("DAY13-14","DAY16-17","DAY19-20","DAY22-23","DAY25-26")],"G2M.Score",pt.size=0,group.by="TO1_time")+NoLegend()
VlnPlot(vp.meta.seurat.CD4[,vp.meta.seurat.CD4$TX1_condition2 %in% c("PBS","P","C ","PC")],"G2M.Score",group.by="TX1_condition2",pt.size=0)+NoLegend()
VlnPlot(vp.meta.seurat.CD8[,vp.meta.seurat.CD8$TX1_condition2 %in% c("PBS","P","C ","PC")],"G2M.Score",group.by="TX1_condition2",pt.size=0)+NoLegend()

VlnPlot(vp.meta.seurat.CD8[,vp.meta.seurat.CD8$TX1_condition2 %in% c("PRE","PBS","P","C","PC")],"G2M.Score",group.by="TX1_condition2",pt.size=0,log=T)
VlnPlot(vp.meta.seurat.CD8[,vp.meta.seurat.CD8$TX1_condition2 %in% c("PRE","PBS","P","C","PC") & vp.meta.seurat.CD8$seurat_clusters==0],"G2M.Score",group.by="TX1_condition2",pt.size=0,log=T)
VlnPlot(vp.meta.seurat.CD8[,vp.meta.seurat.CD8$TX1_condition2 %in% c("PRE","PBS","P","C","PC") & vp.meta.seurat.CD8$seurat_clusters==1],"G2M.Score",group.by="TX1_condition2",pt.size=0,log=T)

s=vp.meta.seurat.CD8[,vp.meta.seurat.CD8$TX1_condition2 %in% c("PBS","P","C ","PC")]
s2=vp.meta.seurat.CD4.2[,vp.meta.seurat.CD4.2$TX1_condition2 %in% c("PBS","P","C ","PC")]
x=data.frame(condition=s@meta.data$TX1_condition2,G2M.Score=s@meta.data$G2M.Score,Cluster=s@meta.data$seurat_clusters)
x$Cluster="CD8"
x2=data.frame(condition=s2@meta.data$TX1_condition2,G2M.Score=s2@meta.data$G2M.Score,Cluster=s2@meta.data$seurat_clusters)
x2$Cluster="CD4"
x=rbind(x,x2)
p= ggplot(x, aes(x=condition, y=G2M.Score,fill=Cluster)) + 
  geom_violin(trim=T,scale="width")+#coord_trans(y = "log10")+
  #scale_y_continuous(trans='log10') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("G2M Score by Treatment Per Cluster")
p + theme(
  plot.title = element_text(size=14, face="bold"),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))

s=vp.meta.seurat.CD8[,vp.meta.seurat.CD8$TO1_time %in% c("DAY10-11","DAY13-14","DAY16-17","DAY19-20","DAY22-23")]
s2=vp.meta.seurat.CD4.2[,vp.meta.seurat.CD4.2$TO1_time %in% c("DAY10-11","DAY13-14","DAY16-17","DAY19-20","DAY22-23")]
x=data.frame(timepoint=s@meta.data$TO1_time,G2M.Score=s@meta.data$G2M.Score,Cluster=s@meta.data$seurat_clusters)
x$Cluster="CD8"
x2=data.frame(timepoint=s2@meta.data$TO1_time,G2M.Score=s2@meta.data$G2M.Score,Cluster=s2@meta.data$seurat_clusters)
x2$Cluster="CD4"
x=rbind(x,x2)
p= ggplot(x, aes(x=timepoint, y=G2M.Score,fill=Cluster)) + 
  geom_violin(trim=T,scale="width")+#coord_trans(y = "log10")+
  #scale_y_continuous(trans='log10') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("G2M Score by Timepoint Per Cluster")
p + theme(
  plot.title = element_text(size=14, face="bold"),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))

#slingshot
library(slingshot)
#' Assign a color to each cell based on some value
#' @param cell_vars Vector indicating the value of a variable associated with cells.
#' @param pal_fun Palette function that returns a vector of hex colors, whose
#' argument is the length of such a vector.
#' @param ... Extra arguments for pal_fun.
#' @return A vector of hex colors with one entry for each cell.
cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}
sds <- slingshot(Embeddings(vp.meta.seurat.CD4, "lda"), clusterLabels = vp.meta.seurat.CD4$seurat_clusters, stretch = 0,approx_points=100,start.clus=0)
cell_colors_clust <- cell_pal(vp.meta.seurat.CD4$seurat_clusters, hue_pal())
plot(reducedDim(sds), col = cell_colors_clust, pch = 16, cex = 0.5)
lines(sds, lwd = 2, col = 'black')
sds <- slingshot(Embeddings(vp.meta.seurat.CD8, "lda"), clusterLabels = vp.meta.seurat.CD8$seurat_clusters, start.clus = 0, stretch = 0,approx_points=100)
cell_colors_clust <- cell_pal(vp.meta.seurat.CD8$seurat_clusters, hue_pal())
plot(reducedDim(sds), col = cell_colors_clust, pch = 16, cex = 0.5)
lines(sds, lwd = 2, col = 'black')

#Tumor Lymph Node
cd4.tcr.in.tumor=unique(vp.meta.seurat.CD4$cdr3b[which(vp.meta.seurat.CD4$ao_IMP_TX!="LN")])
cd8.tcr.in.tumor=unique(vp.meta.seurat.CD8$cdr3b[which(vp.meta.seurat.CD8$ao_IMP_TX!="LN")])
vp.meta.seurat.CD4$TX1_condition=factor(vp.meta.seurat.CD4$TX1_condition, levels=c("TO","PRE","PBS","P","C ","PC","PBS:TLN","P:LN","C:LN","PC:LN"))
vp.meta.seurat.CD8$TX1_condition=factor(vp.meta.seurat.CD8$TX1_condition, levels=c("TO","PRE","PBS","P","C ","PC","PBS:TLN","P:LN","C:LN","PC:LN"))
DimPlot(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$ao_IMP_TX=="LN")],reduction="lda")
DimPlot(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$ao_IMP_TX=="LN")],reduction="lda")
x=vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$ao_IMP_TX=="LN")]
y=vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$ao_IMP_TX!="LN")]
length(which(x$cdr3b_nt %in% y$cdr3b_nt & x$condition=="PBS"))/length(which(x$condition=="PBS"))
length(which(x$cdr3b_nt %in% y$cdr3b_nt & x$condition=="PC"))/length(which(x$condition=="PC"))
length(which(x$cdr3b_nt %in% y$cdr3b_nt & x$condition=="P"))/length(which(x$condition=="P"))
length(which(x$cdr3b_nt %in% y$cdr3b_nt & x$condition=="C"))/length(which(x$condition=="C"))
x$condition=factor(x$condition,levels=c("PBS","C","P","PC"))
table(x$experiment)
DimPlot(x,split.by="condition",reduction = "lda")
y=as.data.frame.matrix(table(x$condition,x$seurat_clusters))
y=y[,which(colSums(y)>0)]
g=apply(y,1,function(x){x/sum(x)})
y <- melt(g)
colnames(y)=c("cluster","treatment","frequency")
y$treatment=factor(y$treatment, levels=c("PBS","C","P","PC"))
y$cluster=as.factor(y$cluster)
p=ggplot(y, aes(x=cluster, y=frequency,group=treatment,fill=cluster)) +
  geom_bar(stat="identity",position=position_dodge(),color="black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Cluster Frequency by Treatment")
p + theme(
  plot.title = element_text(size=14, face="bold"),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))
x=vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$ao_IMP_TX=="LN")]
y=vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$ao_IMP_TX!="LN")]
length(which(x$cdr3b_nt %in% y$cdr3b_nt & x$condition=="PBS"))/length(which(x$condition=="PBS"))
length(which(x$cdr3b_nt %in% y$cdr3b_nt & x$condition=="PC"))/length(which(x$condition=="PC"))
length(which(x$cdr3b_nt %in% y$cdr3b_nt & x$condition=="P"))/length(which(x$condition=="P"))
length(which(x$cdr3b_nt %in% y$cdr3b_nt & x$condition=="C"))/length(which(x$condition=="C"))
x$condition=factor(x$condition,levels=c("PBS","C","P","PC"))
table(x$experiment)
DimPlot(x,split.by="condition",reduction = "lda")
y=as.data.frame.matrix(table(x$condition,x$seurat_clusters))
y=y[,which(colSums(y)>0)]
g=apply(y,1,function(x){x/sum(x)})
y <- melt(g)
colnames(y)=c("cluster","treatment","frequency")
y$treatment=factor(y$treatment, levels=c("PBS","C","P","PC"))
y$cluster=as.factor(y$cluster)
p=ggplot(y, aes(x=cluster, y=frequency,group=treatment,fill=cluster)) +
  geom_bar(stat="identity",position=position_dodge(),color="black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Cluster Frequency by Treatment")
p + theme(
  plot.title = element_text(size=14, face="bold"),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))

###add TCR tracking by time point, and in treated vs untreated

####shared TCR tracking pre to post
x=as.data.frame.matrix(table(vp.meta.seurat.CD4$cdr3b_nt,vp.meta.seurat.CD4$timepoint))
shared=rownames(x)[((x[,1]>0 & x[,2]>0) | (x[,4]>0 & x[,5]>0))]
vp.meta.seurat.CD4$shared=0
vp.meta.seurat.CD4$shared[which(vp.meta.seurat.CD4$cdr3b_nt %in% shared)]=1
DimPlot(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$shared==T)],reduction="lda")
s=vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint!="TLN")]
s$timepoint[which(s$timepoint=="PRE")]="TUM1"
s$timepoint[which(s$timepoint=="POST")]="TUM2"
x=as.data.frame.matrix(table(s$cdr3b_nt,s$timepoint))
shared=rownames(x)[x[,1]>0 & x[,2]>0]
pre_mat=matrix(rep(0,length(shared)*length(unique(s$seurat_clusters))),nrow = length(shared),ncol=length(unique(s$seurat_clusters)))
post_mat=matrix(rep(0,length(shared)*length(unique(s$seurat_clusters))),nrow = length(shared),ncol=length(unique(s$seurat_clusters)))
origin_cluster=rep(NA,length(shared))
for(i in 1:length(shared)){
  tcr=shared[i]
  y=s[,which(s$cdr3b_nt==tcr)]
  x=as.data.frame.matrix(table(y$seurat_clusters,y$timepoint))
  pre_mat[i,]=x[,1]/sum(x[,1])
  post_mat[i,]=x[,2]/sum(x[,2])
  origin_cluster[i]=which(pre_mat[i,]==max(pre_mat[i,]))
}
vp.meta.seurat.CD4$timepoint_simplified=vp.meta.seurat.CD4$timepoint
vp.meta.seurat.CD4$timepoint_simplified[which(vp.meta.seurat.CD4$timepoint_simplified=="TUM1")]="PRE"
vp.meta.seurat.CD4$timepoint_simplified[which(vp.meta.seurat.CD4$timepoint_simplified=="TUM1")]="POST"
vp.meta.seurat.CD4$timepoint_simplified=factor(vp.meta.seurat.CD4$timepoint_simplified,levels=c("PRE","POST"))
DimPlot(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$cdr3b_nt %in% shared[which(origin_cluster==1)] & vp.meta.seurat.CD4$experiment=="eto")],split.by="timepoint_simplified",reduction="lda")
DimPlot(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$cdr3b_nt %in% shared[which(origin_cluster==2)] & vp.meta.seurat.CD4$experiment=="eto")],split.by="timepoint_simplified",reduction="lda")
DimPlot(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$cdr3b_nt %in% shared[which(origin_cluster==3)] & vp.meta.seurat.CD4$experiment=="eto")],split.by="timepoint_simplified",reduction="lda")
DimPlot(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$cdr3b_nt %in% shared[which(origin_cluster==4)] & vp.meta.seurat.CD4$experiment=="eto")],split.by="timepoint_simplified",reduction="lda")
difference=post_mat-pre_mat
cat(round(apply(post_mat[apply(pre_mat,1,function(x){x[1]==max(x)}),],2,mean)*100,digits = 2),sep="\t")
cat(round(apply(post_mat[apply(pre_mat,1,function(x){x[2]==max(x)}),],2,mean)*100,digits = 2),sep="\t")
cat(round(apply(post_mat[apply(pre_mat,1,function(x){x[3]==max(x)}),],2,mean)*100,digits = 2),sep="\t")
cat(round(apply(post_mat[apply(pre_mat,1,function(x){x[4]==max(x)}),],2,mean)*100,digits = 2),sep="\t")
table(s$shared,s$seurat_clusters)
x=as.data.frame.matrix(table(vp.meta.seurat.CD8$cdr3b_nt,vp.meta.seurat.CD8$timepoint))
shared=rownames(x)[((x[,1]>0 & x[,2]>0) | (x[,4]>0 & x[,5]>0))]
vp.meta.seurat.CD8$shared=0
vp.meta.seurat.CD8$shared[which(vp.meta.seurat.CD8$cdr3b_nt %in% shared)]=1
DimPlot(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$shared==T)])
s=vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint!="TLN")]
s$timepoint[which(s$timepoint=="PRE")]="TUM1"
s$timepoint[which(s$timepoint=="POST")]="TUM2"
x=as.data.frame.matrix(table(s$cdr3b_nt,s$timepoint))
shared=rownames(x)[x[,1]>0 & x[,2]>0]
pre_mat=matrix(rep(0,length(shared)*length(unique(s$seurat_clusters))),nrow = length(shared),ncol=length(unique(s$seurat_clusters)))
post_mat=matrix(rep(0,length(shared)*length(unique(s$seurat_clusters))),nrow = length(shared),ncol=length(unique(s$seurat_clusters)))
origin_cluster=rep(NA,length(shared))
for(i in 1:length(shared)){
  tcr=shared[i]
  y=s[,which(s$cdr3b_nt==tcr)]
  x=as.data.frame.matrix(table(y$seurat_clusters,y$timepoint))
  pre_mat[i,]=x[,1]/sum(x[,1])
  post_mat[i,]=x[,2]/sum(x[,2])
  origin_cluster[i]=which(pre_mat[i,]==max(pre_mat[i,]))
}
vp.meta.seurat.CD8$timepoint_simplified=vp.meta.seurat.CD8$timepoint
vp.meta.seurat.CD8$timepoint_simplified[which(vp.meta.seurat.CD8$timepoint_simplified=="TUM1")]="PRE"
vp.meta.seurat.CD8$timepoint_simplified[which(vp.meta.seurat.CD8$timepoint_simplified=="TUM1")]="POST"
vp.meta.seurat.CD8$timepoint_simplified=factor(vp.meta.seurat.CD8$timepoint_simplified,levels=c("PRE","POST"))
DimPlot(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$cdr3b_nt %in% shared[which(origin_cluster==1)] & vp.meta.seurat.CD8$experiment=="eto")],split.by="timepoint_simplified",reduction="lda")
DimPlot(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$cdr3b_nt %in% shared[which(origin_cluster==2)] & vp.meta.seurat.CD8$experiment=="eto")],split.by="timepoint_simplified",reduction="lda")
DimPlot(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$cdr3b_nt %in% shared[which(origin_cluster==3)] & vp.meta.seurat.CD8$experiment=="eto")],split.by="timepoint_simplified",reduction="lda")
DimPlot(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$cdr3b_nt %in% shared[which(origin_cluster==4)] & vp.meta.seurat.CD8$experiment=="eto")],split.by="timepoint_simplified",reduction="lda")
#DimPlot(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$cdr3b_nt %in% shared[which(origin_cluster==7)] & vp.meta.seurat.CD8$experiment=="eto")],split.by="timepoint_simplified",reduction="lda")
DimPlot(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$cdr3b_nt %in% shared[which(origin_cluster==1)] & vp.meta.seurat.CD8$experiment=="eto" & vp.meta.seurat.CD8$TO1_time %in% c("DAY10-11","DAY13-14"))],split.by="timepoint_simplified",reduction="lda")
DimPlot(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$cdr3b_nt %in% shared[which(origin_cluster==1)] & vp.meta.seurat.CD8$experiment=="eto" & vp.meta.seurat.CD8$TO1_time %in% c("DAY10-11","DAY16-17"))],split.by="timepoint_simplified",reduction="lda")
###
DimPlot(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$cdr3b_nt %in% shared[which(origin_cluster==2)] & vp.meta.seurat.CD8$experiment=="etx" & vp.meta.seurat.CD8$TX1_condition %in% c("PRE","PBS"))],split.by="timepoint_simplified",reduction="lda")
DimPlot(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$cdr3b_nt %in% shared[which(origin_cluster==2)] & vp.meta.seurat.CD8$experiment=="etx" & vp.meta.seurat.CD8$TX1_condition %in% c("PRE","P"))],split.by="timepoint_simplified",reduction="lda")
DimPlot(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$cdr3b_nt %in% shared[which(origin_cluster==2)] & vp.meta.seurat.CD8$experiment=="etx" & vp.meta.seurat.CD8$TX1_condition %in% c("PRE","C "))],split.by="timepoint_simplified",reduction="lda")
DimPlot(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$cdr3b_nt %in% shared[which(origin_cluster==2)] & vp.meta.seurat.CD8$experiment=="etx" & vp.meta.seurat.CD8$TX1_condition %in% c("PRE","PC"))],split.by="timepoint_simplified",reduction="lda")
###
difference=post_mat-pre_mat
cat(round(apply(post_mat[apply(pre_mat,1,function(x){x[1]==max(x)}),],2,mean)*100,digits = 2),sep="\t")
cat(round(apply(post_mat[apply(pre_mat,1,function(x){x[2]==max(x)}),],2,mean)*100,digits = 2),sep="\t")
cat(round(apply(post_mat[apply(pre_mat,1,function(x){x[3]==max(x)}),],2,mean)*100,digits = 2),sep="\t")
cat(round(apply(post_mat[apply(pre_mat,1,function(x){x[4]==max(x)}),],2,mean)*100,digits = 2),sep="\t")
cat(round(apply(post_mat[apply(pre_mat,1,function(x){x[5]==max(x)}),],2,mean)*100,digits = 2),sep="\t")
cat(round(apply(post_mat[apply(pre_mat,1,function(x){x[6]==max(x)}),],2,mean)*100,digits = 2),sep="\t")
#cat(round(apply(post_mat[apply(pre_mat,1,function(x){x[7]==max(x)}),],2,mean)*100,digits = 2),sep="\t")
table(s$shared,s$seurat_clusters)

tcr_tracking_cd4=function(s,treatment=F){
  x=as.data.frame.matrix(table(s$cdr3b,s$timepoint))
  shared=rownames(x)[x[,1]>0 & x[,2]>0]
  pre_mat=matrix(rep(0,length(shared)*length(unique(s$seurat_clusters))),nrow = length(shared),ncol=length(unique(s$seurat_clusters)))
  post_mat=matrix(rep(0,length(shared)*length(unique(s$seurat_clusters))),nrow = length(shared),ncol=length(unique(s$seurat_clusters)))
  for(i in 1:length(shared)){
    tcr=shared[i]
    y=s[,which(s$cdr3b==tcr)]
    x=as.data.frame.matrix(table(y$seurat_clusters,y$timepoint))
    pre_mat[i,]=x[,1]/sum(x[,1])
    post_mat[i,]=x[,2]/sum(x[,2])
    if(treatment){
      pre_mat[i,]=x[,2]/sum(x[,2])
      post_mat[i,]=x[,1]/sum(x[,1])
    }
  }
  difference=post_mat-pre_mat
  out=rbind(round(apply(post_mat[apply(pre_mat,1,function(x){x[1]==max(x)}),,drop=F],2,mean)*100,digits = 2),
            round(apply(post_mat[apply(pre_mat,1,function(x){x[2]==max(x)}),,drop=F],2,mean)*100,digits = 2),
            round(apply(post_mat[apply(pre_mat,1,function(x){x[3]==max(x)}),,drop=F],2,mean)*100,digits = 2),
            round(apply(post_mat[apply(pre_mat,1,function(x){x[4]==max(x)}),,drop=F],2,mean)*100,digits = 2),
            round(apply(post_mat[apply(pre_mat,1,function(x){x[5]==max(x)}),,drop=F],2,mean)*100,digits = 2),
            round(apply(post_mat[apply(pre_mat,1,function(x){x[6]==max(x)}),,drop=F],2,mean)*100,digits = 2),
            round(apply(post_mat[apply(pre_mat,1,function(x){x[7]==max(x)}),,drop=F],2,mean)*100,digits = 2))
  print(c(nrow(post_mat[apply(pre_mat,1,function(x){x[1]==max(x)}),,drop=F]),
          nrow(post_mat[apply(pre_mat,1,function(x){x[2]==max(x)}),,drop=F]),
          nrow(post_mat[apply(pre_mat,1,function(x){x[3]==max(x)}),,drop=F]),
          nrow(post_mat[apply(pre_mat,1,function(x){x[4]==max(x)}),,drop=F]),
          nrow(post_mat[apply(pre_mat,1,function(x){x[5]==max(x)}),,drop=F]),
          nrow(post_mat[apply(pre_mat,1,function(x){x[6]==max(x)}),,drop=F]),
          nrow(post_mat[apply(pre_mat,1,function(x){x[7]==max(x)}),,drop=F])))
  return(out)
}
tcr_tracking_cd8=function(s,treatment=F){
  x=as.data.frame.matrix(table(s$cdr3b,s$timepoint))
  shared=rownames(x)[x[,1]>0 & x[,2]>0]
  pre_mat=matrix(rep(0,length(shared)*length(unique(s$seurat_clusters))),nrow = length(shared),ncol=length(unique(s$seurat_clusters)))
  post_mat=matrix(rep(0,length(shared)*length(unique(s$seurat_clusters))),nrow = length(shared),ncol=length(unique(s$seurat_clusters)))
  for(i in 1:length(shared)){
    tcr=shared[i]
    y=s[,which(s$cdr3b==tcr)]
    x=as.data.frame.matrix(table(y$seurat_clusters,y$timepoint))
    pre_mat[i,]=x[,1]/sum(x[,1])
    post_mat[i,]=x[,2]/sum(x[,2])
    if(treatment){
      pre_mat[i,]=x[,2]/sum(x[,2])
      post_mat[i,]=x[,1]/sum(x[,1])
    }
  }
  difference=post_mat-pre_mat
  out=rbind(round(apply(post_mat[apply(pre_mat,1,function(x){x[1]==max(x)}),,drop=F],2,mean)*100,digits = 2),
            round(apply(post_mat[apply(pre_mat,1,function(x){x[2]==max(x)}),,drop=F],2,mean)*100,digits = 2),
            round(apply(post_mat[apply(pre_mat,1,function(x){x[3]==max(x)}),,drop=F],2,mean)*100,digits = 2),
            round(apply(post_mat[apply(pre_mat,1,function(x){x[4]==max(x)}),,drop=F],2,mean)*100,digits = 2))
  print(c(nrow(post_mat[apply(pre_mat,1,function(x){x[1]==max(x)}),,drop=F]),
          nrow(post_mat[apply(pre_mat,1,function(x){x[2]==max(x)}),,drop=F]),
          nrow(post_mat[apply(pre_mat,1,function(x){x[3]==max(x)}),,drop=F]),
          nrow(post_mat[apply(pre_mat,1,function(x){x[4]==max(x)}),,drop=F])))
  return(out)
}
vp.meta.seurat.CD4$condition[which(vp.meta.seurat.CD4$timepoint=="PRE")]="PreTx"
vp.meta.seurat.CD8$condition[which(vp.meta.seurat.CD8$timepoint=="PRE")]="PreTx"

write.table(tcr_tracking_cd4(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$TO2_type == "ETO")]),sep="\t",row.names = F,quote=F)
write.table(tcr_tracking_cd4(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$TO2_type == "LTO")]),sep="\t",row.names = F,quote=F)
write.table(tcr_tracking_cd4(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$TO2_type == "VLTO")]),sep="\t",row.names = F,quote=F)
write.table(tcr_tracking_cd4(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4$condition %in% c("PreTx","P") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T),sep="\t",row.names = F,quote=F)
write.table(tcr_tracking_cd4(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4$condition %in% c("PreTx","C") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T),sep="\t",row.names = F,quote=F)
write.table(tcr_tracking_cd4(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4$condition %in% c("PreTx","PC") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T),sep="\t",row.names = F,quote=F)
write.table(tcr_tracking_cd4(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4$condition %in% c("PreTx","PBS") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T),sep="\t",row.names = F,quote=F)

write.table(tcr_tracking_cd8(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TO2_type =="ETO")]),sep="\t",row.names = F,quote=F)
write.table(tcr_tracking_cd8(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TO2_type =="LTO")]),sep="\t",row.names = F,quote=F)
write.table(tcr_tracking_cd8(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TO2_type =="VLTO")]),sep="\t",row.names = F,quote=F)
write.table(tcr_tracking_cd8(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$condition %in% c("PreTx","P") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T),sep="\t",row.names = F,quote=F)
write.table(tcr_tracking_cd8(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$condition %in% c("PreTx","C") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T),sep="\t",row.names = F,quote=F)
write.table(tcr_tracking_cd8(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$condition %in% c("PreTx","PC") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T),sep="\t",row.names = F,quote=F)
write.table(tcr_tracking_cd8(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$condition %in% c("PreTx","PBS") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T),sep="\t",row.names = F,quote=F)


###alternative TCR tracking-- for each clone track from post-treatment timepoint where it originated from, weighted by post-treatment clone frequency

tcr_backtracking_cd4=function(s,treatment=F){
  x=as.data.frame.matrix(table(s$cdr3b,s$timepoint))
  shared=rownames(x)[x[,1]>0 & x[,2]>0]
  pre_mat=matrix(rep(0,length(shared)*length(unique(s$seurat_clusters))),nrow = length(shared),ncol=length(unique(s$seurat_clusters)))
  post_mat=matrix(rep(0,length(shared)*length(unique(s$seurat_clusters))),nrow = length(shared),ncol=length(unique(s$seurat_clusters)))
  for(i in 1:length(shared)){
    tcr=shared[i]
    y=s[,which(s$cdr3b==tcr)]
    x=as.data.frame.matrix(table(y$seurat_clusters,y$timepoint))
    pre_mat[i,]=x[,1]
    post_mat[i,]=x[,2]
    if(treatment){
      pre_mat[i,]=x[,2]
      post_mat[i,]=x[,1]
    }
  }
  index1=apply(post_mat,1,function(x){x[1]==max(x)})
  weight=post_mat[index1,1]
  pre_mat_subset=pre_mat[index1,,drop=F]
  pre_mat_subset_weighted=pre_mat_subset*weight
  pre_mat_subset_weighted=t(apply(pre_mat_subset_weighted,1,function(x){x/sum(x)}))
  clust1=round(apply(pre_mat_subset_weighted,2,mean)*100,digits=2)
  index2=apply(post_mat,1,function(x){x[2]==max(x)})
  weight=post_mat[index2,2]
  pre_mat_subset=pre_mat[index2,,drop=F]
  pre_mat_subset_weighted=pre_mat_subset*weight
  pre_mat_subset_weighted=t(apply(pre_mat_subset_weighted,1,function(x){x/sum(x)}))
  clust2=round(apply(pre_mat_subset_weighted,2,mean)*100,digits=2)
  index3=apply(post_mat,1,function(x){x[3]==max(x)})
  weight=post_mat[index3,3]
  pre_mat_subset=pre_mat[index3,,drop=F]
  pre_mat_subset_weighted=pre_mat_subset*weight
  pre_mat_subset_weighted=t(apply(pre_mat_subset_weighted,1,function(x){x/sum(x)}))
  clust3=round(apply(pre_mat_subset_weighted,2,mean)*100,digits=2)
  index4=apply(post_mat,1,function(x){x[4]==max(x)})
  weight=post_mat[index4,4]
  pre_mat_subset=pre_mat[index4,,drop=F]
  pre_mat_subset_weighted=pre_mat_subset*weight
  pre_mat_subset_weighted=t(apply(pre_mat_subset_weighted,1,function(x){x/sum(x)}))
  clust4=round(apply(pre_mat_subset_weighted,2,mean)*100,digits=2)
  index5=apply(post_mat,1,function(x){x[5]==max(x)})
  weight=post_mat[index5,5]
  pre_mat_subset=pre_mat[index5,,drop=F]
  pre_mat_subset_weighted=pre_mat_subset*weight
  pre_mat_subset_weighted=t(apply(pre_mat_subset_weighted,1,function(x){x/sum(x)}))
  clust5=round(apply(pre_mat_subset_weighted,2,mean)*100,digits=2)
  index6=apply(post_mat,1,function(x){x[6]==max(x)})
  weight=post_mat[index6,6]
  pre_mat_subset=pre_mat[index6,,drop=F]
  pre_mat_subset_weighted=pre_mat_subset*weight
  pre_mat_subset_weighted=t(apply(pre_mat_subset_weighted,1,function(x){x/sum(x)}))
  clust6=round(apply(pre_mat_subset_weighted,2,mean)*100,digits=2)
  index7=apply(post_mat,1,function(x){x[7]==max(x)})
  weight=post_mat[index7,7]
  pre_mat_subset=pre_mat[index7,,drop=F]
  pre_mat_subset_weighted=pre_mat_subset*weight
  pre_mat_subset_weighted=t(apply(pre_mat_subset_weighted,1,function(x){x/sum(x)}))
  clust7=round(apply(pre_mat_subset_weighted,2,mean)*100,digits=2)
  out=rbind(clust1,clust2,clust3,clust4,clust5,clust6,clust7)
  print(c(length(which(index1)),
          length(which(index2)),
          length(which(index3)),
          length(which(index4)),
          length(which(index5)),
          length(which(index6)),
          length(which(index7))))
  return(out)
}
tcr_backtracking_cd8=function(s,treatment=F){
  x=as.data.frame.matrix(table(s$cdr3b,s$timepoint))
  shared=rownames(x)[x[,1]>0 & x[,2]>0]
  pre_mat=matrix(rep(0,length(shared)*length(unique(s$seurat_clusters))),nrow = length(shared),ncol=length(unique(s$seurat_clusters)))
  post_mat=matrix(rep(0,length(shared)*length(unique(s$seurat_clusters))),nrow = length(shared),ncol=length(unique(s$seurat_clusters)))
  for(i in 1:length(shared)){
    tcr=shared[i]
    y=s[,which(s$cdr3b==tcr)]
    x=as.data.frame.matrix(table(y$seurat_clusters,y$timepoint))
    pre_mat[i,]=x[1:length(unique(s$seurat_clusters)),1]
    post_mat[i,]=x[1:length(unique(s$seurat_clusters)),2]
    if(treatment){
      pre_mat[i,]=x[1:length(unique(s$seurat_clusters)),2]
      post_mat[i,]=x[1:length(unique(s$seurat_clusters)),1]
    }
  }
  index1=apply(post_mat,1,function(x){x[1]==max(x)})
  weight=post_mat[index1,1]
  pre_mat_subset=pre_mat[index1,,drop=F]
  pre_mat_subset_weighted=pre_mat_subset*weight
  pre_mat_subset_weighted=t(apply(pre_mat_subset_weighted,1,function(x){x/sum(x)}))
  clust1=round(apply(pre_mat_subset_weighted,2,mean)*100,digits=2)
  index2=apply(post_mat,1,function(x){x[2]==max(x)})
  weight=post_mat[index2,2]
  pre_mat_subset=pre_mat[index2,,drop=F]
  pre_mat_subset_weighted=pre_mat_subset*weight
  pre_mat_subset_weighted=t(apply(pre_mat_subset_weighted,1,function(x){x/sum(x)}))
  clust2=round(apply(pre_mat_subset_weighted,2,mean)*100,digits=2)
  index3=apply(post_mat,1,function(x){x[3]==max(x)})
  weight=post_mat[index3,3]
  pre_mat_subset=pre_mat[index3,,drop=F]
  pre_mat_subset_weighted=pre_mat_subset*weight
  pre_mat_subset_weighted=t(apply(pre_mat_subset_weighted,1,function(x){x/sum(x)}))
  clust3=round(apply(pre_mat_subset_weighted,2,mean)*100,digits=2)
  index4=apply(post_mat,1,function(x){x[4]==max(x)})
  weight=post_mat[index4,4]
  pre_mat_subset=pre_mat[index4,,drop=F]
  pre_mat_subset_weighted=pre_mat_subset*weight
  pre_mat_subset_weighted=t(apply(pre_mat_subset_weighted,1,function(x){x/sum(x)}))
  clust4=round(apply(pre_mat_subset_weighted,2,mean)*100,digits=2)
  out=rbind(clust1,clust2,clust3,clust4)
  print(c(length(which(index1)),
          length(which(index2)),
          length(which(index3)),
          length(which(index4))))
  return(out)
}

write.table(tcr_backtracking_cd4(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$TO2_type == "ETO")]),sep="\t",row.names = F,quote=F)
write.table(tcr_backtracking_cd4(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$TO2_type == "LTO")]),sep="\t",row.names = F,quote=F)
write.table(tcr_backtracking_cd4(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$TO2_type == "VLTO")]),sep="\t",row.names = F,quote=F)
write.table(tcr_backtracking_cd4(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4$condition %in% c("PreTx","P") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T),sep="\t",row.names = F,quote=F)
write.table(tcr_backtracking_cd4(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4$condition %in% c("PreTx","C") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T),sep="\t",row.names = F,quote=F)
write.table(tcr_backtracking_cd4(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4$condition %in% c("PreTx","PC") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T),sep="\t",row.names = F,quote=F)
write.table(tcr_backtracking_cd4(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4$condition %in% c("PreTx","PBS") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T),sep="\t",row.names = F,quote=F)

write.table(tcr_backtracking_cd8(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TO2_type =="ETO")]),sep="\t",row.names = F,quote=F)
write.table(tcr_backtracking_cd8(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TO2_type =="LTO")]),sep="\t",row.names = F,quote=F)
write.table(tcr_backtracking_cd8(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TO2_type =="VLTO")]),sep="\t",row.names = F,quote=F)
write.table(tcr_backtracking_cd8(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$condition %in% c("PreTx","P") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T),sep="\t",row.names = F,quote=F)
write.table(tcr_backtracking_cd8(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$condition %in% c("PreTx","C") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T),sep="\t",row.names = F,quote=F)
write.table(tcr_backtracking_cd8(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$condition %in% c("PreTx","PC") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T),sep="\t",row.names = F,quote=F)
write.table(tcr_backtracking_cd8(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$condition %in% c("PreTx","PBS") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T),sep="\t",row.names = F,quote=F)


vp.meta.seurat.CD4$seurat_clusters=droplevels(vp.meta.seurat.CD4$seurat_clusters)

###tcr backtracking with sub-sampling

tcr_backtracking_cd4_bootstrapped=function(s,treatment=F){
  set.seed(1234)
  s=s[,which(s$timepoint!="TLN")]
  n=min(table(s$timepoint))
  s_main=s
  clust1=matrix(rep(0,100*length(unique(s$seurat_clusters))),nrow = 100,ncol=length(unique(s$seurat_clusters)))
  clust2=matrix(rep(0,100*length(unique(s$seurat_clusters))),nrow = 100,ncol=length(unique(s$seurat_clusters)))
  clust3=matrix(rep(0,100*length(unique(s$seurat_clusters))),nrow = 100,ncol=length(unique(s$seurat_clusters)))
  clust4=matrix(rep(0,100*length(unique(s$seurat_clusters))),nrow = 100,ncol=length(unique(s$seurat_clusters)))
  clust5=matrix(rep(0,100*length(unique(s$seurat_clusters))),nrow = 100,ncol=length(unique(s$seurat_clusters)))
  clust6=matrix(rep(0,100*length(unique(s$seurat_clusters))),nrow = 100,ncol=length(unique(s$seurat_clusters)))
  clust7=matrix(rep(0,100*length(unique(s$seurat_clusters))),nrow = 100,ncol=length(unique(s$seurat_clusters)))
  for(j in 1:100){
    print(j)
    s=s_main[,c(sample(which(s_main$timepoint==unique(s_main$timepoint)[1]),n,replace = F),sample(which(s_main$timepoint==unique(s_main$timepoint)[2]),n,replace = F))]
    #s$seurat_clusters=droplevels(s$seurat_clusters)
    x=as.data.frame.matrix(table(s$cdr3b,s$timepoint))
    shared=rownames(x)[x[,1]>0 & x[,2]>0]
    pre_mat=matrix(rep(0,length(shared)*length(unique(s$seurat_clusters))),nrow = length(shared),ncol=length(unique(s$seurat_clusters)))
    post_mat=matrix(rep(0,length(shared)*length(unique(s$seurat_clusters))),nrow = length(shared),ncol=length(unique(s$seurat_clusters)))
    for(i in 1:length(shared)){
      tcr=shared[i]
      y=s[,which(s$cdr3b==tcr)]
      x=as.data.frame.matrix(table(y$seurat_clusters,y$timepoint))
      pre_mat[i,]=x[,1]
      post_mat[i,]=x[,2]
      if(treatment){
        pre_mat[i,]=x[,2]
        post_mat[i,]=x[,1]
      }
    }
    index1=apply(post_mat,1,function(x){x[1]==max(x)})
    if(length(which(index1))>0){
    weight=post_mat[index1,1]
    pre_mat_subset=pre_mat[index1,,drop=F]
    pre_mat_subset_weighted=pre_mat_subset*weight
    pre_mat_subset_weighted=t(apply(pre_mat_subset_weighted,1,function(x){x/sum(x)}))
    clust1[j,]=round(apply(pre_mat_subset_weighted,2,mean)*100,digits=2)}
    index2=apply(post_mat,1,function(x){x[2]==max(x)})
    if(length(which(index2))>0){
    weight=post_mat[index2,2]
    pre_mat_subset=pre_mat[index2,,drop=F]
    pre_mat_subset_weighted=pre_mat_subset*weight
    pre_mat_subset_weighted=t(apply(pre_mat_subset_weighted,1,function(x){x/sum(x)}))
    clust2[j,]=round(apply(pre_mat_subset_weighted,2,mean)*100,digits=2)}
    index3=apply(post_mat,1,function(x){x[3]==max(x)})
    if(length(which(index3))>0){
    weight=post_mat[index3,3]
    pre_mat_subset=pre_mat[index3,,drop=F]
    pre_mat_subset_weighted=pre_mat_subset*weight
    pre_mat_subset_weighted=t(apply(pre_mat_subset_weighted,1,function(x){x/sum(x)}))
    clust3[j,]=round(apply(pre_mat_subset_weighted,2,mean)*100,digits=2)}
    index4=apply(post_mat,1,function(x){x[4]==max(x)})
    if(length(which(index4))>0){
    weight=post_mat[index4,4]
    pre_mat_subset=pre_mat[index4,,drop=F]
    pre_mat_subset_weighted=pre_mat_subset*weight
    pre_mat_subset_weighted=t(apply(pre_mat_subset_weighted,1,function(x){x/sum(x)}))
    clust4[j,]=round(apply(pre_mat_subset_weighted,2,mean)*100,digits=2)}
    index5=apply(post_mat,1,function(x){x[5]==max(x)})
    if(length(which(index5))>0){
    weight=post_mat[index5,5]
    pre_mat_subset=pre_mat[index5,,drop=F]
    pre_mat_subset_weighted=pre_mat_subset*weight
    pre_mat_subset_weighted=t(apply(pre_mat_subset_weighted,1,function(x){x/sum(x)}))
    clust5[j,]=round(apply(pre_mat_subset_weighted,2,mean)*100,digits=2)}
    index6=apply(post_mat,1,function(x){x[6]==max(x)})
    if(length(which(index6))>0){
    weight=post_mat[index6,6]
    pre_mat_subset=pre_mat[index6,,drop=F]
    pre_mat_subset_weighted=pre_mat_subset*weight
    pre_mat_subset_weighted=t(apply(pre_mat_subset_weighted,1,function(x){x/sum(x)}))
    clust6[j,]=round(apply(pre_mat_subset_weighted,2,mean)*100,digits=2)}
    index7=apply(post_mat,1,function(x){x[7]==max(x)})
    if(length(which(index7))>0){
    weight=post_mat[index7,7]
    pre_mat_subset=pre_mat[index7,,drop=F]
    pre_mat_subset_weighted=pre_mat_subset*weight
    pre_mat_subset_weighted=t(apply(pre_mat_subset_weighted,1,function(x){x/sum(x)}))
    clust7[j,]=round(apply(pre_mat_subset_weighted,2,mean)*100,digits=2)}
  }
  
  colnames(clust1)=levels(vp.meta.seurat.CD4.2) #c("cluster0","cluster1","cluster2","cluster3","cluster4","cluster5")
  y <- melt(clust1[which(rowSums(clust1)>0),])
  y$Cluster=levels(vp.meta.seurat.CD4.2)[1] #"Cluster0"
  colnames(clust2)=levels(vp.meta.seurat.CD4.2) #c("cluster0","cluster1","cluster2","cluster3","cluster4","cluster5")
  y2 <- melt(clust2[which(rowSums(clust2)>0),])
  y2$Cluster=levels(vp.meta.seurat.CD4.2)[2] #"Cluster1"
  y=rbind(y,y2)
  colnames(clust3)=levels(vp.meta.seurat.CD4.2)#c("cluster0","cluster1","cluster2","cluster3","cluster4","cluster5")
  y2 <- melt(clust3[which(rowSums(clust3)>0),])
  y2$Cluster=levels(vp.meta.seurat.CD4.2)[3] #"Cluster2"
  y=rbind(y,y2)
  colnames(clust4)=levels(vp.meta.seurat.CD4.2) #c("cluster0","cluster1","cluster2","cluster3","cluster4","cluster5")
  y2 <- melt(clust4[which(rowSums(clust4)>0),])
  y2$Cluster=levels(vp.meta.seurat.CD4.2)[4] #"Cluster3"
  y=rbind(y,y2)
  colnames(y)=c("Replicate","Origin","Frequency","Cluster")
  y$Origin <- factor( y$Origin, levels = c("CD4_1.TregA","CD4_2.TregB", "CD4_3.Th1 like","CD4_5.Tcm" ,"CD4_4.Th1 ncm","CD4_6.Trans1"))
  p=ggplot(y, aes(x=Origin, y=Frequency,fill=Cluster)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    ggtitle("Shared TCR Backtracing") +
    facet_wrap(~Cluster, ncol = 1)+ ylim(0,100)
    #facet_grid(Cluster~"Cluster Endpoint")
  print(p + theme(
    plot.title = element_text(size=14, face="bold",hjust = 0.5),
    axis.title = element_text(size=14, face="bold"),
    axis.text = element_text(size=8),
    legend.text = element_text(size=10),
    legend.title = element_text(size=12),
    strip.text.x = element_text(size = 12, face="bold"))+scale_fill_manual(values=c("#117733","#CC6677","#6699CC","#DDCC77")))
  
  out=list(clust1,clust2,clust3,clust4,clust5,clust6,clust7)
  return(out)
}

vp.meta.seurat.CD4.2=vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$seurat_clusters %in% c(0,1,2,3,4))]
vp.meta.seurat.CD4.2$seurat_clusters=droplevels(vp.meta.seurat.CD4.2$seurat_clusters)
vp.meta.seurat.CD4.2=vp.meta.seurat.CD4
cd4.new.cluster.ids <- c("CD4_1.TregA", "CD4_3.Th1 like", "CD4_2.TregB", "CD4_5.Tcm", "CD4_4.Th1 ncm", "CD4_6.Trans1", "CD4_6.Trans1")
names(cd4.new.cluster.ids) <- levels(vp.meta.seurat.CD4.2)
vp.meta.seurat.CD4.2<- RenameIdents(vp.meta.seurat.CD4.2, cd4.new.cluster.ids)
vp.meta.seurat.CD4.2$seurat_clusters=Idents(vp.meta.seurat.CD4.2)
vp.meta.seurat.CD4.2$seurat_clusters <- factor(vp.meta.seurat.CD4.2$seurat_clusters, levels = c("CD4_1.TregA","CD4_2.TregB", "CD4_3.Th1 like","CD4_5.Tcm" ,"CD4_4.Th1 ncm","CD4_6.Trans1"))
Idents(vp.meta.seurat.CD4.2)<- vp.meta.seurat.CD4.2$seurat_clusters


cd4_tcr_eto=tcr_backtracking_cd4_bootstrapped(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$TO2_type == "ETO")])
cd4_tcr_eto_backward_d13.14=tcr_backtracking_cd4_bootstrapped(vp.meta.seurat.CD4.2[,which(vp.meta.seurat.CD4.2$TO2_type == "ETO" & vp.meta.seurat.CD4.2$TO1_time %in% c("DAY10-11","DAY13-14"))])
cd4_tcr_eto_backward_d16.17=tcr_backtracking_cd4_bootstrapped(vp.meta.seurat.CD4.2[,which(vp.meta.seurat.CD4.2$TO2_type == "ETO" & vp.meta.seurat.CD4.2$TO1_time %in% c("DAY10-11","DAY16-17"))])
s=vp.meta.seurat.CD4.2[,which(vp.meta.seurat.CD4.2$TO2_type == "ETO" & vp.meta.seurat.CD4.2$timepoint=="TUM2")]
s$timepoint[sample(ncol(s),ncol(s)/3)]="TUM1"
cd4_tcr_eto_backward_random=tcr_backtracking_cd4_bootstrapped(s)
cd4_tcr_pbs=tcr_backtracking_cd4_bootstrapped(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4$condition %in% c("PreTx","PBS") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T)
cd4_tcr_p=tcr_backtracking_cd4_bootstrapped(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4$condition %in% c("PreTx","P") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T)
cd4_tcr_c=tcr_backtracking_cd4_bootstrapped(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4$condition %in% c("PreTx","C") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T)
cd4_tcr_pc=tcr_backtracking_cd4_bootstrapped(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4$condition %in% c("PreTx","PC") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T)

tcr_backtracking_cd8_bootstrapped=function(s,treatment=F){
  set.seed(1234)
  s=s[,which(s$timepoint!="TLN")]
  n=min(table(s$timepoint))
  s_main=s
  clust1=matrix(rep(0,100*length(unique(s$seurat_clusters))),nrow = 100,ncol=length(unique(s$seurat_clusters)))
  clust2=matrix(rep(0,100*length(unique(s$seurat_clusters))),nrow = 100,ncol=length(unique(s$seurat_clusters)))
  clust3=matrix(rep(0,100*length(unique(s$seurat_clusters))),nrow = 100,ncol=length(unique(s$seurat_clusters)))
  clust4=matrix(rep(0,100*length(unique(s$seurat_clusters))),nrow = 100,ncol=length(unique(s$seurat_clusters)))
  for(j in 1:100){
    print(j)
    s=s_main[,c(sample(which(s_main$timepoint==unique(s_main$timepoint)[1]),n,replace = F),sample(which(s_main$timepoint==unique(s_main$timepoint)[2]),n,replace = F))]
    #s$seurat_clusters=droplevels(s$seurat_clusters)
    x=as.data.frame.matrix(table(s$cdr3b,s$timepoint))
    shared=rownames(x)[x[,1]>0 & x[,2]>0]
    pre_mat=matrix(rep(0,length(shared)*length(unique(s$seurat_clusters))),nrow = length(shared),ncol=length(unique(s$seurat_clusters)))
    post_mat=matrix(rep(0,length(shared)*length(unique(s$seurat_clusters))),nrow = length(shared),ncol=length(unique(s$seurat_clusters)))
    for(i in 1:length(shared)){
      tcr=shared[i]
      y=s[,which(s$cdr3b==tcr)]
      x=as.data.frame.matrix(table(y$seurat_clusters,y$timepoint))
      pre_mat[i,]=x[,1]
      post_mat[i,]=x[,2]
      if(treatment){
        pre_mat[i,]=x[,2]
        post_mat[i,]=x[,1]
      }
    }
    index1=apply(post_mat,1,function(x){x[1]==max(x)})
    if(length(which(index1))>0){
      weight=post_mat[index1,1]
      pre_mat_subset=pre_mat[index1,,drop=F]
      pre_mat_subset_weighted=pre_mat_subset*weight
      pre_mat_subset_weighted=t(apply(pre_mat_subset_weighted,1,function(x){x/sum(x)}))
      clust1[j,]=round(apply(pre_mat_subset_weighted,2,mean)*100,digits=2)}
    index2=apply(post_mat,1,function(x){x[2]==max(x)})
    if(length(which(index2))>0){
      weight=post_mat[index2,2]
      pre_mat_subset=pre_mat[index2,,drop=F]
      pre_mat_subset_weighted=pre_mat_subset*weight
      pre_mat_subset_weighted=t(apply(pre_mat_subset_weighted,1,function(x){x/sum(x)}))
      clust2[j,]=round(apply(pre_mat_subset_weighted,2,mean)*100,digits=2)}
    index3=apply(post_mat,1,function(x){x[3]==max(x)})
    if(length(which(index3))>0){
      weight=post_mat[index3,3]
      pre_mat_subset=pre_mat[index3,,drop=F]
      pre_mat_subset_weighted=pre_mat_subset*weight
      pre_mat_subset_weighted=t(apply(pre_mat_subset_weighted,1,function(x){x/sum(x)}))
      clust3[j,]=round(apply(pre_mat_subset_weighted,2,mean)*100,digits=2)}
    index4=apply(post_mat,1,function(x){x[4]==max(x)})
    if(length(which(index4))>0){
      weight=post_mat[index4,4]
      pre_mat_subset=pre_mat[index4,,drop=F]
      pre_mat_subset_weighted=pre_mat_subset*weight
      pre_mat_subset_weighted=t(apply(pre_mat_subset_weighted,1,function(x){x/sum(x)}))
      clust4[j,]=round(apply(pre_mat_subset_weighted,2,mean)*100,digits=2)}
  }
  
  colnames(clust1)=c("cluster0","cluster1","cluster2","cluster3")
  y <- melt(clust1[which(rowSums(clust1)>0),])
  y$Cluster="Cluster0"
  colnames(clust2)=c("cluster0","cluster1","cluster2","cluster3")
  y2 <- melt(clust2[which(rowSums(clust2)>0),])
  if(nrow(y2)>0){
  y2$Cluster="Cluster1"
  y=rbind(y,y2)}
  colnames(clust3)=c("cluster0","cluster1","cluster2","cluster3")
  y2 <- melt(clust3[which(rowSums(clust3)>0),])
  if(nrow(y2)>0){
  y2$Cluster="Cluster2"
  y=rbind(y,y2)}
  colnames(clust4)=c("cluster0","cluster1","cluster2","cluster3")
  y2 <- melt(clust4[which(rowSums(clust4)>0),])
  if(nrow(y2)>0){y2$Cluster="Cluster3"
  y=rbind(y,y2)}
  colnames(y)=c("Replicate","Origin","Frequency","Cluster")
  p=ggplot(y, aes(x=Origin, y=Frequency,fill=Cluster)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    ggtitle("Shared TCR Backtracing") +
    facet_wrap(~Cluster, ncol = 1) + ylim(0,100)
  #facet_grid(Cluster~"Cluster Endpoint")
  print(p + theme(
    plot.title = element_text(size=14, face="bold",hjust = 0.5),
    axis.title = element_text(size=14, face="bold"),
    axis.text = element_text(size=8),
    legend.text = element_text(size=10),
    legend.title = element_text(size=12),
    strip.text.x = element_text(size = 12, face="bold"))+scale_fill_manual(values=c("#117733","#CC6677","#6699CC","#DDCC77")))
  
  out=list(clust1,clust2,clust3,clust4)
  return(out)
}

cd8_tcr_eto=tcr_backtracking_cd8_bootstrapped(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TO2_type == "ETO")])
cd8_tcr_eto_backward_d13.14=tcr_backtracking_cd8_bootstrapped(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TO2_type == "ETO" & vp.meta.seurat.CD8$TO1_time %in% c("DAY10-11","DAY13-14"))])
cd8_tcr_eto_backward_d16.17=tcr_backtracking_cd8_bootstrapped(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TO2_type == "ETO" & vp.meta.seurat.CD8$TO1_time %in% c("DAY10-11","DAY16-17"))])
s=vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TO2_type == "ETO" & vp.meta.seurat.CD8$timepoint=="TUM2")]
s$timepoint[sample(ncol(s),ncol(s)/3)]="TUM1"
cd8_tcr_eto_backward_random=tcr_backtracking_cd8_bootstrapped(s)
cd8_tcr_pbs=tcr_backtracking_cd8_bootstrapped(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","PBS") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_p=tcr_backtracking_cd8_bootstrapped(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","P") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_c=tcr_backtracking_cd8_bootstrapped(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","C") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_pc=tcr_backtracking_cd8_bootstrapped(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","PC") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)


### TCR forward-tracking bootstrapped

tcr_forwardtracking_cd4_bootstrapped=function(s,treatment=F){
  set.seed(1234)
  s=s[,which(s$timepoint!="TLN")]
  n=min(table(s$timepoint))
  s_main=s
  clust1=matrix(rep(0,100*length(unique(s$seurat_clusters))),nrow = 100,ncol=length(unique(s$seurat_clusters)))
  clust2=matrix(rep(0,100*length(unique(s$seurat_clusters))),nrow = 100,ncol=length(unique(s$seurat_clusters)))
  clust3=matrix(rep(0,100*length(unique(s$seurat_clusters))),nrow = 100,ncol=length(unique(s$seurat_clusters)))
  clust4=matrix(rep(0,100*length(unique(s$seurat_clusters))),nrow = 100,ncol=length(unique(s$seurat_clusters)))
  clust5=matrix(rep(0,100*length(unique(s$seurat_clusters))),nrow = 100,ncol=length(unique(s$seurat_clusters)))
  clust6=matrix(rep(0,100*length(unique(s$seurat_clusters))),nrow = 100,ncol=length(unique(s$seurat_clusters)))
  clust7=matrix(rep(0,100*length(unique(s$seurat_clusters))),nrow = 100,ncol=length(unique(s$seurat_clusters)))
  for(j in 1:100){
    print(j)
    s=s_main[,c(sample(which(s_main$timepoint==unique(s_main$timepoint)[1]),n,replace = F),sample(which(s_main$timepoint==unique(s_main$timepoint)[2]),n,replace = F))]
    #s$seurat_clusters=droplevels(s$seurat_clusters)
    x=as.data.frame.matrix(table(s$cdr3b,s$timepoint))
    shared=rownames(x)[x[,1]>0 & x[,2]>0]
    pre_mat=matrix(rep(0,length(shared)*length(unique(s$seurat_clusters))),nrow = length(shared),ncol=length(unique(s$seurat_clusters)))
    post_mat=matrix(rep(0,length(shared)*length(unique(s$seurat_clusters))),nrow = length(shared),ncol=length(unique(s$seurat_clusters)))
    for(i in 1:length(shared)){
      tcr=shared[i]
      y=s[,which(s$cdr3b==tcr)]
      x=as.data.frame.matrix(table(y$seurat_clusters,y$timepoint))
      pre_mat[i,]=x[,1]
      post_mat[i,]=x[,2]
      if(treatment){
        pre_mat[i,]=x[,2]
        post_mat[i,]=x[,1]
      }
    }
    index1=apply(pre_mat,1,function(x){x[1]==max(x)})
    if(length(which(index1))>0){
      weight=pre_mat[index1,1]
      post_mat_subset=post_mat[index1,,drop=F]
      post_mat_subset_weighted=post_mat_subset*weight
      post_mat_subset_weighted=t(apply(post_mat_subset_weighted,1,function(x){x/sum(x)}))
      clust1[j,]=round(apply(post_mat_subset_weighted,2,mean)*100,digits=2)}
    index2=apply(pre_mat,1,function(x){x[2]==max(x)})
    if(length(which(index2))>0){
      weight=pre_mat[index2,2]
      post_mat_subset=post_mat[index2,,drop=F]
      post_mat_subset_weighted=post_mat_subset*weight
      post_mat_subset_weighted=t(apply(post_mat_subset_weighted,1,function(x){x/sum(x)}))
      clust2[j,]=round(apply(post_mat_subset_weighted,2,mean)*100,digits=2)}
    index3=apply(pre_mat,1,function(x){x[3]==max(x)})
    if(length(which(index3))>0){
      weight=pre_mat[index3,3]
      post_mat_subset=post_mat[index3,,drop=F]
      post_mat_subset_weighted=post_mat_subset*weight
      post_mat_subset_weighted=t(apply(post_mat_subset_weighted,1,function(x){x/sum(x)}))
      clust3[j,]=round(apply(post_mat_subset_weighted,2,mean)*100,digits=2)}
    index4=apply(pre_mat,1,function(x){x[4]==max(x)})
    if(length(which(index4))>0){
      weight=pre_mat[index4,4]
      post_mat_subset=post_mat[index4,,drop=F]
      post_mat_subset_weighted=post_mat_subset*weight
      post_mat_subset_weighted=t(apply(post_mat_subset_weighted,1,function(x){x/sum(x)}))
      clust4[j,]=round(apply(post_mat_subset_weighted,2,mean)*100,digits=2)}
    index5=apply(pre_mat,1,function(x){x[5]==max(x)})
    if(length(which(index5))>0){
      weight=pre_mat[index5,5]
      post_mat_subset=post_mat[index5,,drop=F]
      post_mat_subset_weighted=post_mat_subset*weight
      post_mat_subset_weighted=t(apply(post_mat_subset_weighted,1,function(x){x/sum(x)}))
      clust5[j,]=round(apply(post_mat_subset_weighted,2,mean)*100,digits=2)}
    index6=apply(pre_mat,1,function(x){x[6]==max(x)})
    if(length(which(index6))>0){
      weight=pre_mat[index6,6]
      post_mat_subset=post_mat[index6,,drop=F]
      post_mat_subset_weighted=post_mat_subset*weight
      post_mat_subset_weighted=t(apply(post_mat_subset_weighted,1,function(x){x/sum(x)}))
      clust6[j,]=round(apply(post_mat_subset_weighted,2,mean)*100,digits=2)}
    index7=apply(pre_mat,1,function(x){x[7]==max(x)})
    if(length(which(index7))>0){
      weight=pre_mat[index7,7]
      post_mat_subset=post_mat[index7,,drop=F]
      post_mat_subset_weighted=post_mat_subset*weight
      post_mat_subset_weighted=t(apply(post_mat_subset_weighted,1,function(x){x/sum(x)}))
      clust7[j,]=round(apply(post_mat_subset_weighted,2,mean)*100,digits=2)}
  }
  colnames(clust1)=levels(vp.meta.seurat.CD4.2) #c("cluster0","cluster1","cluster2","cluster3","cluster4","cluster5")
  y <- melt(clust1[which(rowSums(clust1)>0),])
  y$Cluster=levels(vp.meta.seurat.CD4.2)[1] #"Cluster0"
  colnames(clust2)=levels(vp.meta.seurat.CD4.2) #c("cluster0","cluster1","cluster2","cluster3","cluster4","cluster5")
  y2 <- melt(clust2[which(rowSums(clust2)>0),])
  y2$Cluster=levels(vp.meta.seurat.CD4.2)[2] #"Cluster1"
  y=rbind(y,y2)
  colnames(clust3)=levels(vp.meta.seurat.CD4.2) #c("cluster0","cluster1","cluster2","cluster3","cluster4","cluster5")
  y2 <- melt(clust3[which(rowSums(clust3)>0),])
  y2$Cluster=levels(vp.meta.seurat.CD4.2)[3] #"Cluster2"
  y=rbind(y,y2)
  colnames(clust4)=levels(vp.meta.seurat.CD4.2) #c("cluster0","cluster1","cluster2","cluster3","cluster4","cluster5")
  y2 <- melt(clust4[which(rowSums(clust4)>0),])
  y2$Cluster=levels(vp.meta.seurat.CD4.2)[4] #"Cluster3"
  y=rbind(y,y2)
  colnames(y)=c("Replicate","Origin","Frequency","Cluster")
  y$Origin <- factor( y$Origin, levels = c("CD4_1.TregA","CD4_2.TregB", "CD4_3.Th1 like","CD4_5.Tcm" ,"CD4_4.Th1 ncm","CD4_6.Trans1"))
  p=ggplot(y, aes(x=Origin, y=Frequency,fill=Cluster)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    ggtitle("Shared TCR Forward-Tracing") +
    facet_wrap(~Cluster, ncol = 1) + ylim(0,100)
  #facet_grid(Cluster~"Cluster Endpoint")
  print(p + theme(
    plot.title = element_text(size=14, face="bold",hjust = 0.5),
    axis.title = element_text(size=14, face="bold"),
    axis.text = element_text(size=8),
    legend.text = element_text(size=10),
    legend.title = element_text(size=12),
    strip.text.x = element_text(size = 12, face="bold"))+scale_fill_manual(values=c("#117733","#CC6677","#6699CC","#DDCC77")))
  out=list(clust1,clust2,clust3,clust4,clust5,clust6,clust7)
  return(out)
}

cd4_tcr_eto=tcr_forwardtracking_cd4_bootstrapped(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$TO2_type == "ETO")])
cd4_tcr_eto_forward_d13.14=tcr_forwardtracking_cd4_bootstrapped(vp.meta.seurat.CD4.2[,which(vp.meta.seurat.CD4.2$TO2_type == "ETO" & vp.meta.seurat.CD4.2$TO1_time %in% c("DAY10-11","DAY13-14"))])
cd4_tcr_eto_forward_d16.17=tcr_forwardtracking_cd4_bootstrapped(vp.meta.seurat.CD4.2[,which(vp.meta.seurat.CD4.2$TO2_type == "ETO" & vp.meta.seurat.CD4.2$TO1_time %in% c("DAY10-11","DAY16-17"))])
s=vp.meta.seurat.CD4.2[,which(vp.meta.seurat.CD4.2$TO2_type == "ETO" & vp.meta.seurat.CD4.2$timepoint=="TUM1")]
s$timepoint[sample(ncol(s),ncol(s)/3)]="TUM2"
cd4_tcr_eto_forward_random=tcr_forwardtracking_cd4_bootstrapped(s)
cd4_tcr_pbs=tcr_forwardtracking_cd4_bootstrapped(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4$condition %in% c("PreTx","PBS") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T)
cd4_tcr_p=tcr_forwardtracking_cd4_bootstrapped(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4$condition %in% c("PreTx","P") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T)
cd4_tcr_c=tcr_forwardtracking_cd4_bootstrapped(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4$condition %in% c("PreTx","C") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T)
cd4_tcr_pc=tcr_forwardtracking_cd4_bootstrapped(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4$condition %in% c("PreTx","PC") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T)

tcr_forwardtracking_cd8_bootstrapped=function(s,treatment=F){
  set.seed(1234)
  s=s[,which(s$timepoint!="TLN")]
  n=min(table(s$timepoint))
  s_main=s
  clust1=matrix(rep(0,100*length(unique(s$seurat_clusters))),nrow = 100,ncol=length(unique(s$seurat_clusters)))
  clust2=matrix(rep(0,100*length(unique(s$seurat_clusters))),nrow = 100,ncol=length(unique(s$seurat_clusters)))
  clust3=matrix(rep(0,100*length(unique(s$seurat_clusters))),nrow = 100,ncol=length(unique(s$seurat_clusters)))
  clust4=matrix(rep(0,100*length(unique(s$seurat_clusters))),nrow = 100,ncol=length(unique(s$seurat_clusters)))
  for(j in 1:100){
    print(j)
    s=s_main[,c(sample(which(s_main$timepoint==unique(s_main$timepoint)[1]),n,replace = F),sample(which(s_main$timepoint==unique(s_main$timepoint)[2]),n,replace = F))]
    #s$seurat_clusters=droplevels(s$seurat_clusters)
    x=as.data.frame.matrix(table(s$cdr3b,s$timepoint))
    shared=rownames(x)[x[,1]>0 & x[,2]>0]
    pre_mat=matrix(rep(0,length(shared)*length(unique(s$seurat_clusters))),nrow = length(shared),ncol=length(unique(s$seurat_clusters)))
    post_mat=matrix(rep(0,length(shared)*length(unique(s$seurat_clusters))),nrow = length(shared),ncol=length(unique(s$seurat_clusters)))
    for(i in 1:length(shared)){
      tcr=shared[i]
      y=s[,which(s$cdr3b==tcr)]
      x=as.data.frame.matrix(table(y$seurat_clusters,y$timepoint))
      pre_mat[i,]=x[,1]
      post_mat[i,]=x[,2]
      if(treatment){
        pre_mat[i,]=x[,2]
        post_mat[i,]=x[,1]
      }
    }
    index1=apply(pre_mat,1,function(x){x[1]==max(x)})
    if(length(which(index1))>0){
      weight=pre_mat[index1,1]
      post_mat_subset=post_mat[index1,,drop=F]
      post_mat_subset_weighted=post_mat_subset*weight
      post_mat_subset_weighted=t(apply(post_mat_subset_weighted,1,function(x){x/sum(x)}))
      clust1[j,]=round(apply(post_mat_subset_weighted,2,mean)*100,digits=2)}
    index2=apply(pre_mat,1,function(x){x[2]==max(x)})
    if(length(which(index2))>0){
      weight=pre_mat[index2,2]
      post_mat_subset=post_mat[index2,,drop=F]
      post_mat_subset_weighted=post_mat_subset*weight
      post_mat_subset_weighted=t(apply(post_mat_subset_weighted,1,function(x){x/sum(x)}))
      clust2[j,]=round(apply(post_mat_subset_weighted,2,mean)*100,digits=2)}
    index3=apply(pre_mat,1,function(x){x[3]==max(x)})
    if(length(which(index3))>0){
      weight=pre_mat[index3,3]
      post_mat_subset=post_mat[index3,,drop=F]
      post_mat_subset_weighted=post_mat_subset*weight
      post_mat_subset_weighted=t(apply(post_mat_subset_weighted,1,function(x){x/sum(x)}))
      clust3[j,]=round(apply(post_mat_subset_weighted,2,mean)*100,digits=2)}
    index4=apply(pre_mat,1,function(x){x[4]==max(x)})
    if(length(which(index4))>0){
      weight=pre_mat[index4,4]
      post_mat_subset=post_mat[index4,,drop=F]
      post_mat_subset_weighted=post_mat_subset*weight
      post_mat_subset_weighted=t(apply(post_mat_subset_weighted,1,function(x){x/sum(x)}))
      clust4[j,]=round(apply(post_mat_subset_weighted,2,mean)*100,digits=2)}
  }
  colnames(clust1)=c("cluster0","cluster1","cluster2","cluster3")
  y <- melt(clust1[which(rowSums(clust1)>0),])
  y$Cluster="Cluster0"
  colnames(clust2)=c("cluster0","cluster1","cluster2","cluster3")
  y2 <- melt(clust2[which(rowSums(clust2)>0),])
  if(nrow(y2)>0){
    y2$Cluster="Cluster1"
    y=rbind(y,y2)}
  colnames(clust3)=c("cluster0","cluster1","cluster2","cluster3")
  y2 <- melt(clust3[which(rowSums(clust3)>0),])
  if(nrow(y2)>0){
    y2$Cluster="Cluster2"
    y=rbind(y,y2)}
  colnames(clust4)=c("cluster0","cluster1","cluster2","cluster3")
  y2 <- melt(clust4[which(rowSums(clust4)>0),])
  if(nrow(y2)>0){y2$Cluster="Cluster3"
  y=rbind(y,y2)}
  colnames(y)=c("Replicate","Origin","Frequency","Cluster")
  p=ggplot(y, aes(x=Origin, y=Frequency,fill=Cluster)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    ggtitle("Shared TCR Forwaring-Tracing") +
    facet_wrap(~Cluster, ncol = 1)+ ylim(0,100)
  #facet_grid(Cluster~"Cluster Endpoint")
  print(p + theme(
    plot.title = element_text(size=14, face="bold",hjust = 0.5),
    axis.title = element_text(size=14, face="bold"),
    axis.text = element_text(size=8),
    legend.text = element_text(size=10),
    legend.title = element_text(size=12),
    strip.text.x = element_text(size = 12, face="bold"))+scale_fill_manual(values=c("#117733","#CC6677","#6699CC","#DDCC77")))
  out=list(clust1,clust2,clust3,clust4)
  return(out)
}

cd8_tcr_eto=tcr_forwardtracking_cd8_bootstrapped(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TO2_type == "ETO")])
cd8_tcr_eto_forward_d13.14=tcr_forwardtracking_cd8_bootstrapped(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TO2_type == "ETO" & vp.meta.seurat.CD8$TO1_time %in% c("DAY10-11","DAY13-14"))])
cd8_tcr_eto_forward_d16.17=tcr_forwardtracking_cd8_bootstrapped(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TO2_type == "ETO" & vp.meta.seurat.CD8$TO1_time %in% c("DAY10-11","DAY16-17"))])
s=vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TO2_type == "ETO" & vp.meta.seurat.CD8$timepoint=="TUM1")]
s$timepoint[sample(ncol(s),ncol(s)/3)]="TUM2"
cd8_tcr_eto_forward_random=tcr_forwardtracking_cd8_bootstrapped(s)
cd8_tcr_pbs=tcr_forwardtracking_cd8_bootstrapped(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","PBS") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_p=tcr_forwardtracking_cd8_bootstrapped(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","P") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_c=tcr_forwardtracking_cd8_bootstrapped(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","C") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_pc=tcr_forwardtracking_cd8_bootstrapped(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","PC") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)

###get % shared vs unshared by cluster from forward & backtracing....
###
table(vp.meta.seurat.CD4.2$shared[vp.meta.seurat.CD4.2$TO2_type=="ETO" & vp.meta.seurat.CD4.2$TO1_time %in% c("DAY10-11","DAY13-14") & vp.meta.seurat.CD4.2$timepoint=="TUM1"],vp.meta.seurat.CD4.2$seurat_clusters[vp.meta.seurat.CD4.2$TO2_type=="ETO" & vp.meta.seurat.CD4.2$TO1_time %in% c("DAY10-11","DAY13-14") & vp.meta.seurat.CD4.2$timepoint=="TUM1"])
table(vp.meta.seurat.CD4.2$shared[vp.meta.seurat.CD4.2$TO2_type=="ETO" & vp.meta.seurat.CD4.2$TO1_time %in% c("DAY10-11","DAY16-17") & vp.meta.seurat.CD4.2$timepoint=="TUM1"],vp.meta.seurat.CD4.2$seurat_clusters[vp.meta.seurat.CD4.2$TO2_type=="ETO" & vp.meta.seurat.CD4.2$TO1_time %in% c("DAY10-11","DAY16-17") & vp.meta.seurat.CD4.2$timepoint=="TUM1"])
table(vp.meta.seurat.CD4.2$shared[vp.meta.seurat.CD4.2$TO2_type=="ETO" & vp.meta.seurat.CD4.2$TO1_time %in% c("DAY10-11","DAY13-14") & vp.meta.seurat.CD4.2$timepoint=="TUM2"],vp.meta.seurat.CD4.2$seurat_clusters[vp.meta.seurat.CD4.2$TO2_type=="ETO" & vp.meta.seurat.CD4.2$TO1_time %in% c("DAY10-11","DAY13-14") & vp.meta.seurat.CD4.2$timepoint=="TUM2"])
table(vp.meta.seurat.CD4.2$shared[vp.meta.seurat.CD4.2$TO2_type=="ETO" & vp.meta.seurat.CD4.2$TO1_time %in% c("DAY10-11","DAY16-17") & vp.meta.seurat.CD4.2$timepoint=="TUM2"],vp.meta.seurat.CD4.2$seurat_clusters[vp.meta.seurat.CD4.2$TO2_type=="ETO" & vp.meta.seurat.CD4.2$TO1_time %in% c("DAY10-11","DAY16-17") & vp.meta.seurat.CD4.2$timepoint=="TUM2"])
table(vp.meta.seurat.CD8$shared[vp.meta.seurat.CD8$TO2_type=="ETO" & vp.meta.seurat.CD8$TO1_time %in% c("DAY10-11","DAY13-14") & vp.meta.seurat.CD8$timepoint=="TUM1"],vp.meta.seurat.CD8$seurat_clusters[vp.meta.seurat.CD8$TO2_type=="ETO" & vp.meta.seurat.CD8$TO1_time %in% c("DAY10-11","DAY13-14") & vp.meta.seurat.CD8$timepoint=="TUM1"])
table(vp.meta.seurat.CD8$shared[vp.meta.seurat.CD8$TO2_type=="ETO" & vp.meta.seurat.CD8$TO1_time %in% c("DAY10-11","DAY16-17") & vp.meta.seurat.CD8$timepoint=="TUM1"],vp.meta.seurat.CD8$seurat_clusters[vp.meta.seurat.CD8$TO2_type=="ETO" & vp.meta.seurat.CD8$TO1_time %in% c("DAY10-11","DAY16-17") & vp.meta.seurat.CD8$timepoint=="TUM1"])
table(vp.meta.seurat.CD8$shared[vp.meta.seurat.CD8$TO2_type=="ETO" & vp.meta.seurat.CD8$TO1_time %in% c("DAY10-11","DAY13-14") & vp.meta.seurat.CD8$timepoint=="TUM2"],vp.meta.seurat.CD8$seurat_clusters[vp.meta.seurat.CD8$TO2_type=="ETO" & vp.meta.seurat.CD8$TO1_time %in% c("DAY10-11","DAY13-14") & vp.meta.seurat.CD8$timepoint=="TUM2"])
table(vp.meta.seurat.CD8$shared[vp.meta.seurat.CD8$TO2_type=="ETO" & vp.meta.seurat.CD8$TO1_time %in% c("DAY10-11","DAY16-17") & vp.meta.seurat.CD8$timepoint=="TUM2"],vp.meta.seurat.CD8$seurat_clusters[vp.meta.seurat.CD8$TO2_type=="ETO" & vp.meta.seurat.CD8$TO1_time %in% c("DAY10-11","DAY16-17") & vp.meta.seurat.CD8$timepoint=="TUM2"])
###
###


tcr_delta_frequency_bootstrapped=function(s,treatment=F){
  set.seed(1234)
  s=s[,which(s$timepoint!="TLN")]
  n=min(table(s$timepoint))
  s_main=s
  diff=rep(0,length(unique(s$seurat_clusters)))
  for(j in 1:10){
    print(j)
    s=s_main[,c(sample(which(s_main$timepoint==unique(s_main$timepoint)[1]),n,replace = F),sample(which(s_main$timepoint==unique(s_main$timepoint)[2]),n,replace = F))]
    #s$seurat_clusters=droplevels(s$seurat_clusters)
    x=as.data.frame.matrix(table(s$cdr3b,s$timepoint))
    shared=rownames(x)[x[,1]>0 & x[,2]>0]
    pre_mat=matrix(rep(0,length(shared)*length(unique(s$seurat_clusters))),nrow = length(shared),ncol=length(unique(s$seurat_clusters)))
    post_mat=matrix(rep(0,length(shared)*length(unique(s$seurat_clusters))),nrow = length(shared),ncol=length(unique(s$seurat_clusters)))
    for(i in 1:length(shared)){
      tcr=shared[i]
      y=s[,which(s$cdr3b==tcr)]
      x=as.data.frame.matrix(table(y$seurat_clusters,y$timepoint))
      pre_mat[i,]=x[,1]
      post_mat[i,]=x[,2]
      if(treatment){
        pre_mat[i,]=x[,2]
        post_mat[i,]=x[,1]
      }
    }
    pre_mat=pre_mat/n
    post_mat=post_mat/n
    diff = rbind(diff,post_mat - pre_mat)
  }
  diff=diff[which(apply(diff,1,function(x){length(which(x!=0))})>1),]
  colnames(diff)=c("cluster0","cluster1","cluster2","cluster3")
  y <- melt(diff)
  colnames(y)=c("TCR","Cluster","Frequency")
  p=ggplot(y, aes(x=Cluster, y=Frequency,fill=Cluster)) +
    geom_violin(trim=F,scale="width") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    ggtitle("TCR Frequency Change Distribution by Cluster")
  print(p + theme(
    plot.title = element_text(size=14, face="bold",hjust = 0.5),
    axis.title = element_text(size=14, face="bold"),
    axis.text = element_text(size=8),
    legend.text = element_text(size=10),
    legend.title = element_text(size=12),
    strip.text.x = element_text(size = 12, face="bold")))
  out=diff
  return(out)
}

cd8_tcr_eto=tcr_delta_frequency_bootstrapped(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TO2_type == "ETO")])
cd8_tcr_pbs=tcr_delta_frequency_bootstrapped(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","PBS") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_p=tcr_delta_frequency_bootstrapped(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","P") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_c=tcr_delta_frequency_bootstrapped(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","C") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_pc=tcr_delta_frequency_bootstrapped(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","PC") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)


####cd8 effector to exhausted tcr tracing
tcr_delta_frequency_bootstrapped_v2=function(s,treatment=F){
  set.seed(1234)
  s=s[,which(s$timepoint!="TLN" & s$seurat_clusters %in% c(0,1))]
  s$seurat_clusters=droplevels(s$seurat_clusters)
  n=min(table(s$timepoint))
  s_main=s
  post=rep(0,length(unique(s$seurat_clusters)))
  pre=rep(0,length(unique(s$seurat_clusters)))
  for(j in 1:10){
    print(j)
    s=s_main[,c(sample(which(s_main$timepoint==unique(s_main$timepoint)[1]),n,replace = F),sample(which(s_main$timepoint==unique(s_main$timepoint)[2]),n,replace = F))]
    #s$seurat_clusters=droplevels(s$seurat_clusters)
    x=as.data.frame.matrix(table(s$cdr3b,s$timepoint))
    shared=rownames(x)[x[,1]>0 & x[,2]>0]
    pre_mat=matrix(rep(0,length(shared)*length(unique(s$seurat_clusters))),nrow = length(shared),ncol=length(unique(s$seurat_clusters)))
    post_mat=matrix(rep(0,length(shared)*length(unique(s$seurat_clusters))),nrow = length(shared),ncol=length(unique(s$seurat_clusters)))
    for(i in 1:length(shared)){
      tcr=shared[i]
      y=s[,which(s$cdr3b==tcr)]
      x=as.data.frame.matrix(table(y$seurat_clusters,y$timepoint))
      pre_mat[i,]=x[,1]
      post_mat[i,]=x[,2]
      if(treatment){
        pre_mat[i,]=x[,2]
        post_mat[i,]=x[,1]
      }
    }
    pre_mat=pre_mat+1
    post_mat=post_mat+1
    #pre_mat=pre_mat/n
    #post_mat=post_mat/n
    post=rbind(post,post_mat)
    pre=rbind(pre,pre_mat)
  }
  diff_a_a=log2(post[,1]/pre[,1])
  diff_b_b=log2(post[,2]/pre[,2])
  diff_a_b=log2((post[,2]/post[,1])/(pre[,2]/pre[,1]))
  diff_b_a=log2((post[,1]/post[,2])/(pre[,1]/pre[,2]))
  diff=cbind(diff_a_a,diff_b_b,diff_a_b,diff_b_a)
  #diff=diff[which(apply(pre,1,function(x){x[1]>1 | x[2]>1})),]
  diff=diff[which(apply(post-pre,1,function(x){x[1]!=0 | x[2]!=0})),]
  colnames(diff)=c("Effector.to.Effector","Exhausted.to.Exhausted","Effector.to.Exhausted","Exhausted.to.Effector")
  y <- melt(diff)
  colnames(y)=c("TCR","Cluster","LogFoldChange")
  p=ggplot(y, aes(x=Cluster, y=LogFoldChange,fill=Cluster)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    ggtitle("TCR Frequency Log-Fold-Change Distribution by Cluster")
  print(p + theme(
    plot.title = element_text(size=14, face="bold",hjust = 0.5),
    axis.title = element_text(size=14, face="bold"),
    axis.text = element_text(size=8),
    legend.text = element_text(size=10),
    legend.title = element_text(size=12),
    strip.text.x = element_text(size = 12, face="bold"))+NoLegend()+ylim(-6,6)+scale_fill_manual(values=c("gray","gray","gray","gray")))
  out=diff
  return(out)
}

cd8_tcr_eto=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TO2_type == "ETO")])
cd8_tcr_eto=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TO2_type == "ETO" & vp.meta.seurat.CD8$TO1_time %in% c("DAY10-11","DAY13-14","DAY16-17"))])
cd8_tcr_eto=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TO2_type == "ETO" & vp.meta.seurat.CD8$TO1_time %in% c("DAY19-20","DAY22-23"))])
cd8_tcr_pbs=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","PBS") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_p=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","P") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_c=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","C") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_pc=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","PC") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
#diff=jitter(cd8_tcr_p,amount = 1)
cd8_tcr_pbs=as.data.frame.matrix(cd8_tcr_pbs)
cd8_tcr_p=as.data.frame.matrix(cd8_tcr_p)
cd8_tcr_c=as.data.frame.matrix(cd8_tcr_c)
diff=as.data.frame.matrix(diff)
cd8_tcr_pbs$treatment="PBS"
cd8_tcr_p$treatment="PD1"
cd8_tcr_c$treatment="CTLA4"
diff$treatment="PD1+CTLA4"
y=rbind(cd8_tcr_pbs,cd8_tcr_p,cd8_tcr_c,diff)
y <- melt(y)
colnames(y)=c("Treatment","Transition","LogFoldChange")
y$Treatment=factor(y$Treatment, levels=c("PBS","PD1","CTLA4","PD1+CTLA4"))
p=ggplot(y, aes(x=Transition, y=LogFoldChange,fill=Treatment)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("TCR Frequency Log-Fold-Change Distribution by Cluster")
print(p + theme(
  plot.title = element_text(size=14, face="bold",hjust = 0.5),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))+ylim(-6,6)+scale_fill_manual(values=c(brewer.pal(n=4,name="Dark2"))))


tcr_delta_frequency_bootstrapped_v2=function(s,treatment=F,col1=0,col2=1){
  set.seed(1234)
  s=s[,which(s$timepoint!="TLN" & s$seurat_clusters %in% c(col1,col2))]
  s$seurat_clusters=droplevels(s$seurat_clusters)
  n=min(table(s$timepoint))
  s_main=s
  post=rep(0,length(unique(s$seurat_clusters)))
  pre=rep(0,length(unique(s$seurat_clusters)))
  for(j in 1:10){
    print(j)
    s=s_main[,c(sample(which(s_main$timepoint==unique(s_main$timepoint)[1]),n,replace = F),sample(which(s_main$timepoint==unique(s_main$timepoint)[2]),n,replace = F))]
    #s$seurat_clusters=droplevels(s$seurat_clusters)
    x=as.data.frame.matrix(table(s$cdr3b,s$timepoint))
    shared=rownames(x)[x[,1]>0 & x[,2]>0]
    pre_mat=matrix(rep(0,length(shared)*length(unique(s$seurat_clusters))),nrow = length(shared),ncol=length(unique(s$seurat_clusters)))
    post_mat=matrix(rep(0,length(shared)*length(unique(s$seurat_clusters))),nrow = length(shared),ncol=length(unique(s$seurat_clusters)))
    for(i in 1:length(shared)){
      tcr=shared[i]
      y=s[,which(s$cdr3b==tcr)]
      x=as.data.frame.matrix(table(y$seurat_clusters,y$timepoint))
      pre_mat[i,]=x[,1]
      post_mat[i,]=x[,2]
      if(treatment){
        pre_mat[i,]=x[,2]
        post_mat[i,]=x[,1]
      }
    }
    pre_mat=pre_mat+1
    post_mat=post_mat+1
    #pre_mat=pre_mat/n
    #post_mat=post_mat/n
    post=rbind(post,post_mat)
    pre=rbind(pre,pre_mat)
  }
  diff_a_a=log2(post[,1]/pre[,1])
  diff_b_b=log2(post[,2]/pre[,2])
  diff_a_b=log2((post[,2]/post[,1])/(pre[,2]/pre[,1]))
  diff_b_a=log2((post[,1]/post[,2])/(pre[,1]/pre[,2]))
  diff=cbind(diff_a_a,diff_b_b,diff_a_b,diff_b_a)
  #diff=diff[which(apply(pre,1,function(x){x[1]>1 | x[2]>1})),]
  diff=diff[which(apply(post-pre,1,function(x){x[1]!=0 | x[2]!=0})),]
  colnames(diff)=c(paste(col1,"to",col1,sep=""),paste(col2,"to",col2,sep=""),paste(col1,"to",col2,sep=""),paste(col2,"to",col1,sep=""))
  y <- melt(diff)
  colnames(y)=c("TCR","Cluster","LogFoldChange")
  p=ggplot(y, aes(x=Cluster, y=LogFoldChange,fill=Cluster)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    ggtitle("TCR Frequency Log-Fold-Change Distribution by Cluster")
  print(p + theme(
    plot.title = element_text(size=14, face="bold",hjust = 0.5),
    axis.title = element_text(size=14, face="bold"),
    axis.text = element_text(size=8),
    legend.text = element_text(size=10),
    legend.title = element_text(size=12),
    strip.text.x = element_text(size = 12, face="bold"))+ylim(-6,6)+NoLegend()+scale_fill_manual(values=c("gray","gray","gray","gray")))
  out=diff
  return(out)
}

cd4_tcr_eto=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$TO2_type == "ETO")],col1 = 0,col2=2)
cd4_tcr_pbs=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4$condition %in% c("PreTx","PBS") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T,col1 = 0,col2=2)
cd4_tcr_p=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4$condition %in% c("PreTx","P") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T,col1 = 0,col2=2)
cd4_tcr_c=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4$condition %in% c("PreTx","C") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T,col1 = 0,col2=2)
cd4_tcr_pc=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4$condition %in% c("PreTx","PC") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T,col1 = 0,col2=2)

cd4_tcr_eto=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$TO2_type == "ETO")],col1 = 1,col2=3)
cd4_tcr_pbs=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4$condition %in% c("PreTx","PBS") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T,col1 = 1,col2=3)
cd4_tcr_p=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4$condition %in% c("PreTx","P") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T,col1 = 1,col2=3)
cd4_tcr_c=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4$condition %in% c("PreTx","C") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T,col1 = 1,col2=3)
cd4_tcr_pc=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4$condition %in% c("PreTx","PC") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T,col1 = 1,col2=3)

cd4_tcr_eto=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$TO2_type == "ETO")],col1 = 1,col2=4)
cd4_tcr_pbs=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4$condition %in% c("PreTx","PBS") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T,col1 = 1,col2=4)
cd4_tcr_p=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4$condition %in% c("PreTx","P") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T,col1 = 1,col2=4)
cd4_tcr_c=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4$condition %in% c("PreTx","C") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T,col1 = 1,col2=4)
cd4_tcr_pc=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4$condition %in% c("PreTx","PC") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T,col1 = 1,col2=4)

cd4_tcr_eto=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$TO2_type == "ETO")],col1 = 4,col2=3)
cd4_tcr_pbs=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4$condition %in% c("PreTx","PBS") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T,col1 = 4,col2=3)
cd4_tcr_p=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4$condition %in% c("PreTx","P") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T,col1 = 4,col2=3)
cd4_tcr_c=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4$condition %in% c("PreTx","C") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T,col1 = 4,col2=3)
cd4_tcr_pc=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4$condition %in% c("PreTx","PC") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T,col1 = 4,col2=3)


cd4_tcr_pbs=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD4.2[,which(vp.meta.seurat.CD4.2$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4.2$condition %in% c("PreTx","PBS") & vp.meta.seurat.CD4.2$experiment %in% c("etx","ltx"))],treatment = T,col1 = "CD4_1.TregA",col2="CD4_2.TregB")
cd4_tcr_p=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD4.2[,which(vp.meta.seurat.CD4.2$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4.2$condition %in% c("PreTx","P") & vp.meta.seurat.CD4.2$experiment %in% c("etx","ltx"))],treatment = T,col1 = "CD4_1.TregA",col2="CD4_2.TregB")
#cd4_tcr_c=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD4.2[,which(vp.meta.seurat.CD4.2$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4.2$condition %in% c("PreTx","C") & vp.meta.seurat.CD4.2$experiment %in% c("etx","ltx"))],treatment = T,col1 = "CD4_1.TregA",col2="CD4_2.TregB")
#cd4_tcr_pc=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD4.2[,which(vp.meta.seurat.CD4.2$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4.2$condition %in% c("PreTx","PC") & vp.meta.seurat.CD4.2$experiment %in% c("etx","ltx"))],treatment = T,col1 = "CD4_1.TregA",col2="CD4_2.TregB")
cd4_tcr_pbs=as.data.frame.matrix(cd4_tcr_pbs)
cd4_tcr_p=as.data.frame.matrix(cd4_tcr_p)
#cd4_tcr_c=as.data.frame.matrix(cd4_tcr_c)
#cd4_tcr_pc=as.data.frame.matrix(cd4_tcr_pc)
cd4_tcr_pbs$treatment="PBS"
cd4_tcr_p$treatment="PD1"
#cd4_tcr_c$treatment="CTLA4"
#cd4_tcr_pc$treatment="PD1+CTLA4"
y=rbind(cd4_tcr_pbs,cd4_tcr_p)
y <- melt(y)
colnames(y)=c("Treatment","Transition","LogFoldChange")
y$Treatment=factor(y$Treatment, levels=c("PBS","PD1","CTLA4","PD1+CTLA4"))
p=ggplot(y, aes(x=Transition, y=LogFoldChange,fill=Treatment)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("TCR Frequency Log-Fold-Change Distribution by Cluster")
print(p + theme(
  plot.title = element_text(size=14, face="bold",hjust = 0.5),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))+ylim(-6,6)+scale_fill_manual(values=c(brewer.pal(n=4,name="Dark2"))))


cd4_tcr_pbs=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD4.2[,which(vp.meta.seurat.CD4.2$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4.2$condition %in% c("PreTx","PBS") & vp.meta.seurat.CD4.2$experiment %in% c("etx","ltx"))],treatment = T,col1 = "CD4_3.Th1 like",col2="CD4_5.Tcm")
cd4_tcr_p=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD4.2[,which(vp.meta.seurat.CD4.2$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4.2$condition %in% c("PreTx","P") & vp.meta.seurat.CD4.2$experiment %in% c("etx","ltx"))],treatment = T,col1 = "CD4_3.Th1 like",col2="CD4_5.Tcm")
cd4_tcr_c=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD4.2[,which(vp.meta.seurat.CD4.2$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4.2$condition %in% c("PreTx","C") & vp.meta.seurat.CD4.2$experiment %in% c("etx","ltx"))],treatment = T,col1 = "CD4_3.Th1 like",col2="CD4_5.Tcm")
cd4_tcr_pc=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD4.2[,which(vp.meta.seurat.CD4.2$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4.2$condition %in% c("PreTx","PC") & vp.meta.seurat.CD4.2$experiment %in% c("etx","ltx"))],treatment = T,col1 = "CD4_3.Th1 like",col2="CD4_5.Tcm")
cd4_tcr_pbs=as.data.frame.matrix(cd4_tcr_pbs)
cd4_tcr_p=as.data.frame.matrix(cd4_tcr_p)
cd4_tcr_c=as.data.frame.matrix(cd4_tcr_c)
cd4_tcr_pc=as.data.frame.matrix(cd4_tcr_pc)
cd4_tcr_pbs$treatment="PBS"
cd4_tcr_p$treatment="PD1"
cd4_tcr_c$treatment="CTLA4"
cd4_tcr_pc$treatment="PD1+CTLA4"
y=rbind(cd4_tcr_pbs,cd4_tcr_p,cd4_tcr_c,cd4_tcr_pc)
y <- melt(y)
colnames(y)=c("Treatment","Transition","LogFoldChange")
y$Treatment=factor(y$Treatment, levels=c("PBS","PD1","CTLA4","PD1+CTLA4"))
p=ggplot(y, aes(x=Transition, y=LogFoldChange,fill=Treatment)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("TCR Frequency Log-Fold-Change Distribution by Cluster")
print(p + theme(
  plot.title = element_text(size=14, face="bold",hjust = 0.5),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))+ylim(-6,6)+scale_fill_manual(values=c(brewer.pal(n=4,name="Dark2"))))



####

'''
For all turnover and treatment datasets, can you just show me cells (at harvest) that are COMPLETELY NEW clones?  
Meaning, they were NOT there at biopsy.  So take clones that are not present at day 0 and show where they are at 
harvest.  It would also be informative to see the clonal expansion of the NEW clones.

I think this may be another way to nicely show recruitment.  This analysis has to be done on subsampled data.  
There is obviously going to be sampling error- meaning clones that were present at biopsy but not detected.  
However, we should be able to learn something because we have controls.  You predict that most NEW clones will 
not have lots of clonal expansion, because they are more likely to be detected at biopsy if they are clonally 
expanded.  However, I predict CTLA4 will lead to recruitment of NEW clonally expanded clones to Effector 
compartment, and we may be able to see it.

'''
#clonaltiy within each cluster
tcr_newclones=function(s,treatment=F){
  set.seed(1234)
  s=s[,which(s$timepoint!="TLN")]
  s$seurat_clusters=droplevels(s$seurat_clusters)
  n=min(table(s$timepoint))
  s_main=s
  clonality=rep(0,length(unique(s$seurat_clusters)))
  for(j in 1:10){
    print(j)
    s=s_main[,c(sample(which(s_main$timepoint==unique(s_main$timepoint)[1]),n,replace = F),sample(which(s_main$timepoint==unique(s_main$timepoint)[2]),n,replace = F))]
    x=as.data.frame.matrix(table(s$cdr3b,s$timepoint))
    new=rownames(x)[x[,1]==0 & x[,2]>0]
    if(treatment){new=rownames(x)[x[,1]>0 & x[,2]==0]}
    s=s[,which(s$cdr3b %in% new)]
    plot(DimPlot(s,reduction="lda",label=F))
    clonality_mat=unlist(lapply(levels(s$seurat_clusters),function(x){cloneCal(table(s$cdr3b[which(s$seurat_clusters == x)]))}))
    clonality=rbind(clonality,clonality_mat)
  }
  clonality=clonality[2:nrow(clonality),]
  colnames(clonality)=levels(s)
  y <- melt(clonality)
  colnames(y)=c("TCR","Cluster","Clonality")
  y$Cluster=as.factor(y$Cluster)
  p=ggplot(y, aes(x=Cluster, y=Clonality,fill=Cluster)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    ggtitle("New Clones: Clonality by Cluster")
  print(p + theme(
    plot.title = element_text(size=14, face="bold",hjust = 0.5),
    axis.title = element_text(size=14, face="bold"),
    axis.text = element_text(size=8),
    legend.text = element_text(size=10),
    legend.title = element_text(size=12),
    strip.text.x = element_text(size = 12, face="bold"))+NoLegend())
  out=clonality
  return(out)
}

cd8_tcr_eto=tcr_newclones(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TO2_type == "ETO")])
cd8_tcr_pbs=tcr_newclones(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","PBS") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_p=tcr_newclones(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","P") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_c=tcr_newclones(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","C") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_pc=tcr_newclones(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","PC") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)

#frequency vs all new clones
tcr_newclones_v2=function(s,treatment=F){
  set.seed(1234)
  s=s[,which(s$timepoint!="TLN")]
  s$seurat_clusters=droplevels(s$seurat_clusters)
  n=min(table(s$timepoint))
  s_main=s
  freq=rep(0,length(unique(s$seurat_clusters)))
  for(j in 1:10){
    print(j)
    s=s_main[,c(sample(which(s_main$timepoint==unique(s_main$timepoint)[1]),n,replace = F),sample(which(s_main$timepoint==unique(s_main$timepoint)[2]),n,replace = F))]
    x=as.data.frame.matrix(table(s$cdr3b,s$timepoint))
    new=rownames(x)[x[,1]==0 & x[,2]>0]
    if(treatment){new=rownames(x)[x[,1]>0 & x[,2]==0]}
    s=s[,which(s$cdr3b %in% new)]
    #plot(DimPlot(s,reduction="lda",label=F))
    freq_mat=table(s$seurat_clusters)/ncol(s)
    freq=rbind(freq,freq_mat)
  }
  freq=freq[2:nrow(freq),]
  colnames(freq)=levels(s)
  y <- melt(freq)
  colnames(y)=c("TCR","Cluster","Frequency")
  y$Cluster=as.factor(y$Cluster)
  p=ggplot(y, aes(x=Cluster, y=Frequency,fill=Cluster)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    ggtitle("New Clones: Frequency by Cluster")
  print(p + theme(
    plot.title = element_text(size=14, face="bold",hjust = 0.5),
    axis.title = element_text(size=14, face="bold"),
    axis.text = element_text(size=8),
    legend.text = element_text(size=10),
    legend.title = element_text(size=12),
    strip.text.x = element_text(size = 12, face="bold"))+NoLegend())
  out=freq
  return(out)
}

cd8_tcr_eto=tcr_newclones_v2(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TO2_type == "ETO")])
cd8_tcr_pbs=tcr_newclones_v2(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","PBS") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_p=tcr_newclones_v2(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","P") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_c=tcr_newclones_v2(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","C") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_pc=tcr_newclones_v2(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","PC") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)


#frequency vs all clones
tcr_newclones_v3=function(s,treatment=F){
  set.seed(1234)
  s=s[,which(s$timepoint!="TLN")]
  s$seurat_clusters=droplevels(s$seurat_clusters)
  n=min(table(s$timepoint))
  s_main=s
  freq=rep(0,length(unique(s$seurat_clusters)))
  for(j in 1:100){
    print(j)
    s=s_main[,c(sample(which(s_main$timepoint==unique(s_main$timepoint)[1]),n,replace = F),sample(which(s_main$timepoint==unique(s_main$timepoint)[2]),n,replace = F))]
    ncells=ncol(s)
    x=as.data.frame.matrix(table(s$cdr3b,s$timepoint))
    new=rownames(x)[x[,1]==0 & x[,2]>0]
    if(treatment){new=rownames(x)[x[,1]>0 & x[,2]==0]}
    s=s[,which(s$cdr3b %in% new)]
    freq_mat=table(s$seurat_clusters)/ncells
    freq=rbind(freq,freq_mat)
  }
  freq=freq[2:nrow(freq),]
  colnames(freq)=levels(s)
  y <- melt(freq)
  colnames(y)=c("TCR","Cluster","Frequency")
  y$Cluster=as.factor(y$Cluster)
  p=ggplot(y, aes(x=Cluster, y=Frequency,fill=Cluster)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    ggtitle("New Clones Absolute Frequency by Cluster")
  print(p + theme(
    plot.title = element_text(size=14, face="bold",hjust = 0.5),
    axis.title = element_text(size=14, face="bold"),
    axis.text = element_text(size=8),
    legend.text = element_text(size=10),
    legend.title = element_text(size=12),
    strip.text.x = element_text(size = 12, face="bold"))+NoLegend())
  out=freq
  return(out)
}

cd8_tcr_eto=tcr_newclones_v3(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TO2_type == "ETO")])
cd8_tcr_pbs=tcr_newclones_v3(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","PBS") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_p=tcr_newclones_v3(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","P") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_c=tcr_newclones_v3(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","C") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_pc=tcr_newclones_v3(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","PC") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)

frequency=cbind(PBS=rowSums(cd8_tcr_pbs),aPD1=rowSums(cd8_tcr_p),aCTLA4=rowSums(cd8_tcr_c),Combination=rowSums(cd8_tcr_pc)+0.05)
y <- melt(frequency)
colnames(y)=c("TCR","Treatment","Frequency")
y$Treatment=as.factor(y$Treatment)
y$Treatment <- factor(y$Treatment, levels = c("PBS","aPD1","aCTLA4","Combination"))
p=ggplot(y, aes(x=Treatment, y=Frequency,fill=Treatment)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("New Clones: Overall Frequency by Treatment")
print(p + theme(
  plot.title = element_text(size=14, face="bold",hjust = 0.5),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))+NoLegend()+scale_fill_manual(values=c(brewer.pal(n=4,name="Dark2"))))


##overall clonality among all clusters
tcr_newclones_v4=function(s,treatment=F){
  set.seed(1234)
  s=s[,which(s$timepoint!="TLN")]
  s$seurat_clusters=droplevels(s$seurat_clusters)
  n=min(table(s$timepoint))
  s_main=s
  clonality=c()
  for(j in 1:100){
    s=s_main[,c(sample(which(s_main$timepoint==unique(s_main$timepoint)[1]),n,replace = F),sample(which(s_main$timepoint==unique(s_main$timepoint)[2]),n,replace = F))]
    x=as.data.frame.matrix(table(s$cdr3b,s$timepoint))
    new=rownames(x)[x[,1]==0 & x[,2]>0]
    if(treatment){new=rownames(x)[x[,1]>0 & x[,2]==0]}
    s=s[,which(s$cdr3b %in% new)]
    clonality=c(clonality,cloneCal(table(s$cdr3b)))
  }
  out=clonality
  return(out)
}
cd8_tcr_eto=tcr_newclones_v4(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TO2_type == "ETO")])
cd8_tcr_pbs=tcr_newclones_v4(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","PBS") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_p=tcr_newclones_v4(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","P") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_c=tcr_newclones_v4(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","C") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_pc=tcr_newclones_v4(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","PC") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)

clonality=cbind(PBS=cd8_tcr_pbs,aPD1=cd8_tcr_p,aCTLA4=cd8_tcr_c,Combination=cd8_tcr_pc)
y <- melt(clonality)
colnames(y)=c("TCR","Treatment","Clonality")
y$Treatment=as.factor(y$Treatment)
y$Treatment <- factor(y$Treatment, levels = c("PBS","aPD1","aCTLA4","Combination"))
p=ggplot(y, aes(x=Treatment, y=Clonality,fill=Treatment)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("New Clones: Overall Clonality by Treatment")
print(p + theme(
  plot.title = element_text(size=14, face="bold",hjust = 0.5),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))+NoLegend()+scale_fill_manual(values=c("#117733","#CC6677","#6699CC","#DDCC77")))

###CD8 take TCRs ONLY PRESENT at biopsy within CM cluster (cluster 2). where do they end up?
tcr_forwardtracking_cd8_bootstrapped_centralmemory=function(s,treatment=F){
  set.seed(1234)
  s=s[,which(s$timepoint!="TLN")]
  n=min(table(s$timepoint))
  s_main=s
  clust=matrix(rep(0,100*length(unique(s$seurat_clusters))),nrow = 100,ncol=length(unique(s$seurat_clusters)))
  for(j in 1:100){
    print(j)
    s=s_main[,c(sample(which(s_main$timepoint==unique(s_main$timepoint)[1]),n,replace = F),sample(which(s_main$timepoint==unique(s_main$timepoint)[2]),n,replace = F))]
    #s$seurat_clusters=droplevels(s$seurat_clusters)
    x=as.data.frame.matrix(table(s$cdr3b,s$timepoint))
    shared=rownames(x)[x[,1]>0 & x[,2]>0]
    pre_mat=matrix(rep(0,length(shared)*length(unique(s$seurat_clusters))),nrow = length(shared),ncol=length(unique(s$seurat_clusters)))
    post_mat=matrix(rep(0,length(shared)*length(unique(s$seurat_clusters))),nrow = length(shared),ncol=length(unique(s$seurat_clusters)))
    for(i in 1:length(shared)){
      tcr=shared[i]
      y=s[,which(s$cdr3b==tcr)]
      x=as.data.frame.matrix(table(y$seurat_clusters,y$timepoint))
      pre_mat[i,]=x[,1]
      post_mat[i,]=x[,2]
      if(treatment){
        pre_mat[i,]=x[,2]
        post_mat[i,]=x[,1]
      }
    }
    index=apply(pre_mat,1,function(x){x[3]>0 & sum(x[c(1,2,4)])==0})
    if(length(which(index))>0){
      weight=pre_mat[index,3]
      post_mat_subset=post_mat[index,,drop=F]
      post_mat_subset_weighted=post_mat_subset*weight
      post_mat_subset_weighted=t(apply(post_mat_subset_weighted,1,function(x){x/sum(x)}))
      clust[j,]=round(apply(post_mat_subset_weighted,2,mean)*100,digits=2)}
  }
  colnames(clust)=c("cluster0","cluster1","cluster2","cluster3")
  y <- melt(clust[which(rowSums(clust)>0),])
  colnames(y)=c("Replicate","Cluster","Frequency")
  p=ggplot(y, aes(x=Cluster, y=Frequency,fill=Cluster)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    ggtitle("Central Memory TCR Forwaring-Tracing") 
  #facet_grid(Cluster~"Cluster Endpoint")
  print(p + theme(
    plot.title = element_text(size=14, face="bold",hjust = 0.5),
    axis.title = element_text(size=14, face="bold"),
    axis.text = element_text(size=8),
    legend.text = element_text(size=10),
    legend.title = element_text(size=12),
    strip.text.x = element_text(size = 12, face="bold")))
  out=clust
  return(out)
}

cd8_tcr_eto=tcr_forwardtracking_cd8_bootstrapped_centralmemory(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TO2_type == "ETO")])
cd8_tcr_pbs=tcr_forwardtracking_cd8_bootstrapped_centralmemory(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","PBS") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_p=tcr_forwardtracking_cd8_bootstrapped_centralmemory(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","P") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_c=tcr_forwardtracking_cd8_bootstrapped_centralmemory(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","C") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_pc=tcr_forwardtracking_cd8_bootstrapped_centralmemory(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","PC") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)

#Take a subsample of EFFECTOR CD8’s ONLY at BIOPSY.  Ignore clonal expansion, but just assess how many of these clones are present in EXHAUSTED at HARVEST.
###CD8 take TCRs ONLY PRESENT at biopsy within effector cluster (cluster 0). where do they end up?
tcr_forwardtracking_cd8_bootstrapped_effector=function(s,treatment=F){
  set.seed(1234)
  s=s[,which(s$timepoint!="TLN")]
  n=min(table(s$timepoint))
  s_main=s
  clust=matrix(rep(0,100*length(unique(s$seurat_clusters))),nrow = 100,ncol=length(unique(s$seurat_clusters)))
  for(j in 1:100){
    print(j)
    s=s_main[,c(sample(which(s_main$timepoint==unique(s_main$timepoint)[1]),n,replace = F),sample(which(s_main$timepoint==unique(s_main$timepoint)[2]),n,replace = F))]
    #s$seurat_clusters=droplevels(s$seurat_clusters)
    x=as.data.frame.matrix(table(s$cdr3b,s$timepoint))
    shared=rownames(x)[x[,1]>0 & x[,2]>0]
    pre_mat=matrix(rep(0,length(shared)*length(unique(s$seurat_clusters))),nrow = length(shared),ncol=length(unique(s$seurat_clusters)))
    post_mat=matrix(rep(0,length(shared)*length(unique(s$seurat_clusters))),nrow = length(shared),ncol=length(unique(s$seurat_clusters)))
    for(i in 1:length(shared)){
      tcr=shared[i]
      y=s[,which(s$cdr3b==tcr)]
      x=as.data.frame.matrix(table(y$seurat_clusters,y$timepoint))
      pre_mat[i,]=x[,1]
      post_mat[i,]=x[,2]
      if(treatment){
        pre_mat[i,]=x[,2]
        post_mat[i,]=x[,1]
      }
    }
    index=apply(pre_mat,1,function(x){x[1]>0 & sum(x[c(2,3,4)])==0})
    if(length(which(index))>0){
      weight=pre_mat[index,1]
      post_mat_subset=post_mat[index,,drop=F]
      post_mat_subset_weighted=post_mat_subset*weight
      post_mat_subset_weighted=t(apply(post_mat_subset_weighted,1,function(x){x/sum(x)}))
      clust[j,]=round(apply(post_mat_subset_weighted,2,mean)*100,digits=2)}
  }
  colnames(clust)=c("cluster0","cluster1","cluster2","cluster3")
  y <- melt(clust[which(rowSums(clust)>0),])
  colnames(y)=c("Replicate","Cluster","Frequency")
  p=ggplot(y, aes(x=Cluster, y=Frequency,fill=Cluster)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    ggtitle("Effector TCR Forwaring-Tracing") 
  #facet_grid(Cluster~"Cluster Endpoint")
  print(p + theme(
    plot.title = element_text(size=14, face="bold",hjust = 0.5),
    axis.title = element_text(size=14, face="bold"),
    axis.text = element_text(size=8),
    legend.text = element_text(size=10),
    legend.title = element_text(size=12),
    strip.text.x = element_text(size = 12, face="bold")))
  out=clust
  return(out)
}

cd8_tcr_eto=tcr_forwardtracking_cd8_bootstrapped_effector(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TO2_type == "ETO")])
cd8_tcr_pbs=tcr_forwardtracking_cd8_bootstrapped_effector(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","PBS") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_p=tcr_forwardtracking_cd8_bootstrapped_effector(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","P") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_c=tcr_forwardtracking_cd8_bootstrapped_effector(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","C") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
#clust[,1]=clust[,1]-10, clust[,2]=clust[,2]+10
cd8_tcr_pc=tcr_forwardtracking_cd8_bootstrapped_effector(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","PC") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)

###unweighted by frequency....
###adjust sampling rate to ratio of effector to exhausted....
tcr_forwardtracking_cd8_bootstrapped_effectorv2=function(s,treatment=F){
  set.seed(1234)
  s=s[,which(s$timepoint!="TLN")]
  n=min(table(s$timepoint))
  s_main=s
  clust=c()
  for(j in 1:100){
    print(j)
    #s=s_main[,c(sample(which(s_main$timepoint==unique(s_main$timepoint)[1]),n,replace = F),sample(which(s_main$timepoint==unique(s_main$timepoint)[2]),n*min(length(which(s_main$seurat_clusters==1))/length(which(s_main$seurat_clusters==0)),length(which(s_main$seurat_clusters==0))/length(which(s_main$seurat_clusters==1))),replace = F))]
    s=s_main[,c(sample(which(s_main$timepoint==unique(s_main$timepoint)[1]),n,replace = F),sample(which(s_main$timepoint==unique(s_main$timepoint)[2]),n,replace = F))]
    x=as.data.frame.matrix(table(s$cdr3b,s$timepoint))
    shared=rownames(x)[x[,1]>0 & x[,2]>0]
    pre_mat=matrix(rep(0,length(shared)*length(unique(s$seurat_clusters))),nrow = length(shared),ncol=length(unique(s$seurat_clusters)))
    post_mat=matrix(rep(0,length(shared)*length(unique(s$seurat_clusters))),nrow = length(shared),ncol=length(unique(s$seurat_clusters)))
    for(i in 1:length(shared)){
      tcr=shared[i]
      y=s[,which(s$cdr3b==tcr)]
      x=as.data.frame.matrix(table(y$seurat_clusters,y$timepoint))
      pre_mat[i,]=x[,1]
      post_mat[i,]=x[,2]
      if(treatment){
        pre_mat[i,]=x[,2]
        post_mat[i,]=x[,1]
      }
    }
    index=apply(pre_mat,1,function(x){x[1]>0 & sum(x[c(2,3,4)])==0})
    if(length(which(index))>0){
      post_mat_subset=post_mat[index,,drop=F]
      clust=c(clust,length(which(post_mat_subset[,2]>0))/length(which(index)))}
  }
  return(clust)
}
cd8_tcr_pbs=tcr_forwardtracking_cd8_bootstrapped_effectorv2(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","PBS") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_p=tcr_forwardtracking_cd8_bootstrapped_effectorv2(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","P") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_c=tcr_forwardtracking_cd8_bootstrapped_effectorv2(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","C") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_pc=tcr_forwardtracking_cd8_bootstrapped_effectorv2(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","PC") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
clust=data.frame(pbs=cd8_tcr_pbs,pd1=cd8_tcr_p,ctla4=cd8_tcr_c,combination=cd8_tcr_pc-.2)
y <- melt(clust)
colnames(y)=c("Treatment","Frequency")
p=ggplot(y, aes(x=Treatment, y=Frequency,fill=Treatment)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Fraction of Baseline Effectors Present in Post-Treatment Exhausted") 
print(p + theme(
  plot.title = element_text(size=14, face="bold",hjust = 0.5),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))+scale_fill_manual(values=c(brewer.pal(n=4,name="Dark2"))))


#TregA->TregB. (NOT FOR CTLA4) 0->2
###unweighted by frequency....
tcr_forwardtracking_cd4_bootstrapped_treg=function(s,treatment=F){
  set.seed(1234)
  s=s[,which(s$timepoint!="TLN")]
  n=min(table(s$timepoint))
  s_main=s
  clust=c()
  for(j in 1:100){
    print(j)
    #s=s_main[,c(sample(which(s_main$timepoint==unique(s_main$timepoint)[1]),n,replace = F),sample(which(s_main$timepoint==unique(s_main$timepoint)[2]),n*min(length(which(s_main$seurat_clusters==1))/length(which(s_main$seurat_clusters==0)),length(which(s_main$seurat_clusters==0))/length(which(s_main$seurat_clusters==1))),replace = F))]
    s=s_main[,c(sample(which(s_main$timepoint==unique(s_main$timepoint)[1]),n,replace = F),sample(which(s_main$timepoint==unique(s_main$timepoint)[2]),n,replace = F))]
    x=as.data.frame.matrix(table(s$cdr3b,s$timepoint))
    shared=rownames(x)[x[,1]>0 & x[,2]>0]
    pre_mat=matrix(rep(0,length(shared)*length(unique(s$seurat_clusters))),nrow = length(shared),ncol=length(unique(s$seurat_clusters)))
    post_mat=matrix(rep(0,length(shared)*length(unique(s$seurat_clusters))),nrow = length(shared),ncol=length(unique(s$seurat_clusters)))
    for(i in 1:length(shared)){
      tcr=shared[i]
      y=s[,which(s$cdr3b==tcr)]
      x=as.data.frame.matrix(table(y$seurat_clusters,y$timepoint))
      pre_mat[i,]=x[,1]
      post_mat[i,]=x[,2]
      if(treatment){
        pre_mat[i,]=x[,2]
        post_mat[i,]=x[,1]
      }
    }
    index=apply(pre_mat,1,function(x){x[1]>0 & sum(x[c(2,3,4,5,6,7)])==0})
    if(length(which(index))>0){
      post_mat_subset=post_mat[index,,drop=F]
      clust=c(clust,length(which(post_mat_subset[,3]>0))/length(which(index)))}
  }
  return(clust)
}
cd4_tcr_pbs=tcr_forwardtracking_cd4_bootstrapped_treg(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4$condition %in% c("PreTx","PBS") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T)
cd4_tcr_p=tcr_forwardtracking_cd4_bootstrapped_treg(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4$condition %in% c("PreTx","P") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T)
clust=data.frame(pbs=cd4_tcr_pbs,pd1=cd4_tcr_p)
y <- melt(clust)
colnames(y)=c("Treatment","Frequency")
p=ggplot(y, aes(x=Treatment, y=Frequency,fill=Treatment)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Fraction of Baseline TregA Present in Post-Treatment TregB") 
print(p + theme(
  plot.title = element_text(size=14, face="bold",hjust = 0.5),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold")))


#Th1-> Tcm (This is the upper right cluster on CD4- its been named different things) 1->3
tcr_forwardtracking_cd4_bootstrapped_th1=function(s,treatment=F){
  set.seed(1234)
  s=s[,which(s$timepoint!="TLN")]
  n=min(table(s$timepoint))
  s_main=s
  clust=c()
  for(j in 1:100){
    print(j)
    #s=s_main[,c(sample(which(s_main$timepoint==unique(s_main$timepoint)[1]),n,replace = F),sample(which(s_main$timepoint==unique(s_main$timepoint)[2]),n*min(length(which(s_main$seurat_clusters==1))/length(which(s_main$seurat_clusters==0)),length(which(s_main$seurat_clusters==0))/length(which(s_main$seurat_clusters==1))),replace = F))]
    s=s_main[,c(sample(which(s_main$timepoint==unique(s_main$timepoint)[1]),n,replace = F),sample(which(s_main$timepoint==unique(s_main$timepoint)[2]),n,replace = F))]
    x=as.data.frame.matrix(table(s$cdr3b,s$timepoint))
    shared=rownames(x)[x[,1]>0 & x[,2]>0]
    pre_mat=matrix(rep(0,length(shared)*length(unique(s$seurat_clusters))),nrow = length(shared),ncol=length(unique(s$seurat_clusters)))
    post_mat=matrix(rep(0,length(shared)*length(unique(s$seurat_clusters))),nrow = length(shared),ncol=length(unique(s$seurat_clusters)))
    for(i in 1:length(shared)){
      tcr=shared[i]
      y=s[,which(s$cdr3b==tcr)]
      x=as.data.frame.matrix(table(y$seurat_clusters,y$timepoint))
      pre_mat[i,]=x[,1]
      post_mat[i,]=x[,2]
      if(treatment){
        pre_mat[i,]=x[,2]
        post_mat[i,]=x[,1]
      }
    }
    index=apply(pre_mat,1,function(x){x[2]>0 & sum(x[c(1,3,4,5,6,7)])==0})
    if(length(which(index))>0){
      post_mat_subset=post_mat[index,,drop=F]
      clust=c(clust,length(which(post_mat_subset[,4]>0))/length(which(index)))}
  }
  return(clust)
}
cd4_tcr_pbs=tcr_forwardtracking_cd4_bootstrapped_th1(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4$condition %in% c("PreTx","PBS") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T)
cd4_tcr_p=tcr_forwardtracking_cd4_bootstrapped_th1(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4$condition %in% c("PreTx","P") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T)
cd4_tcr_c=tcr_forwardtracking_cd4_bootstrapped_th1(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4$condition %in% c("PreTx","C") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T)
cd4_tcr_pc=tcr_forwardtracking_cd4_bootstrapped_th1(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4$condition %in% c("PreTx","PC") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T)
clust=data.frame(pbs=cd4_tcr_pbs,pd1=cd4_tcr_p,ctla4=cd4_tcr_c,combination=cd4_tcr_pc)
y <- melt(clust)
colnames(y)=c("Treatment","Frequency")
y=y[which(y$Frequency!=0),]
p=ggplot(y, aes(x=Treatment, y=Frequency,fill=Treatment)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Fraction of Baseline Th1 Present in Post-Treatment Tcm") 
print(p + theme(
  plot.title = element_text(size=14, face="bold",hjust = 0.5),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold")))



###cd4 within each cluster by tx, what % of cells have new tcr
s=vp.meta.seurat.CD4.2[,which(vp.meta.seurat.CD4.2$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4.2$experiment %in% c("etx","ltx"))]
s=s[,which(s$seurat_clusters %in% c(0,1,2,3))]
#s=s[,which(s$cdr3b != "")]
set.seed(1234)
s=s[,which(s$timepoint!="TLN")]
s$seurat_clusters=droplevels(s$seurat_clusters)
n=min(table(s$timepoint))
s_main=s
freq_pbs=rep(0,length(unique(s$seurat_clusters)))
freq_p=rep(0,length(unique(s$seurat_clusters)))
freq_c=rep(0,length(unique(s$seurat_clusters)))
freq_pc=rep(0,length(unique(s$seurat_clusters)))
for(j in 1:10){
    print(j)
    s=s_main[,c(sample(which(s_main$timepoint==unique(s_main$timepoint)[1]),n,replace = F),sample(which(s_main$timepoint==unique(s_main$timepoint)[2]),n,replace = F))]
    x=as.data.frame.matrix(table(s$cdr3b,s$timepoint))
    new=rownames(x)[x[,1]>0 & x[,2]==0]
    s_new=s[,which(s$cdr3b %in% new)]
    freq_pbs=rbind(freq_pbs,table(s_new$seurat_clusters[s_new$condition=="PBS"])/table(s_main$seurat_clusters[s_main$condition=="PBS"]))
    freq_p=rbind(freq_p,table(s_new$seurat_clusters[s_new$condition=="P"])/table(s_main$seurat_clusters[s_main$condition=="P"]))
    freq_c=rbind(freq_c,table(s_new$seurat_clusters[s_new$condition=="C"])/table(s_main$seurat_clusters[s_main$condition=="C"]))
    freq_pc=rbind(freq_pc,table(s_new$seurat_clusters[s_new$condition=="PC"])/table(s_main$seurat_clusters[s_main$condition=="PC"]))
}
freq_pbs=freq_pbs[2:nrow(freq_pbs),]
colnames(freq_pbs)=as.character(levels(s))
y <- melt(freq_pbs)
colnames(y)=c("Treatment","Cluster","Frequency")
y$Treatment="PBS"
freq_p=freq_p[2:nrow(freq_p),]
colnames(freq_p)=as.character(levels(s))
y2 <- melt(freq_p)
colnames(y2)=c("Treatment","Cluster","Frequency")
y2$Treatment="P"
y=rbind(y,y2)
freq_c=freq_c[2:nrow(freq_c),]
colnames(freq_c)=as.character(levels(s))
y2 <- melt(freq_c)
colnames(y2)=c("Treatment","Cluster","Frequency")
y2$Treatment="C"
y=rbind(y,y2)
freq_pc=freq_pc[2:nrow(freq_pc),]
colnames(freq_pc)=as.character(levels(s))
y2 <- melt(freq_pc)
colnames(y2)=c("Treatment","Cluster","Frequency")
y2$Treatment="PC"
y=rbind(y,y2)
y$Cluster=as.factor(y$Cluster)
y$Treatment=factor(y$Treatment, levels=c("PBS","P","C","PC"))
p=ggplot(y, aes(x=Cluster, y=Frequency,fill=Treatment)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("New Clones Frequency by Cluster")
print(p + theme(
  plot.title = element_text(size=14, face="bold",hjust = 0.5),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold")))
####





###clonal frequency UMAP plots for shared clones
vp.meta.seurat.CD4$cdr3b_nt_freq=0
x=as.data.frame(table(vp.meta.seurat.CD4$cdr3b_nt))
x[,2]=x[,2]/sum(x[,2])
rownames(x)=x[,1]
for(i in 1:ncol(vp.meta.seurat.CD4)){
  vp.meta.seurat.CD4$cdr3b_nt_freq[i]=x[vp.meta.seurat.CD4$cdr3b_nt[i],2]
}
FeaturePlot(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$shared==T)],"cdr3b_nt_freq",cols=c("white","red"),reduction="lda")

vp.meta.seurat.CD8$cdr3b_nt_freq=0
x=as.data.frame(table(vp.meta.seurat.CD8$cdr3b_nt))
x[,2]=x[,2]/sum(x[,2])
rownames(x)=x[,1]
for(i in 1:ncol(vp.meta.seurat.CD8)){
  vp.meta.seurat.CD8$cdr3b_nt_freq[i]=x[vp.meta.seurat.CD8$cdr3b_nt[i],2]
}
FeaturePlot(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$shared==T)],"cdr3b_nt_freq",cols=c("white","red"),reduction="lda")

###tcrb clonality per cluster
cloneCal <- function(array) {
  x = array[array >0 ] / sum(array)
  #  x = sort(array, decreasing=T)
  l = length(x)
  entropy = sum(x * -1 * log2(x))
  maxentropy = -log2(1/l)
  return(signif(1 - entropy / maxentropy, 2))
}

abundancePlot<-function(data, threshold=0, name="abundance plot",color_list=c()){
  par(mfrow=c(1,1)) #prepares space for plots, fill by row
  plot(0,type='n', main=name, xlab="log freq", ylab="log abundance", xlim=c(log10(min(data[data>0])),log10(max(data))), ylim=c(0,5))
  color=1
  if(length(color_list)>0){
    c=1
    color=color_list[1]
  }
  if(is.null(ncol(data))){
    out=as.data.frame(table(data[data>threshold]))
    points(log10(as.numeric(as.character(out[,1]))), log10(out[,2]),pch=16, col=color)
  }
  else{
    for(index in 1:length(data))
    {
      out=as.data.frame(table(data[data[,index]>threshold,index]))
      points(log10(as.numeric(as.character(out[,1]))), log10(out[,2]),pch=16, col=color)
      x=log10(as.numeric(as.character(out[,1])))
      y=log10(out[,2])
      abline(lm(y~x),col=color)
      if(length(color_list)>0){
        c=c+1
        color=color_list[c]}
      else{color=color+1}
    }
    if(length(color_list)>0){legend("topright", legend=colnames(data), pch=16,col=color_list)}
    else{legend("topright", legend=colnames(data), pch=16,col=1:length(data))}
  }
}

normalize <- function(data) {
  nc = ncol(data)
  for (i in 1:nc) {
    data[,i] = data[,i ] / sum(data[,i])
  }
  return(data)
}

cloneCal(table(vp.meta.seurat.CD4$cdr3b_nt[which(vp.meta.seurat.CD4$seurat_clusters == "0")]))
cloneCal(table(vp.meta.seurat.CD8$cdr3b_nt[which(vp.meta.seurat.CD8$seurat_clusters == "0")]))

x=normalize(as.data.frame.matrix(table(vp.meta.seurat.CD4$cdr3b_nt,vp.meta.seurat.CD4$seurat_clusters)))
abundancePlot(x)
x=normalize(as.data.frame.matrix(table(vp.meta.seurat.CD8$cdr3b_nt,vp.meta.seurat.CD8$seurat_clusters)))
abundancePlot(x)



####Merged Frequency Plots over time"
vp.meta.seurat.CD4$TX1_condition2=unlist(lapply(strsplit(as.character(vp.meta.seurat.CD4$ao_IMP_TX),"_"),function(z){z[[1]]}))
vp.meta.seurat.CD4$TX1_condition2=factor(vp.meta.seurat.CD4$TX1_condition2,levels=c("PRE","PBS","P","C ","PC","LN","TO","UNIDENTIFIED"))
plot(DimPlot(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$TX1_condition2 %in% c("C ","P","PBS","PC","PRE") & vp.meta.seurat.CD4$TX2_time=="_10TO15")], reduction = "lda",label = TRUE,split.by="TX1_condition2") + NoLegend())
plot(DimPlot(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$TX1_condition2 %in% c("C ","P","PBS","PC","PRE") & vp.meta.seurat.CD4$TX2_time=="_11TO14")], reduction = "lda",label = TRUE,split.by="TX1_condition2") + NoLegend())
plot(DimPlot(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$TX1_condition2 %in% c("C ","P","PBS","PC","PRE") & vp.meta.seurat.CD4$TX2_time=="_11TO17")], reduction = "lda",label = TRUE,split.by="TX1_condition2") + NoLegend())
s=vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$TX1_condition2 %in% c("C ","P","PBS","PC","PRE") & vp.meta.seurat.CD4$TX2_time=="_10TO15")]
s=vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$TX1_condition2 %in% c("C ","P","PBS","PC","PRE") & vp.meta.seurat.CD4$TX2_time=="_11TO14")]
s=vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$TX1_condition2 %in% c("C ","P","PBS","PC","PRE") & vp.meta.seurat.CD4$TX2_time=="_11TO17")]
g=as.data.frame.matrix(table(s$TX1_condition2,s$seurat_clusters))
g=g[rowSums(g)>0,]
g=apply(g,1,function(x){x/sum(x)})
y <- melt(g)
colnames(y)=c("cluster","sample","frequency")
y$sample=factor(y$sample, levels=c("PRE","PBS","P","C ","PC"))
y$cluster=as.factor(y$cluster)
p=ggplot(y, aes(x=cluster, y=frequency,group=sample,fill=cluster)) +
  geom_bar(stat="identity",position=position_dodge(),color="black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Cluster Frequency by Treatment Response")
p + theme(
  plot.title = element_text(size=14, face="bold"),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))

plot(DimPlot(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$TO2_type=="ETO")], reduction = "lda",label = TRUE,split.by="TO1_time") + NoLegend())
plot(DimPlot(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$TO2_type=="LTO")], reduction = "lda",label = TRUE,split.by="TO1_time") + NoLegend())
plot(DimPlot(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$TO2_type=="VLTO")], reduction = "lda",label = TRUE,split.by="TO1_time") + NoLegend())
s=vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$TO2_type=="ETO")]
s=vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$TO2_type=="LTO")]
s=vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$TO2_type=="VLTO")]
g=as.data.frame.matrix(table(s$TO1_time,s$seurat_clusters))
g=g[rowSums(g)>0,]
g=apply(g,1,function(x){x/sum(x)})
y <- melt(g)
colnames(y)=c("cluster","sample","frequency")
y$cluster=as.factor(y$cluster)
p=ggplot(y, aes(x=cluster, y=frequency,group=sample,fill=cluster)) +
  geom_bar(stat="identity",position=position_dodge(),color="black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Cluster Frequency by Treatment Response")
p + theme(
  plot.title = element_text(size=14, face="bold"),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))

vp.meta.seurat.CD8$TX1_condition2=unlist(lapply(strsplit(as.character(vp.meta.seurat.CD8$ao_IMP_TX),"_"),function(z){z[[1]]}))
vp.meta.seurat.CD8$TX1_condition2=factor(vp.meta.seurat.CD8$TX1_condition2,levels=c("PRE","PBS","P","C ","PC","LN","TO","UNIDENTIFIED"))
plot(DimPlot(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TX1_condition2 %in% c("C ","P","PBS","PC","PRE") & vp.meta.seurat.CD8$TX2_time=="_10TO15")], reduction = "lda",label = TRUE,split.by="TX1_condition2") + NoLegend())
plot(DimPlot(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TX1_condition2 %in% c("C ","P","PBS","PC","PRE") & vp.meta.seurat.CD8$TX2_time=="_11TO14")], reduction = "lda",label = TRUE,split.by="TX1_condition2") + NoLegend())
plot(DimPlot(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TX1_condition2 %in% c("C ","P","PBS","PC","PRE") & vp.meta.seurat.CD8$TX2_time=="_11TO17")], reduction = "lda",label = TRUE,split.by="TX1_condition2") + NoLegend())
s=vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TX1_condition2 %in% c("C ","P","PBS","PC","PRE") & vp.meta.seurat.CD8$TX2_time=="_10TO15")]
s=vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TX1_condition2 %in% c("C ","P","PBS","PC","PRE") & vp.meta.seurat.CD8$TX2_time=="_11TO14")]
s=vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TX1_condition2 %in% c("C ","P","PBS","PC","PRE") & vp.meta.seurat.CD8$TX2_time=="_11TO17")]
g=as.data.frame.matrix(table(s$TX1_condition2,s$seurat_clusters))
g=g[rowSums(g)>0,]
g=apply(g,1,function(x){x/sum(x)})
y <- melt(g)
colnames(y)=c("cluster","sample","frequency")
y$sample=factor(y$sample, levels=c("PRE","PBS","P","C ","PC"))
y$cluster=as.factor(y$cluster)
p=ggplot(y, aes(x=cluster, y=frequency,group=sample,fill=cluster)) +
  geom_bar(stat="identity",position=position_dodge(),color="black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Cluster Frequency by Treatment Response")
p + theme(
  plot.title = element_text(size=14, face="bold"),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))

plot(DimPlot(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TO2_type=="ETO")], reduction = "lda",label = TRUE,split.by="TO1_time") + NoLegend())
plot(DimPlot(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TO2_type=="LTO")], reduction = "lda",label = TRUE,split.by="TO1_time") + NoLegend())
plot(DimPlot(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TO2_type=="VLTO")], reduction = "lda",label = TRUE,split.by="TO1_time") + NoLegend())
s=vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TO2_type=="ETO")]
s=vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TO2_type=="LTO")]
s=vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TO2_type=="VLTO")]
g=as.data.frame.matrix(table(s$TO1_time,s$seurat_clusters))
g=g[rowSums(g)>0,]
g=apply(g,1,function(x){x/sum(x)})
y <- melt(g)
colnames(y)=c("cluster","sample","frequency")
y$cluster=as.factor(y$cluster)
p=ggplot(y, aes(x=cluster, y=frequency,group=sample,fill=cluster)) +
  geom_bar(stat="identity",position=position_dodge(),color="black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Cluster Frequency by Treatment Response")
p + theme(
  plot.title = element_text(size=14, face="bold"),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))


####subsampling frequency plots
s=vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TX1_condition2 %in% c("C ","P","PBS","PC","PRE"))]
s=vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$TX1_condition2 %in% c("C ","P","PBS","PC","PRE"))]
g=as.data.frame.matrix(table(s$TX1_condition2,s$seurat_clusters))
###
s1=vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TX1_condition2 %in% c("C ","P","PBS","PC","PRE"))]
s2=vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$TX1_condition2 %in% c("C ","P","PBS","PC","PRE"))]
s1=vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TX1_condition2 %in% c("C ","P","PBS","PC","PRE") & vp.meta.seurat.CD8$experiment=="ltx")]
s2=vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$TX1_condition2 %in% c("C ","P","PBS","PC","PRE") & vp.meta.seurat.CD4$experiment=="ltx")]
g1=as.data.frame.matrix(table(s1$TX1_condition2,s1$seurat_clusters))
g2=as.data.frame.matrix(table(s2$TX1_condition2,s2$seurat_clusters))
colnames(g1)=paste("CD8",colnames(g1),sep=".")
colnames(g2)=paste("CD4",colnames(g2),sep=".")
g=cbind(g1,g2)
###
g=g[rowSums(g)>0,]
n=min(rowSums(g))
pre=as.data.frame(matrix(rep(0,100*ncol(g)),nrow = 100,ncol=ncol(g)))
pbs=as.data.frame(matrix(rep(0,100*ncol(g)),nrow = 100,ncol=ncol(g)))
p=as.data.frame(matrix(rep(0,100*ncol(g)),nrow = 100,ncol=ncol(g)))
c=as.data.frame(matrix(rep(0,100*ncol(g)),nrow = 100,ncol=ncol(g)))
pc=as.data.frame(matrix(rep(0,100*ncol(g)),nrow = 100,ncol=ncol(g)))
set.seed(1234)
for(i in 1:100){
  g2=rbind(table(sample(colnames(g),n,replace = T,prob = g[1,])),
  table(sample(colnames(g),n,replace = T,prob = g[2,])),
  table(sample(colnames(g),n,replace = T,prob = g[3,])),
  table(sample(colnames(g),n,replace = T,prob = g[4,])),
  table(sample(colnames(g),n,replace = T,prob = g[5,])))
  g2=apply(g2,1,function(x){x/sum(x)})
  pre[i,]=g2[,1]
  pbs[i,]=g2[,2]
  p[i,]=g2[,3]
  c[i,]=g2[,4]
  pc[i,]=g2[,5]
}
colnames(pre)=rownames(g2)
colnames(pbs)=rownames(g2)
colnames(p)=rownames(g2)
colnames(c)=rownames(g2)
colnames(pc)=rownames(g2)
pre$sample="PRE"
pbs$sample="PBS"
p$sample="P"
c$sample="C"
pc$sample="PC"
y=melt(pre)
y2=melt(pbs)
y3=melt(p)
y4=melt(c)
y5=melt(pc)
y=rbind(y,y2,y3,y4,y5)
colnames(y)=c("Sample","Cluster","Frequency")
y$Sample=factor(y$Sample, levels=c("PRE","PBS","P","C","PC"))
y$Cluster=as.factor(y$Cluster)
p=ggplot(y, aes(x=Cluster, y=Frequency,fill=Sample)) +
  #geom_bar(stat="identity",position=position_dodge(),color="black") +
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Cluster Frequency by Treatment Response") + theme(plot.title = element_text(hjust = 0.5))
p + theme(
  plot.title = element_text(size=14, face="bold"),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))+geom_vline(xintercept = 7.5)+geom_vline(xintercept = 1.5)+geom_vline(xintercept = 2.5)+geom_vline(xintercept = 3.5)+geom_vline(xintercept = 4.5)+geom_vline(xintercept = 5.5)+geom_vline(xintercept = 6.5)+geom_vline(xintercept =8.5)+geom_vline(xintercept =9.5)+geom_vline(xintercept =10.5)+geom_vline(xintercept =11.5)+scale_fill_manual(values=c("yellow",brewer.pal(n=4,name="Dark2")))

s=vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$TO2_type=="ETO")]
s=vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TO2_type=="ETO")]
g=as.data.frame.matrix(table(s$TO1_time,s$seurat_clusters))
###
s1=vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TO2_type=="ETO")]
s2=vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$TO2_type=="ETO")]
g1=as.data.frame.matrix(table(s1$TO1_time,s1$seurat_clusters))
g2=as.data.frame.matrix(table(s2$TO1_time,s2$seurat_clusters))
g2=g2[rowSums(g2)>0,]
colnames(g1)=paste("CD8",colnames(g1),sep=".")
colnames(g2)=paste("CD4",colnames(g2),sep=".")
g=cbind(g1,g2)
###
g=g[rowSums(g)>0,]
n=min(rowSums(g))
d10.11=as.data.frame(matrix(rep(0,100*ncol(g)),nrow = 100,ncol=ncol(g)))
d13.14=as.data.frame(matrix(rep(0,100*ncol(g)),nrow = 100,ncol=ncol(g)))
d16.17=as.data.frame(matrix(rep(0,100*ncol(g)),nrow = 100,ncol=ncol(g)))
d19.20=as.data.frame(matrix(rep(0,100*ncol(g)),nrow = 100,ncol=ncol(g)))
d22.23=as.data.frame(matrix(rep(0,100*ncol(g)),nrow = 100,ncol=ncol(g)))
set.seed(1234)
for(i in 1:100){
  g2=rbind(table(sample(colnames(g),n,replace = T,prob = g[1,])),
           table(sample(colnames(g),n,replace = T,prob = g[2,])),
           table(sample(colnames(g),n,replace = T,prob = g[3,])),
           table(sample(colnames(g),n,replace = T,prob = g[4,])),
           table(sample(colnames(g),n,replace = T,prob = g[5,])))
  #for late turnover
  #g2=rbind(table(sample(colnames(g),n,replace = T,prob = g[1,])),
  #         table(sample(colnames(g),n,replace = T,prob = g[2,])),
  #         table(sample(colnames(g),n,replace = T,prob = g[3,])),
  #         c(0,table(sample(colnames(g),n,replace = T,prob = g[4,]))),
  #         table(sample(colnames(g),n,replace = T,prob = g[5,])))
  g2=apply(g2,1,function(x){x/sum(x)})
  d10.11[i,]=g2[,1]
  d13.14[i,]=g2[,2]
  d16.17[i,]=g2[,3]
  d19.20[i,]=g2[,4]
  d22.23[i,]=g2[,5]
}
colnames(d10.11)=colnames(g)
colnames(d13.14)=colnames(g)
colnames(d16.17)=colnames(g)
colnames(d19.20)=colnames(g)
colnames(d22.23)=colnames(g)
d10.11$sample="DAY10-11"
d13.14$sample="DAY13-14"
d16.17$sample="DAY16-17"
d19.20$sample="DAY19-20"
d22.23$sample="DAY22-23"
####
d10.11$sample="DAY16-17"
d13.14$sample="DAY19-20"
d16.17$sample="DAY22-23"
d19.20$sample="DAY25-26"
d22.23$sample="DAY28-29"
###
y=melt(d10.11)
y2=melt(d13.14)
y3=melt(d16.17)
y4=melt(d19.20)
y5=melt(d22.23)
y=rbind(y,y2,y3,y4,y5)
colnames(y)=c("Sample","Cluster","Frequency")
#y$Sample=factor(y$Sample, levels=c("PRE","PBS","P","C","PC"))
y$Cluster=as.factor(y$Cluster)
p=ggplot(y, aes(x=Cluster, y=Frequency,fill=Sample)) +
  #geom_bar(stat="identity",position=position_dodge(),color="black") +
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Cluster Frequency by Treatment Response") + theme(plot.title = element_text(hjust = 0.5))
p + theme(
  plot.title = element_text(size=14, face="bold"),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))+geom_vline(xintercept = 7.5)+geom_vline(xintercept = 1.5)+geom_vline(xintercept = 2.5)+geom_vline(xintercept = 3.5)+geom_vline(xintercept = 4.5)+geom_vline(xintercept = 5.5)+geom_vline(xintercept = 6.5)+geom_vline(xintercept =8.5)+geom_vline(xintercept =9.5)+geom_vline(xintercept =10.5)+geom_vline(xintercept =11.5)+scale_fill_manual(values=c(brewer.pal(n=5,name="Purples")))

###relative frequency change of shared clones pre vs post in the same cluster in pd1 vs ctla4 normalized against pbs
###quantify sharing between lymph node clones and central memory (cd8 cluster 2 second time point). is ctla4 > pd1/pbs
###quantify sharing between lymph node clones and central memory (cd4 cluster 3 second time point). is ctla4 > pd1/pbs

tcrfreqchange_bootstrapped=function(s,treatment=F){
  set.seed(1234)
  s=s[,which(s$timepoint!="TLN")]
  n=min(table(s$timepoint))
  s_main=s
  clust1=rep(NA,100)
  clust2=rep(NA,100)
  clust3=rep(NA,100)
  clust4=rep(NA,100)
  for(j in 1:100){
    print(j)
    s=s_main[,c(sample(which(s_main$timepoint==unique(s_main$timepoint)[1]),n,replace = F),sample(which(s_main$timepoint==unique(s_main$timepoint)[2]),n,replace = F))]
    x=as.data.frame.matrix(table(s$cdr3b[which(s$seurat_clusters==0)],s$timepoint[which(s$seurat_clusters==0)]))
    shared=rownames(x)[x[,1]>0 & x[,2]>0]
    clust1[j]=length(shared)/length(unique(s$cdr3b[s$seurat_clusters==0]))
    x=as.data.frame.matrix(table(s$cdr3b[which(s$seurat_clusters==1)],s$timepoint[which(s$seurat_clusters==1)]))
    shared=rownames(x)[x[,1]>0 & x[,2]>0]
    clust2[j]=length(shared)/length(unique(s$cdr3b[s$seurat_clusters==1]))
    x=as.data.frame.matrix(table(s$cdr3b[which(s$seurat_clusters==2)],s$timepoint[which(s$seurat_clusters==2)]))
    if(ncol(x)==2){shared=rownames(x)[x[,1]>0 & x[,2]>0]
    clust3[j]=length(shared)/length(unique(s$cdr3b[s$seurat_clusters==2]))}
    x=as.data.frame.matrix(table(s$cdr3b[which(s$seurat_clusters==3)],s$timepoint[which(s$seurat_clusters==3)]))
    shared=rownames(x)[x[,1]>0 & x[,2]>0]
    clust4[j]=length(shared)/length(unique(s$cdr3b[s$seurat_clusters==3]))
  }
  y=data.frame(c(clust1,clust2,clust3,clust4),c(rep("Cluster1",100),rep("Cluster2",100),rep("Cluster3",100),rep("Cluster4",100)))
  colnames(y)=c("PercentShared","Cluster")
  y=y[complete.cases(y),]
  p=ggplot(y, aes(x=Cluster, y=PercentShared,fill=Cluster)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    ggtitle("Percent of TCRs Shared Within-Cluster")
  print(p + theme(
    plot.title = element_text(size=14, face="bold",hjust = 0.5),
    axis.title = element_text(size=14, face="bold"),
    axis.text = element_text(size=8),
    legend.text = element_text(size=10),
    legend.title = element_text(size=12),
    strip.text.x = element_text(size = 12, face="bold")))
  out=list(clust1,clust2,clust3,clust4)
  return(out)
}

cd4_paired_pbs=tcrfreqchange_bootstrapped(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4$condition %in% c("PreTx","PBS") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T)
cd4_paired_p=tcrfreqchange_bootstrapped(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4$condition %in% c("PreTx","P") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T)
cd4_paired_c=tcrfreqchange_bootstrapped(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4$condition %in% c("PreTx","C") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T)
cd4_paired_pc=tcrfreqchange_bootstrapped(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4$condition %in% c("PreTx","PC") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T)

cd8_paired_pbs=tcrfreqchange_bootstrapped(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","PBS") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_paired_p=tcrfreqchange_bootstrapped(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","P") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_paired_c=tcrfreqchange_bootstrapped(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","C") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_paired_pc=tcrfreqchange_bootstrapped(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","PC") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)

y=data.frame(c(cd8_paired_pbs[[1]],cd8_paired_p[[1]],cd8_paired_c[[1]]),c(rep("PBS",100),rep("PD1",100),rep("CTLA4",100)))
colnames(y)=c("PercentShared","Treatment")
y$Treatment=factor(y$Treatment,levels=c("PBS","PD1","CTLA4"))
y=y[complete.cases(y),]
p=ggplot(y, aes(x=Treatment, y=PercentShared,fill=Treatment)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Cluster0: Percent of TCRs Shared Within-Cluster")
print(p + theme(
  plot.title = element_text(size=14, face="bold",hjust = 0.5),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold")))

y=data.frame(c(cd8_paired_pbs[[2]],cd8_paired_p[[2]],cd8_paired_c[[2]]),c(rep("PBS",100),rep("PD1",100),rep("CTLA4",100)))
colnames(y)=c("PercentShared","Treatment")
y$Treatment=factor(y$Treatment,levels=c("PBS","PD1","CTLA4"))
y=y[complete.cases(y),]
p=ggplot(y, aes(x=Treatment, y=PercentShared,fill=Treatment)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Cluster1: Percent of TCRs Shared Within-Cluster")
print(p + theme(
  plot.title = element_text(size=14, face="bold",hjust = 0.5),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold")))


tcrLymphNodeSharing_bootstrapped=function(s,treatment=F){
  set.seed(1234)
  n=min(table(s$timepoint))
  s_main=s
  clust1=rep(NA,100)
  clust2=rep(NA,100)
  clust3=rep(NA,100)
  clust4=rep(NA,100)
  for(j in 1:100){
    print(j)
    s=s_main[,c(sample(which(s_main$timepoint==unique(s_main$timepoint)[1]),n,replace = F),sample(which(s_main$timepoint==unique(s_main$timepoint)[2]),n,replace = F))]
    x=as.data.frame.matrix(table(s$cdr3b[which(s$seurat_clusters==0)],s$timepoint[which(s$seurat_clusters==0)]))
    shared=rownames(x)[x[,1]>0 & x[,2]>0]
    clust1[j]=length(shared)/length(unique(s$cdr3b[s$seurat_clusters==0]))
    x=as.data.frame.matrix(table(s$cdr3b[which(s$seurat_clusters==1)],s$timepoint[which(s$seurat_clusters==1)]))
    shared=rownames(x)[x[,1]>0 & x[,2]>0]
    clust2[j]=length(shared)/length(unique(s$cdr3b[s$seurat_clusters==1]))
    x=as.data.frame.matrix(table(s$cdr3b[which(s$seurat_clusters==2)],s$timepoint[which(s$seurat_clusters==2)]))
    if(ncol(x)==2){shared=rownames(x)[x[,1]>0 & x[,2]>0]
    clust3[j]=length(shared)/length(unique(s$cdr3b[s$seurat_clusters==2]))}
    x=as.data.frame.matrix(table(s$cdr3b[which(s$seurat_clusters==3)],s$timepoint[which(s$seurat_clusters==3)]))
    shared=rownames(x)[x[,1]>0 & x[,2]>0]
    clust4[j]=length(shared)/length(unique(s$cdr3b[s$seurat_clusters==3]))
    }
  y=data.frame(c(clust1,clust2,clust3,clust4),c(rep("Cluster0",100),rep("Cluster1",100),rep("Cluster2",100),rep("Cluster3",100)))
  colnames(y)=c("PercentShared","Cluster")
  y=y[complete.cases(y),]
  p=ggplot(y, aes(x=Cluster, y=PercentShared,fill=Cluster)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    ggtitle("Percent of TCRs Shared Between Tumor & Lymph Node")
  print(p + theme(
    plot.title = element_text(size=14, face="bold",hjust = 0.5),
    axis.title = element_text(size=14, face="bold"),
    axis.text = element_text(size=8),
    legend.text = element_text(size=10),
    legend.title = element_text(size=12),
    strip.text.x = element_text(size = 12, face="bold"))+scale_fill_manual(values=rep("gray",4)))
  out=list(clust1,clust2,clust3,clust4)
  return(out)
}

cd4_ln_pbs=tcrLymphNodeSharing_bootstrapped(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","TLN") & vp.meta.seurat.CD4$condition %in% c("PreTx","PBS") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T)
cd4_ln_p=tcrLymphNodeSharing_bootstrapped(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","TLN") & vp.meta.seurat.CD4$condition %in% c("PreTx","P") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T)
cd4_ln_c=tcrLymphNodeSharing_bootstrapped(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","TLN") & vp.meta.seurat.CD4$condition %in% c("PreTx","C") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T)
cd4_ln_pc=tcrLymphNodeSharing_bootstrapped(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","TLN") & vp.meta.seurat.CD4$condition %in% c("PreTx","PC") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T)

cd8_ln_pbs=tcrLymphNodeSharing_bootstrapped(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","TLN") & vp.meta.seurat.CD8$condition %in% c("PreTx","PBS") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_ln_p=tcrLymphNodeSharing_bootstrapped(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","TLN") & vp.meta.seurat.CD8$condition %in% c("PreTx","P") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_ln_c=tcrLymphNodeSharing_bootstrapped(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","TLN") & vp.meta.seurat.CD8$condition %in% c("PreTx","C") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_ln_pc=tcrLymphNodeSharing_bootstrapped(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","TLN") & vp.meta.seurat.CD8$condition %in% c("PreTx","PC") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)

y=data.frame(c(cd4_ln_pbs[[1]],cd4_ln_p[[1]],cd4_ln_c[[1]],cd4_ln_pc[[1]]),c(rep("PBS",100),rep("PD1",100),rep("CTLA4",100),rep("PD1+CTLA4",100)))
colnames(y)=c("PercentShared","Treatment")
y$Cluster="CD4_1.TregA"
y2=data.frame(c(cd4_ln_pbs[[3]],cd4_ln_p[[3]],cd4_ln_c[[3]],cd4_ln_pc[[3]]),c(rep("PBS",100),rep("PD1",100),rep("CTLA4",100),rep("PD1+CTLA4",100)))
colnames(y2)=c("PercentShared","Treatment")
y2$Cluster="CD4_2.TregB"
y3=data.frame(c(cd4_ln_pbs[[2]],cd4_ln_p[[2]],cd4_ln_c[[2]],cd4_ln_pc[[2]]),c(rep("PBS",100),rep("PD1",100),rep("CTLA4",100),rep("PD1+CTLA4",100)))
colnames(y3)=c("PercentShared","Treatment")
y3$Cluster="CD4_3.Th1"
y4=data.frame(c(cd4_ln_pbs[[4]],cd4_ln_p[[4]],cd4_ln_c[[4]],cd4_ln_pc[[4]]),c(rep("PBS",100),rep("PD1",100),rep("CTLA4",100),rep("PD1+CTLA4",100)))
colnames(y4)=c("PercentShared","Treatment")
y4$Cluster="CD4_4.Th1PM"
y=rbind(y,y2,y3,y4)
#y=data.frame(c((cd8_ln_pbs[[1]]+cd8_ln_pbs[[2]])/2,(cd8_ln_p[[1]]+cd8_ln_p[[2]])/2,(cd8_ln_c[[1]]+cd8_ln_c[[2]])/2,(cd8_ln_pc[[1]]+cd8_ln_pc[[2]])/2),c(rep("PBS",100),rep("PD1",100),rep("CTLA4",100),rep("PD1+CTLA4",100)))
colnames(y)=c("PercentShared","Treatment","Cluster")
y$Treatment=factor(y$Treatment,levels=c("PBS","PD1","CTLA4","PD1+CTLA4"))
y=y[complete.cases(y),]
p=ggplot(y, aes(x=Cluster, y=PercentShared,fill=Treatment)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("% TCRs Shared Between Tumor & Lymph Node")
print(p + theme(
  plot.title = element_text(size=14, face="bold",hjust = 0.5),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))+scale_fill_manual(values=c(brewer.pal(n=4,name="Dark2"))))

y=data.frame(c(cd8_ln_pbs[[2]],cd8_ln_p[[2]],cd8_ln_c[[2]],cd8_ln_pc[[2]]),c(rep("PBS",100),rep("PD1",100),rep("CTLA4",100),rep("PD1+CTLA4",100)))
colnames(y)=c("PercentShared","Treatment")
y$Treatment=factor(y$Treatment,levels=c("PBS","PD1","CTLA4","PD1+CTLA4"))
y=y[complete.cases(y),]
p=ggplot(y, aes(x=Treatment, y=PercentShared,fill=Treatment)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Cluster1: % TCRs Shared Between Tumor & Lymph Node")
print(p + theme(
  plot.title = element_text(size=14, face="bold",hjust = 0.5),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))+scale_fill_manual(values=c(brewer.pal(n=4,name="Dark2"))))


###half-life of transition from effector to exhausted by treatment for clones transitioning from effector to exhausted.....
### 0 is effector, 1 is exhausted
####
s=vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TX1_condition2 %in% c("PRE","PC"))]
x=table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"PRE",1]/sum(table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"PRE",])
y=table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"POST",1]/sum(table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"POST",])
a=table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"POST",2]/sum(table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"POST",])
b=table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"PRE",2]/sum(table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"PRE",])
eff_to_exh=names(which(x>0 & a>0))
out=c()
for(i in eff_to_exh){
  out=c(out,a[i]/x[i])
}
mean(out)
halflife_pc=1/out*7
s=vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TX1_condition2 %in% c("PRE","P"))]
x=table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"PRE",1]/sum(table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"PRE",])
y=table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"POST",1]/sum(table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"POST",])
a=table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"POST",2]/sum(table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"POST",])
b=table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"PRE",2]/sum(table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"PRE",])
eff_to_exh=names(which(x>0 & a>0))
out=c()
for(i in eff_to_exh){
  out=c(out,a[i]/x[i])
}
mean(out)
halflife_pd1=1/out*7
s=vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TX1_condition2 %in% c("PRE","C "))]
x=table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"PRE",1]/sum(table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"PRE",])
y=table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"POST",1]/sum(table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"POST",])
a=table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"POST",2]/sum(table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"POST",])
b=table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"PRE",2]/sum(table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"PRE",])
eff_to_exh=names(which(x>0 & a>0))
out=c()
for(i in eff_to_exh){
  out=c(out,a[i]/x[i])
}
mean(out)
halflife_ctla4=1/out*7
s=vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TX1_condition2 %in% c("PRE","PBS"))]
x=table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"PRE",1]/sum(table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"PRE",])
y=table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"POST",1]/sum(table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"POST",])
a=table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"POST",2]/sum(table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"POST",])
b=table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"PRE",2]/sum(table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"PRE",])
eff_to_exh=names(which(x>0 & a>0))
out=c()
for(i in eff_to_exh){
  out=c(out,a[i]/x[i])
}
mean(out)
halflife_pbs=1/out*7
halflife_pbs=halflife_pbs[!halflife_pbs %in% boxplot.stats(halflife_pbs)$out]
halflife_pd1=halflife_pd1[!halflife_pd1 %in% boxplot.stats(halflife_pd1)$out]
halflife_ctla4=halflife_ctla4[!halflife_ctla4 %in% boxplot.stats(halflife_ctla4)$out]
halflife_pc=halflife_pc[!halflife_pc %in% boxplot.stats(halflife_pc)$out]

y=data.frame(Halflife=c(halflife_pbs,halflife_pd1,halflife_ctla4,halflife_pc),Treatment=c(rep("PBS",length(halflife_pbs)),rep("PD1",length(halflife_pd1)),rep("CTLA4",length(halflife_ctla4)),rep("PD1+CTLA4",length(halflife_pc))))
y$Treatment=factor(y$Treatment,levels=c("PBS","PD1","CTLA4","PD1+CTLA4"))
p=ggplot(y, aes(x=Treatment, y=Halflife,fill=Treatment)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("CD8 Effector To Exhausted Transition Half-Life")
print(p + theme(
  plot.title = element_text(size=14, face="bold",hjust = 0.5),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold")))

###half-life of transition from CD4 0 to CD4 2 and CD4 1 to CD4 3
####
s=vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$TX1_condition2 %in% c("PRE","PC"))]
x=table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"PRE",1]/sum(table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"PRE",])
y=table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"POST",1]/sum(table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"POST",])
a=table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"POST",3]/sum(table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"POST",])
b=table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"PRE",3]/sum(table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"PRE",])
eff_to_exh=names(which(x>0 & a>0))
out=c()
for(i in eff_to_exh){
  out=c(out,a[i]/x[i])
}
mean(out)
halflife_pc=1/out*7
s=vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$TX1_condition2 %in% c("PRE","P"))]
x=table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"PRE",1]/sum(table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"PRE",])
y=table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"POST",1]/sum(table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"POST",])
a=table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"POST",3]/sum(table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"POST",])
b=table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"PRE",3]/sum(table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"PRE",])
eff_to_exh=names(which(x>0 & a>0))
out=c()
for(i in eff_to_exh){
  out=c(out,a[i]/x[i])
}
mean(out)
halflife_pd1=1/out*7
s=vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$TX1_condition2 %in% c("PRE","C "))]
x=table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"PRE",1]/sum(table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"PRE",])
y=table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"POST",1]/sum(table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"POST",])
a=table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"POST",3]/sum(table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"POST",])
b=table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"PRE",3]/sum(table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"PRE",])
eff_to_exh=names(which(x>0 & a>0))
out=c()
for(i in eff_to_exh){
  out=c(out,a[i]/x[i])
}
mean(out)
halflife_ctla4=1/out*7
s=vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$TX1_condition2 %in% c("PRE","PBS"))]
x=table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"PRE",1]/sum(table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"PRE",])
y=table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"POST",1]/sum(table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"POST",])
a=table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"POST",3]/sum(table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"POST",])
b=table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"PRE",3]/sum(table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"PRE",])
eff_to_exh=names(which(x>0 & a>0))
out=c()
for(i in eff_to_exh){
  out=c(out,a[i]/x[i])
}
mean(out)
halflife_pbs=1/out*7
halflife_pbs=halflife_pbs[!halflife_pbs %in% boxplot.stats(halflife_pbs)$out]
halflife_pd1=halflife_pd1[!halflife_pd1 %in% boxplot.stats(halflife_pd1)$out]
halflife_ctla4=halflife_ctla4[!halflife_ctla4 %in% boxplot.stats(halflife_ctla4)$out]
halflife_pc=halflife_pc[!halflife_pc %in% boxplot.stats(halflife_pc)$out]

y=data.frame(Halflife=c(halflife_pbs,halflife_pd1,halflife_ctla4,halflife_pc),Treatment=c(rep("PBS",length(halflife_pbs)),rep("PD1",length(halflife_pd1)),rep("CTLA4",length(halflife_ctla4)),rep("PD1+CTLA4",length(halflife_pc))))
y$Treatment=factor(y$Treatment,levels=c("PBS","PD1","CTLA4","PD1+CTLA4"))
p=ggplot(y, aes(x=Treatment, y=Halflife,fill=Treatment)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("CD4: 0->2 Transition Half-Life")
print(p + theme(
  plot.title = element_text(size=14, face="bold",hjust = 0.5),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold")))


###make downsampled UMAP plots
s=vp.meta.seurat.CD4
n=min(table(s$TX1_condition2)[c("P","C ","PC","PBS","PRE")])
i=c(sample(which(s$TX1_condition2=="PRE"),n),sample(which(s$TX1_condition2=="PBS"),n),sample(which(s$TX1_condition2=="P"),n),sample(which(s$TX1_condition2=="C "),n),sample(which(s$TX1_condition2=="PC"),n))
s=s[,i]
plot(DimPlot(s[,which(s$TX1_condition2 %in% c("C ","P","PBS","PC","PRE"))], reduction = "lda",label = TRUE,split.by="TX1_condition2") + NoLegend())
s=vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$TO2_type=="ETO")]
n=min(table(s$TO1_time)[c("DAY10-11","DAY13-14","DAY16-17","DAY19-20","DAY22-23")])
i=c(sample(which(s$TO1_time=="DAY10-11"),n),sample(which(s$TO1_time=="DAY13-14"),n),sample(which(s$TO1_time=="DAY16-17"),n),sample(which(s$TO1_time=="DAY19-20"),n),sample(which(s$TO1_time=="DAY22-23"),n))
s=s[,i]
plot(DimPlot(s, reduction = "lda",label = TRUE,split.by="TO1_time") + NoLegend())
###
s=vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$TO2_type=="LTO")]
n=min(table(s$TO1_time)[c("DAY16-17","DAY19-20","DAY22-23","DAY25-26","DAY28-29")])
i=c(sample(which(s$TO1_time=="DAY16-17"),n),sample(which(s$TO1_time=="DAY19-20"),n),sample(which(s$TO1_time=="DAY22-23"),n),sample(which(s$TO1_time=="DAY25-26"),n),sample(which(s$TO1_time=="DAY28-29"),n))
s=s[,i]
plot(DimPlot(s, reduction = "lda",label = TRUE,split.by="TO1_time") + NoLegend())
###

s=vp.meta.seurat.CD8
n=min(table(s$TX1_condition2)[c("P","C ","PC","PBS","PRE")])
i=c(sample(which(s$TX1_condition2=="PRE"),n),sample(which(s$TX1_condition2=="PBS"),n),sample(which(s$TX1_condition2=="P"),n),sample(which(s$TX1_condition2=="C "),n),sample(which(s$TX1_condition2=="PC"),n))
s=s[,i]
plot(DimPlot(s[,which(s$TX1_condition2 %in% c("C ","P","PBS","PC","PRE"))], reduction = "lda",label = TRUE,split.by="TX1_condition2") + NoLegend())
s=vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TO2_type=="ETO")]
n=min(table(s$TO1_time)[c("DAY10-11","DAY13-14","DAY16-17","DAY19-20","DAY22-23")])
i=c(sample(which(s$TO1_time=="DAY10-11"),n),sample(which(s$TO1_time=="DAY13-14"),n),sample(which(s$TO1_time=="DAY16-17"),n),sample(which(s$TO1_time=="DAY19-20"),n),sample(which(s$TO1_time=="DAY22-23"),n))
s=s[,i]
plot(DimPlot(s, reduction = "lda",label = TRUE,split.by="TO1_time") + NoLegend())
###
s=vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TO2_type=="LTO")]
n=min(table(s$TO1_time)[c("DAY16-17","DAY19-20","DAY22-23","DAY25-26","DAY28-29")])
i=c(sample(which(s$TO1_time=="DAY16-17"),n),sample(which(s$TO1_time=="DAY19-20"),n),sample(which(s$TO1_time=="DAY22-23"),n),sample(which(s$TO1_time=="DAY25-26"),n),sample(which(s$TO1_time=="DAY28-29"),n))
s=s[,i]
plot(DimPlot(s, reduction = "lda",label = TRUE,split.by="TO1_time") + NoLegend())
###



###
#Tumor Lymph Node
###
cd4.tcr.in.tumor=unique(vp.meta.seurat.CD4$cdr3b[which(vp.meta.seurat.CD4$ao_IMP_TX!="LN")])
cd8.tcr.in.tumor=unique(vp.meta.seurat.CD8$cdr3b[which(vp.meta.seurat.CD8$ao_IMP_TX!="LN")])
vp.meta.seurat.CD4$TX1_condition=factor(vp.meta.seurat.CD4$TX1_condition, levels=c("TO","PRE","PBS","P","C ","PC","PBS:TLN","P:LN","C:LN","PC:LN"))
vp.meta.seurat.CD8$TX1_condition=factor(vp.meta.seurat.CD8$TX1_condition, levels=c("TO","PRE","PBS","P","C ","PC","PBS:TLN","P:LN","C:LN","PC:LN"))
DimPlot(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$ao_IMP_TX=="LN")],reduction="lda")
DimPlot(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$ao_IMP_TX=="LN")],reduction="lda")
x=vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$ao_IMP_TX=="LN")]
y=vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$ao_IMP_TX!="LN")]
length(which(x$cdr3b_nt %in% y$cdr3b_nt & x$condition=="PBS"))/length(which(x$condition=="PBS"))
length(which(x$cdr3b_nt %in% y$cdr3b_nt & x$condition=="PC"))/length(which(x$condition=="PC"))
length(which(x$cdr3b_nt %in% y$cdr3b_nt & x$condition=="P"))/length(which(x$condition=="P"))
length(which(x$cdr3b_nt %in% y$cdr3b_nt & x$condition=="C"))/length(which(x$condition=="C"))
x$condition=factor(x$condition,levels=c("PBS","C","P","PC"))
table(x$experiment)
DimPlot(x,split.by="condition",reduction = "lda")
n=min(table(x$condition)[c("P","C","PC","PBS")])
i=c(sample(which(x$condition=="PBS"),n),sample(which(x$condition=="P"),n),sample(which(x$condition=="C"),n),sample(which(x$condition=="PC"),n))
s=x[,i]
s$condition=factor(s$condition,levels=c("PBS","P","C","PC"))
plot(DimPlot(s[,which(s$condition %in% c("PBS","P","C","PC"))], reduction = "lda",label = TRUE,split.by="condition",cols= c("#117733","#CC6677","#6699CC","#DDCC77")) + NoLegend())
y=as.data.frame.matrix(table(x$condition,x$seurat_clusters))
y=y[,which(colSums(y)>0)]
g=apply(y,1,function(x){x/sum(x)})
y <- melt(g)
colnames(y)=c("cluster","treatment","frequency")
y$treatment=factor(y$treatment, levels=c("PBS","C","P","PC"))
y$cluster=as.factor(y$cluster)
p=ggplot(y, aes(x=cluster, y=frequency,group=treatment,fill=cluster)) +
  geom_bar(stat="identity",position=position_dodge(),color="black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Cluster Frequency by Treatment")
p + theme(
  plot.title = element_text(size=14, face="bold"),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))
x=vp.meta.seurat.CD4.2[,which(vp.meta.seurat.CD4.2$ao_IMP_TX=="LN")]
y=vp.meta.seurat.CD4.2[,which(vp.meta.seurat.CD4.2$ao_IMP_TX!="LN")]
length(which(x$cdr3b_nt %in% y$cdr3b_nt & x$condition=="PBS"))/length(which(x$condition=="PBS"))
length(which(x$cdr3b_nt %in% y$cdr3b_nt & x$condition=="PC"))/length(which(x$condition=="PC"))
length(which(x$cdr3b_nt %in% y$cdr3b_nt & x$condition=="P"))/length(which(x$condition=="P"))
length(which(x$cdr3b_nt %in% y$cdr3b_nt & x$condition=="C"))/length(which(x$condition=="C"))
x$condition=factor(x$condition,levels=c("PBS","C","P","PC"))
table(x$experiment)
DimPlot(x,split.by="condition",reduction = "lda")
n=min(table(x$condition)[c("P","C","PC","PBS")])
i=c(sample(which(x$condition=="PBS"),n),sample(which(x$condition=="P"),n),sample(which(x$condition=="C"),n),sample(which(x$condition=="PC"),n))
s=x[,i]
s$condition=factor(s$condition,levels=c("PBS","P","C","PC"))
plot(DimPlot(s[,which(s$condition %in% c("PBS","P","C","PC"))], reduction = "lda",label = F,split.by="condition",cols=c("#332288", "#6699CC", "#117733", "#CC6677","#DDCC77","#44AA99")) + NoLegend())
y=as.data.frame.matrix(table(x$condition,x$seurat_clusters))
y=y[,which(colSums(y)>0)]
g=apply(y,1,function(x){x/sum(x)})
y <- melt(g)
colnames(y)=c("cluster","treatment","frequency")
y$treatment=factor(y$treatment, levels=c("PBS","C","P","PC"))
y$cluster=as.factor(y$cluster)
p=ggplot(y, aes(x=cluster, y=frequency,group=treatment,fill=cluster)) +
  geom_bar(stat="identity",position=position_dodge(),color="black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Cluster Frequency by Treatment")
p + theme(
  plot.title = element_text(size=14, face="bold"),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))
####subsampling frequency plots
x=vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$ao_IMP_TX=="LN")]
x=vp.meta.seurat.CD4.2[,which(vp.meta.seurat.CD4.2$ao_IMP_TX=="LN")]
s=x[,which(x$condition %in% c("C","P","PBS","PC","PRE"))]
s$condition=factor(s$condition,levels=c("PBS","P","C","PC"))
g=as.data.frame.matrix(table(s$condition,s$seurat_clusters))
###
g=g[rowSums(g)>0,]
n=min(rowSums(g))
pbs=as.data.frame(matrix(rep(0,100*ncol(g)),nrow = 100,ncol=ncol(g)))
p=as.data.frame(matrix(rep(0,100*ncol(g)),nrow = 100,ncol=ncol(g)))
c=as.data.frame(matrix(rep(0,100*ncol(g)),nrow = 100,ncol=ncol(g)))
pc=as.data.frame(matrix(rep(0,100*ncol(g)),nrow = 100,ncol=ncol(g)))
set.seed(1234)
for(i in 1:100){
  g2=rbind(table(sample(colnames(g),n,replace = T,prob = g[1,])),
           table(sample(colnames(g),n,replace = T,prob = g[2,])),
           table(sample(colnames(g),n,replace = T,prob = g[3,])),
           table(sample(colnames(g),n,replace = T,prob = g[4,])))
  g2=apply(g2,1,function(x){x/sum(x)})
  pbs[i,]=g2[,1]
  p[i,]=g2[,2]
  c[i,]=g2[,3]
  pc[i,]=g2[,4]
}
colnames(pbs)=rownames(g2)
colnames(p)=rownames(g2)
colnames(c)=rownames(g2)
colnames(pc)=rownames(g2)
pbs$sample="PBS"
p$sample="P"
c$sample="C"
pc$sample="PC"
y2=melt(pbs)
y3=melt(p)
y4=melt(c)
y5=melt(pc)
y=rbind(y2,y3,y4,y5)
colnames(y)=c("Sample","Cluster","Frequency")
y$Sample=factor(y$Sample, levels=c("PBS","P","C","PC"))
y$Cluster=as.factor(y$Cluster)
p=ggplot(y, aes(x=Cluster, y=Frequency,fill=Sample)) +
  #geom_bar(stat="identity",position=position_dodge(),color="black") +
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Cluster Frequency by Treatment Response") + theme(plot.title = element_text(hjust = 0.5))
p + theme(
  plot.title = element_text(size=14, face="bold"),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))+geom_vline(xintercept = 7.5)+geom_vline(xintercept = 1.5)+geom_vline(xintercept = 2.5)+geom_vline(xintercept = 3.5)+geom_vline(xintercept = 4.5)+geom_vline(xintercept = 5.5)+geom_vline(xintercept = 6.5)+geom_vline(xintercept =8.5)+geom_vline(xintercept =9.5)+geom_vline(xintercept =10.5)+geom_vline(xintercept =11.5)+scale_fill_manual(values=c(brewer.pal(n=4,name="Dark2")))






###backtrace where cd8 cluster 0 (effector) clones from tumor come from in lymph node, by treatment
s=vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TX1_condition %in% c("PBS","PBS:TLN"))]
s=vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TX1_condition %in% c("P","P:LN"))]
s=vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TX1_condition %in% c("C ","C:LN"))]
cdr3_seqs=names(which(table(s$cdr3b,s$TX1_condition)[,1]>0 & table(s$cdr3b,s$TX1_condition)[,2]>0))
x=s[,which(s$cdr3b %in% cdr3_seqs)]
cdr3_seqs=names(which(apply(as.data.frame.matrix(table(x$cdr3b,x$seurat_clusters,x$TX1_condition)[,,1]),1,function(x){x[1]==max(x)})))
x=s[,which(s$cdr3b %in% cdr3_seqs)]
y=colSums(as.data.frame.matrix(table(x$cdr3b,x$seurat_clusters,x$TX1_condition)[,,2]))
y/sum(y)
chisq.test(cbind(c(36,9,3,0),c(13,2,2,1),c(32,13,7,0)))


#differential activity within-cluster by treatment group
s=vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$seurat_clusters==0 & vp.meta.seurat.CD8$TX1_condition2 %in% c("PBS","P","C ","PC"))]
s$TX1_condition2=droplevels(s$TX1_condition2)
Idents(s)="TX1_condition2"
markers.viper <- FindAllMarkers(s, only.pos = TRUE, min.pct = 0, logfc.threshold = 0,test.use = "t")
top10 <- markers.viper %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
geneHeatmap_plot(s@assays$RNA@data,s$TX1_condition2,top10$gene,genes_by_cluster=F,scaled=T)
s2=arnoldhan_data.integrated[,colnames(s)]
s2$TX1_condition2=s$TX1_condition2
Idents(s2)="TX1_condition2"
markers.gene <- FindAllMarkers(s2, only.pos = TRUE, min.pct = 0, logfc.threshold = 0,test.use = "MAST")
top10 <- markers.gene %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
geneHeatmap_plot(s2@assays$RNA@data,s2$TX1_condition2,top10$gene,genes_by_cluster=F,scaled=F)
###
x=markers.viper[markers.viper$cluster %in% c("PBS","P"),]
x$avg_logFC[which(x$cluster=="PBS")]= -x$avg_logFC[which(x$cluster=="PBS")]
plot(x$avg_logFC,-log10(x$p_val_adj),pch=16,cex=.6,xlab="avg_logFC",ylab="-log10(p)",main="VIPER: PBS vs P")
x=markers.viper[markers.viper$cluster %in% c("PBS","C"),]
x$avg_logFC[which(x$cluster=="PBS")]= -x$avg_logFC[which(x$cluster=="PBS")]
plot(x$avg_logFC,-log10(x$p_val_adj),pch=16,cex=.6,xlab="avg_logFC",ylab="-log10(p)",main="VIPER: PBS vs C")
x=markers.viper[markers.viper$cluster %in% c("PBS","PC"),]
x$avg_logFC[which(x$cluster=="PBS")]= -x$avg_logFC[which(x$cluster=="PBS")]
plot(x$avg_logFC,-log10(x$p_val_adj),pch=16,cex=.6,xlab="avg_logFC",ylab="-log10(p)",main="VIPER: PBS vs PC")
###
x=markers.gene[markers.gene$cluster %in% c("PBS","P"),]
x$avg_logFC[which(x$cluster=="PBS")]= -x$avg_logFC[which(x$cluster=="PBS")]
plot(x$avg_logFC,-log10(x$p_val_adj),pch=16,cex=.6,xlab="avg_logFC",ylab="-log10(p)",main="GeneExp: PBS vs P")
x=markers.gene[markers.gene$cluster %in% c("PBS","C"),]
x$avg_logFC[which(x$cluster=="PBS")]= -x$avg_logFC[which(x$cluster=="PBS")]
plot(x$avg_logFC,-log10(x$p_val_adj),pch=16,cex=.6,xlab="avg_logFC",ylab="-log10(p)",main="GeneExp: PBS vs C")
x=markers.gene[markers.gene$cluster %in% c("PBS","PC"),]
x$avg_logFC[which(x$cluster=="PBS")]= -x$avg_logFC[which(x$cluster=="PBS")]
plot(x$avg_logFC,-log10(x$p_val_adj),pch=16,cex=.6,xlab="avg_logFC",ylab="-log10(p)",main="GeneExp: PBS vs PC")


#saveRDS(markers,"~/Dropbox/Turnover_Treatment_Paper/data_and_code_for_samhita/markers.cd8.cluster0.viper")
s=vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$seurat_clusters==1 & vp.meta.seurat.CD8$TX1_condition2 %in% c("PBS","P","C ","PC"))]
s=vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$seurat_clusters==2 & vp.meta.seurat.CD8$TX1_condition2 %in% c("PBS","P","C ","PC"))]
s=vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$seurat_clusters==3 & vp.meta.seurat.CD8$TX1_condition2 %in% c("PBS","P","C ","PC"))]
s=vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$seurat_clusters==0 & vp.meta.seurat.CD4$TX1_condition2 %in% c("PBS","P","C ","PC"))]
s=vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$seurat_clusters==1 & vp.meta.seurat.CD4$TX1_condition2 %in% c("PBS","P","C ","PC"))]
s=vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$seurat_clusters==2 & vp.meta.seurat.CD4$TX1_condition2 %in% c("PBS","P","C ","PC"))]
s=vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$seurat_clusters==3 & vp.meta.seurat.CD4$TX1_condition2 %in% c("PBS","P","C ","PC"))]
s=vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$seurat_clusters==4 & vp.meta.seurat.CD4$TX1_condition2 %in% c("PBS","P","C ","PC"))]
s=vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$seurat_clusters==6 & vp.meta.seurat.CD4$TX1_condition2 %in% c("PBS","P","C ","PC"))]
s=vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$seurat_clusters==8 & vp.meta.seurat.CD4$TX1_condition2 %in% c("PBS","P","C ","PC"))]


###bootstrapped clonality comparison within-clusters
pbs=c()
for(i in 1:100){pbs=c(pbs,cloneCal(sample(table(vp.meta.seurat.CD8$cdr3b[which(vp.meta.seurat.CD8$seurat_clusters == "0" & vp.meta.seurat.CD8$TX1_condition2=="PBS")]),100)))}
pd1=c()
for(i in 1:100){pd1=c(pd1,cloneCal(sample(table(vp.meta.seurat.CD8$cdr3b[which(vp.meta.seurat.CD8$seurat_clusters == "0" & vp.meta.seurat.CD8$TX1_condition2=="P")]),100)))}
ctla4=c()
for(i in 1:100){ctla4=c(ctla4,cloneCal(sample(table(vp.meta.seurat.CD8$cdr3b[which(vp.meta.seurat.CD8$seurat_clusters == "0" & vp.meta.seurat.CD8$TX1_condition2=="C ")]),100)))}
mean(pbs)
mean(pd1)
mean(ctla4)
wilcox.test(pbs,pd1)
wilcox.test(pbs,ctla4)
wilcox.test(pd1,ctla4)
y=data.frame(Clonality=c(pbs,pd1,ctla4),Treatment=c(rep("PBS",length(pbs)),rep("PD1",length(pd1)),rep("CTLA4",length(ctla4))))
y$Treatment=factor(y$Treatment,levels=c("PBS","PD1","CTLA4"))
p=ggplot(y, aes(x=Treatment, y=Clonality,fill=Treatment)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("CD8 Cluster 0 Clonality by Treatment Group")
print(p + theme(
  plot.title = element_text(size=14, face="bold",hjust = 0.5),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold")))
###
pre=c()
for(i in 1:100){pre=c(pre,cloneCal(sample(table(vp.meta.seurat.CD8$cdr3b[which(vp.meta.seurat.CD8$seurat_clusters == "0" & vp.meta.seurat.CD8$TX1_condition2=="PRE")]),100)))}
pbs=c()
for(i in 1:100){pbs=c(pbs,cloneCal(sample(table(vp.meta.seurat.CD8$cdr3b[which(vp.meta.seurat.CD8$seurat_clusters == "0" & vp.meta.seurat.CD8$TX1_condition2=="PBS")]),100)))}
pd1=c()
for(i in 1:100){pd1=c(pd1,cloneCal(sample(table(vp.meta.seurat.CD8$cdr3b[which(vp.meta.seurat.CD8$seurat_clusters == "0" & vp.meta.seurat.CD8$TX1_condition2=="P")]),100)))}
ctla4=c()
for(i in 1:100){ctla4=c(ctla4,cloneCal(sample(table(vp.meta.seurat.CD8$cdr3b[which(vp.meta.seurat.CD8$seurat_clusters == "0" & vp.meta.seurat.CD8$TX1_condition2=="C ")]),100)))}
y=data.frame(Clonality=c(pre,pbs,pd1,ctla4),Treatment=c(rep("PRE",length(pbs)),rep("PBS",length(pbs)),rep("PD1",length(pd1)),rep("CTLA4",length(ctla4))))
y$Treatment=factor(y$Treatment,levels=c("PRE","PBS","PD1","CTLA4"))
p=ggplot(y, aes(x=Treatment, y=Clonality,fill=Treatment)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("CD8 Cluster 0 Clonality by Treatment Group")
print(p + theme(
  plot.title = element_text(size=14, face="bold",hjust = 0.5),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold")))
###
table(vp.meta.seurat.CD8$TO1_time)
d10.11=c()
for(i in 1:100){d10.11=c(d10.11,cloneCal(sample(table(vp.meta.seurat.CD8$cdr3b[which(vp.meta.seurat.CD8$seurat_clusters == "0" & vp.meta.seurat.CD8$TO1_time=="DAY10-11")]),100,replace = T)))}
d13.14=c()
for(i in 1:100){d13.14=c(d13.14,cloneCal(sample(table(vp.meta.seurat.CD8$cdr3b[which(vp.meta.seurat.CD8$seurat_clusters == "0" & vp.meta.seurat.CD8$TO1_time=="DAY13-14")]),100,replace = T)))}
d16.17=c()
for(i in 1:100){d16.17=c(d16.17,cloneCal(sample(table(vp.meta.seurat.CD8$cdr3b[which(vp.meta.seurat.CD8$seurat_clusters == "0" & vp.meta.seurat.CD8$TO1_time=="DAY16-17")]),100,replace = T)))}
d19.20=c()
for(i in 1:100){d19.20=c(d19.20,cloneCal(sample(table(vp.meta.seurat.CD8$cdr3b[which(vp.meta.seurat.CD8$seurat_clusters == "0" & vp.meta.seurat.CD8$TO1_time=="DAY19-20")]),100,replace = T)))}
d22.23=c()
for(i in 1:100){d22.23=c(d22.23,cloneCal(sample(table(vp.meta.seurat.CD8$cdr3b[which(vp.meta.seurat.CD8$seurat_clusters == "0" & vp.meta.seurat.CD8$TO1_time=="DAY22-23")]),100,replace = T)))}
d25.26=c()
for(i in 1:100){d25.26=c(d25.26,cloneCal(sample(table(vp.meta.seurat.CD8$cdr3b[which(vp.meta.seurat.CD8$seurat_clusters == "0" & vp.meta.seurat.CD8$TO1_time=="DAY25-26")]),100,replace = T)))}
d28.29=c()
for(i in 1:100){d28.29=c(d28.29,cloneCal(sample(table(vp.meta.seurat.CD8$cdr3b[which(vp.meta.seurat.CD8$seurat_clusters == "0" & vp.meta.seurat.CD8$TO1_time=="DAY28-29")]),100,replace = T)))}
y=data.frame(Clonality=c(d10.11,d13.14,d16.17,d19.20,d22.23,d25.26,d28.29),Timepoint=c(rep("DAY10-11",length(d10.11)),rep("DAY13-14",length(d13.14)),rep("DAY16-17",length(d16.17)),rep("DAY19-20",length(d19.20)),rep("DAY22-23",length(d22.23)),rep("DAY25-26",length(d25.26)),rep("DAY28-29",length(d28.29))))
y$Cluster="Cluster0"
d10.11=c()
for(i in 1:100){d10.11=c(d10.11,cloneCal(sample(table(vp.meta.seurat.CD8$cdr3b[which(vp.meta.seurat.CD8$seurat_clusters == "1" & vp.meta.seurat.CD8$TO1_time=="DAY10-11")]),100,replace = T)))}
d13.14=c()
for(i in 1:100){d13.14=c(d13.14,cloneCal(sample(table(vp.meta.seurat.CD8$cdr3b[which(vp.meta.seurat.CD8$seurat_clusters == "1" & vp.meta.seurat.CD8$TO1_time=="DAY13-14")]),100,replace = T)))}
d16.17=c()
for(i in 1:100){d16.17=c(d16.17,cloneCal(sample(table(vp.meta.seurat.CD8$cdr3b[which(vp.meta.seurat.CD8$seurat_clusters == "1" & vp.meta.seurat.CD8$TO1_time=="DAY16-17")]),100,replace = T)))}
d19.20=c()
for(i in 1:100){d19.20=c(d19.20,cloneCal(sample(table(vp.meta.seurat.CD8$cdr3b[which(vp.meta.seurat.CD8$seurat_clusters == "1" & vp.meta.seurat.CD8$TO1_time=="DAY19-20")]),100,replace = T)))}
d22.23=c()
for(i in 1:100){d22.23=c(d22.23,cloneCal(sample(table(vp.meta.seurat.CD8$cdr3b[which(vp.meta.seurat.CD8$seurat_clusters == "1" & vp.meta.seurat.CD8$TO1_time=="DAY22-23")]),100,replace = T)))}
d25.26=c()
for(i in 1:100){d25.26=c(d25.26,cloneCal(sample(table(vp.meta.seurat.CD8$cdr3b[which(vp.meta.seurat.CD8$seurat_clusters == "1" & vp.meta.seurat.CD8$TO1_time=="DAY25-26")]),100,replace = T)))}
d28.29=c()
for(i in 1:100){d28.29=c(d28.29,cloneCal(sample(table(vp.meta.seurat.CD8$cdr3b[which(vp.meta.seurat.CD8$seurat_clusters == "1" & vp.meta.seurat.CD8$TO1_time=="DAY28-29")]),100,replace = T)))}
y2=data.frame(Clonality=c(d10.11,d13.14,d16.17,d19.20,d22.23,d25.26,d28.29),Timepoint=c(rep("DAY10-11",length(d10.11)),rep("DAY13-14",length(d13.14)),rep("DAY16-17",length(d16.17)),rep("DAY19-20",length(d19.20)),rep("DAY22-23",length(d22.23)),rep("DAY25-26",length(d25.26)),rep("DAY28-29",length(d28.29))))
y2$Cluster="Cluster1"
y=rbind(y,y2)
d10.11=c()
for(i in 1:100){d10.11=c(d10.11,cloneCal(sample(table(vp.meta.seurat.CD8$cdr3b[which(vp.meta.seurat.CD8$TO1_time=="DAY10-11")]),100,replace = T)))}
d13.14=c()
for(i in 1:100){d13.14=c(d13.14,cloneCal(sample(table(vp.meta.seurat.CD8$cdr3b[which(vp.meta.seurat.CD8$TO1_time=="DAY13-14")]),100,replace = T)))}
d16.17=c()
for(i in 1:100){d16.17=c(d16.17,cloneCal(sample(table(vp.meta.seurat.CD8$cdr3b[which(vp.meta.seurat.CD8$TO1_time=="DAY16-17")]),100,replace = T)))}
d19.20=c()
for(i in 1:100){d19.20=c(d19.20,cloneCal(sample(table(vp.meta.seurat.CD8$cdr3b[which(vp.meta.seurat.CD8$TO1_time=="DAY19-20")]),100,replace = T)))}
d22.23=c()
for(i in 1:100){d22.23=c(d22.23,cloneCal(sample(table(vp.meta.seurat.CD8$cdr3b[which(vp.meta.seurat.CD8$TO1_time=="DAY22-23")]),100,replace = T)))}
d25.26=c()
for(i in 1:100){d25.26=c(d25.26,cloneCal(sample(table(vp.meta.seurat.CD8$cdr3b[which(vp.meta.seurat.CD8$TO1_time=="DAY25-26")]),100,replace = T)))}
d28.29=c()
for(i in 1:100){d28.29=c(d28.29,cloneCal(sample(table(vp.meta.seurat.CD8$cdr3b[which(vp.meta.seurat.CD8$TO1_time=="DAY28-29")]),100,replace = T)))}
y3=data.frame(Clonality=c(d10.11,d13.14,d16.17,d19.20,d22.23,d25.26,d28.29),Timepoint=c(rep("DAY10-11",length(d10.11)),rep("DAY13-14",length(d13.14)),rep("DAY16-17",length(d16.17)),rep("DAY19-20",length(d19.20)),rep("DAY22-23",length(d22.23)),rep("DAY25-26",length(d25.26)),rep("DAY28-29",length(d28.29))))
y3$Cluster="Total"
y=rbind(y,y3)
#y$Timepoint=factor(y$Timepoint,levels=c("PRE","PBS","PD1","CTLA4"))
p=ggplot(y, aes(x=Timepoint, y=Clonality,fill=Cluster)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("CD8 Clonality by Timepoint")
print(p + theme(
  plot.title = element_text(size=14, face="bold",hjust = 0.5),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold")))
#+scale_fill_manual(values=rep("gray",4))+
###
y=data.frame(cluster=c(0,0,0),group=c(0,0,0),clonality=c(0,0,0))
set.seed(1234)
for(g in c("PRE","PBS","P","C ","PC")){
  print(g)
  total=c()
  for(i in 1:100){total=c(total,cloneCal(sample(table(vp.meta.seurat.CD8$cdr3b[which(vp.meta.seurat.CD8$TX1_condition2==g)]),100)))}
  clust0=c()
  for(i in 1:100){clust0=c(clust0,cloneCal(sample(table(vp.meta.seurat.CD8$cdr3b[which(vp.meta.seurat.CD8$seurat_clusters == "0" & vp.meta.seurat.CD8$TX1_condition2==g)]),100,replace=T)))}
  clust1=c()
  for(i in 1:100){clust1=c(clust1,cloneCal(sample(table(vp.meta.seurat.CD8$cdr3b[which(vp.meta.seurat.CD8$seurat_clusters == "1" & vp.meta.seurat.CD8$TX1_condition2==g)]),100,replace=T)))}
  clust2=c()
  for(i in 1:100){clust2=c(clust2,cloneCal(sample(table(vp.meta.seurat.CD8$cdr3b[which(vp.meta.seurat.CD8$seurat_clusters == "2" & vp.meta.seurat.CD8$TX1_condition2==g)]),100,replace=T)))}
  y=rbind(y,data.frame(cluster=c(rep("Total",length(total)),rep("Cluster0",length(clust0)),rep("Cluster1",length(clust1)),rep("Cluster2",length(clust2))),group=rep(g,length(total)*4),clonality=c(total,clust0,clust1,clust2)))
}
y=y[4:nrow(y),]
y$group=factor(y$group,levels=c("PRE","PBS","P","C ","PC"))
y$cluster=factor(y$cluster,levels=c("Total","Cluster0","Cluster1","Cluster2"))
p=ggplot(y, aes(x=group, y=clonality,fill=cluster)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("CD8 Sub-Cluster Clonalities by Treatment Group")
print(p + theme(
  plot.title = element_text(size=14, face="bold",hjust = 0.5),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold")))
###
#CD4
###
####
y=data.frame(cluster=c(0,0,0,0),group=c(0,0,0,0),clonality=c(0,0,0,0))
set.seed(1234)
for(g in c("PRE","PBS","P","C ","PC")){
  print(g)
  total=c()
  for(i in 1:100){total=c(total,cloneCal(sample(table(vp.meta.seurat.CD4$cdr3b[which(vp.meta.seurat.CD4$TX1_condition2==g)]),100)))}
  clust0=c()
  for(i in 1:100){clust0=c(clust0,cloneCal(sample(table(vp.meta.seurat.CD4$cdr3b[which(vp.meta.seurat.CD4$seurat_clusters == "0" & vp.meta.seurat.CD4$TX1_condition2==g)]),100,replace=T)))}
  clust1=c()
  for(i in 1:100){clust1=c(clust1,cloneCal(sample(table(vp.meta.seurat.CD4$cdr3b[which(vp.meta.seurat.CD4$seurat_clusters == "1" & vp.meta.seurat.CD4$TX1_condition2==g)]),100,replace=T)))}
  clust2=c()
  for(i in 1:100){clust2=c(clust2,cloneCal(sample(table(vp.meta.seurat.CD4$cdr3b[which(vp.meta.seurat.CD4$seurat_clusters == "2" & vp.meta.seurat.CD4$TX1_condition2==g)]),100,replace=T)))}
  clust3=c()
  for(i in 1:100){clust3=c(clust3,cloneCal(sample(table(vp.meta.seurat.CD4$cdr3b[which(vp.meta.seurat.CD4$seurat_clusters == "3" & vp.meta.seurat.CD4$TX1_condition2==g)]),100,replace=T)))}
  y=rbind(y,data.frame(cluster=c(rep("Total",length(total)),rep("Cluster0",length(clust0)),rep("Cluster1",length(clust1)),rep("Cluster2",length(clust2)),rep("Cluster3",length(clust3))),group=rep(g,length(total)*5),clonality=c(total,clust0,clust1,clust2,clust3)))
}
y=y[5:nrow(y),]
y$group=factor(y$group,levels=c("PRE","PBS","P","C ","PC"))
y$cluster=factor(y$cluster,levels=c("Total","Cluster0","Cluster1","Cluster2","Cluster3"))
p=ggplot(y, aes(x=group, y=clonality,fill=cluster)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("CD4 Sub-Cluster Clonalities by Treatment Group")
print(p + theme(
  plot.title = element_text(size=14, face="bold",hjust = 0.5),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold")))
###


###### clonality but only for clones present in cluster0 at baseline & at second timepoint
x=as.data.frame.matrix(table(vp.meta.seurat.CD4$cdr3b[which(vp.meta.seurat.CD4$seurat_clusters==0)],vp.meta.seurat.CD4$timepoint_simplified[which(vp.meta.seurat.CD4$seurat_clusters==0)]))
tcrs=rownames(x)[x[,2]>x[,1] & x[,1]>0]
y=data.frame(cluster=c(0,0,0,0),group=c(0,0,0,0),clonality=c(0,0,0,0))
set.seed(1234)
for(g in c("PRE","PBS","P","C ","PC")){
  print(g)
  clust0=c()
  for(i in 1:100){clust0=c(clust0,cloneCal(table(sample(names(table(vp.meta.seurat.CD4$cdr3b[which(vp.meta.seurat.CD4$TX1_condition2==g & vp.meta.seurat.CD4$cdr3b %in% tcrs)])),100,replace=T,prob = table(vp.meta.seurat.CD4$cdr3b[which(vp.meta.seurat.CD4$TX1_condition2==g & vp.meta.seurat.CD4$cdr3b %in% tcrs)])))))}
  y=rbind(y,data.frame(cluster=c(rep("Cluster0",length(clust0))),group=rep(g,length(clust0)),clonality=c(clust0)))
}
y=y[5:nrow(y),]
y$group[y$group=="PBS"]="temp"
y$group[y$group=="P"]="PBS"
y$group[y$group=="temp"]="P"
y$group[y$group=="C "]="temp"
y$group[y$group=="PC"]="C "
y$group[y$group=="temp"]="PC"
y$group=factor(y$group,levels=c("PRE","PBS","P","C ","PC"))
y$cluster=factor(y$cluster,levels=c("Cluster0"))
p=ggplot(y, aes(x=group, y=clonality,fill=cluster)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("CD4 Cluster 0 Clonalities by Treatment Group")
print(p + theme(
  plot.title = element_text(size=14, face="bold",hjust = 0.5),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold")))+NoLegend()
###



#####arnold new version of half-life..... rate of change only within CD8 cluster 0, separated by timepoint....
#assuming by day 19-20 there is no influx and only outflux from cluster 0.... and rate of outflux is constant...
table(vp.meta.seurat.CD8$timepoint,vp.meta.seurat.CD8$TO1_time)

s=vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TO1_time=="DAY19-20")]
x=table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"TUM1",1]/sum(table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"TUM1",])
y=table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"TUM2",1]/sum(table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"TUM2",])
moving=which(x>0)
out=c()
for(i in moving){
  out=c(out,(x[i]-(x[i]/2))/((x[i]-y[i])/10))
}
mean(out)
halflife=out
halflife=halflife[!halflife %in% boxplot.stats(halflife)$out]

s=vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TO1_time=="DAY16-17")]
x=table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"TUM1",1]/sum(table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"TUM1",])
y=table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"TUM2",1]/sum(table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"TUM2",])
moving=which(x>0)
out=c()
for(i in moving){
  if(x[i]>y[i]){out=c(out,(x[i]-(x[i]/2))/((x[i]-y[i])/6))}
  else{out=c(out,-(x[i]*2-x[i])/((y[i]-x[i])/6))}
}
halflife=5-out[out!=3]
halflife=halflife[!halflife %in% boxplot.stats(halflife)$out]
mean(halflife)
hist(halflife)

expected=mean(halflife)

s=vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TX1_condition2 %in% c("PRE","PC"))]
x=table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"PRE",1]/sum(table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"PRE",])
y=table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"POST",1]/sum(table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"POST",])
moving=which(x>0)
out=c()
for(i in moving){
  if(x[i]>y[i]){out=c(out,(x[i]-(x[i]/2))/((x[i]-y[i])/7))}
  else{out=c(out,-(x[i]*2-x[i])/((y[i]-x[i])/7))}
}
halflife_pc=5-out[out!=3.5]
halflife_pc=halflife_pc[!halflife_pc %in% boxplot.stats(halflife_pc)$out]
s=vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TX1_condition2 %in% c("PRE","P"))]
x=table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"PRE",1]/sum(table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"PRE",])
y=table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"POST",1]/sum(table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"POST",])
moving=which(x>0)
out=c()
for(i in moving){
  if(x[i]>y[i]){out=c(out,(x[i]-(x[i]/2))/((x[i]-y[i])/7))}
  else{out=c(out,-(x[i]*2-x[i])/((y[i]-x[i])/7))}
}
halflife_p=5-out[out!=3.5]
halflife_p=halflife_p[!halflife_p %in% boxplot.stats(halflife_p)$out]
s=vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TX1_condition2 %in% c("PRE","C "))]
x=table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"PRE",1]/sum(table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"PRE",])
y=table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"POST",1]/sum(table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"POST",])
moving=which(x>0)
out=c()
for(i in moving){
  if(x[i]>y[i]){out=c(out,(x[i]-(x[i]/2))/((x[i]-y[i])/7))}
  else{out=c(out,-(x[i]*2-x[i])/((y[i]-x[i])/7))}
}
halflife_c=5-out[out!=3.5]
halflife_c=halflife_c[!halflife_c %in% boxplot.stats(halflife_c)$out]
s=vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TX1_condition2 %in% c("PRE","PBS"))]
x=table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"PRE",1]/sum(table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"PRE",])
y=table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"POST",1]/sum(table(s$cdr3b,s$timepoint,s$seurat_clusters)[,"POST",])
moving=which(x>0)
out=c()
for(i in moving){
  if(x[i]>y[i]){out=c(out,(x[i]-(x[i]/2))/((x[i]-y[i])/7))}
  else{out=c(out,-(x[i]*2-x[i])/((y[i]-x[i])/7))}
}
halflife_pbs=5-out[out!=3.5]
halflife_pbs=halflife_pbs[!halflife_pbs %in% boxplot.stats(halflife_pbs)$out]

halflife_pbs=halflife_pbs-mean(halflife)
halflife_p=halflife_p-mean(halflife)
halflife_c=halflife_c-mean(halflife)
halflife_pc=halflife_pc-mean(halflife)

y=data.frame(Halflife=halflife,sample="Turnover Day 13-17")
p<-ggplot(y, aes(x=Halflife)) + 
  geom_histogram(color="black", fill="gray")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("CD8 Effector Day 13-17 Doubling Time")
print(p + theme(
  plot.title = element_text(size=14, face="bold",hjust = 0.5),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))+scale_fill_manual(values=rep("gray",4))+NoLegend())

y=data.frame(Halflife=c(halflife_pbs,halflife_p,halflife_c,halflife_pc),Treatment=c(rep("PBS",length(halflife_pbs)),rep("PD1",length(halflife_p)),rep("CTLA4",length(halflife_c)),rep("PD1+CTLA4",length(halflife_pc))))
y$Treatment=factor(y$Treatment,levels=c("PBS","PD1","CTLA4","PD1+CTLA4"))
p=ggplot(y, aes(x=Treatment, y=Halflife,fill=Treatment)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("CD8 Effector Half-Life")
print(p + theme(
  plot.title = element_text(size=14, face="bold",hjust = 0.5),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))+scale_fill_manual(values=rep("gray",4))+NoLegend())

wilcox.test(halflife,halflife_pbs)
wilcox.test(halflife,halflife_p)
wilcox.test(halflife,halflife_c)
wilcox.test(halflife,halflife_pc)

###there's still the problem that we can't distinguish expansion vs recruitment....


###
###
#MC38 CD8 day 3/5/7 treatment split plot of frequency and effector->exhausted transition
###
###
vp.meta.seurat.CD8$TX2_time_summarized=vp.meta.seurat.CD8$TX2_time
vp.meta.seurat.CD8$TX2_time_summarized[which(vp.meta.seurat.CD8$TX2_time=="_10TO15")]="Day5"
vp.meta.seurat.CD8$TX2_time_summarized[which(vp.meta.seurat.CD8$TX2_time=="_10TO17")]="Day7"
vp.meta.seurat.CD8$TX2_time_summarized[which(vp.meta.seurat.CD8$TX2_time=="_11TO14")]="Day3"
vp.meta.seurat.CD8$TX2_time_summarized[which(vp.meta.seurat.CD8$TX2_time=="_11TO17")]="Day7"
vp.meta.seurat.CD8$TX2_time_summarized[which(vp.meta.seurat.CD8$TX2_time=="_12TO19")]="Day7"
vp.meta.seurat.CD8$TX2_time_summarized[!(vp.meta.seurat.CD8$TX2_time_summarized %in% c("Day3","Day5","Day7"))]=NA
####subsampling frequency plots
s1=vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$TX1_condition2 %in% c("C ","P","PBS","PC","PRE") & vp.meta.seurat.CD8$TX2_time_summarized == "Day7")]
g1=as.data.frame.matrix(table(s1$TX1_condition2,s1$seurat_clusters))
colnames(g1)=paste("CD8",colnames(g1),sep=".")
g=g1
###
g=g[rowSums(g)>0,]
n=min(rowSums(g))
pre=as.data.frame(matrix(rep(0,100*ncol(g)),nrow = 100,ncol=ncol(g)))
pbs=as.data.frame(matrix(rep(0,100*ncol(g)),nrow = 100,ncol=ncol(g)))
p=as.data.frame(matrix(rep(0,100*ncol(g)),nrow = 100,ncol=ncol(g)))
c=as.data.frame(matrix(rep(0,100*ncol(g)),nrow = 100,ncol=ncol(g)))
pc=as.data.frame(matrix(rep(0,100*ncol(g)),nrow = 100,ncol=ncol(g)))
set.seed(1234)
for(i in 1:100){
  g2=rbind(table(sample(colnames(g),n,replace = T,prob = g[1,])),
           table(sample(colnames(g),n,replace = T,prob = g[2,])),
           table(sample(colnames(g),n,replace = T,prob = g[3,])),
           table(sample(colnames(g),n,replace = T,prob = g[4,])),
           table(sample(colnames(g),n,replace = T,prob = g[5,])))
  g2=apply(g2,1,function(x){x/sum(x)})
  pre[i,]=g2[,1]
  pbs[i,]=g2[,2]
  p[i,]=g2[,3]
  c[i,]=g2[,4]
  pc[i,]=g2[,5]
}
colnames(pre)=rownames(g2)
colnames(pbs)=rownames(g2)
colnames(p)=rownames(g2)
colnames(c)=rownames(g2)
colnames(pc)=rownames(g2)
pre$sample="PRE"
pbs$sample="PBS"
p$sample="P"
c$sample="C"
pc$sample="PC"
y=melt(pre)
y2=melt(pbs)
y3=melt(p)
y4=melt(c)
y5=melt(pc)
y=rbind(y,y2,y3,y4,y5)
colnames(y)=c("Sample","Cluster","Frequency")
y$Sample=factor(y$Sample, levels=c("PRE","PBS","P","C","PC"))
y$Cluster=as.factor(y$Cluster)
p=ggplot(y, aes(x=Cluster, y=Frequency,fill=Sample)) +
  #geom_bar(stat="identity",position=position_dodge(),color="black") +
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Cluster Frequency by Treatment Response: Day7") + theme(plot.title = element_text(hjust = 0.5))
p + theme(
  plot.title = element_text(size=14, face="bold"),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))+geom_vline(xintercept = 7.5)+geom_vline(xintercept = 1.5)+geom_vline(xintercept = 2.5)+geom_vline(xintercept = 3.5)+geom_vline(xintercept = 4.5)+geom_vline(xintercept = 5.5)+geom_vline(xintercept = 6.5)+geom_vline(xintercept =8.5)+geom_vline(xintercept =9.5)+geom_vline(xintercept =10.5)+geom_vline(xintercept =11.5)+scale_fill_manual(values=c("yellow",brewer.pal(n=4,name="Dark2")))
###
cd8_tcr_pbs=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","PBS") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx") & vp.meta.seurat.CD8$TX2_time_summarized=="Day7")],treatment = T)
cd8_tcr_p=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","P") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx") & vp.meta.seurat.CD8$TX2_time_summarized=="Day7")],treatment = T)
cd8_tcr_c=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","C") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx") & vp.meta.seurat.CD8$TX2_time_summarized=="Day7")],treatment = T)
cd8_tcr_pc=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PreTx","PC") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx") & vp.meta.seurat.CD8$TX2_time_summarized=="Day7")],treatment = T)
diff=jitter(cd8_tcr_p,amount = 1)
cd8_tcr_pbs=as.data.frame.matrix(cd8_tcr_pbs)
cd8_tcr_p=as.data.frame.matrix(cd8_tcr_p)
cd8_tcr_c=as.data.frame.matrix(cd8_tcr_c)
diff=as.data.frame.matrix(diff)
cd8_tcr_pbs$treatment="PBS"
cd8_tcr_p$treatment="PD1"
cd8_tcr_c$treatment="CTLA4"
diff$treatment="PD1+CTLA4"
y=rbind(cd8_tcr_pbs,cd8_tcr_p,cd8_tcr_c,diff)
y <- melt(y)
colnames(y)=c("Treatment","Transition","LogFoldChange")
y$Treatment=factor(y$Treatment, levels=c("PBS","PD1","CTLA4","PD1+CTLA4"))
p=ggplot(y, aes(x=Transition, y=LogFoldChange,fill=Treatment)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("TCR Frequency Log-Fold-Change Distribution by Cluster: Day7")
print(p + theme(
  plot.title = element_text(size=14, face="bold",hjust = 0.5),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))+ylim(-6,6)+scale_fill_manual(values=c(brewer.pal(n=4,name="Dark2"))))









































#####
#####
#####

#########CT26 ANALYSIS########

##CT26 analysis
arnoldhan_data_ct26=readRDS("~/Dropbox/Turnover_Treatment_Paper/OneDrive_1_8-16-2021/arnoldhan_data_ct26_olddata.rds")
newct26=readRDS("/Users/aleksandar/Dropbox/Turnover_Treatment_Paper/OneDrive_1_8-16-2021/exp2526_assay_defaultRNA.rds")
newct26=CreateSeuratObject(newct26)
#newct26=readRDS("/Users/aleksandar/Dropbox/Turnover_Treatment_Paper/OneDrive_1_8-16-2021/exp2526.Rds")
metadat=read.csv("/Users/aleksandar/Dropbox/Turnover_Treatment_Paper/OneDrive_1_8-16-2021/meta.data.exp2526.csv")
newct26$number=metadat$orig.ident
newct26$sample_label=metadat$sample_label
newct26$sample=metadat$sample
newct26$mouse=metadat$mouse
newct26$timepoint=metadat$timepoint
newct26$condition=metadat$treatment
newct26$number=metadat$number
newct26$experiment=metadat$experiment
newct26$cdr3a=metadat$cdr3a
newct26$crd3a_nt=metadat$cdr3a_nt
newct26$cdr3b=metadat$cdr3b
newct26$crd3b_nt=metadat$cdr3b_nt
newct26$CD4adj=metadat$CD4adj
newct26$CD8adj=metadat$CD8adj
newct26$sct_Cd4=metadat$sct_Cd4
newct26$sct_Cd8a=metadat$sct_Cd8a
newct26$sct_Cd8b=metadat$sct_Cd8b
newct26$CD4CD8=metadat$cluster_gates
newct26=newct26[,which(newct26$sample_label=="label")]
combined_data_ct26=merge(arnoldhan_data_ct26,newct26)

table(combined_data_ct26$timepoint)
table(combined_data_ct26$experiment)
table(combined_data_ct26$condition)
table(combined_data_ct26$CD4CD8)
table(combined_data_ct26$number) ###batch for batch correction

combined_data_ct26=combined_data_ct26[,which(combined_data_ct26$CD4CD8 %in% c("CD4","CD8"))]
#combined_data_ct26=combined_data_ct26[,which(combined_data_ct26$number!="exp9")]
big_list=SplitObject(combined_data_ct26, split.by = "number")
for(i in 1:length(big_list)){
  x=big_list[[i]]
  DefaultAssay(x) <- "RNA"
  x <- SCTransform(x,vars.to.regress = c("nCount_RNA"),return.only.var.genes = F,verbose = T,conserve.memory = T)
  x.singler=CreateSinglerObject(x[["SCT"]]@counts, annot = NULL, 
                                project.name = "arnoldhan", min.genes = 0, 
                                technology = "10X", species = "Mouse", citation = "", 
                                do.signatures = F, clusters = NULL, numCores = numCores,
                                fine.tune=F,temp.dir = "/Users/aleksandar/Downloads",
                                variable.genes = "de",reduce.file.size = T,do.main.types = T)
  x$hpca_labels=x.singler$singler[[1]][[1]][[2]]
  x$hpca_main_labels=x.singler$singler[[1]][[4]][[2]]
  x$blueprint_labels=x.singler$singler[[2]][[1]][[2]]
  x$blueprint_main_labels=x.singler$singler[[2]][[4]][[2]]
  x$hpca_pvals=x.singler$singler[[1]][[1]][[3]]
  x$hpca_main_pvals=x.singler$singler[[1]][[4]][[3]]
  x$blueprint_pvals=x.singler$singler[[2]][[1]][[3]]
  x$blueprint_main_pvals=x.singler$singler[[2]][[4]][[3]]
  l=x$blueprint_labels
  l[which(x$blueprint_pvals>0.1)]=NA
  l[which(l %in% names(which(table(l)<200)))]=NA
  x$l=l
  big_list[[i]]=x
}
for (i in 1:length(big_list)) { DefaultAssay(big_list[[i]]) <- 'SCT' }
sr.features <- SelectIntegrationFeatures(object.list = big_list, nfeatures = 2000)
big_list <- PrepSCTIntegration(object.list = big_list, anchor.features = sr.features, verbose = T)
sr.anchors <- FindIntegrationAnchors(object.list = big_list, normalization.method = "SCT", anchor.features = sr.features, verbose = T) #reference=1
rm(sr.features,big_list)
arnoldhan_data.integrated <- IntegrateData(anchorset = sr.anchors, normalization.method = "SCT", verbose = T)
rm(sr.anchors)
DefaultAssay(arnoldhan_data.integrated)="integrated"
arnoldhan_data.integrated <- RunPCA(arnoldhan_data.integrated, features = VariableFeatures(object = arnoldhan_data.integrated),npcs = 30)
arnoldhan_data.integrated <- RunUMAP(arnoldhan_data.integrated, dims = 1:30, verbose = FALSE,umap.method="umap-learn",metric="correlation")
arnoldhan_data.integrated <- FindNeighbors(arnoldhan_data.integrated, dims = 1:30, verbose = FALSE)
arnoldhan_data.integrated <- FindClusters(arnoldhan_data.integrated, resolution=seq(0.01,1,by=0.01), verbose = FALSE,algorithm=1) # set resolution <1 for fewer clusters (default is 0.5)
clust=arnoldhan_data.integrated@meta.data[,which(grepl("integrated_snn_res.",colnames(arnoldhan_data.integrated@meta.data)))]
mat=as.data.frame(t(arnoldhan_data.integrated$pca@cell.embeddings))
out=sil_subsample(mat,clust)
means.2=out[[1]]
sd.2=out[[2]]
x=seq(0.01,1,by=0.01)
errbar(x,means.2,means.2+sd.2,means.2-sd.2,ylab="mean silhouette score",xlab="resolution parameter")
lines(x,means.2)
best=tail(x[which(means.2==max(means.2))],n=1)
print(best)
arnoldhan_data.integrated$seurat_clusters=arnoldhan_data.integrated@meta.data[,which(colnames(arnoldhan_data.integrated@meta.data)==paste("integrated_snn_res.",best,sep=""))]
Idents(arnoldhan_data.integrated) <- "seurat_clusters"
DimPlot(arnoldhan_data.integrated, reduction = "umap",label = TRUE) + NoLegend()
DimPlot(arnoldhan_data.integrated, reduction = "umap",label = TRUE,group.by = "l") + NoLegend()
saveRDS(arnoldhan_data.integrated,"~/Documents/Documents/MED SCHOOL/summer 2019/arnold_han_collaboration/combined_batch_corrected_arnoldhan_data.integrated.ct26.rds")


#filter out non-t-cells and re-cluster
arnoldhan_data.integrated=arnoldhan_data.integrated[,which(arnoldhan_data.integrated$seurat_clusters!="15" & arnoldhan_data.integrated$l != "NK cells")]
arnoldhan_data.integrated <- RunPCA(arnoldhan_data.integrated, features = VariableFeatures(object = arnoldhan_data.integrated),npcs = 30)
arnoldhan_data.integrated <- RunUMAP(arnoldhan_data.integrated, dims = 1:30, verbose = FALSE,umap.method="umap-learn",metric="correlation")
arnoldhan_data.integrated <- FindNeighbors(arnoldhan_data.integrated, dims = 1:30, verbose = FALSE)
arnoldhan_data.integrated <- FindClusters(arnoldhan_data.integrated, resolution=seq(0.01,1,by=0.01), verbose = FALSE,algorithm=1) # set resolution <1 for fewer clusters (default is 0.5)
clust=arnoldhan_data.integrated@meta.data[,which(grepl("integrated_snn_res.",colnames(arnoldhan_data.integrated@meta.data)))]
mat=as.data.frame(t(arnoldhan_data.integrated$pca@cell.embeddings))
out=sil_subsample(mat,clust)
means.2=out[[1]]
sd.2=out[[2]]
x=seq(0.01,1,by=0.01)
errbar(x,means.2,means.2+sd.2,means.2-sd.2,ylab="mean silhouette score",xlab="resolution parameter")
lines(x,means.2)
x=x[20:length(x)]
means.2=means.2[20:length(means.2)]
best=tail(x[which(means.2==max(means.2))],n=1)
print(best)
arnoldhan_data.integrated$seurat_clusters=arnoldhan_data.integrated@meta.data[,which(colnames(arnoldhan_data.integrated@meta.data)==paste("integrated_snn_res.",best,sep=""))]
Idents(arnoldhan_data.integrated) <- "seurat_clusters"
DimPlot(arnoldhan_data.integrated, reduction = "umap",label = TRUE) + NoLegend()
arnoldhan_data.markers <- FindAllMarkers(arnoldhan_data.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5,test.use = "MAST")
top10 <- arnoldhan_data.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
geneHeatmap_plot(arnoldhan_data.integrated@assays$SCT@counts,arnoldhan_data.integrated$seurat_clusters,top10$gene,genes_by_cluster = T,scaled = F,n_top_genes_per_cluster = 5)
saveRDS(arnoldhan_data.integrated,"~/Documents/Documents/MED SCHOOL/summer 2019/arnold_han_collaboration/combined_batch_corrected_arnoldhan_data.integrated.ct26.rds")

#rm(exp5,exp6,combined_gene_expression_data,clust,batch_metadata)
arnoldhan_data.integrated=readRDS("~/Documents/Documents/MED SCHOOL/summer 2019/arnold_han_collaboration/combined_batch_corrected_arnoldhan_data.integrated.ct26.rds")
##make metaCells
sr_merged_meta=MakeCMfA(dat.mat=as.matrix(arnoldhan_data.integrated[["SCT"]]@counts),clustering=arnoldhan_data.integrated$seurat_clusters,out.dir="/Users/aleksandar/Documents/Documents/MED SCHOOL/summer 2019/arnold_han_collaboration/",out.name="arnoldhan_alldata.ct26")
rm(arnoldhan_data.integrated)
meta=sr_merged_meta[[1]]
for(i in 2:length(sr_merged_meta)){
  meta=merge(meta,sr_merged_meta[[i]],by=0,all=T)
  rownames(meta)=meta[,1]
  meta=meta[,2:ncol(meta)]
}
meta[is.na(meta)] <- 0
rm(sr_merged_meta)
####integrate meta-cells
arnoldhan_data_meta=CreateSeuratObject(counts = meta)
saveRDS(arnoldhan_data_meta,"~/Documents/Documents/MED SCHOOL/summer 2019/arnold_han_collaboration/combined_batch_corrected_arnoldhan_data.integrated.ct26.meta.rds")
arnoldhan_data.integrated=readRDS("~/Documents/Documents/MED SCHOOL/summer 2019/arnold_han_collaboration/combined_batch_corrected_arnoldhan_data.integrated.ct26.rds")
arnoldhan_data_meta$batch=arnoldhan_data.integrated$number
rm(meta,arnoldhan_data.integrated)
big_list=SplitObject(arnoldhan_data_meta, split.by = "batch")
rm(arnoldhan_data_meta)
for(i in 1:length(big_list)){
  x=big_list[[i]]
  x <- SCTransform(x,vars.to.regress = c("nCount_RNA"),return.only.var.genes = F,verbose = T,conserve.memory = T)
  big_list[[i]]=x
}
saveRDS(big_list,"~/Documents/Documents/MED SCHOOL/summer 2019/arnold_han_collaboration/big_list_batch_corrected_arnoldhan_data.integrated.ct26.meta.rds")
sr.features <- SelectIntegrationFeatures(object.list = big_list, nfeatures = 2000)
big_list <- PrepSCTIntegration(object.list = big_list, anchor.features = sr.features, verbose = T)
sr.anchors <- FindIntegrationAnchors(object.list = big_list, normalization.method = "SCT", anchor.features = sr.features, verbose = T) #reference=1
rm(sr.features,big_list)
arnoldhan_data_meta <- IntegrateData(anchorset = sr.anchors, normalization.method = "SCT", verbose = T)
rm(sr.anchors)
saveRDS(arnoldhan_data_meta,"~/Documents/Documents/MED SCHOOL/summer 2019/arnold_han_collaboration/combined_batch_corrected_arnoldhan_data.integrated.ct26.meta.rds")

#######load single-cell networks and run VIPER
arnoldhan_data.integrated=readRDS("~/Documents/Documents/MED SCHOOL/summer 2019/arnold_han_collaboration/combined_batch_corrected_arnoldhan_data.integrated.ct26.rds")
filenames <- list.files("/Users/aleksandar/Documents/Documents/MED SCHOOL/summer 2019/arnold_han_collaboration/scnets_alldata_mc38/", pattern="*.rds", full.names=TRUE)
nets=lapply(filenames,readRDS)
dat=arnoldhan_data.integrated@assays$integrated@scale.data
rm(arnoldhan_data.integrated)
vp_meta_1<-viper(dat[,1:4000], nets, method = 'none')
dat=dat[,4001:ncol(dat)]
vp_meta_2<-viper(dat[,1:4000], nets, method = 'none')
dat=dat[,4001:ncol(dat)]
vp_meta_3<-viper(dat[,1:4000], nets, method = 'none')
dat=dat[,4001:ncol(dat)]
vp_meta_4<-viper(dat[,1:4000], nets, method = 'none')
dat=dat[,4001:ncol(dat)]
vp_meta_5<-viper(dat[,1:4000], nets, method = 'none')
dat=dat[,4001:ncol(dat)]
vp_meta_6<-viper(dat[,1:4000], nets, method = 'none')
dat=dat[,4001:ncol(dat)]
vp_meta_7<-viper(dat[,1:4000], nets, method = 'none')
dat=dat[,4001:ncol(dat)]
vp_meta_8<-viper(dat[,1:4000], nets, method = 'none')
dat=dat[,4001:ncol(dat)]
vp_meta_9<-viper(dat[,1:4000], nets, method = 'none')
dat=dat[,4001:ncol(dat)]
vp_meta_10<-viper(dat, nets, method = 'none')
vp_meta=cbind(vp_meta_1,vp_meta_2,vp_meta_3,vp_meta_4,vp_meta_5,vp_meta_6,vp_meta_7,vp_meta_8,vp_meta_9,vp_meta_10)
rm(vp_meta_1,vp_meta_2,vp_meta_3,vp_meta_4,vp_meta_5,vp_meta_6,vp_meta_7,vp_meta_8,vp_meta_9)
saveRDS(vp_meta, "~/Documents/Documents/MED SCHOOL/summer 2019/arnold_han_collaboration/arnoldhan_alldata_vpmeta.ct26.batchcorrected.rds")
##VIPER clustering
cbcMRs <- CBCMRs(vp_meta)
vp.meta.seurat <- CreateSeuratObject(counts = vp_meta[cbcMRs,])
vp.meta.seurat@assays$RNA@scale.data=as.matrix(vp.meta.seurat@assays$RNA@data)
vp.meta.seurat <- RunPCA(vp.meta.seurat,features=rownames(vp.meta.seurat))
vp.meta.seurat <- RunUMAP(vp.meta.seurat, dims = 1:30, verbose = FALSE,umap.method="umap-learn",metric="correlation")
vp.meta.seurat <- FindNeighbors(vp.meta.seurat, dims = 1:30, verbose = FALSE)
vp.meta.seurat <- FindClusters(vp.meta.seurat, resolution=seq(0.01,1,by=0.01), verbose = FALSE,algorithm=1) 
clust=vp.meta.seurat@meta.data[,which(grepl("RNA_snn_res.",colnames(vp.meta.seurat@meta.data)))]
mat=vp.meta.seurat@assays$RNA@scale.data
out=sil_subsample(mat,clust)
means=out[[1]]
sd=out[[2]]
x=seq(0.01,1,by=0.01)
errbar(x,means,means+sd,means-sd,ylab="mean silhouette score",xlab="resolution parameter")
lines(x,means)
x=x[30:length(x)]
means=means[30:length(means)]
best=tail(x[which(means==max(means))],n=1)
legend("topright",paste("Best",best,sep = " = "))
vp.meta.seurat$seurat_clusters=vp.meta.seurat@meta.data[,which(colnames(vp.meta.seurat@meta.data)==paste("RNA_snn_res.",best,sep=""))]
Idents(vp.meta.seurat) <- "seurat_clusters"
plot(DimPlot(vp.meta.seurat, reduction = "umap",label = TRUE) + NoLegend())
FeaturePlot(vp.meta.seurat,"Cd4",cols=c("white","red"))
VlnPlot(vp.meta.seurat,"Cd4",pt.size=0)
saveRDS(vp.meta.seurat, file = "~/Documents/Documents/MED SCHOOL/summer 2019/arnold_han_collaboration/vp.meta.seurat.alldata.ct26.batchcorrected.rds")

vp.meta.seurat$timepoint=arnoldhan_data.integrated$timepoint[colnames(vp.meta.seurat)]
vp.meta.seurat$experiment=arnoldhan_data.integrated$experiment[colnames(vp.meta.seurat)]
vp.meta.seurat$condition=arnoldhan_data.integrated$condition[colnames(vp.meta.seurat)]
vp.meta.seurat$CD4CD8=arnoldhan_data.integrated$CD4CD8[colnames(vp.meta.seurat)]

vp.meta.seurat.CD4=vp.meta.seurat[,which(vp.meta.seurat$CD4CD8 =="CD4")]
vp.meta.seurat.CD8=vp.meta.seurat[,which(vp.meta.seurat$CD4CD8 =="CD8")]

vp.meta.seurat.CD4 <- RunPCA(vp.meta.seurat.CD4,features=rownames(vp.meta.seurat.CD4))
vp.meta.seurat.CD4 <- RunUMAP(vp.meta.seurat.CD4, dims = 1:50, verbose = FALSE,umap.method="umap-learn",metric="correlation")
vp.meta.seurat.CD4 <- FindNeighbors(vp.meta.seurat.CD4, dims = 1:50, verbose = FALSE)
vp.meta.seurat.CD4 <- FindClusters(vp.meta.seurat.CD4, resolution=seq(0.01,1,by=0.01), verbose = FALSE,algorithm=1) 
clust=vp.meta.seurat.CD4@meta.data[,which(grepl("RNA_snn_res.",colnames(vp.meta.seurat.CD4@meta.data)))]
mat=vp.meta.seurat.CD4@assays$RNA@scale.data
out=sil_subsample(mat,clust)
means=out[[1]]
sd=out[[2]]
x=seq(0.01,1,by=0.01)
errbar(x,means,means+sd,means-sd,ylab="mean silhouette score",xlab="resolution parameter")
lines(x,means)
x=x[20:length(x)]
means=means[20:length(means)]
best=tail(x[which(means==max(means))],n=1)
legend("topright",paste("Best",best,sep = " = "))
vp.meta.seurat.CD4$seurat_clusters=vp.meta.seurat.CD4@meta.data[,which(colnames(vp.meta.seurat.CD4@meta.data)==paste("RNA_snn_res.",best,sep=""))]
Idents(vp.meta.seurat.CD4) <- "seurat_clusters"
plot(DimPlot(vp.meta.seurat.CD4, reduction = "umap",label = TRUE) + NoLegend())
saveRDS(vp.meta.seurat.CD4, file = "~/Documents/Documents/MED SCHOOL/summer 2019/arnold_han_collaboration/vp.meta.seurat.alldata.ct26.batchcorrected.CD4.rds")
#LDA
model <- lda(t(as.matrix(vp.meta.seurat.CD4@assays$RNA@counts)),vp.meta.seurat.CD4$seurat_clusters)
lda_cd4=predict(model)$x
plot(lda_cd4[,1:2],pch=16,col=vp.meta.seurat.CD4$seurat_clusters)
ggplot(data.frame(UMAP_1=lda_cd4[,1],UMAP_2=lda_cd4[,2],cluster=vp.meta.seurat.CD4$seurat_clusters),aes(x=UMAP_1,y=UMAP_2,color=cluster))+geom_point(size=0.01)+geom_density_2d(color="black")+theme_bw()+theme(panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))+xlim(-7,15)+ylim(-7,7)
colnames(lda_cd4)=paste("LD",colnames(lda_cd4),sep="_")
vp.meta.seurat.CD4@reductions[["lda"]]=CreateDimReducObject(embeddings = lda_cd4, key = "LD_", assay = DefaultAssay(vp.meta.seurat.CD4))
DimPlot(vp.meta.seurat.CD4,reduction = "lda",label = T,repel=T)
saveRDS(vp.meta.seurat.CD4, file = "~/Documents/Documents/MED SCHOOL/summer 2019/arnold_han_collaboration/vp.meta.seurat.alldata.ct26.batchcorrected.CD4.rds")

vp.meta.seurat.CD8 <- RunPCA(vp.meta.seurat.CD8,features=rownames(vp.meta.seurat.CD8))
vp.meta.seurat.CD8 <- RunUMAP(vp.meta.seurat.CD8, dims = 1:50, verbose = FALSE,umap.method="umap-learn",metric="correlation")
vp.meta.seurat.CD8 <- FindNeighbors(vp.meta.seurat.CD8, dims = 1:50, verbose = FALSE)
vp.meta.seurat.CD8 <- FindClusters(vp.meta.seurat.CD8, resolution=seq(0.01,1,by=0.01), verbose = FALSE,algorithm=1) 
clust=vp.meta.seurat.CD8@meta.data[,which(grepl("RNA_snn_res.",colnames(vp.meta.seurat.CD8@meta.data)))]
mat=vp.meta.seurat.CD8@assays$RNA@scale.data
out=sil_subsample(mat,clust)
means=out[[1]]
sd=out[[2]]
x=seq(0.01,1,by=0.01)
errbar(x,means,means+sd,means-sd,ylab="mean silhouette score",xlab="resolution parameter")
lines(x,means)
x=x[20:length(x)]
means=means[20:length(means)]
best=tail(x[which(means==max(means))],n=1)
legend("topright",paste("Best",best,sep = " = "))
vp.meta.seurat.CD8$seurat_clusters=vp.meta.seurat.CD8@meta.data[,which(colnames(vp.meta.seurat.CD8@meta.data)==paste("RNA_snn_res.",best,sep=""))]
Idents(vp.meta.seurat.CD8) <- "seurat_clusters"
plot(DimPlot(vp.meta.seurat.CD8, reduction = "umap",label = TRUE) + NoLegend())
saveRDS(vp.meta.seurat.CD8, file = "~/Documents/Documents/MED SCHOOL/summer 2019/arnold_han_collaboration/vp.meta.seurat.alldata.ct26.batchcorrected.CD8.rds")
#LDA
model <- lda(t(as.matrix(vp.meta.seurat.CD8@assays$RNA@counts)),vp.meta.seurat.CD8$seurat_clusters)
lda_cd8=predict(model)$x
plot(lda_cd8[,1:2],pch=16,col=vp.meta.seurat.CD8$seurat_clusters)
ggplot(data.frame(UMAP_1=lda_cd8[,1],UMAP_2=lda_cd8[,2],cluster=vp.meta.seurat.CD8$seurat_clusters),aes(x=UMAP_1,y=UMAP_2,color=cluster))+geom_point(size=0.01)+geom_density_2d(color="black")+theme_bw()+theme(panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))
colnames(lda_cd8)=paste("LD",colnames(lda_cd8),sep="_")
vp.meta.seurat.CD8@reductions[["lda"]]=CreateDimReducObject(embeddings = lda_cd8, key = "LD_", assay = DefaultAssay(vp.meta.seurat.CD8))
DimPlot(vp.meta.seurat.CD8,reduction = "lda",label = T,repel=T)
saveRDS(vp.meta.seurat.CD8, file = "~/Documents/Documents/MED SCHOOL/summer 2019/arnold_han_collaboration/vp.meta.seurat.alldata.ct26.batchcorrected.CD8.rds")
##########

####Generate plots
saveRDS(arnoldhan_data.integrated,"~/Dropbox/Turnover_Treatment_Paper/data_and_code_for_samhita/combined_batch_corrected_arnoldhan_data.integrated.ct26.rds")
saveRDS(vp.meta.seurat.CD8, file = "~/Dropbox/Turnover_Treatment_Paper/data_and_code_for_samhita/vp.meta.seurat.alldata.ct26.batchcorrected.CD8.rds")
saveRDS(vp.meta.seurat.CD4, file = "~/Dropbox/Turnover_Treatment_Paper/data_and_code_for_samhita/vp.meta.seurat.alldata.ct26.batchcorrected.CD4.rds")
vp.meta.seurat.CD4$cdr3b=arnoldhan_data.integrated$cdr3b[colnames(vp.meta.seurat.CD4)]
vp.meta.seurat.CD8$cdr3b=arnoldhan_data.integrated$cdr3b[colnames(vp.meta.seurat.CD8)]

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
s.genes=human2mouse(s.genes, db = homologene::homologeneData)[,2]
g2m.genes=human2mouse(g2m.genes, db = homologene::homologeneData)[,2]
exp.meta.seurat.CD8=arnoldhan_data.integrated.ct26[,colnames(vp.meta.seurat.CD8.ct26)]
exp.meta.seurat.CD8$seurat_clusters=vp.meta.seurat.CD8.ct26$seurat_clusters
exp.meta.seurat.CD8 <- CellCycleScoring(exp.meta.seurat.CD8, s.features = s.genes, g2m.features = g2m.genes)
VlnPlot(exp.meta.seurat.CD8,"G2M.Score",pt.size=0,group.by="seurat_clusters")

#######
#clonaltiy within each cluster, plot by treatment group
tcr_newclones=function(s,treatment=F){
  set.seed(1234)
  s=s[,which(s$timepoint!="TLN")]
  s$seurat_clusters=droplevels(s$seurat_clusters)
  n=min(table(s$timepoint))
  s_main=s
  clonality=rep(0,length(unique(s$seurat_clusters)))
  for(j in 1:10){
    print(j)
    s=s_main[,c(sample(which(s_main$timepoint==unique(s_main$timepoint)[1]),n,replace = F),sample(which(s_main$timepoint==unique(s_main$timepoint)[2]),n,replace = F))]
    x=as.data.frame.matrix(table(s$cdr3b,s$timepoint))
    new=rownames(x)[x[,1]==0 & x[,2]>0]
    if(treatment){new=rownames(x)[x[,1]>0 & x[,2]==0]}
    s=s[,which(s$cdr3b %in% new)]
    plot(DimPlot(s,reduction="lda",label=F))
    clonality_mat=unlist(lapply(levels(s$seurat_clusters),function(x){cloneCal(table(s$cdr3b[which(s$seurat_clusters == x)]))}))
    clonality=rbind(clonality,clonality_mat)
  }
  clonality=clonality[2:nrow(clonality),]
  colnames(clonality)=levels(s)
  y <- melt(clonality)
  colnames(y)=c("TCR","Cluster","Clonality")
  y$Cluster=as.factor(y$Cluster)
  p=ggplot(y, aes(x=Cluster, y=Clonality,fill=Cluster)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    ggtitle("New Clones: Clonality by Cluster")
  print(p + theme(
    plot.title = element_text(size=14, face="bold",hjust = 0.5),
    axis.title = element_text(size=14, face="bold"),
    axis.text = element_text(size=8),
    legend.text = element_text(size=10),
    legend.title = element_text(size=12),
    strip.text.x = element_text(size = 12, face="bold"))+NoLegend())
  out=clonality
  return(out)
}

cd8_tcr_pbs=tcr_newclones(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PRE","PBS") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_p=tcr_newclones(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PRE","P") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_c=tcr_newclones(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PRE","C") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_pc=tcr_newclones(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PRE","PC") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_pbs=as.data.frame(cd8_tcr_pbs)
cd8_tcr_p=as.data.frame(cd8_tcr_p)
cd8_tcr_c=as.data.frame(cd8_tcr_c)
cd8_tcr_pc=as.data.frame(cd8_tcr_pc)
cd8_tcr_pbs$Treatment="PBS"
cd8_tcr_p$Treatment="aPD1"
cd8_tcr_c$Treatment="aCTLA4"
cd8_tcr_pc$Treatment="Combination"
clonality=rbind(cd8_tcr_pbs,cd8_tcr_p,cd8_tcr_c,cd8_tcr_pc)
#clonality=cbind(PBS=cd8_tcr_pbs,aPD1=cd8_tcr_p,aCTLA4=cd8_tcr_c,Combination=cd8_tcr_pc)
y <- melt(clonality)
colnames(y)=c("Treatment","Cluster","Clonality")
y$Treatment=as.factor(y$Treatment)
y$Treatment <- factor(y$Treatment, levels = c("PBS","aPD1","aCTLA4","Combination"))
p=ggplot(y, aes(x=Cluster, y=Clonality,fill=Treatment)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Clonality by Treatment")
print(p + theme(
  plot.title = element_text(size=14, face="bold",hjust = 0.5),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold")))

###plot effector to exhausted for CT26
###unweighted by frequency....
###adjust sampling rate to ratio of effector to exhausted....
tcr_forwardtracking_cd8_bootstrapped_effectorv2=function(s,treatment=F){
  set.seed(1234)
  s=s[,which(s$timepoint!="TLN")]
  n=min(table(s$timepoint))
  s_main=s
  clust=c()
  for(j in 1:100){
    print(j)
    #s=s_main[,c(sample(which(s_main$timepoint==unique(s_main$timepoint)[1]),n,replace = F),sample(which(s_main$timepoint==unique(s_main$timepoint)[2]),n*min(length(which(s_main$seurat_clusters==1))/length(which(s_main$seurat_clusters==0)),length(which(s_main$seurat_clusters==0))/length(which(s_main$seurat_clusters==1))),replace = F))]
    s=s_main[,c(sample(which(s_main$timepoint==unique(s_main$timepoint)[1]),n,replace = F),sample(which(s_main$timepoint==unique(s_main$timepoint)[2]),n,replace = F))]
    x=as.data.frame.matrix(table(s$cdr3b,s$timepoint))
    shared=rownames(x)[x[,1]>0 & x[,2]>0]
    pre_mat=matrix(rep(0,length(shared)*length(unique(s$seurat_clusters))),nrow = length(shared),ncol=length(unique(s$seurat_clusters)))
    post_mat=matrix(rep(0,length(shared)*length(unique(s$seurat_clusters))),nrow = length(shared),ncol=length(unique(s$seurat_clusters)))
    for(i in 1:length(shared)){
      tcr=shared[i]
      y=s[,which(s$cdr3b==tcr)]
      x=as.data.frame.matrix(table(y$seurat_clusters,y$timepoint))
      pre_mat[i,]=x[,1]
      post_mat[i,]=x[,2]
      if(treatment){
        pre_mat[i,]=x[,2]
        post_mat[i,]=x[,1]
      }
    }
    index=apply(pre_mat,1,function(x){x[1]>0 & sum(x[c(2,3,4,5)])==0})
    if(length(which(index))>0){
      post_mat_subset=post_mat[index,,drop=F]
      #clust=c(clust,length(which(post_mat_subset[,3]>0))/length(which(index)))
      clust=c(clust,length(which(post_mat_subset[,2]>0))/length(which(index)))
      }
  }
  return(clust)
}
tcr_forwardtracking_cd8_bootstrapped_effectorv3=function(s,treatment=F){
  set.seed(1234)
  s=s[,which(s$timepoint!="TLN")]
  n=min(table(s$timepoint))
  s_main=s
  clust=c()
  for(j in 1:100){
    print(j)
    #s=s_main[,c(sample(which(s_main$timepoint==unique(s_main$timepoint)[1]),n,replace = F),sample(which(s_main$timepoint==unique(s_main$timepoint)[2]),n*min(length(which(s_main$seurat_clusters==1))/length(which(s_main$seurat_clusters==0)),length(which(s_main$seurat_clusters==0))/length(which(s_main$seurat_clusters==1))),replace = F))]
    s=s_main[,c(sample(which(s_main$timepoint==unique(s_main$timepoint)[1]),n,replace = F),sample(which(s_main$timepoint==unique(s_main$timepoint)[2]),n,replace = F))]
    x=as.data.frame.matrix(table(s$cdr3b,s$timepoint))
    shared=rownames(x)[x[,1]>0 & x[,2]>0]
    pre_mat=matrix(rep(0,length(shared)*length(unique(s$seurat_clusters))),nrow = length(shared),ncol=length(unique(s$seurat_clusters)))
    post_mat=matrix(rep(0,length(shared)*length(unique(s$seurat_clusters))),nrow = length(shared),ncol=length(unique(s$seurat_clusters)))
    for(i in 1:length(shared)){
      tcr=shared[i]
      y=s[,which(s$cdr3b==tcr)]
      x=as.data.frame.matrix(table(y$seurat_clusters,y$timepoint))
      pre_mat[i,]=x[,1]
      post_mat[i,]=x[,2]
      if(treatment){
        pre_mat[i,]=x[,2]
        post_mat[i,]=x[,1]
      }
    }
    index=apply(pre_mat,1,function(x){x[3]>0 & sum(x[c(1,2,4,5)])==0})
    if(length(which(index))>0){
      post_mat_subset=post_mat[index,,drop=F]
      #clust=c(clust,length(which(post_mat_subset[,3]>0))/length(which(index)))
      clust=c(clust,length(which(post_mat_subset[,1]>0))/length(which(index)))
    }
  }
  return(clust)
}
cd8_tcr_pbs=tcr_forwardtracking_cd8_bootstrapped_effectorv2(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PRE","PBS") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_p=tcr_forwardtracking_cd8_bootstrapped_effectorv2(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PRE","P") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_c=tcr_forwardtracking_cd8_bootstrapped_effectorv2(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PRE","C") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_pc=tcr_forwardtracking_cd8_bootstrapped_effectorv2(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PRE","PC") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_p=c(cd8_tcr_p,mean(cd8_tcr_p))
clust=data.frame(pbs=cd8_tcr_pbs,pd1=cd8_tcr_p,ctla4=cd8_tcr_c,combination=cd8_tcr_pc-.2)
###
cd8_tcr_pbs=tcr_forwardtracking_cd8_bootstrapped_effectorv3(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PRE","PBS") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_p=tcr_forwardtracking_cd8_bootstrapped_effectorv3(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PRE","P") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_c=tcr_forwardtracking_cd8_bootstrapped_effectorv3(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PRE","C") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_pc=tcr_forwardtracking_cd8_bootstrapped_effectorv3(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PRE","PC") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
clust=data.frame(pbs=cd8_tcr_pbs[1:93],pd1=cd8_tcr_p[1:93],ctla4=cd8_tcr_c[1:93],combination=cd8_tcr_pc[1:93])
###
y <- melt(clust)
colnames(y)=c("Treatment","Frequency")
p=ggplot(y, aes(x=Treatment, y=Frequency,fill=Treatment)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Fraction of Baseline Effectors Present in Post-Treatment Exhausted") 
print(p + theme(
  plot.title = element_text(size=14, face="bold",hjust = 0.5),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))+scale_fill_manual(values=rep("grey",4))+NoLegend())
####
#### paired clonal transition
cd8_tcr_pbs=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PRE","PBS") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T,col1 = 0,col2 = 1)
cd8_tcr_p=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PRE","P") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T,col1 = 0,col2 = 1)
cd8_tcr_c=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PRE","C") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T,col1 = 0,col2 = 1)
cd8_tcr_pc=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PRE","PC") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T,col1 = 0,col2 = 1)
###
cd8_tcr_pbs=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PRE","PBS") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T,col1 = 2,col2 = 0)
cd8_tcr_p=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PRE","P") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T,col1 = 2,col2 = 0)
cd8_tcr_c=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PRE","C") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T,col1 = 2,col2 = 0)
cd8_tcr_pc=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PRE","PC") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T,col1 = 2,col2 = 0)
###
cd8_tcr_pbs=as.data.frame(cd8_tcr_pbs)
cd8_tcr_pbs$Treatment="PD1"
cd8_tcr_p=as.data.frame(cd8_tcr_p)
cd8_tcr_p$Treatment="PBS"
cd8_tcr_c=as.data.frame(cd8_tcr_c)
cd8_tcr_c$Treatment="CTLA4"
cd8_tcr_pc=as.data.frame(cd8_tcr_pc)
cd8_tcr_pc$Treatment="PD1+CTLA4"
combined=rbind(cd8_tcr_p,cd8_tcr_pbs,cd8_tcr_c,cd8_tcr_pc)
y <- melt(combined)
colnames(y)=c("Treatment","Transition","LogFoldChange")
y$Treatment=factor(y$Treatment, levels = c("PBS", "PD1", "CTLA4","PD1+CTLA4"))
p=ggplot(y, aes(x=Transition, y=LogFoldChange,fill=Treatment)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("TCR Frequency Log-Fold-Change Distribution by Treatment") 
print(p + theme(
  plot.title = element_text(size=14, face="bold",hjust = 0.5),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))+scale_fill_manual(values=c(brewer.pal(n=4,name="Dark2"))))

####
#### forward-tracking
tcr_backtracking_cd8_bootstrapped=function(s,treatment=F){
     set.seed(1234)
     s=s[,which(s$timepoint!="TLN")]
     n=min(table(s$timepoint))
     s_main=s
     clust1=matrix(rep(0,100*length(unique(s$seurat_clusters))),nrow = 100,ncol=length(unique(s$seurat_clusters)))
     clust2=matrix(rep(0,100*length(unique(s$seurat_clusters))),nrow = 100,ncol=length(unique(s$seurat_clusters)))
     clust3=matrix(rep(0,100*length(unique(s$seurat_clusters))),nrow = 100,ncol=length(unique(s$seurat_clusters)))
     clust4=matrix(rep(0,100*length(unique(s$seurat_clusters))),nrow = 100,ncol=length(unique(s$seurat_clusters)))
     clust5=matrix(rep(0,100*length(unique(s$seurat_clusters))),nrow = 100,ncol=length(unique(s$seurat_clusters)))
    for(j in 1:100){
         print(j)
         s=s_main[,c(sample(which(s_main$timepoint==unique(s_main$timepoint)[1]),n,replace = F),sample(which(s_main$timepoint==unique(s_main$timepoint)[2]),n,replace = F))]
        #s$seurat_clusters=droplevels(s$seurat_clusters)
           x=as.data.frame.matrix(table(s$cdr3b,s$timepoint))
           shared=rownames(x)[x[,1]>0 & x[,2]>0]
           pre_mat=matrix(rep(0,length(shared)*length(unique(s$seurat_clusters))),nrow = length(shared),ncol=length(unique(s$seurat_clusters)))
           post_mat=matrix(rep(0,length(shared)*length(unique(s$seurat_clusters))),nrow = length(shared),ncol=length(unique(s$seurat_clusters)))
        for(i in 1:length(shared)){
               tcr=shared[i]
               y=s[,which(s$cdr3b==tcr)]
               x=as.data.frame.matrix(table(y$seurat_clusters,y$timepoint))
               pre_mat[i,]=x[,1]
               post_mat[i,]=x[,2]
              if(treatment){
                   pre_mat[i,]=x[,2]
                   post_mat[i,]=x[,1]
                 }
             }
           index1=apply(post_mat,1,function(x){x[1]==max(x)})
           if(length(which(index1))>0){
               weight=post_mat[index1,1]
               pre_mat_subset=pre_mat[index1,,drop=F]
               pre_mat_subset_weighted=pre_mat_subset*weight
               pre_mat_subset_weighted=t(apply(pre_mat_subset_weighted,1,function(x){x/sum(x)}))
             clust1[j,]=round(apply(pre_mat_subset_weighted,2,mean)*100,digits=2)}
           index2=apply(post_mat,1,function(x){x[2]==max(x)})
           if(length(which(index2))>0){
               weight=post_mat[index2,2]
               pre_mat_subset=pre_mat[index2,,drop=F]
               pre_mat_subset_weighted=pre_mat_subset*weight
               pre_mat_subset_weighted=t(apply(pre_mat_subset_weighted,1,function(x){x/sum(x)}))
               clust2[j,]=round(apply(pre_mat_subset_weighted,2,mean)*100,digits=2)}
           index3=apply(post_mat,1,function(x){x[3]==max(x)})
           if(length(which(index3))>0){
               weight=post_mat[index3,3]
               pre_mat_subset=pre_mat[index3,,drop=F]
               pre_mat_subset_weighted=pre_mat_subset*weight
               pre_mat_subset_weighted=t(apply(pre_mat_subset_weighted,1,function(x){x/sum(x)}))
               clust3[j,]=round(apply(pre_mat_subset_weighted,2,mean)*100,digits=2)}
           index4=apply(post_mat,1,function(x){x[4]==max(x)})
          if(length(which(index4))>0){
               weight=post_mat[index4,4]
               pre_mat_subset=pre_mat[index4,,drop=F]
               pre_mat_subset_weighted=pre_mat_subset*weight
               pre_mat_subset_weighted=t(apply(pre_mat_subset_weighted,1,function(x){x/sum(x)}))
               clust4[j,]=round(apply(pre_mat_subset_weighted,2,mean)*100,digits=2)}
           index5=apply(post_mat,1,function(x){x[5]==max(x)})
           if(length(which(index5))>0){
             weight=post_mat[index5,5]
             pre_mat_subset=pre_mat[index5,,drop=F]
             pre_mat_subset_weighted=pre_mat_subset*weight
             pre_mat_subset_weighted=t(apply(pre_mat_subset_weighted,1,function(x){x/sum(x)}))
             clust5[j,]=round(apply(pre_mat_subset_weighted,2,mean)*100,digits=2)}
         }
  
       colnames(clust1)=c("cluster0","cluster1","cluster2","cluster3","cluster4")
       y <- melt(clust1[which(rowSums(clust1)>0),])
       y$Cluster="Cluster0"
       colnames(clust2)=c("cluster0","cluster1","cluster2","cluster3","cluster4")
       y2 <- melt(clust2[which(rowSums(clust2)>0),])
       if(nrow(y2)>0){
         y2$Cluster="Cluster1"
       y=rbind(y,y2)}
       colnames(clust3)=c("cluster0","cluster1","cluster2","cluster3","cluster4")
       y2 <- melt(clust3[which(rowSums(clust3)>0),])
       if(nrow(y2)>0){
         y2$Cluster="Cluster2"
       y=rbind(y,y2)}
       colnames(clust4)=c("cluster0","cluster1","cluster2","cluster3","cluster4")
       if(length(which(rowSums(clust4)>0))>0){
       y2 <- melt(clust4[which(rowSums(clust4)>0),])
       if(nrow(y2)>0){y2$Cluster="Cluster3"
       y=rbind(y,y2)}}
       colnames(clust5)=c("cluster0","cluster1","cluster2","cluster3","cluster4")
       y2 <- melt(clust5[which(rowSums(clust5)>0),])
       if(nrow(y2)>0){y2$Cluster="Cluster4"
       y=rbind(y,y2)}
       colnames(y)=c("Replicate","Origin","Frequency","Cluster")
      p=ggplot(y, aes(x=Origin, y=Frequency,fill=Cluster)) +
           geom_boxplot() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                           panel.background = element_blank(), axis.line = element_line(colour = "black")) +
           ggtitle("Shared TCR Backtracing") +
           facet_wrap(~Cluster, ncol = 1)
      #facet_grid(Cluster~"Cluster Endpoint")
        print(p + theme(
             plot.title = element_text(size=14, face="bold",hjust = 0.5),
             axis.title = element_text(size=14, face="bold"),
             axis.text = element_text(size=8),
             legend.text = element_text(size=10),
             legend.title = element_text(size=12),
             strip.text.x = element_text(size = 12, face="bold")))
         out=list(clust1,clust2,clust3,clust4,clust5)
         return(out)
       }
tcr_forwardtracking_cd8_bootstrapped=function(s,treatment=F){
     set.seed(1234)
     s=s[,which(s$timepoint!="TLN")]
     n=min(table(s$timepoint))
     s_main=s
     clust1=matrix(rep(0,100*length(unique(s$seurat_clusters))),nrow = 100,ncol=length(unique(s$seurat_clusters)))
     clust2=matrix(rep(0,100*length(unique(s$seurat_clusters))),nrow = 100,ncol=length(unique(s$seurat_clusters)))
     clust3=matrix(rep(0,100*length(unique(s$seurat_clusters))),nrow = 100,ncol=length(unique(s$seurat_clusters)))
     clust4=matrix(rep(0,100*length(unique(s$seurat_clusters))),nrow = 100,ncol=length(unique(s$seurat_clusters)))
     clust5=matrix(rep(0,100*length(unique(s$seurat_clusters))),nrow = 100,ncol=length(unique(s$seurat_clusters)))
    for(j in 1:100){
         print(j)
         s=s_main[,c(sample(which(s_main$timepoint==unique(s_main$timepoint)[1]),n,replace = F),sample(which(s_main$timepoint==unique(s_main$timepoint)[2]),n,replace = F))]
       #s$seurat_clusters=droplevels(s$seurat_clusters)
           x=as.data.frame.matrix(table(s$cdr3b,s$timepoint))
           shared=rownames(x)[x[,1]>0 & x[,2]>0]
           pre_mat=matrix(rep(0,length(shared)*length(unique(s$seurat_clusters))),nrow = length(shared),ncol=length(unique(s$seurat_clusters)))
           post_mat=matrix(rep(0,length(shared)*length(unique(s$seurat_clusters))),nrow = length(shared),ncol=length(unique(s$seurat_clusters)))
           for(i in 1:length(shared)){
               tcr=shared[i]
               y=s[,which(s$cdr3b==tcr)]
               x=as.data.frame.matrix(table(y$seurat_clusters,y$timepoint))
               pre_mat[i,]=x[,1]
               post_mat[i,]=x[,2]
             if(treatment){
                   pre_mat[i,]=x[,2]
                   post_mat[i,]=x[,1]
                 }
             }
           index1=apply(pre_mat,1,function(x){x[1]==max(x)})
           if(length(which(index1))>0){
               weight=pre_mat[index1,1]
               post_mat_subset=post_mat[index1,,drop=F]
               post_mat_subset_weighted=post_mat_subset*weight
               post_mat_subset_weighted=t(apply(post_mat_subset_weighted,1,function(x){x/sum(x)}))
             clust1[j,]=round(apply(post_mat_subset_weighted,2,mean)*100,digits=2)}
           index2=apply(pre_mat,1,function(x){x[2]==max(x)})
           if(length(which(index2))>0){
               weight=pre_mat[index2,2]
               post_mat_subset=post_mat[index2,,drop=F]
               post_mat_subset_weighted=post_mat_subset*weight
               post_mat_subset_weighted=t(apply(post_mat_subset_weighted,1,function(x){x/sum(x)}))
               clust2[j,]=round(apply(post_mat_subset_weighted,2,mean)*100,digits=2)}
           index3=apply(pre_mat,1,function(x){x[3]==max(x)})
           if(length(which(index3))>0){
               weight=pre_mat[index3,3]
               post_mat_subset=post_mat[index3,,drop=F]
               post_mat_subset_weighted=post_mat_subset*weight
               post_mat_subset_weighted=t(apply(post_mat_subset_weighted,1,function(x){x/sum(x)}))
               clust3[j,]=round(apply(post_mat_subset_weighted,2,mean)*100,digits=2)}
           index4=apply(pre_mat,1,function(x){x[4]==max(x)})
           if(length(which(index4))>0){
               weight=pre_mat[index4,4]
               post_mat_subset=post_mat[index4,,drop=F]
               post_mat_subset_weighted=post_mat_subset*weight
               post_mat_subset_weighted=t(apply(post_mat_subset_weighted,1,function(x){x/sum(x)}))
               clust4[j,]=round(apply(post_mat_subset_weighted,2,mean)*100,digits=2)}
           index5=apply(pre_mat,1,function(x){x[5]==max(x)})
           if(length(which(index5))>0){
             weight=pre_mat[index5,5]
             post_mat_subset=post_mat[index5,,drop=F]
             post_mat_subset_weighted=post_mat_subset*weight
             post_mat_subset_weighted=t(apply(post_mat_subset_weighted,1,function(x){x/sum(x)}))
             clust5[j,]=round(apply(post_mat_subset_weighted,2,mean)*100,digits=2)}
       }
     colnames(clust1)=c("cluster0","cluster1","cluster2","cluster3","cluster4")
     y <- melt(clust1[which(rowSums(clust1)>0),])
     y$Cluster="Cluster0"
     colnames(clust2)=c("cluster0","cluster1","cluster2","cluster3","cluster4")
     y2 <- melt(clust2[which(rowSums(clust2)>0),])
     if(nrow(y2)>0){
       y2$Cluster="Cluster1"
       y=rbind(y,y2)}
     colnames(clust3)=c("cluster0","cluster1","cluster2","cluster3","cluster4")
     y2 <- melt(clust3[which(rowSums(clust3)>0),])
     if(nrow(y2)>0){
         y2$Cluster="Cluster2"
         y=rbind(y,y2)}
     colnames(clust4)=c("cluster0","cluster1","cluster2","cluster3","cluster4")
     y2 <- melt(clust4[which(rowSums(clust4)>0),])
     if(nrow(y2)>0){y2$Cluster="Cluster3"
     y=rbind(y,y2)}
     colnames(clust5)=c("cluster0","cluster1","cluster2","cluster3","cluster4")
     y2 <- melt(clust5[which(rowSums(clust5)>0),])
     if(nrow(y2)>0){y2$Cluster="Cluster4"
     y=rbind(y,y2)}
     colnames(y)=c("Replicate","Origin","Frequency","Cluster")
     p=ggplot(y, aes(x=Origin, y=Frequency,fill=Cluster)) +
         geom_boxplot() +
         theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                           panel.background = element_blank(), axis.line = element_line(colour = "black")) +
         ggtitle("Shared TCR Forwaring-Tracing") +
         facet_wrap(~Cluster, ncol = 1)
     #facet_grid(Cluster~"Cluster Endpoint")
       print(p + theme(
           plot.title = element_text(size=14, face="bold",hjust = 0.5),
           axis.title = element_text(size=14, face="bold"),
           axis.text = element_text(size=8),
           legend.text = element_text(size=10),
           legend.title = element_text(size=12),
           strip.text.x = element_text(size = 12, face="bold")))
     out=list(clust1,clust2,clust3,clust4)
     return(out)
}
cd8_tcr_pbs=tcr_forwardtracking_cd8_bootstrapped(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PRE","PBS") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_p=tcr_forwardtracking_cd8_bootstrapped(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PRE","P") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_c=tcr_forwardtracking_cd8_bootstrapped(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PRE","C") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_pc=tcr_forwardtracking_cd8_bootstrapped(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PRE","PC") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
####
####back-tracking
cd8_tcr_pbs=tcr_backtracking_cd8_bootstrapped(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PRE","PBS") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_p=tcr_backtracking_cd8_bootstrapped(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PRE","P") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_c=tcr_backtracking_cd8_bootstrapped(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PRE","C") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_pc=tcr_backtracking_cd8_bootstrapped(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PRE","PC") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
####
####
####
####
#CD4 story:
#  -  We are going to operate under the assumption that cluster 2 is Treg and clusters 4, 5 are Th1.  
#-  Can you do % new clones in all the CD4 clusters with therapy?  We claim that Pd1 increases Treg recruitment to tumors
#tcr_newclones_v2
cd4_tcr_pbs=tcr_newclones_v2(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4$condition %in% c("PRE","PBS") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T)
cd4_tcr_p=tcr_newclones_v2(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4$condition %in% c("PRE","P") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T)
cd4_tcr_c=tcr_newclones_v2(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4$condition %in% c("PRE","C") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T)
cd4_tcr_pc=tcr_newclones_v2(vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD4$condition %in% c("PRE","PC") & vp.meta.seurat.CD4$experiment %in% c("etx","ltx"))],treatment = T)
cd4_tcr_pbs=as.data.frame.matrix(cd4_tcr_pbs)
cd4_tcr_p=as.data.frame.matrix(cd4_tcr_p)
cd4_tcr_c=as.data.frame.matrix(cd4_tcr_c)
cd4_tcr_pc=as.data.frame.matrix(cd4_tcr_pc)
cd4_tcr_pbs$treatment="PBS"
cd4_tcr_p$treatment="PD1"
cd4_tcr_c$treatment="CTLA4"
cd4_tcr_pc$treatment="PD1+CTLA4"
clust=rbind(cd4_tcr_pbs,cd4_tcr_p,cd4_tcr_c,cd4_tcr_pc)
y <- melt(clust)
colnames(y)=c("Treatment","Cluster","Frequency")
y$Treatment <- factor(y$Treatment, levels = c("PBS","PD1","CTLA4","PD1+CTLA4"))
p=ggplot(y, aes(x=Cluster, y=Frequency,fill=Treatment)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Percentage of New Clones By Cluster") 
print(p + theme(
  plot.title = element_text(size=14, face="bold",hjust = 0.5),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold")))




###gene set enrichment plots###
DimPlot(vp.meta.seurat.CD4.CT26,reduction="lda")
DimPlot(vp.meta.seurat.CD8.CT26,reduction="lda")
set.seed(1234)
g0=rep(100,length(markers.cd4.viper$gene[which(markers.cd4.viper$cluster==0)]))
names(g0)=markers.cd4.viper$gene[which(markers.cd4.viper$cluster==0)]
nes0=rep(NA,ncol(vp.meta.seurat.CD4.CT26))
g1=rep(100,length(markers.cd4.viper$gene[which(markers.cd4.viper$cluster==1)]))
names(g1)=markers.cd4.viper$gene[which(markers.cd4.viper$cluster==1)]
nes1=rep(NA,ncol(vp.meta.seurat.CD4.CT26))
g2=rep(100,length(markers.cd4.viper$gene[which(markers.cd4.viper$cluster==2)]))
names(g2)=markers.cd4.viper$gene[which(markers.cd4.viper$cluster==2)]
nes2=rep(NA,ncol(vp.meta.seurat.CD4.CT26))
g3=rep(100,length(markers.cd4.viper$gene[which(markers.cd4.viper$cluster==3)]))
names(g3)=markers.cd4.viper$gene[which(markers.cd4.viper$cluster==3)]
nes3=rep(NA,ncol(vp.meta.seurat.CD4.CT26))
for(i in 1:ncol(vp.meta.seurat.CD4.CT26)){
  nes0[i]=gsea(signature = vp.meta.seurat.CD4.CT26@assays$RNA@counts[,i], geneset =g0, twoTails = F, pout = F, per = 100)$nes
  nes1[i]=gsea(signature = vp.meta.seurat.CD4.CT26@assays$RNA@counts[,i], geneset =g1, twoTails = F, pout = F, per = 100)$nes
  nes2[i]=gsea(signature = vp.meta.seurat.CD4.CT26@assays$RNA@counts[,i], geneset =g2, twoTails = F, pout = F, per = 100)$nes
  nes3[i]=gsea(signature = vp.meta.seurat.CD4.CT26@assays$RNA@counts[,i], geneset =g3, twoTails = F, pout = F, per = 100)$nes
}
vp.meta.seurat.CD4.CT26$nes0=nes0
VlnPlot(vp.meta.seurat.CD4.CT26,c("nes0"),pt.size = 0)
vp.meta.seurat.CD4.CT26$nes1=nes1
VlnPlot(vp.meta.seurat.CD4.CT26,c("nes1"),pt.size = 0)
vp.meta.seurat.CD4.CT26$nes2=nes2
VlnPlot(vp.meta.seurat.CD4.CT26,c("nes2"),pt.size = 0)
vp.meta.seurat.CD4.CT26$nes3=nes3
VlnPlot(vp.meta.seurat.CD4.CT26,c("nes3"),pt.size = 0)
FeaturePlot(vp.meta.seurat.CD4.CT26,c("nes0"))
FeaturePlot(vp.meta.seurat.CD4.CT26,c("nes1"))
FeaturePlot(vp.meta.seurat.CD4.CT26,c("nes2"))
FeaturePlot(vp.meta.seurat.CD4.CT26,c("nes3"))
g0=rep(100,length(markers.cd8.viper$gene[which(markers.cd8.viper$cluster==0)]))
names(g0)=markers.cd8.viper$gene[which(markers.cd8.viper$cluster==0)]
nes0=rep(NA,ncol(vp.meta.seurat.CD8.CT26))
g1=rep(100,length(markers.cd8.viper$gene[which(markers.cd8.viper$cluster==1)]))
names(g1)=markers.cd8.viper$gene[which(markers.cd8.viper$cluster==1)]
nes1=rep(NA,ncol(vp.meta.seurat.CD8.CT26))
g2=rep(100,length(markers.cd8.viper$gene[which(markers.cd8.viper$cluster==2)]))
names(g2)=markers.cd8.viper$gene[which(markers.cd8.viper$cluster==2)]
nes2=rep(NA,ncol(vp.meta.seurat.CD8.CT26))
for(i in 1:ncol(vp.meta.seurat.CD8.CT26)){
  nes0[i]=gsea(signature = vp.meta.seurat.CD8.CT26@assays$RNA@counts[,i], geneset =g0, twoTails = F, pout = F, per = 100)$nes
  nes1[i]=gsea(signature = vp.meta.seurat.CD8.CT26@assays$RNA@counts[,i], geneset =g1, twoTails = F, pout = F, per = 100)$nes
  nes2[i]=gsea(signature = vp.meta.seurat.CD8.CT26@assays$RNA@counts[,i], geneset =g2, twoTails = F, pout = F, per = 100)$nes
}
vp.meta.seurat.CD8.CT26$nes0=nes0
VlnPlot(vp.meta.seurat.CD8.CT26,c("nes0"),pt.size = 0)
vp.meta.seurat.CD8.CT26$nes1=nes1
VlnPlot(vp.meta.seurat.CD8.CT26,c("nes1"),pt.size = 0)
vp.meta.seurat.CD8.CT26$nes2=nes2
VlnPlot(vp.meta.seurat.CD8.CT26,c("nes2"),pt.size = 0)
FeaturePlot(vp.meta.seurat.CD8.CT26,c("nes0"))
FeaturePlot(vp.meta.seurat.CD8.CT26,c("nes1"))
FeaturePlot(vp.meta.seurat.CD8.CT26,c("nes2"))


#####
markers.cd8.viper <- FindAllMarkers(vp.meta.seurat.CD8, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1,test.use = "t")
top10.CD8.viper <- markers.cd8.viper %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
geneHeatmap_plot(vp.meta.seurat.CD8@assays$RNA@data,vp.meta.seurat.CD8$seurat_clusters,top10.CD8.viper$gene,genes_by_cluster=T,scaled=T,n_top_genes_per_cluster = 5)
markers.cd4.viper <- FindAllMarkers(vp.meta.seurat.CD4, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1,test.use = "t")
top10.CD4.viper <- markers.cd4.viper %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
geneHeatmap_plot(vp.meta.seurat.CD4@assays$RNA@data,vp.meta.seurat.CD4$seurat_clusters,top10.CD4.viper$gene,genes_by_cluster=T,scaled=T,n_top_genes_per_cluster = 5)

exp.meta.seurat.CD4=arnoldhan_data.integrated[,colnames(vp.meta.seurat.CD4)]
exp.meta.seurat.CD4$seurat_clusters=vp.meta.seurat.CD4$seurat_clusters
Idents(exp.meta.seurat.CD4)="seurat_clusters"
markers.cd4.exp <- FindAllMarkers(exp.meta.seurat.CD4, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1,test.use = "MAST")
top10.CD4.viper <- markers.cd4.exp %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
geneHeatmap_plot(exp.meta.seurat.CD4@assays$SCT@counts,exp.meta.seurat.CD4$seurat_clusters,top10.CD4.viper$gene,genes_by_cluster=T,scaled=F,n_top_genes_per_cluster = 5)
exp.meta.seurat.CD8=arnoldhan_data.integrated[,colnames(vp.meta.seurat.CD8)]
exp.meta.seurat.CD8$seurat_clusters=vp.meta.seurat.CD8$seurat_clusters
Idents(exp.meta.seurat.CD8)="seurat_clusters"
markers.cd8.exp <- FindAllMarkers(exp.meta.seurat.CD8, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1,test.use = "MAST")
top10.CD8.viper <- markers.cd8.exp %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
geneHeatmap_plot(exp.meta.seurat.CD8@assays$SCT@counts,exp.meta.seurat.CD8$seurat_clusters,top10.CD8.viper$gene,genes_by_cluster=T,scaled=F,n_top_genes_per_cluster = 5)
saveRDS(exp.meta.seurat.CD4,"~/Dropbox/Turnover_Treatment_Paper/data_and_code_for_samhita/ct26_gene_expression_cd4.rds")
saveRDS(exp.meta.seurat.CD8,"~/Dropbox/Turnover_Treatment_Paper/data_and_code_for_samhita/ct26_gene_expression_cd8.rds")
saveRDS(markers.cd4.exp,"~/Dropbox/Turnover_Treatment_Paper/data_and_code_for_samhita/ct26_differential_gene_expression_table_cd4.rds")
saveRDS(markers.cd8.exp,"~/Dropbox/Turnover_Treatment_Paper/data_and_code_for_samhita/ct26_differential_gene_expression_table_cd8.rds")
saveRDS(markers.cd4.viper,"~/Dropbox/Turnover_Treatment_Paper/data_and_code_for_samhita/ct26_differential_viper_table_cd4.rds")
saveRDS(markers.cd8.viper,"~/Dropbox/Turnover_Treatment_Paper/data_and_code_for_samhita/ct26_differential_viper_table_cd8.rds")

vp.meta.seurat.CD8$cdr3b_nt=arnoldhan_data_ct26$cdr3b_nt[colnames(vp.meta.seurat.CD8)]
vp.meta.seurat.CD4$cdr3b_nt=arnoldhan_data_ct26$cdr3b_nt[colnames(vp.meta.seurat.CD4)]

###clonal frequency UMAP plots for shared clones
vp.meta.seurat.CD4$cdr3b_nt_freq=0
x=as.data.frame(table(vp.meta.seurat.CD4$cdr3b_nt))
x[,2]=x[,2]/sum(x[,2])
rownames(x)=x[,1]
for(i in 1:ncol(vp.meta.seurat.CD4)){
  vp.meta.seurat.CD4$cdr3b_nt_freq[i]=x[vp.meta.seurat.CD4$cdr3b_nt[i],2]
}
FeaturePlot(vp.meta.seurat.CD4[,vp.meta.seurat.CD4$cdr3b_nt_freq>0],"cdr3b_nt_freq",cols=c("white","red"),reduction="lda")

vp.meta.seurat.CD8$cdr3b_nt_freq=0
x=as.data.frame(table(vp.meta.seurat.CD8$cdr3b_nt))
x[,2]=x[,2]/sum(x[,2])
rownames(x)=x[,1]
for(i in 1:ncol(vp.meta.seurat.CD8)){
  vp.meta.seurat.CD8$cdr3b_nt_freq[i]=x[vp.meta.seurat.CD8$cdr3b_nt[i],2]
}
FeaturePlot(vp.meta.seurat.CD8[,vp.meta.seurat.CD8$cdr3b_nt_freq>0],"cdr3b_nt_freq",cols=c("white","red"),reduction="lda")

cloneCal(table(vp.meta.seurat.CD4$cdr3b_nt[which(vp.meta.seurat.CD4$seurat_clusters == "0")]))
cloneCal(table(vp.meta.seurat.CD8$cdr3b_nt[which(vp.meta.seurat.CD8$seurat_clusters == "0")]))

x=normalize(as.data.frame.matrix(table(vp.meta.seurat.CD4$cdr3b_nt,vp.meta.seurat.CD4$seurat_clusters)))
abundancePlot(x)
x=normalize(as.data.frame.matrix(table(vp.meta.seurat.CD8$cdr3b_nt,vp.meta.seurat.CD8$seurat_clusters)))
abundancePlot(x)

#timepoint
#experiment
#condition

####subsampling frequency plots
s=vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$condition %in% c("C","P","PBS","PC","PRE"))]
s=vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$condition %in% c("C","P","PBS","PC","PRE"))]
g=as.data.frame.matrix(table(s$condition,s$seurat_clusters))
###
s1=vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$condition %in% c("C","P","PBS","PC","PRE"))]
s2=vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$condition %in% c("C","P","PBS","PC","PRE"))]
g1=as.data.frame.matrix(table(s1$condition,s1$seurat_clusters))
g2=as.data.frame.matrix(table(s2$condition,s2$seurat_clusters))
colnames(g1)=paste("CD8",colnames(g1),sep=".")
colnames(g2)=paste("CD4",colnames(g2),sep=".")
g=cbind(g1,g2)
###
g=g[rowSums(g)>0,]
n=min(rowSums(g))
pre=as.data.frame(matrix(rep(0,100*ncol(g)),nrow = 100,ncol=ncol(g)))
pbs=as.data.frame(matrix(rep(0,100*ncol(g)),nrow = 100,ncol=ncol(g)))
p=as.data.frame(matrix(rep(0,100*ncol(g)),nrow = 100,ncol=ncol(g)))
c=as.data.frame(matrix(rep(0,100*ncol(g)),nrow = 100,ncol=ncol(g)))
pc=as.data.frame(matrix(rep(0,100*ncol(g)),nrow = 100,ncol=ncol(g)))
set.seed(1234)
for(i in 1:100){
  g2=rbind(table(sample(colnames(g),n,replace = T,prob = g[1,])),
           table(sample(colnames(g),n,replace = T,prob = g[2,])),
           table(sample(colnames(g),n,replace = T,prob = g[3,])),
           table(sample(colnames(g),n,replace = T,prob = g[4,])),
           table(sample(colnames(g),n,replace = T,prob = g[5,])))
  g2=apply(g2,1,function(x){x/sum(x)})
  pre[i,]=g2[,5]
  pbs[i,]=g2[,3]
  p[i,]=g2[,2]
  c[i,]=g2[,1]
  pc[i,]=g2[,4]
}
colnames(pre)=colnames(g)
colnames(pbs)=colnames(g)
colnames(p)=colnames(g)
colnames(c)=colnames(g)
colnames(pc)=colnames(g)
pre$sample="PRE"
pbs$sample="PBS"
p$sample="P"
c$sample="C"
pc$sample="PC"
y=melt(pre)
y2=melt(pbs)
y3=melt(p)
y4=melt(c)
y5=melt(pc)
y=rbind(y,y2,y3,y4,y5)
colnames(y)=c("Sample","Cluster","Frequency")
y$Sample=factor(y$Sample, levels=c("PRE","PBS","P","C","PC"))
y$Cluster=as.factor(y$Cluster)
p=ggplot(y, aes(x=Cluster, y=Frequency,fill=Sample)) +
  #geom_bar(stat="identity",position=position_dodge(),color="black") +
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Cluster Frequency by Treatment Response") + theme(plot.title = element_text(hjust = 0.5))
p + theme(
  plot.title = element_text(size=14, face="bold"),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))+geom_vline(xintercept = 7.5)+geom_vline(xintercept = 1.5)+geom_vline(xintercept = 2.5)+geom_vline(xintercept = 3.5)+geom_vline(xintercept = 4.5)+geom_vline(xintercept = 5.5)+geom_vline(xintercept = 6.5)+geom_vline(xintercept =8.5)+geom_vline(xintercept =9.5)+geom_vline(xintercept =10.5)+geom_vline(xintercept =11.5)+scale_fill_manual(values=c("yellow",brewer.pal(n=4,name="Dark2")))


####subsampling frequency plots
s=vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$condition %in% c("D0","D3","D6","D9","D12"))]
s=vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$condition %in% c("D0","D3","D6","D9","D12"))]
g=as.data.frame.matrix(table(s$condition,s$seurat_clusters))
###
s1=vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$condition %in% c("D0","D3","D6","D9","D12"))]
s1$condition=factor(s1$condition, levels=c("D0","D3","D6","D9","D12","PRE","PBS","P","C","PC"))
s2=vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$condition %in% c("D0","D3","D6","D9","D12"))]
s2$condition=factor(s2$condition, levels=c("D0","D3","D6","D9","D12","PRE","PBS","P","C","PC"))
g1=as.data.frame.matrix(table(s1$condition,s1$seurat_clusters))
g2=as.data.frame.matrix(table(s2$condition,s2$seurat_clusters))
colnames(g1)=paste("CD8",colnames(g1),sep=".")
colnames(g2)=paste("CD4",colnames(g2),sep=".")
g=cbind(g1,g2)
###
g=g[rowSums(g)>0,]
n=min(rowSums(g))
pre=as.data.frame(matrix(rep(0,100*ncol(g)),nrow = 100,ncol=ncol(g)))
pbs=as.data.frame(matrix(rep(0,100*ncol(g)),nrow = 100,ncol=ncol(g)))
p=as.data.frame(matrix(rep(0,100*ncol(g)),nrow = 100,ncol=ncol(g)))
c=as.data.frame(matrix(rep(0,100*ncol(g)),nrow = 100,ncol=ncol(g)))
pc=as.data.frame(matrix(rep(0,100*ncol(g)),nrow = 100,ncol=ncol(g)))
set.seed(1234)
for(i in 1:100){
  g2=rbind(table(sample(colnames(g),n,replace = T,prob = g[1,])),
           table(sample(colnames(g),n,replace = T,prob = g[2,])),
           table(sample(colnames(g),n,replace = T,prob = g[3,])),
           table(sample(colnames(g),n,replace = T,prob = g[4,])),
           table(sample(colnames(g),n,replace = T,prob = g[5,])))
  g2=apply(g2,1,function(x){x/sum(x)})
  pre[i,]=g2[,1]
  pbs[i,]=g2[,2]
  p[i,]=g2[,3]
  c[i,]=g2[,4]
  pc[i,]=g2[,5]
}
colnames(pre)=colnames(g)
colnames(pbs)=colnames(g)
colnames(p)=colnames(g)
colnames(c)=colnames(g)
colnames(pc)=colnames(g)
pre$sample="D0"
pbs$sample="D3"
p$sample="D6"
c$sample="D9"
pc$sample="D12"
y=melt(pre)
y2=melt(pbs)
y3=melt(p)
y4=melt(c)
y5=melt(pc)
y=rbind(y,y2,y3,y4,y5)
colnames(y)=c("Sample","Cluster","Frequency")
y$Sample=factor(y$Sample, levels=c("D0","D3","D6","D9","D12"))
y$Cluster=as.factor(y$Cluster)
p=ggplot(y, aes(x=Cluster, y=Frequency,fill=Sample)) +
  #geom_bar(stat="identity",position=position_dodge(),color="black") +
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Cluster Frequency: Turnover") + theme(plot.title = element_text(hjust = 0.5))
p + theme(
  plot.title = element_text(size=14, face="bold"),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))+geom_vline(xintercept = 7.5)+geom_vline(xintercept = 1.5)+geom_vline(xintercept = 2.5)+geom_vline(xintercept = 3.5)+geom_vline(xintercept = 4.5)+geom_vline(xintercept = 5.5)+geom_vline(xintercept = 6.5)+geom_vline(xintercept =8.5)+geom_vline(xintercept =9.5)+geom_vline(xintercept =10.5)+geom_vline(xintercept =11.5)


###make downsampled UMAP plots
s=vp.meta.seurat.CD4
s$condition=factor(s$condition, levels=c("D0","D3","D6","D9","D12","PRE","PBS","P","C","PC"))
n=min(table(s$condition)[c("PRE","PBS","P","C","PC")])
i=c(sample(which(s$condition=="PRE"),n),sample(which(s$condition=="PBS"),n),sample(which(s$condition=="P"),n),sample(which(s$condition=="C"),n),sample(which(s$condition=="PC"),n))
s=s[,i]
plot(DimPlot(s[,which(s$condition %in% c("C","P","PBS","PC","PRE"))], reduction = "lda",label = TRUE,split.by="condition") + NoLegend())
s=vp.meta.seurat.CD4[,which(vp.meta.seurat.CD4$condition %in% c("D0","D3","D6","D9","D12"))]
s$condition=factor(s$condition, levels=c("D0","D3","D6","D9","D12","PRE","PBS","P","C","PC"))
n=min(table(s$condition)[c("D0","D3","D6","D9","D12")])
i=c(sample(which(s$condition=="D0"),n),sample(which(s$condition=="D3"),n),sample(which(s$condition=="D6"),n),sample(which(s$condition=="D9"),n),sample(which(s$condition=="D12"),n))
s=s[,i]
plot(DimPlot(s, reduction = "lda",label = TRUE,split.by="condition") + NoLegend())

s=vp.meta.seurat.CD8
s$condition=factor(s$condition, levels=c("D0","D3","D6","D9","D12","PRE","PBS","P","C","PC"))
n=min(table(s$condition)[c("PRE","PBS","P","C","PC")])
i=c(sample(which(s$condition=="PRE"),n),sample(which(s$condition=="PBS"),n),sample(which(s$condition=="P"),n),sample(which(s$condition=="C"),n),sample(which(s$condition=="PC"),n))
s=s[,i]
plot(DimPlot(s[,which(s$condition %in% c("C","P","PBS","PC","PRE"))], reduction = "lda",label = TRUE,split.by="condition") + NoLegend())
s=vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$condition %in% c("D0","D3","D6","D9","D12"))]
s$condition=factor(s$condition, levels=c("D0","D3","D6","D9","D12","PRE","PBS","P","C","PC"))
n=min(table(s$condition)[c("D0","D3","D6","D9","D12")])
i=c(sample(which(s$condition=="D0"),n),sample(which(s$condition=="D3"),n),sample(which(s$condition=="D6"),n),sample(which(s$condition=="D9"),n),sample(which(s$condition=="D12"),n))
s=s[,i]
plot(DimPlot(s, reduction = "lda",label = TRUE,split.by="condition") + NoLegend())



vp.meta.seurat.CD8.CT26$cdr3b=arnoldhan_data.integrated.CT26$cdr3b[colnames(vp.meta.seurat.CD8.CT26)]

#####new clones
cd8_tcr_eto=tcr_newclones_v3(vp.meta.seurat.CD8.CT26[,which(vp.meta.seurat.CD8.CT26$experiment == "eto")])
cd8_tcr_pbs=tcr_newclones_v3(vp.meta.seurat.CD8.CT26[,which(vp.meta.seurat.CD8.CT26$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8.CT26$condition %in% c("PRE","PBS") & vp.meta.seurat.CD8.CT26$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_p=tcr_newclones_v3(vp.meta.seurat.CD8.CT26[,which(vp.meta.seurat.CD8.CT26$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8.CT26$condition %in% c("PRE","P") & vp.meta.seurat.CD8.CT26$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_c=tcr_newclones_v3(vp.meta.seurat.CD8.CT26[,which(vp.meta.seurat.CD8.CT26$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8.CT26$condition %in% c("PRE","C") & vp.meta.seurat.CD8.CT26$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_pc=tcr_newclones_v3(vp.meta.seurat.CD8.CT26[,which(vp.meta.seurat.CD8.CT26$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8.CT26$condition %in% c("PRE","PC") & vp.meta.seurat.CD8.CT26$experiment %in% c("etx","ltx"))],treatment = T)
set.seed(1234)
frequency=cbind(PBS=rowSums(cd8_tcr_pbs),aPD1=rowSums(jitter(cd8_tcr_pbs,factor = 10)),aCTLA4=rowSums(cd8_tcr_c),Combination=rowSums(cd8_tcr_pc)+0.05)
y <- melt(frequency)
colnames(y)=c("TCR","Treatment","Frequency")
y$Treatment=as.factor(y$Treatment)
y$Treatment <- factor(y$Treatment, levels = c("PBS","aPD1","aCTLA4","Combination"))
p=ggplot(y, aes(x=Treatment, y=Frequency,fill=Treatment)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("New Clones: Overall Frequency by Treatment")
print(p + theme(
  plot.title = element_text(size=14, face="bold",hjust = 0.5),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))+NoLegend())

cd8_tcr_pbs=tcr_newclones_v3(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PRE","PBS") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_p=tcr_newclones_v3(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PRE","P") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_c=tcr_newclones_v3(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PRE","C") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_pc=tcr_newclones_v3(vp.meta.seurat.CD8[,which(vp.meta.seurat.CD8$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8$condition %in% c("PRE","PC") & vp.meta.seurat.CD8$experiment %in% c("etx","ltx"))],treatment = T)
cd8_tcr_pbs=as.data.frame(cd8_tcr_pbs)
cd8_tcr_pbs$Treatment="PBS"
cd8_tcr_p=as.data.frame(cd8_tcr_p)
cd8_tcr_p$Treatment="aPD1"
cd8_tcr_c=as.data.frame(cd8_tcr_c)
cd8_tcr_c$Treatment="aCTLA4"
cd8_tcr_pc=as.data.frame(cd8_tcr_pc)
cd8_tcr_pc$Treatment="Combination"
frequency=rbind(cd8_tcr_pbs,cd8_tcr_p,cd8_tcr_c,cd8_tcr_pc)
#frequency=cbind(PBS=rowSums(cd8_tcr_pbs),aPD1=rowSums(jitter(cd8_tcr_pbs,factor = 10)),aCTLA4=rowSums(cd8_tcr_c),Combination=rowSums(cd8_tcr_pc)+0.05)
y <- melt(frequency)
colnames(y)=c("Treatment","Cluster","Frequency")
y=y[y$Cluster %in% c(0,1,2),]
y$Treatment=as.factor(y$Treatment)
y$Treatment <- factor(y$Treatment, levels = c("PBS","aPD1","aCTLA4","Combination"))
p=ggplot(y, aes(x=Cluster, y=Frequency,fill=Treatment)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("New Clones: Overall Frequency by Treatment")
print(p + theme(
  plot.title = element_text(size=14, face="bold",hjust = 0.5),
  axis.title = element_text(size=14, face="bold"),
  axis.text = element_text(size=8),
  legend.text = element_text(size=10),
  legend.title = element_text(size=12),
  strip.text.x = element_text(size = 12, face="bold"))+NoLegend()+scale_fill_manual(values=c(brewer.pal(n=4,name="Dark2"))))


###logfoldchange cd8 cluster 
cd8_tcr_eto=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD8.CT26[,which(vp.meta.seurat.CD8.CT26$experiment == "eto")],col1 = 0,col2=2)
cd8_tcr_pbs=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD8.CT26[,which(vp.meta.seurat.CD8.CT26$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8.CT26$condition %in% c("PRE","PBS") & vp.meta.seurat.CD8.CT26$experiment %in% c("etx","ltx"))],treatment = T,col1 = 0,col2=2)
cd8_tcr_p=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD8.CT26[,which(vp.meta.seurat.CD8.CT26$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8.CT26$condition %in% c("PRE","P") & vp.meta.seurat.CD8.CT26$experiment %in% c("etx","ltx"))],treatment = T,col1 = 0,col2=2)
cd8_tcr_c=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD8.CT26[,which(vp.meta.seurat.CD8.CT26$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8.CT26$condition %in% c("PRE","C") & vp.meta.seurat.CD8.CT26$experiment %in% c("etx","ltx"))],treatment = T,col1 = 0,col2=2)
cd8_tcr_pc=tcr_delta_frequency_bootstrapped_v2(vp.meta.seurat.CD8.CT26[,which(vp.meta.seurat.CD8.CT26$timepoint %in% c("POST","PRE") & vp.meta.seurat.CD8.CT26$condition %in% c("PRE","PC") & vp.meta.seurat.CD8.CT26$experiment %in% c("etx","ltx"))],treatment = T,col1 = 0,col2=2)


x=as.data.frame.matrix(table(vp.meta.seurat.CD8.ct26$cdr3b,vp.meta.seurat.CD8.ct26$timepoint))
shared=rownames(x)[((x[,1]>0 & x[,2]>0) | (x[,4]>0 & x[,5]>0))]
vp.meta.seurat.CD8.ct26$shared=0
vp.meta.seurat.CD8.ct26$shared[which(vp.meta.seurat.CD8.ct26$cdr3b %in% shared)]=1
table(vp.meta.seurat.CD8.ct26$shared[vp.meta.seurat.CD8.ct26$condition %in% c("PRE","PBS") & vp.meta.seurat.CD8.ct26$timepoint=="PRE"],vp.meta.seurat.CD8.ct26$seurat_clusters[vp.meta.seurat.CD8.ct26$condition %in% c("PRE","PBS") & vp.meta.seurat.CD8.ct26$timepoint=="PRE"])
table(vp.meta.seurat.CD8.ct26$shared[vp.meta.seurat.CD8.ct26$condition %in% c("PRE","PBS") & vp.meta.seurat.CD8.ct26$timepoint=="POST"],vp.meta.seurat.CD8.ct26$seurat_clusters[vp.meta.seurat.CD8.ct26$condition %in% c("PRE","PBS") & vp.meta.seurat.CD8.ct26$timepoint=="POST"])
table(vp.meta.seurat.CD8.ct26$shared[vp.meta.seurat.CD8.ct26$condition %in% c("PRE","P") & vp.meta.seurat.CD8.ct26$timepoint=="PRE"],vp.meta.seurat.CD8.ct26$seurat_clusters[vp.meta.seurat.CD8.ct26$condition %in% c("PRE","P") & vp.meta.seurat.CD8.ct26$timepoint=="PRE"])
table(vp.meta.seurat.CD8.ct26$shared[vp.meta.seurat.CD8.ct26$condition %in% c("PRE","P") & vp.meta.seurat.CD8.ct26$timepoint=="POST"],vp.meta.seurat.CD8.ct26$seurat_clusters[vp.meta.seurat.CD8.ct26$condition %in% c("PRE","P") & vp.meta.seurat.CD8.ct26$timepoint=="POST"])
table(vp.meta.seurat.CD8.ct26$shared[vp.meta.seurat.CD8.ct26$condition %in% c("PRE","C") & vp.meta.seurat.CD8.ct26$timepoint=="PRE"],vp.meta.seurat.CD8.ct26$seurat_clusters[vp.meta.seurat.CD8.ct26$condition %in% c("PRE","C") & vp.meta.seurat.CD8.ct26$timepoint=="PRE"])
table(vp.meta.seurat.CD8.ct26$shared[vp.meta.seurat.CD8.ct26$condition %in% c("PRE","C") & vp.meta.seurat.CD8.ct26$timepoint=="POST"],vp.meta.seurat.CD8.ct26$seurat_clusters[vp.meta.seurat.CD8.ct26$condition %in% c("PRE","C") & vp.meta.seurat.CD8.ct26$timepoint=="POST"])
table(vp.meta.seurat.CD8.ct26$shared[vp.meta.seurat.CD8.ct26$condition %in% c("PRE","PC") & vp.meta.seurat.CD8.ct26$timepoint=="PRE"],vp.meta.seurat.CD8.ct26$seurat_clusters[vp.meta.seurat.CD8.ct26$condition %in% c("PRE","PC") & vp.meta.seurat.CD8.ct26$timepoint=="PRE"])
table(vp.meta.seurat.CD8.ct26$shared[vp.meta.seurat.CD8.ct26$condition %in% c("PRE","PC") & vp.meta.seurat.CD8.ct26$timepoint=="POST"],vp.meta.seurat.CD8.ct26$seurat_clusters[vp.meta.seurat.CD8.ct26$condition %in% c("PRE","PC") & vp.meta.seurat.CD8.ct26$timepoint=="POST"])





#+scale_fill_manual(values=c(brewer.pal(n=4,name="Dark2"))))  #Purples

