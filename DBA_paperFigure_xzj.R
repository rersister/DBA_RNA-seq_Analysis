
####
@Author: xzj
@Date: 20240118
@Time: 
@Description: Figures in paper
####
###fig1F
rna <- read_csv("H:/dba/cellline/rna.rpkm.csv")
ery_gene=c('TF','TFRC','GYPA','GATA1','KLF1','LMO2','NFE2','BCL11A',
           'AHSP','ALAD','ALAS1','ALAS2','BACH1','CPOX','FECH','FLVCR1',
           'GLRX5','IBA57','ISCA1','ISCU','PPOX','UROD','UROS',
           'HBA2','HBB','HBG1','HBG2')
rna=rna[ery_gene,]
pheatmap(rna,cluster_cols = T)
pro <- read_csv("H:/dba/cellline/pro.count.csv")
ery_gene=c('TF','TFRC','GYPA','GATA1','KLF1','LMO2','NFE2','BCL11A',
           'AHSP','ALAD','ALAS1','ALAS2','BACH1','CPOX','FECH','FLVCR1',
           'GLRX5','IBA57','ISCA1','ISCU','PPOX','UROD','UROS',
           'HBA2','HBB','HBG1','HBG2')
pro=pro[ery_gene,]
pheatmap(pro,cluster_cols = T)


#########fig3a bone marrow single cell map
library(Seurat)
load("H:/dba/scrna/bm.RData")
DimPlot(bm,label = T)
#########fig3b gsva
library(testSctpa)
sub <- bm[, (sample(colnames(bm), size =10000, replace=F))]
DefaultAssay(sub)='RNA'
sub= cal_PAS(seurat_object = sub,
             tool = 'AUCell',   ## GSVA, ssGSEA, plage, zscore or Vision
             normalize = 'log',
             species = 'human', 
             pathway='GO.bp')
sub[['go']]=sub[['PAS']]

sub1=subset(sub,features=c('GO-HEMOGLOBIN-METABOLIC-PROCESS','GO-RIBOSOME-BIOGENESIS','GO-RIBOSOME-BIOGENESIS','GO-ERYTHROCYTE-HOMEOSTASIS',
                           'GO-TRANSLATIONAL-INITIATION','GO-OXIDATIVE-PHOSPHORYLATION','GO-RIBOSOME-ASSEMBLY','GO-RRNA-TRANSCRIPTION',
                           'GO-RIBOSOMAL-LARGE-SUBUNIT-ASSEMBLY','GO-RIBOSOMAL-SMALL-SUBUNIT-BIOGENESIS','GO-MATURATION-OF-LSU-RRNA',
                           'GO-GAS-TRANSPORT','GO-ELECTRON-TRANSPORT-CHAIN','GO-OXYGEN-TRANSPORT','GO-HEME-BIOSYNTHETIC-PROCESS',
                           'GO-RIBOSOMAL-LARGE-SUBUNIT-BIOGENESIS','GO-RIBOSOMAL-SMALL-SUBUNIT-ASSEMBLY',
                           'GO-ERYTHROCYTE-DEVELOPMENT','GO-RIBOSOMAL-LARGE-SUBUNIT-EXPORT-FROM-NUCLEUS',
                           'GO-RIBOSOMAL-SMALL-SUBUNIT-EXPORT-FROM-NUCLEUS','GO-MITOCHONDRIAL-ELECTRON-TRANSPORT-NADH-TO-UBIQUINONE',
                           'GO-PROTOPORPHYRINOGEN-IX-BIOSYNTHETIC-PROCESS','GO-NADH-DEHYDROGENASE-COMPLEX-ASSEMBLY',
                           'GO-ADP-TRANSPORT','GO-ATP-TRANSPORT','GO-RESPIRATORY-ELECTRON-TRANSPORT-CHAIN',
                           'GO-ATP-SYNTHESIS-COUPLED-ELECTRON-TRANSPORT','GO-PEPTIDE-BIOSYNTHETIC-PROCESS',
                           'GO-POSITIVE-REGULATION-OF-RIBOSOME-BIOGENESIS','GO-SEQUESTERING-OF-IRON-ION'))
avg <- AverageExpression(sub1,return.seurat =F)
avg.mat=avg[["PAS"]][,c('T.cell','B.cell','NK','DC','Basophil','Eosinophil','Neutrophil','Macrophage','Monocyte','Megakaryocyte','GMP','CMP.LMPP','HSC',
                       'MEP','BFU-E','CFU-E','ProE','BasoE','PolyE','OrthoE')]
pheatmap::pheatmap(avg.mat,scale = 'row',cluster_cols = F,cluster_rows = F)
#########fig3c module score
DefaultAssay(sub)='RNA'
gene.list=list(RP=c('RPL21','RPL17','RPS10','RPS29','RPL36A','RPS17','RPS28','RPS26','RPL5','RPL13A','RPL11','RPS3A','RPS27A','RPS7','RPL9',
               'RPS4X','RPL38','RPS27','RPLP0','RPL15','RPS24','RPL37A','RPS21','RPL22','RPL23','RPL37','RPL34','RPL35A','RPL27A','RPL31','RPL36',
               'RPL41','RPL39','RPLP1','RPL27','RPS15','RPS9','RPL8','RPS19','RPL28','RPLP2','RACK1','RPL12','RPS18','RPL18A','RPL18','RPL13','RPS16','RPL30',
               'RPS11','RPS14','RPS15A','RPL19','RPL32','RPL23A','RPS25','RPL29','RPS12','RPL26','RPL10A','RPS6','RPS5','RPL35','RPS2','RPL14','RPL6','RPS3','RPSA',
               'RPL7A','RPS8','RPL3','RPL4','RPL10','RPL24','RPS20','RPL7','RPS13','RPS23'),
               oxphos=c('NDUFS8','ATP5G2','ATP5H','NDUFB11','UQCRB','NDUFA11','ATP5D','NDUFS7','COX4I1','NDUFA12','CYC1','NDUFS6',
                        'NDUFAB1','COX6C','ATP5O','NDUFA6','NDUFA8','ATP5G1','COX5A','COX7C','NDUFA4','COX6B1','NDUFA1','COX8A',
                        'NDUFB9','NDUFS4','NDUFS5','SDHB','COX5B','NDUFA2','UQCRC1','COX7B','NDUFB4','COX7A2','UQCRH','UQCR10','NDUFB1','COX7A2L','NDUFB3'))
bm <- AddModuleScore(object = bm, features = gene.list)
p1=VlnPlot(bm,features = c('Cluster1'))
colnames(p1$data)=c('score','ct')
p1$data$term='RP'
p2=VlnPlot(bm,features = c('Cluster2'))
colnames(p2$data)=c('score','ct')
p2$data$term='oxphos'
dat=rbind(p1$data,p2$data)
dat1=subset(dat, ct %in% c('HSC','CMP.LMPP','MEP','BFU-E','CFU-E','ProE','BasoE','PolyE','OrthoE'))
term_mean=aggregate(dat1$score, by=list(dat1$ct,dat1$term), FUN=mean)
colnames(term_mean)=c('ct','term','mean')
term_mean$ct=factor(term_mean$ct,levels=c('HSC','CMP.LMPP','MEP','BFU-E','CFU-E','ProE','BasoE','PolyE','OrthoE'))
table(term_mean$ct)
ggplot(data = term_mean,aes(x=ct,y=mean,group = term,color=term,shape=term))+
  geom_point()+
  geom_line()+
  xlab("Year")+#横坐标名称
  ylab("数值")+#纵坐标名称
  theme_bw() +#去掉背景灰色
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),#以上theme中代码用于去除网格线且保留坐标轴边框
        text = element_text(family = "STXihei"),#设置中文字体的显示
        legend.position = c(.075,.915),#更改图例的位置，放至图内部的左上角
        legend.box.background = element_rect(color="black"))#+为图例田间边框线 scale_x_continuous(limits = c(2000,2018),breaks = seq(2000,2018,1))#更改横坐标刻度值

#########fig3d Correlation between oxidative phosphorylation and the ribosome pathway
library(ggplot2)
library(ggpubr)
library(IDPmisc)
library(MASS)
library(viridis)
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
DefaultAssay(sub)='RNA'
sub= cal_PAS(seurat_object = sub,
             tool = 'AUCell',   ## GSVA, ssGSEA, plage, zscore or Vision
             normalize = 'log',
             species = 'human', 
             pathway='kegg')
sub[['kegg']]=sub[['PAS']]
term <- FetchData(sub, vars=c('Ribosome','Oxidative-phosphorylation'))
sub <- AddMetaData(sub, metadata = term)
sub$density <- get_density(sub$Ribosome, sub$Oxidative.phosphorylation, n = 100)
meta=sub@meta.data
ggplot(meta,aes(x=meta$Ribosome, y=meta$Oxidative.phosphorylation,color = meta$density)) + 
  geom_point() + 
  scale_color_viridis()+stat_cor(data=meta[c('Ribosome','Oxidative.phosphorylation')], method = "pearson")+
  stat_smooth(method="lm",se=T)+
  theme_bw()+theme(panel.grid=element_blank())



#########fig4a  figs3b
require(tidyverse)
require(ggplot2)
require(ggsci)
count <- read.csv("H:/dba/patient/rna.rpkm.csv")
colnames(count)=c('gene','DBA1','DBA2','DBA3','DBA4','DBA5','DBA6','DBA7','DBA8','DBA9','DBA10','ITP1','ITP2','ITP3')
count1=as.data.frame(subset(count,gene %in% c('RANGAP1')))
df=melt(count1)
df$group[df$variable %in% c('DBA1','DBA2','DBA3','DBA4','DBA5','DBA6','DBA7','DBA8','DBA9','DBA10')]='DBA'
df$group[df$variable %in% c('ITP1','ITP2','ITP3')]='ITP'

df %>% 
  ggplot(aes(x = gene, y = value, fill = group)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "Count")+
  scale_x_discrete(name = "gene") +
  theme_bw()+ggpubr::stat_compare_means(comparisons = list(c("DBA", "ITP")), label = "p.signif")+
  theme(plot.title = element_text(size = 14, face =  "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11))+
  scale_fill_lancet()




count1=as.data.frame(subset(count,gene %in% c('RAN')))
df=melt(count1)
df$group[df$variable %in% c('DBA1','DBA2','DBA3','DBA4','DBA5','DBA6','DBA7','DBA8','DBA9','DBA10')]='DBA'
df$group[df$variable %in% c('ITP1','ITP2','ITP3')]='ITP'

df %>% 
  ggplot(aes(x = gene, y = value, fill = group)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "Count")+
  scale_x_discrete(name = "gene") +
  theme_bw()+ggpubr::stat_compare_means(comparisons = list(c("DBA", "ITP")), label = "p.signif")+
  theme(plot.title = element_text(size = 14, face =  "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11))+
  scale_fill_lancet()


#########fig6a
setwd('F:/DBA/data')
fig1a=read.csv(file = 'fig1a_keggterm.csv')
paletteFunc <- colorRampPalette(c("red", "blue"));
palette     <- paletteFunc(12)
barplot(foxo3_target_go$`P-value`, col=palette)
a=data.frame(pathway=fig1a$Ingenuity.Canonical.Pathways,FDR=fig1a$X.log.p.value.)
a$pathway <- factor(x =a$pathway, levels = c('Mitochondrial Dysfunction','Oxidative Phosphorylation','EIF2 Signaling'
                                             ,'Regulation of eIF4 and p70S6K Signaling','mTOR Signaling','Heme Biosynthesis II'
                                             ,'Glutathione Redox Reactions I','Nucleotide Excision Repair Pathway','Tetrapyrrole Biosynthesis II',
                                             'Adenine and Adenosine Salvage I','Ascorbate Recycling (Cytosolic)','Oxidized GTP and dGTP Detoxification'))
f <- ggplot(foxo3_target_go, aes(x = log10, y =Name))
f + geom_bar(stat = "identity", aes(fill = Name))+
  theme_set(theme_bw())+scale_fill_manual(values=palette)+
  theme(axis.title.x = element_text(size = 14),axis.title.y = element_text(size = 14),axis.text.x = element_text(size=12),axis.text.y = element_text(size=12))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+###去网格线
  theme(axis.text.x = element_text(size = 14, family = "myFont", face = "bold", vjust = 0.5, hjust = 0.5, angle = 45))+ ###改坐标标题
  theme(axis.line = element_line(colour = "black"))###加坐标轴

####fig6b
count <- read.csv("H:/dba/patient/rna.rpkm.csv")
colnames(count)=c('gene','DBA1','DBA2','DBA3','DBA4','DBA5','DBA6','DBA7','DBA8','DBA9','DBA10','ITP1','ITP2','ITP3')
pheatmap::pheatmap(count[gene.list$oxphos,],scale = "row", clustering_distance_rows = "correlation",
                   clustering_method = "average",color = colorRampPalette(c("#63B8FF", "white", "#FF8C69"))(50))

###fig6c
count <- read.csv("H:/dba/patient/rna.rpkm.csv")
colnames(count)=c('gene','DBA1','DBA2','DBA3','DBA4','DBA5','DBA6','DBA7','DBA8','DBA9','DBA10','ITP1','ITP2','ITP3')
pheatmap::pheatmap(count[gene.list$RP,],scale = "row", clustering_distance_rows = "correlation",
                   clustering_method = "average",color = colorRampPalette(c("#63B8FF", "white", "#FF8C69"))(50))
