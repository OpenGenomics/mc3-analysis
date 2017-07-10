#######################INSTALL LIBRARIES######################
library(UpSetR)
library(data.table)
library(ggplot2)
library(reshape2)
library(plyr)
library(viridis)



#######################MAKE THE UPSETR FIGURE################
data <- fread("mc3.v0.2.8.CallerMatrix_v3.txt")
upset(data, sets = c("INDELOCATOR", "MUSE", "MUTECT", "PINDEL", "RADIA", "SOMATICSNIPER","VARSCANI","VARSCANS"), sets.bar.color = "#BE312D",order.by = "freq", empty.intersections = "on",point.size = 5, name.size = 12,line.size = 2)
dev.copy2pdf(file="maf.v0.2.8.UpSetR.v4.pdf",height=11,width=22)
##############################################################



######################MAKE THE MSI-MMR FIGURE#################
mc3 <- fread("mc3.v0.2.8.selected.dat")
tcga <- unique(mc3$Tumor_Sample_Barcode)

missense <- mc3[which(mc3$Variant_Classification  == "Missense_Mutation"),]
mmrgenes = c("POLE","MLH1","MLH3","MGMT","MSH6","MSH3","MSH2","PMS1","PMS2")

cancers = c("UVM","PCPG","PRAD","LGG","THCA","ACC","MESO","TGCT","UCS","GBM","KICH","HNSC","KIRC","LAML","STAD","KIRP","BLCA","BRCA","SARC","LIHC","THYM","CESC","OV","CHOL","LUAD","UCEC","DLBC","LUSC","ESCA","PAAD","SKCM","READ","COAD")


missense$MSIMMR = ifelse(missense$Hugo_Symbol %in% mmrgenes,1,0)
#Adjust this on be PASS MISMMR only 
#missense$MSIMMR = ifelse(missense$MSIMMR == 1 & missense$FILTER == "PASS",1,0)
msisamps <- missense[which(missense$MSIMMR == 1),]
msamp <- unique(msisamps$Tumor_Sample_Barcode)

yo$MSIMMR = ifelse(yo$x %in% msamp,1,0)

p <- ggplot(yo,aes(x=factor(yo$MSIMMR),y=freq))
p <- p + geom_boxplot()
p <- p + scale_y_log10()
p <- p + ylab("Log10 Mutation Frequency per sample (Missense only)")
p <- p + xlab("MSI/MMR Status      0=No 1=Yes")
p <- p + theme_minimal()
#p

#ALL VARIANTS###
mc3 <- missense
tcga <- unique(mc3$Tumor_Sample_Barcode)
tabdat <- table(mc3$Tumor_Sample_Barcode,mc3$CENTER)
caller <- c("MUTECT","MUSE","VARSCANS","SOMATICSNIPER","RADIA")
tabdat2 <- tabdat[,grepl("MUTECT|MUSE|VARSCANS|SOMATICSNIPER|RADIA", colnames(tabdat))]
tabdat <- tabdat2
sampCallerCnts <- data.frame(row.names=row.names(tabdat))
for(i in caller){
    tmp <- data.frame(callr=rowSums(tabdat[,which(grepl(i,strsplit(colnames(tabdat),"\\|")))]))
    sampCallerCnts <- cbind(sampCallerCnts,tmp)
}
colnames(sampCallerCnts) = caller
dim(sampCallerCnts)
sampCallerCnts$total = rowSums(tabdat)
cantyp <- unique(data.frame(tumor=mc3$Tumor_Sample_Barcode,cancer=mc3$CODE))
sampCallerCnts$CT = cantyp$cancer
new.sampCallerCnts <- sampCallerCnts[order(sampCallerCnts$total),]
new.sampCallerCnts$Tumor_Sample_Barcode = row.names(new.sampCallerCnts)
new.sampCallerCnts$MSIMMR = ifelse(new.sampCallerCnts$Tumor_Sample_Barcode %in% msamp,1,0)

msc <- melt(new.sampCallerCnts,id=c("Tumor_Sample_Barcode","CT","MSIMMR"))

examples <- c("BRCA","COAD","GBM","SKCM")
canmsc <- msc[which(msc$CT %in% examples),]
titlect <- "All Variants"

p <- ggplot(canmsc,aes(x=variable,y=value,color=factor(MSIMMR)))
p <- p + geom_boxplot()
p <- p + scale_y_log10()
p <- p + ylab("Log10 Mutation Frequency per sample (Missense only)")
p <- p + xlab("")
p <- p + theme_minimal()
p <- p + theme(legend.position="bottom",axis.text.x = element_text(angle = 90, hjust = 1))
p <- p + scale_color_discrete(name="MMR_MSI",labels=c("No","Yes"))
p <- p + facet_grid(.~CT)
p

pdf("AllVariants_MSI_MMR.pdf",height = 7,width = 22)
print(p)
dev.off()

#Pass only Variants
mc3 <- missense[which(missense$FILTER == "PASS"),]
tcga <- unique(mc3$Tumor_Sample_Barcode)
tabdat <- table(mc3$Tumor_Sample_Barcode,mc3$CENTER)
caller <- c("MUTECT","MUSE","VARSCANS","SOMATICSNIPER","RADIA")
tabdat2 <- tabdat[,grepl("MUTECT|MUSE|VARSCANS|SOMATICSNIPER|RADIA", colnames(tabdat))]
tabdat <- tabdat2
sampCallerCnts <- data.frame(row.names=row.names(tabdat))
for(i in caller){
    tmp <- data.frame(callr=rowSums(tabdat[,which(grepl(i,strsplit(colnames(tabdat),"\\|")))]))
    sampCallerCnts <- cbind(sampCallerCnts,tmp)
}
colnames(sampCallerCnts) = caller
dim(sampCallerCnts)
sampCallerCnts$total = rowSums(tabdat)
cantyp <- unique(data.frame(tumor=mc3$Tumor_Sample_Barcode,cancer=mc3$CODE))
sampCallerCnts$CT = cantyp$cancer
new.sampCallerCnts <- sampCallerCnts[order(sampCallerCnts$total),]
new.sampCallerCnts$Tumor_Sample_Barcode = row.names(new.sampCallerCnts)
new.sampCallerCnts$MSIMMR = ifelse(new.sampCallerCnts$Tumor_Sample_Barcode %in% msamp,1,0)


msc <- melt(new.sampCallerCnts,id=c("Tumor_Sample_Barcode","CT","MSIMMR"))
examples <- c("BRCA","COAD","GBM","SKCM")
canmsc <- msc[which(msc$CT %in% examples),]
titlect <- "PASS only Variants"

p <- ggplot(canmsc,aes(x=variable,y=value,color=factor(MSIMMR)))
p <- p + geom_boxplot()
p <- p + scale_y_log10()
p <- p + ylab("Log10 Mutation Frequency per sample (Missense only)")
p <- p + xlab("")
p <- p + theme_minimal()
p <- p + theme(legend.position="bottom",axis.text.x = element_text(angle = 90, hjust = 1))
p <- p + scale_color_discrete(name="MMR_MSI",labels=c("No","Yes"))
p <- p + facet_grid(.~CT)
p

pdf("PASSonly_MSI_MMR.pdf",height = 7,width = 22)
print(p)
dev.off()
##############################################################




#######################PURITY FIGURE#########################
dat <- fread("selected_absolute.txt")
mc3 <- unique(dat)

p <- ggplot(dat=mc3,aes(x=reorder(CODE, purity, FUN=median),y=purity))
p <- p + geom_boxplot()
p <- p + xlab("Cancer Type")
p <- p + ylab("ABSOLUTE Purity")
p <- p + coord_flip()
p <- p + theme_minimal()
p

pdf("Absolute_Purity.pdf",height=4,width=6)
print(p)
dev.off()
##############################################################





######################HAMBURGER###############################
###############INDELS#############
mc3 <- fread("mc3.v0.2.8.selected_v2.dat")
indels = c("INS","DEL")
idel <- mc3[which(mc3$Variant_Type %in% indels),] # Limit to SNV
mc3 <- idel

tcga <- unique(mc3$Tumor_Sample_Barcode)
tabdat <- table(mc3$Tumor_Sample_Barcode,mc3$CENTER)
caller <- c("PINDEL","INDELOCATOR","VARSCANI")
tabdat2 <- tabdat[,grepl("PINDEL|INDELOCATOR|VARSCANI", colnames(tabdat))]
tabdat <- tabdat2

#Get the SNV callers only
sampCallerCnts <- data.frame(row.names=row.names(tabdat))
for(i in caller){
    tmp <- data.frame(callr=rowSums(tabdat[,which(grepl(i,strsplit(colnames(tabdat),"\\|")))]))
    sampCallerCnts <- cbind(sampCallerCnts,tmp)
}
colnames(sampCallerCnts) = caller
dim(sampCallerCnts)
sampCallerCnts$total = rowSums(tabdat)
cantyp <- unique(data.frame(tumor=mc3$Tumor_Sample_Barcode,cancer=mc3$CODE))
sampCallerCnts$CT = cantyp$cancer
new.sampCallerCnts <- sampCallerCnts[order(sampCallerCnts$total),]
new.sampCallerCnts$rn = row.names(new.sampCallerCnts)

require(dplyr)
temp <- data.frame(total=new.sampCallerCnts$total,CT=new.sampCallerCnts$CT,rn=row.names(new.sampCallerCnts))

sorted <- temp %>%
          arrange(CT, total) %>%
          group_by(CT) %>%
          mutate(rank=row_number())
s <- data.frame(sorted)


t.new.samps <- merge(new.sampCallerCnts, s, by="rn")
this.df <- t.new.samps[c(1,3,4,2,5,6,9)]
mt <- melt(this.df, id.vars = c("CT.x","rank","rn"))
#DEFINE THE correct colors for this output
VictorianCols=c("#B1A2A7","#C9A784","#8C7851","#5DA5A1","#C45331","#E7960A","#F6E849","#016384","#D8CDB7","#086453","#F7D87B")
INDEL_Cols=c("#B1A2A7","#C9A784","#8C7851","#D8CDB7")
myColors <- INDEL_Cols
SNV_Cols=c("#5DA5A1","#C45331","#E7960A","#F6E849","#016384","#D8CDB7","#086453","#F7D87B")

mt$facet_ord = factor(mt$CT.x, levels=c("UVM","PCPG","PRAD","LGG","THCA","ACC","MESO","TGCT","UCS","GBM","KICH","HNSC","KIRC","LAML","STAD","KIRP","BLCA","BRCA","SARC","LIHC","THYM","CESC","OV","CHOL","LUAD","UCEC","DLBC","LUSC","ESCA","PAAD","SKCM","READ","COAD"))

mt_indel <- mt
#INDELS Gram
p <- ggplot(mt_indel,aes(x=rank,y=value,color=variable))
p <- p + geom_point(size=.25)
p <- p + facet_grid(facet_ord~.,scales="free_y")
p <- p + theme_bw() + theme( strip.background  = element_blank(),
        panel.grid.major = element_line(colour = "grey80"),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y = element_line(color="black"), 
        plot.title = element_text(angle = 100) )
p <- p + scale_y_log10(breaks = c(0,1,10,100,1000,10000,100000,1000000),limits = c(1, 1000000))
p <- p + ylab("Number of 'Called' Mutations per Sample")
p <- p + scale_colour_manual(values=INDEL_Cols)
p <- p + guides(colour = guide_legend(override.aes = list(size=3)))
p <- p + labs(colour = "CALLER")
p <- p + coord_flip()
p
indel_p <- p

#SAVE AS A PDF
pdf("Vertical_INDEL.pdf",width = 5,height=18)
print(p)
dev.off()

#SAVE AS A SVG
ggsave(file="Vertical_INDEL.svg", plot=indel_p,width = 5,height=18,useDingbats=FALSE)


###########SNVS###############3
mc3 <- fread("mc3.v0.2.8.selected_v2.dat")
snp <- mc3[which(mc3$Variant_Type == "SNP"),] # Limit to SNV
mc3 <- snp

tcga <- unique(mc3$Tumor_Sample_Barcode)
tabdat <- table(mc3$Tumor_Sample_Barcode,mc3$CENTER)
caller <- c("MUTECT","MUSE","VARSCANS","SOMATICSNIPER","RADIA")
tabdat2 <- tabdat[,grepl("MUTECT|MUSE|VARSCANS|SOMATICSNIPER|RADIA", colnames(tabdat))]
tabdat <- tabdat2

#Get the SNV callers only
sampCallerCnts <- data.frame(row.names=row.names(tabdat))
for(i in caller){
    tmp <- data.frame(callr=rowSums(tabdat[,which(grepl(i,strsplit(colnames(tabdat),"\\|")))]))
    sampCallerCnts <- cbind(sampCallerCnts,tmp)
}
colnames(sampCallerCnts) = caller
dim(sampCallerCnts)
sampCallerCnts$total = rowSums(tabdat)
cantyp <- unique(data.frame(tumor=mc3$Tumor_Sample_Barcode,cancer=mc3$CODE))
sampCallerCnts$CT = cantyp$cancer
new.sampCallerCnts <- sampCallerCnts[order(sampCallerCnts$total),]
new.sampCallerCnts$rn = row.names(new.sampCallerCnts)

#Rank the total mutations by Sample/Cancertype
require(dplyr)
temp <- data.frame(total=new.sampCallerCnts$total,CT=new.sampCallerCnts$CT,rn=row.names(new.sampCallerCnts))

sorted <- temp %>%
          arrange(CT, total) %>%
          group_by(CT) %>%
          mutate(rank=row_number())
s <- data.frame(sorted)


t.new.samps <- merge(new.sampCallerCnts, s, by="rn")
this.df <- t.new.samps[c(1,6,5,4,3,2,7,8,11)]
mt <- melt(this.df, id.vars = c("CT.x","rank","rn"))

#DEFINE THE correct colors for this output
myColors <- c("#5DA5A1","#C45331","#E7960A","#F6E849","#016384","#D8CDB7")
SNV_Cols=c("#5DA5A1","#C45331","#E7960A","#F6E849","#016384","#D8CDB7","#086453","#F7D87B")
names(myColors) <- factor(unique(plop$call)[5:1])
names(myColors)[6] <- "total.x"

mt$facet_ord = factor(mt$CT.x, levels=c("UVM","PCPG","PRAD","LGG","THCA","ACC","MESO","TGCT","UCS","GBM","KICH","HNSC","KIRC","LAML","STAD","KIRP","BLCA","BRCA","SARC","LIHC","THYM","CESC","OV","CHOL","LUAD","UCEC","DLBC","LUSC","ESCA","PAAD","SKCM","READ","COAD"))

mt_snv = mt


p <- ggplot(mt_snv,aes(x=rank,y=value,color=variable))
p <- p + geom_point(size=.25)
p <- p + facet_grid(facet_ord~.,scales="free_y")
p <- p + theme_bw() + theme( strip.background  = element_blank(),
        panel.grid.major = element_line(colour = "grey80"),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y = element_line(color="black"), 
        plot.title = element_text(angle = 100))
p <- p + scale_y_log10(breaks = c(0,1,10,100,1000,10000,100000,1000000),limits = c(1, 1000000))
p <- p + ylab("Number of 'Called' Mutations per Sample")
p <- p + scale_colour_manual(values=SNV_Cols)
p <- p + guides(colour = guide_legend(override.aes = list(size=3)))
p <- p + labs(colour = "CALLER")
p <- p + coord_flip()
p

snv_p <- p

#SAVE AS A PDF
pdf("Vertical_SNV.pdf",width=5,height=18)
print(p)
dev.off()

#SAVE AS A SVG
ggsave(file="Vertical_SNV.svg", plot=snv_p,width = 5,height=18, useDingbats=FALSE)
#############################################################



###############################Disparate mutation count###############
#TOP50 - KANDOTH GENES
dat <- read.table(file="test.d.targeted.KANDOTH.txt",header=FALSE,sep="\t")
colnames(dat) <- c("GENE","CONTROL","PASS","DIFF","GENE_LENGTH","WEIGHTED_DIFF","TARGETED","DIFF_GTarget","DIFF_TPass")
pd <- c("PASS","DIFF_TPass")
dat_o <- dat[order(dat$DIFF_TPass),]
dat_oe <- tail(dat_o,50) #extreme values
md <- melt(dat_oe)
mymd <- md[which(md$variable %in% pd),]
mymd$porder = (factor(mymd$GENE,dat_oe[order(dat_oe$GENE_LENGTH),]$GENE))

p <- ggplot(mymd,aes(x=porder,y=value,fill=variable))
p <- p + geom_bar(stat="identity",position="stack")
p <- p + theme_minimal()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5),
        axis.text=element_text(size=6) )
p <- p + scale_fill_manual(values=viridis(3),name = "",labels = c("Pass", "Not Pass"))
p <- p + xlab("")
p <- p + ylab("")
p

pdf("PASS_TARGET_Kandoth.pdf",width=10,height=4)
print(p)
dev.off()

#TOP50 Disparate mutations counts 
dat <- read.table(file="test.d.targeted.ALL.txt",header=FALSE,sep="\t")
colnames(dat) <- c("GENE","CONTROL","PASS","DIFF","GENE_LENGTH","WEIGHTED_DIFF","TARGETED","DIFF_GTarget","DIFF_TPass")
pd <- c("PASS","DIFF_TPass")
dat_o <- dat[order(dat$DIFF),]
dat_oc <- dat_o[which(dat_o$GENE_LENGTH != 'NaN'),]
dat_oce <- tail(dat_oc,50) #extreme values
mdo <- melt(dat_oce)
mymdo <- mdo[which(mdo$variable %in% pd),]
mymdo$porder = (factor(mymdo$GENE,dat_oce[order(dat_oce$GENE_LENGTH),]$GENE))


p <- ggplot(mymdo,aes(x=porder,y=value,fill=variable))
p <- p + geom_bar(stat="identity",position="stack")
p <- p + theme_minimal()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5),
        axis.text=element_text(size=6) )
p <- p + scale_fill_manual(values=viridis(3),name = "",labels = c("Pass", "Not Pass"))
p <- p + xlab("")
p <- p + ylab("")
p

pdf("PASS_TARGETED_MostDiparate.pdf",width=10,height=4)
print(p)
dev.off()

######################################################################












