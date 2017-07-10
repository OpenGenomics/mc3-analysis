library(UpSetR)
library(data.table)

data <- fread("Data/mc3.v0.2.8.upsetR.txt")
pdf(file="Figures/maf.v0.2.8.UpSetR.pdf",height=11,width=22)
upset(data, sets = c("INDELOCATOR", "MUSE", "MUTECT", "PINDEL", "RADIA", "SOMATICSNIPER","VARSCANI","VARSCANS"), sets.bar.color = "#BE312D",order.by = "freq", empty.intersections = "on",point.size = 5, line.size = 2)
dev.off()
#dev.copy2pdf(file="Figures/maf.v0.2.8.UpSetR.pdf",height=11,width=22)
