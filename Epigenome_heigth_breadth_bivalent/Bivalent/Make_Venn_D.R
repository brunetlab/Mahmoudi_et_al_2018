setwd('/Volumes/MyBook_5/Salah_Fibro_data/ChIP-seq/GoodvsBad/Bivalent/')

library('Vennerable')

## 2016-06-07
# make bivalent domain with aging venn diagrams based on IDR-like calls

### bivalent differences with age
Venn.bivalent.age <- Venn(SetNames = c("Young (3m)", "Old (29m)"), 
                          Weight = c(`10`=1349-886, `11`=886, `01`=1516-886))

pdf("2016-06-07_Bivalent_domains_EarFibro_aging_VennDiagram.pdf")
plot(Venn.bivalent.age, show = list(Faces = FALSE))
dev.off()


### bivalent prop with age
Venn.bivalent.YOUNG <- Venn(SetNames = c("Young H3K4me3", "Young H3K27me3"), 
                          Weight = c(`10`=16852-1349, `11`=1349, `01`=3666-1349))

pdf("2016-06-07_Bivalent_domains_EarFibro_YOung3m_VennDiagram.pdf")
plot(Venn.bivalent.YOUNG, show = list(Faces = FALSE))
dev.off()

### bivalent prop with age
Venn.bivalent.OLD <- Venn(SetNames = c("Old H3K4me3", "Old H3K27me3"), 
                            Weight = c(`10`=16461-1516, `11`=1516, `01`=3624-1516))

pdf("2016-06-07_Bivalent_domains_EarFibro_Old29m_VennDiagram.pdf")
plot(Venn.bivalent.OLD, show = list(Faces = FALSE))
dev.off()
