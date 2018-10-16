setwd('/Volumes/MyBook_5/Salah_Fibro_data/ChIP-seq/H3K4me3_breadth_old_method/')

library('Vennerable')

## 2016-07-29
# make broad domain with aging venn diagrams based on IDR-like calls

### bivalent differences with age
Venn.bivalent.age <- Venn(SetNames = c("Young (3m)", "Old (29m)"), 
                          Weight = c(`10`=141, `11`=473, `01`=118))

pdf("2016-07-29_Broad_H3K4me3_domains_EarFibro_aging_VennDiagram.pdf")
plot(Venn.bivalent.age, show = list(Faces = FALSE))
dev.off()

