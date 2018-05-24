# 2. 
library(genomation)
library(GenomicRanges)
library(dplyr)
chip.files = list.files('./Data/ChIP', full.names=TRUE)

lchip = lapply(chip.files, readGeneric, zero.based=TRUE)
lnams = basename(chip.files)
lnams = sub('_GRCh38.bed.gz','',lnams)
lnams = sub('K562_','',lnams)
names(lchip) = lnams
lchip = GRangesList(lchip)
lchip = reduce(lchip)

uchip = unlist(lchip)

# intersection
ovlap      = as.data.frame(findOverlaps(uchip, lchip))
ovlap$set1 = names(uchip)[ovlap$queryHits]
ovlap$set2 = names(lchip)[ovlap$subjectHits]
inter      = ovlap %>%
    group_by(set1, set2) %>%
    summarize(intersection = length(unique(queryHits))) 
inter$key = with(inter, paste(set1, set2))

# union
union            = expand.grid(seq(lchip), seq(lchip))
union$set1       = names(lchip)[union$Var1]
union$set2       = names(lchip)[union$Var2]
union$length1    = elementNROWS(lchip)[union$Var1]
union$length2    = elementNROWS(lchip)[union$Var2]
union$sum_length = with(union, (length1 + length2))
union$key = with(union, paste(set1, set2))
union = union[,-c(1:4)]

# jaccard
dist = merge(inter, union, by='key')
dist$perc  = with(dist, round(intersection/length1,2))
dist$union = with(dist, sum_length - intersection)
dist$jacc  = with(dist, intersection/union)

sum(jacc$jacc > 1)
summary(jacc$jacc)
ss = subset(jacc, jacc>1)
# heatmap
mat.perc = data.table::dcast(dist, set1~set2, value.var='perc') %>%
    mutate(set1 = NULL) 
rownames(mat.perc) = colnames(mat.perc)


library(ComplexHeatmap)
pdf('PercentageHeatmap.pdf', width = 60, height=60)
    Heatmap(mat.perc)
dev.off()


