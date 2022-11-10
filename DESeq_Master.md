```R
#load necessary packages and functions 
library(DESeq2)
library(ggplot2)
library(tidyverse)

source("RNAseqfunctions.R")
#functions loaded printed below
cat(paste0(readLines("RNAseqfunctions.R"), collapse="\n"))
```

    Loading required package: S4Vectors
    
    Loading required package: stats4
    
    Loading required package: BiocGenerics
    
    Loading required package: parallel
    
    
    Attaching package: ‘BiocGenerics’
    
    
    The following objects are masked from ‘package:parallel’:
    
        clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
        clusterExport, clusterMap, parApply, parCapply, parLapply,
        parLapplyLB, parRapply, parSapply, parSapplyLB
    
    
    The following objects are masked from ‘package:stats’:
    
        IQR, mad, sd, var, xtabs
    
    
    The following objects are masked from ‘package:base’:
    
        anyDuplicated, append, as.data.frame, basename, cbind, colnames,
        dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
        grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
        order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
        rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
        union, unique, unsplit, which.max, which.min
    
    
    
    Attaching package: ‘S4Vectors’
    
    
    The following objects are masked from ‘package:base’:
    
        expand.grid, I, unname
    
    
    Loading required package: IRanges
    
    Loading required package: GenomicRanges
    
    Loading required package: GenomeInfoDb
    
    Loading required package: SummarizedExperiment
    
    Loading required package: MatrixGenerics
    
    Loading required package: matrixStats
    
    
    Attaching package: ‘MatrixGenerics’
    
    
    The following objects are masked from ‘package:matrixStats’:
    
        colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
        colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
        colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
        colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
        colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
        colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
        colWeightedMeans, colWeightedMedians, colWeightedSds,
        colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
        rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
        rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
        rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
        rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
        rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
        rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
        rowWeightedSds, rowWeightedVars
    
    
    Loading required package: Biobase
    
    Welcome to Bioconductor
    
        Vignettes contain introductory material; view with
        'browseVignettes()'. To cite Bioconductor, see
        'citation("Biobase")', and for packages 'citation("pkgname")'.
    
    
    
    Attaching package: ‘Biobase’
    
    
    The following object is masked from ‘package:MatrixGenerics’:
    
        rowMedians
    
    
    The following objects are masked from ‘package:matrixStats’:
    
        anyMissing, rowMedians
    
    
    ── [1mAttaching packages[22m ─────────────────────────────────────── tidyverse 1.3.1 ──
    
    [32m✔[39m [34mtibble [39m 3.1.2     [32m✔[39m [34mdplyr  [39m 1.0.7
    [32m✔[39m [34mtidyr  [39m 1.1.3     [32m✔[39m [34mstringr[39m 1.4.0
    [32m✔[39m [34mreadr  [39m 1.4.0     [32m✔[39m [34mforcats[39m 0.5.1
    [32m✔[39m [34mpurrr  [39m 0.3.4     
    
    ── [1mConflicts[22m ────────────────────────────────────────── tidyverse_conflicts() ──
    [31m✖[39m [34mdplyr[39m::[32mcollapse()[39m   masks [34mIRanges[39m::collapse()
    [31m✖[39m [34mdplyr[39m::[32mcombine()[39m    masks [34mBiobase[39m::combine(), [34mBiocGenerics[39m::combine()
    [31m✖[39m [34mdplyr[39m::[32mcount()[39m      masks [34mmatrixStats[39m::count()
    [31m✖[39m [34mdplyr[39m::[32mdesc()[39m       masks [34mIRanges[39m::desc()
    [31m✖[39m [34mtidyr[39m::[32mexpand()[39m     masks [34mS4Vectors[39m::expand()
    [31m✖[39m [34mdplyr[39m::[32mfilter()[39m     masks [34mstats[39m::filter()
    [31m✖[39m [34mdplyr[39m::[32mfirst()[39m      masks [34mS4Vectors[39m::first()
    [31m✖[39m [34mdplyr[39m::[32mlag()[39m        masks [34mstats[39m::lag()
    [31m✖[39m [34mggplot2[39m::[32mPosition()[39m masks [34mBiocGenerics[39m::Position(), [34mbase[39m::Position()
    [31m✖[39m [34mpurrr[39m::[32mreduce()[39m     masks [34mGenomicRanges[39m::reduce(), [34mIRanges[39m::reduce()
    [31m✖[39m [34mdplyr[39m::[32mrename()[39m     masks [34mS4Vectors[39m::rename()
    [31m✖[39m [34mdplyr[39m::[32mslice()[39m      masks [34mIRanges[39m::slice()
    
    Warning message in file(filename, "r", encoding = encoding):
    “cannot open file 'RNAseqfunctions.R': No such file or directory”



    Error in file(filename, "r", encoding = encoding): cannot open the connection
    Traceback:


    1. source("RNAseqfunctions.R")

    2. file(filename, "r", encoding = encoding)


load counts and sample table created in mouseexp4/DESeq/CaEf_analysis/loadsamples.ipynb


```R
load("2022_03-31_RNAseqcountsandsamples.RData")

```


```R
#load counts and table as DDS object with Assembly22 names  
dds_22 <- DEseqimport(Ca_counts, Ca_samples, catx2gene= FALSE)
```


```R
#load counts and table as DDS object with orf/gene names 
dds_gene <- DEseqimport(Ca_counts, Ca_samples, catx2gene= TRUE)
```

Generate PCA plots


```R
#normalize/transform samples for PCA analysis 
Ca_vst <- varianceStabilizingTransformation(dds_gene, blind = TRUE, fitType = "local")
```


```R
#generate PCA plots 

#without sample names
PCA(Ca_vst, "Principle Component Analysis of Mouse and in vitro samples \n")  

#with sample names
PCA(Ca_vst, "Principle Component Analysis of Mouse and in vitro samples \n", samplename=TRUE) + xlim(-50,75)
```


```R
#generate plots with specific coloring 
PCAdata <- plotPCA(Ca_vst, intgroup = c("condition", "species2", "environment" ), returnData = TRUE)
percentVar <- round(100 * attr(PCAdata, "percentVar"))
ggplot(PCAdata, aes(x = PC1, y = PC2, color = environment, shape = species2)) + 
        geom_point(size = 3) + xlab(paste0("PC1: ", percentVar[1], 
        "% variance")) + ylab(paste0("PC2: ", percentVar[2], 
        "% variance")) + coord_fixed() 
ggsave("PCAplotNEW.pdf", last_plot())
```


```R
#do DE anaylsis of Ca mouse vs. Ca invitro
#generate DEseq object 
dds_DE <- DESeq(dds_gene)
#do LFC analysis comparing Ca mouse vs. invitro 
dds_DELFC <- lfcShrink(dds_DE, contrast = c("condition", "mouse_Ca","invitro_Ca"), type = 'ashr')
#save in CSV file for outside anaylsis 
write.csv(dds_DELFC, file =("LFCres_Caresponse_mousevinvitro_CGDname.csv"))
```


```R
#generate same objects as above for more readability (gene/orf name) 
ddsgene_DE <- DESeq(dds_gene)
ddsgene_DELFC2 <- lfcShrink(ddsgene_DE, contrast = c("condition", "mouse_Ca","invitro_Ca"), type = 'ashr')
```


```R
#make volcano plot 
#volplot(LFC sample table, title, see inside function for more changes
volplot(ddsgene_DELFC2, "Ca mouse vs. Ca invitro \n Gene Expression")

library(EnhancedVolcano)

EnhancedVolcano(ddsgene_DELFC2,
    lab=rownames(ddsgene_DELFC2),
    x = 'log2FoldChange',
    y = 'padj',
    selectLab = c('REC8','SPO22','orf19.4789','DLH1','orf19.4955','CSK1','SPR3','orf19.5212','orf19.810','SPO11','PHR3','orf19.2713','orf19.3373','SPO1','orf19.7464','KAR4','orf19.553','DIT1','orf19.6405','YBL053','SPO75','orf19.1780','RON1','CEK2','STE4','STE3','MFALPHA','STE2','orf19.1476','PTH2','orf19.3306'),
    xlab = bquote(~Log[2]~ 'fold change'),
    ylab = bquote(~Log[10]~ ~adjusted~ 'P value'),
    pCutoff = 0.01,
    FCcutoff = 2,
    # xlim = c(),
    legendPosition = 'bottom',
    legendLabSize = 10,
    legendIconSize = 3.0,
    title = "CaEf mouse vs. Ca invitro \n Gene Expression in hypoxic conditions",
    labSize = 3.0,
    labCol = 'black',
    boxedLabels = TRUE,
    drawConnectors = TRUE,
    widthConnectors = 0.5,
    colConnectors = 'black', 
    maxoverlapsConnectors = 30, 
    subtitle = "")

```


```R


CaEfCacsv <- lfcShrink(ddsgene_DE, contrast = c("condition", "mouse_CaEf","invitro_Ca"), type = 'ashr') 

save(CaEfCacsv, file = "CaEfCaplot.Rdata")
# dev.copy(pdf,'mouseCaEf_v_invitroCa_volcano2.pdf', width = 3.125, height = 2.15)
# plot
# dev.off()
```

working with in vitro and mouse samples seperately


```R
#subset mouse and in vitro samples
samples_invitro <- filter(Ca_samples, condition == 'invitro_Ca' | condition == 'invitro_CaEf' )
invitrolist <-  samples_invitro$samplename
counts_invitro <- Ca_counts[,invitrolist]

samples_mouse <- filter(Ca_samples, condition == 'mouse_Ca' | condition == 'mouse_CaEf' )
mouselist <-  samples_mouse$samplename
counts_mouse <- Ca_counts[,mouselist]
```


```R
#write csv files saving LFCDEexpression with CGD names
DEseqimport(sampledata = samples_invitro, countdata = counts_invitro, catx2gene = FALSE) %>% 
    DESeq %>% 
    lfcShrink(., contrast = c("condition", "invitro_CaEf","invitro_Ca"), type = 'ashr') %>%
    write.csv(., file =("LFCres_invitro_CaEfvCa_CGDname.csv"))    

DEseqimport(sampledata = samples_mouse, countdata = counts_mouse, catx2gene = FALSE) %>% 
    DESeq %>% 
    lfcShrink(., contrast = c("condition", "mouse_CaEf","mouse_Ca"), type = 'ashr') %>%
    write.csv(., file =("LFCres_mouse_CaEfvCa_CGDname.csv"))    
```


```R

```


```R
#generate LFC with gene/orf names and make volcano plot - in vitro
plot <- DEseqimport(sampledata = samples_invitro, countdata = counts_invitro, catx2gene = TRUE) %>% 
    DESeq %>% 
    lfcShrink(., contrast = c("condition", "invitro_CaEf","invitro_Ca"), type = 'ashr') %>%
    volplot(., "CaEf invitro vs. Ca invitro \n Gene Expression in hypoxic conditions")
dev.copy(pdf, "invitro_CaEf_v_Ca_volcano_NEW.pdf")
plot
dev.off()
```


```R
#generate LFC with gene/orf names and make volcano plot - mouse 
DEseqimport(sampledata = samples_mouse, countdata = counts_mouse, catx2gene = TRUE) %>% 
    DESeq %>% 
    lfcShrink(., contrast = c("condition", "mouse_CaEf","mouse_Ca"), type = 'ashr') %>%
    volplot(., "CaEf mouse vs. Ca mouse \n Gene Expression in hypoxic conditions")
```
