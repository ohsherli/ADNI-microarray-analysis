20230130 cleaning up ADNI raw data
================

``` r
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
    ## ✔ ggplot2 3.4.1     ✔ purrr   1.0.1
    ## ✔ tibble  3.1.8     ✔ dplyr   1.1.0
    ## ✔ tidyr   1.3.0     ✔ stringr 1.5.0
    ## ✔ readr   2.1.3     ✔ forcats 0.5.2
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

``` r
library(ggplot2)
library(readxl)
library(RColorBrewer)

library(pheatmap)
library(polycor)
library(corrplot)
```

    ## corrplot 0.92 loaded

``` r
library(factoextra)
```

    ## Welcome! Want to learn more? See two factoextra-related books at https://goo.gl/ve3WBa

``` r
library(sva)
```

    ## Loading required package: mgcv
    ## Loading required package: nlme
    ## 
    ## Attaching package: 'nlme'
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     collapse
    ## 
    ## This is mgcv 1.8-41. For overview type 'help("mgcv-package")'.
    ## Loading required package: genefilter
    ## 
    ## Attaching package: 'genefilter'
    ## 
    ## The following object is masked from 'package:readr':
    ## 
    ##     spec
    ## 
    ## Loading required package: BiocParallel

``` r
library(randomForest)
```

    ## randomForest 4.7-1.1
    ## Type rfNews() to see new features/changes/bug fixes.
    ## 
    ## Attaching package: 'randomForest'
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     combine
    ## 
    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     margin

``` r
library(caret)
```

    ## Loading required package: lattice
    ## 
    ## Attaching package: 'caret'
    ## 
    ## The following object is masked from 'package:purrr':
    ## 
    ##     lift

``` r
library(tree)
library(Metrics)
```

    ## 
    ## Attaching package: 'Metrics'
    ## 
    ## The following objects are masked from 'package:caret':
    ## 
    ##     precision, recall

``` r
library(e1071)
library(Boruta)

# devtools::install_github("plotly/plotly.R")
library(plotly)
```

    ## 
    ## Attaching package: 'plotly'
    ## 
    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     last_plot
    ## 
    ## The following object is masked from 'package:stats':
    ## 
    ##     filter
    ## 
    ## The following object is masked from 'package:graphics':
    ## 
    ##     layout

``` r
library(cowplot)

library(limma)

# install.packages("here")
library(here)
```

    ## here() starts at /Users/ohsherli/Desktop/NTU/EG_Lab/Jan2023_ADNI_analysis

``` r
knitr::opts_knit$set(root.dir = "~/Desktop/NTU/EG_Lab/Jan2023_ADNI_analysis")
```

Reading in files:

``` r
## main gene expression file
adni_gene_expression_raw <- read.csv("../data/ADNI_raw_files/ADNI_Gene_Expression_Profile/ADNI_Gene_Expression_Profile.csv", header = FALSE)
# adni_gene_expression_raw %>% head
adni_gene_expression_raw[9:15,1:10]
```

    ##               V1        V2        V3    V4    V5    V6    V7    V8    V9   V10
    ## 9       ProbeSet LocusLink    Symbol                                          
    ## 10   11715100_at   LOC8355  HIST1H3G 2.237 2.294  2.14 2.062  2.04 2.439 1.955
    ## 11 11715101_s_at   LOC8355  HIST1H3G 2.624 2.416 2.322   2.5 2.395 2.309 2.451
    ## 12 11715102_x_at   LOC8355  HIST1H3G 1.873 1.884 1.999 1.851  2.08 1.997 1.539
    ## 13 11715103_x_at LOC126282 TNFAIP8L1  2.92 2.668 3.634 3.632 3.278 3.578 3.362
    ## 14 11715104_s_at  LOC92736     OTOP2 2.147 2.156 2.516 2.283 2.251 2.235 1.992
    ## 15   11715105_at LOC284099  C17ORF78 2.268  2.13 1.957 2.347 2.154 2.055 2.323

``` r
### first 8 rows of the dataframe are information about the subjects
adni_gene_expression_subjectInfo <- adni_gene_expression_raw[1:8,c(1,4:ncol(adni_gene_expression_raw))] %>%
  column_to_rownames(var = "V1")
colnames(adni_gene_expression_subjectInfo) <- c(paste0("subject_", 1:ncol(adni_gene_expression_subjectInfo)))
adni_gene_expression_subjectInfo[,(ncol(adni_gene_expression_subjectInfo)-5):ncol(adni_gene_expression_subjectInfo)]
```

    ##                  subject_740 subject_741 subject_742 subject_743 subject_744
    ## Phase                 ADNIGO       ADNI2       ADNI2       ADNI2       ADNI2
    ## Visit                     bl         v03         v03         v03         v06
    ## SubjectID         009_S_2381  053_S_4557  073_S_4300  041_S_4014  007_S_0101
    ## 260/280                 1.87        2.03        2.11        1.94        2.06
    ## 260/230                 1.45        1.33        0.27        1.72        1.35
    ## RIN                      6.6         6.8         6.2         5.8         6.7
    ## Affy Plate                 8           5           3           1           4
    ## YearofCollection        2011        2012        2011        2011        2012
    ##                  subject_745
    ## Phase                       
    ## Visit                       
    ## SubjectID                   
    ## 260/280                     
    ## 260/230                     
    ## RIN                         
    ## Affy Plate                  
    ## YearofCollection

``` r
# last column is empty (used to store gene symbol details in original df)
adni_gene_expression_subjectInfo <- adni_gene_expression_subjectInfo[,1:(ncol(adni_gene_expression_subjectInfo)-1)]
adni_gene_expression_subjectInfo <- adni_gene_expression_subjectInfo %>% t %>% data.frame #transpose for easier cleaning

adni_gene_expression_subjectInfo %>% nrow # 744
```

    ## [1] 744

``` r
adni_gene_expression_subjectInfo$RIN <- as.numeric(adni_gene_expression_subjectInfo$RIN)
```

    ## Warning: NAs introduced by coercion

``` r
adni_gene_expression_subjectInfo$Phase <- factor(adni_gene_expression_subjectInfo$Phase)
adni_gene_expression_subjectInfo$Affy.Plate <- factor(adni_gene_expression_subjectInfo$Affy.Plate)
adni_gene_expression_subjectInfo$YearofCollection <- factor(adni_gene_expression_subjectInfo$YearofCollection)

# adni_gene_expression_subjectInfo %>% head
```

``` r
### collate the rest of the rows of the dataframe to make gene expression matrix
adni_gene_expression_values <- adni_gene_expression_raw[9:nrow(adni_gene_expression_raw),]

# rename columns to match subjects and additional information
colnames(adni_gene_expression_values) <- c("ProbeSet", "LocusLink", "Symbol", 
                                           c(paste0("subject_", 1:(ncol(adni_gene_expression_values)-4))),
                                           "geneInfo")
# adni_gene_expression_values[(ncol(adni_gene_expression_values)-5):ncol(adni_gene_expression_values)]

# remove first row (column headings already added)
adni_gene_expression_values <- adni_gene_expression_values[2:nrow(adni_gene_expression_values), ]

# change microarray values to all numeric
adni_gene_expression_values[,4:(ncol(adni_gene_expression_values)-1)] <- lapply(adni_gene_expression_values[,4:(ncol(adni_gene_expression_values)-1)], as.numeric)
# adni_gene_expression_values %>% head(20)

# fix gene name for SEPT2
which(adni_gene_expression_values$Symbol == "2-Sep")
```

    ## [1]   267   268 34688 49209

``` r
adni_gene_expression_values$Symbol[267] <- "SEPT2"
adni_gene_expression_values$Symbol[268] <- "SEPT2"
adni_gene_expression_values$Symbol[34688] <- "SEPT2"
adni_gene_expression_values$Symbol[49209] <- "SEPT2"
```

``` r
# separate gene expression dataframe into probe/gene information and gene expression values
adni_gene_expression_geneInfo <- adni_gene_expression_values[,c(1:3,ncol(adni_gene_expression_values))]

# keep ProbeSet as rownames (most complete set of gene "names", remove all other geneInfo columns)
adni_gene_expression_values <- adni_gene_expression_values[,c(1, 4:(ncol(adni_gene_expression_values)-1))] 
rownames(adni_gene_expression_values) <- adni_gene_expression_values$ProbeSet
adni_gene_expression_values <- adni_gene_expression_values[,c(2:ncol(adni_gene_expression_values))] %>% 
  t %>%
  as.data.frame()
# adni_gene_expression_values %>% head
adni_gene_expression_values %>% nrow # 744
```

    ## [1] 744

``` r
# save(adni_gene_expression_values, file = "adni_gene_expression_values.RData")
```

``` r
# import diagnosis data
diagnosis_info <- read.csv("~/Desktop/NTU/EG_Lab/data/ADNI_raw_files/Diagnosis/DXSUM_PDXCONV_ADNIALL.csv")
# diagnosis_info %>% head

diagnosis_info.split_phase <- diagnosis_info %>% split(diagnosis_info$Phase)
# diagnosis_info.split_phase$ADNIGO %>% head()
# diagnosis_info.split_phase$ADNI2
```

ADNI2: DXCHANGE  
DXSUM  
Diagnostic Summary 1. Which best describes the participant’s change in
cognitive status from last visit to current visit? At Screening Visit,
indicate initial diagnosis using one of the ‘Stable’ options. N 1  
1=Stable: NL; 2=Stable: MCI; 3=Stable: Dementia; 4=Conversion: NL to
MCI; 5=Conversion: MCI to Dementia; 6=Conversion: NL to Dementia;
7=Reversion: MCI to NL; 8=Reversion: Dementia to MCI; 9=Reversion:
Dementia to NL

ADNIGO: ADNIGO  
DXCHANGE  
DXSUM Diagnostic Summary 1. Which best describes the participant’s
change in cognitive status from last visit to current visit: N 1  
1=Stable: NL to NL; 2=Stable: MCI to MCI; 3=Stable: Dementia to
Dementia; 4=Conversion: NL to MCI; 5=Conversion: MCI to Dementia;
6=Conversion: NL to Dementia; 7=Reversion: MCI to NL; 8=Reversion:
Dementia to MCI; 9=Reversion: Dementia to NL

``` r
# adni_gene_expression_subjectInfo %>% head
adni_gene_expression_subjectInfo.split_phase <- adni_gene_expression_subjectInfo %>% 
  split(adni_gene_expression_subjectInfo$Phase)
# adni_gene_expression_subjectInfo.split_phase
```

### Diagnosis information for adni2

``` r
# list of adni2 diagnosis split by subject
#adni2_diagnosis.split_subject <- split(diagnosis_info.split_phase$ADNI2, diagnosis_info.split_phase$ADNI2$PTID) 
# adni2_diagnosis.split_subject$"002_S_0295"
# which(adni2_diagnosis.split_subject$"002_S_0295"$VISCODE == "v06")

# initialise empty vector to store diagnosis values
gene_expression_adni2_diagnosis_values <- rep(0,nrow(adni_gene_expression_subjectInfo.split_phase$ADNI2))
names(gene_expression_adni2_diagnosis_values) <- rownames(adni_gene_expression_subjectInfo.split_phase$ADNI2)

for (i in 1:nrow(adni_gene_expression_subjectInfo.split_phase$ADNI2)) {
  
  subjectID_i <- adni_gene_expression_subjectInfo.split_phase$ADNI2$SubjectID[i]
  visit_i <- adni_gene_expression_subjectInfo.split_phase$ADNI2$Visit[i]
  
  # find index of subject 
  adni2_diagnosis.subject_index_i <- which(diagnosis_info.split_phase$ADNI2$PTID == subjectID_i)
  # find index of visit
  adni2_diagnosis.visit_index_i <- which(diagnosis_info.split_phase$ADNI2$VISCODE == visit_i)
  
  # select row based on visit
  adni2_diagnosis.subject_visit_i <-
    diagnosis_info.split_phase$ADNI2[intersect(adni2_diagnosis.subject_index_i, adni2_diagnosis.visit_index_i),]
  
  # get DXCHANGE column value
  diagnosis_value_i <- adni2_diagnosis.subject_visit_i$DXCHANGE
  
  if(length(diagnosis_value_i) == 0){
    diagnosis_value_i = NA
  }
  
  gene_expression_adni2_diagnosis_values[i] <- diagnosis_value_i # store diagnosis value into vector
}
```

``` r
# vector of subjects that have no diagnosis values
which(is.na(gene_expression_adni2_diagnosis_values))
```

    ##  subject_73  subject_79 subject_124 subject_210 subject_361 subject_386 
    ##          47          52          76         135         223         243 
    ## subject_540 subject_636 subject_679 
    ##         335         385         407

``` r
 # subject_73  subject_79 subject_124 subject_210 subject_361 subject_386 subject_540 subject_636 subject_679 
 #         47          52          76         135         223         243         335         385         407

# check DXSUM_PDCONV_ADNIALL sheet for individual patient information

adni_gene_expression_subjectInfo.split_phase$ADNI2[47,] # 072_S_4102 v04
```

    ##            Phase Visit  SubjectID X260.280 X260.230 RIN Affy.Plate
    ## subject_73 ADNI2   v04 072_S_4102     2.06     1.64 6.7          8
    ##            YearofCollection
    ## subject_73             2011

``` r
# conversion to dementia at visit 11 but still MCI at visit 3 and visit 5
# assumption of MCI at visit 4
gene_expression_adni2_diagnosis_values[47] <- 2

adni_gene_expression_subjectInfo.split_phase$ADNI2[52,] # 006_S_4192 v04
```

    ##            Phase Visit  SubjectID X260.280 X260.230 RIN Affy.Plate
    ## subject_79 ADNI2   v04 006_S_4192     2.08     0.28 5.9          5
    ##            YearofCollection
    ## subject_79             2011

``` r
# status is 3 at visit 1, 3, 5 and 11
# assumption of status 3 at visit 4
gene_expression_adni2_diagnosis_values[52] <- 3

adni_gene_expression_subjectInfo.split_phase$ADNI2[76,] # 013_S_4268 v04
```

    ##             Phase Visit  SubjectID X260.280 X260.230 RIN Affy.Plate
    ## subject_124 ADNI2   v04 013_S_4268     2.02     0.68 7.3          2
    ##             YearofCollection
    ## subject_124             2011

``` r
# status is 2 at visit 1, 3, 5, conversion to 7 at visit 11 (reversion to NL)
# assumption of status 2 at visit 4
gene_expression_adni2_diagnosis_values[76] <- 2

adni_gene_expression_subjectInfo.split_phase$ADNI2[135,] # 072_S_4103 v04
```

    ##             Phase Visit  SubjectID X260.280 X260.230 RIN Affy.Plate
    ## subject_210 ADNI2   v04 072_S_4103     2.05      1.9 7.4          8
    ##             YearofCollection
    ## subject_210             2011

``` r
# status is 1 at visit 1, 3, 5
# assumption of status 1 at visit 4
gene_expression_adni2_diagnosis_values[135] <- 1

adni_gene_expression_subjectInfo.split_phase$ADNI2[223,] # 009_S_4564 v02
```

    ##             Phase Visit  SubjectID X260.280 X260.230 RIN Affy.Plate
    ## subject_361 ADNI2   v02 009_S_4564     2.04     1.17 7.6          5
    ##             YearofCollection
    ## subject_361             2012

``` r
# status is 3 at visit 1 and 3
# assumption of status 2 at visit 2
gene_expression_adni2_diagnosis_values[223] <- 2

adni_gene_expression_subjectInfo.split_phase$ADNI2[243,] # 022_S_4266 v04
```

    ##             Phase Visit  SubjectID X260.280 X260.230 RIN Affy.Plate
    ## subject_386 ADNI2   v04 022_S_4266     2.06     1.71 7.4          4
    ##             YearofCollection
    ## subject_386             2012

``` r
# status is 1 at visit 1, 3 and 5
# assumption of status 1 at visit 4
gene_expression_adni2_diagnosis_values[243] <- 1

adni_gene_expression_subjectInfo.split_phase$ADNI2[335,] # 114_S_4379 v04
```

    ##             Phase Visit  SubjectID X260.280 X260.230 RIN Affy.Plate
    ## subject_540 ADNI2   v04 114_S_4379     2.05     1.38 5.7          5
    ##             YearofCollection
    ## subject_540             2012

``` r
# status is 3 at visit 1, 3 and 5
# assumption of status 3 at visit 4
gene_expression_adni2_diagnosis_values[335] <- 3

adni_gene_expression_subjectInfo.split_phase$ADNI2[385,] # 094_S_4282 v04
```

    ##             Phase Visit  SubjectID X260.280 X260.230 RIN Affy.Plate
    ## subject_636 ADNI2   v04 094_S_4282     1.94     1.48 6.5          2
    ##             YearofCollection
    ## subject_636             2011

``` r
# status is 3 at visit 1, 3 and 5
# assumption of status 3 at visit 4
gene_expression_adni2_diagnosis_values[385] <- 3

adni_gene_expression_subjectInfo.split_phase$ADNI2[407,] # 037_S_4432 v04
```

    ##             Phase Visit  SubjectID X260.280 X260.230 RIN Affy.Plate
    ## subject_679 ADNI2   v04 037_S_4432     1.87     1.29 7.1          6
    ##             YearofCollection
    ## subject_679             2012

``` r
# status is 2 at visit 1 and 3, status is 5 at visit 5 (conversion to dementia)
# remove subject from later analysis (tag as NA for removal later)
gene_expression_adni2_diagnosis_values[407] <- NA

gene_expression_adni2_diagnosis_values %>% table
```

    ## .
    ##   1   2   3   5 
    ## 171 211  61   5

``` r
#   1   2   3   5 
# 171 212  61   5 

# gene_expression_adni2_diagnosis_values %>% head
```

``` r
# new adni2 subject info including diagnosis
adni2_diagnosis_subjectInfo <- adni_gene_expression_subjectInfo.split_phase$ADNI2 %>%
  mutate(Diagnosis = gene_expression_adni2_diagnosis_values,
         .after = "SubjectID")
# adni2_diagnosis_subjectInfo %>% head
```

## Diagnosis information for ADNIGO

``` r
# initialise empty vector to store diagnosis values
gene_expression_ADNIGO_diagnosis_values <- rep(0,nrow(adni_gene_expression_subjectInfo.split_phase$ADNIGO))
names(gene_expression_ADNIGO_diagnosis_values) <- rownames(adni_gene_expression_subjectInfo.split_phase$ADNIGO)

for (i in 1:nrow(adni_gene_expression_subjectInfo.split_phase$ADNIGO)) {
  
  subjectID_i <- adni_gene_expression_subjectInfo.split_phase$ADNIGO$SubjectID[i]
  visit_i <- adni_gene_expression_subjectInfo.split_phase$ADNIGO$Visit[i]
  
  # find index of subject 
  ADNIGO_diagnosis.subject_index_i <- which(diagnosis_info.split_phase$ADNIGO$PTID == subjectID_i)
  # find index of visit
  ADNIGO_diagnosis.visit_index_i <- which(diagnosis_info.split_phase$ADNIGO$VISCODE == visit_i)
  
  # select row based on visit
  ADNIGO_diagnosis.subject_visit_i <-
    diagnosis_info.split_phase$ADNIGO[intersect(ADNIGO_diagnosis.subject_index_i, ADNIGO_diagnosis.visit_index_i),]
  
  # get DXCHANGE column value
  diagnosis_value_i <- ADNIGO_diagnosis.subject_visit_i$DXCHANGE
  
  if(length(diagnosis_value_i) == 0){
    diagnosis_value_i = NA
  }
  
  gene_expression_ADNIGO_diagnosis_values[i] <- diagnosis_value_i # store diagnosis value into vector
}
```

``` r
which(is.na(gene_expression_ADNIGO_diagnosis_values)) 
```

    ## subject_60 
    ##         22

``` r
# subject_60 
#         22 

# check DXSUM_PDCONV_ADNIALL sheet for individual patient information

adni_gene_expression_subjectInfo.split_phase$ADNIGO[22,] # 073_S_2182 m03
```

    ##             Phase Visit  SubjectID X260.280 X260.230 RIN Affy.Plate
    ## subject_60 ADNIGO   m03 073_S_2182     2.06     1.71 7.6          9
    ##            YearofCollection
    ## subject_60             2011

``` r
# status at bl is 2 and m06 is 2
# assumption of status 2 at m03
gene_expression_ADNIGO_diagnosis_values[22] <- 2

gene_expression_ADNIGO_diagnosis_values %>% table
```

    ## .
    ##   1   2   3   4   5   7   8 
    ##  73 159  37  10  13   2   1

``` r
 #  1   2   3   4   5   7   8 
 # 73 159  37  10  13   2   1 

# gene_expression_ADNIGO_diagnosis_values %>% head
```

``` r
# new adnigo subject info including diagnosis
adnigo_diagnosis_subjectInfo <- adni_gene_expression_subjectInfo.split_phase$ADNIGO %>%
  mutate(Diagnosis = gene_expression_ADNIGO_diagnosis_values,
         .after = "SubjectID")
# adnigo_diagnosis_subjectInfo %>% head
```

### Combine patient info tables from ADNI2 and ADNIGO (with added diagnosis column)

``` r
adni_gene_expression_subjectInfo_Dx <- rbind(adni2_diagnosis_subjectInfo, adnigo_diagnosis_subjectInfo)

# reorder rows
adni_gene_expression_subjectInfo_Dx <- adni_gene_expression_subjectInfo_Dx[
  match(rownames(adni_gene_expression_values), 
        rownames(adni_gene_expression_subjectInfo_Dx)),]

# adni_gene_expression_subjectInfo_Dx %>% head

adni_gene_expression_subjectInfo_Dx$Diagnosis %>% table
```

    ## .
    ##   1   2   3   4   5   7   8 
    ## 244 370  98  10  18   2   1

``` r
## diagnosis dictionary file
adni_diagnosis_metadata_raw <- read.csv("../data/ADNI_raw_files/Diagnosis/DXSUM_PDXCONV_ADNIALL.csv")
# adni_diagnosis_metadata_raw %>% head(10)
```

``` r
write.csv(adni_gene_expression_values, file = "20230201_adni_gene_expression_values.csv")
write.csv(adni_gene_expression_subjectInfo_Dx, file = "20230201_adni_gene_expression_subjectInfo_Dx.csv")
```

## Compiling further patient information (age, sex, MMSE scores)

``` r
ADNIMERGE <- read.csv("/Users/ohsherli/Desktop/NTU/EG_Lab/data/ADNI_raw_files/Data___Database/ADNIMERGE.csv")
APOERES <- read.csv("/Users/ohsherli/Desktop/NTU/EG_Lab/data/ADNI_raw_files/APOERES.csv")
ADNI_mRNA_Age <- read.csv("/Users/ohsherli/Desktop/NTU/EG_Lab/data/ADNI_raw_files/ADNI_mRNA_Age.csv")
```

**From APOERES_DICT:** APGEN1 APOEGO2 ADNI GO & 2 APOE Genotypes -4
Genotype - Allele 1 APGEN2 APOEGO2 ADNI GO & 2 APOE Genotypes -4
Genotype - Allele 2

``` r
# save.image(file = "20230201_workspace_cleaning_adni.RData")
```

``` r
sessionInfo()
```

    ## R version 4.2.2 (2022-10-31)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Big Sur ... 10.16
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] here_1.0.1           limma_3.54.1         cowplot_1.1.1       
    ##  [4] plotly_4.10.1        Boruta_8.0.0         e1071_1.7-13        
    ##  [7] Metrics_0.1.4        tree_1.0-42          caret_6.0-93        
    ## [10] lattice_0.20-45      randomForest_4.7-1.1 sva_3.46.0          
    ## [13] BiocParallel_1.32.5  genefilter_1.80.0    mgcv_1.8-41         
    ## [16] nlme_3.1-160         factoextra_1.0.7     corrplot_0.92       
    ## [19] polycor_0.8-1        pheatmap_1.0.12      RColorBrewer_1.1-3  
    ## [22] readxl_1.4.1         forcats_0.5.2        stringr_1.5.0       
    ## [25] dplyr_1.1.0          purrr_1.0.1          readr_2.1.3         
    ## [28] tidyr_1.3.0          tibble_3.1.8         ggplot2_3.4.1       
    ## [31] tidyverse_1.3.2     
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] backports_1.4.1        plyr_1.8.8             lazyeval_0.2.2        
    ##   [4] splines_4.2.2          listenv_0.9.0          GenomeInfoDb_1.34.9   
    ##   [7] digest_0.6.31          foreach_1.5.2          htmltools_0.5.4       
    ##  [10] fansi_1.0.4            magrittr_2.0.3         memoise_2.0.1         
    ##  [13] googlesheets4_1.0.1    tzdb_0.3.0             recipes_1.0.3         
    ##  [16] globals_0.16.2         Biostrings_2.66.0      annotate_1.76.0       
    ##  [19] modelr_0.1.10          gower_1.0.0            matrixStats_0.63.0    
    ##  [22] hardhat_1.2.0          timechange_0.1.1       colorspace_2.1-0      
    ##  [25] blob_1.2.3             rvest_1.0.3            ggrepel_0.9.3         
    ##  [28] haven_2.5.1            xfun_0.37              crayon_1.5.2          
    ##  [31] RCurl_1.98-1.10        jsonlite_1.8.4         survival_3.5-5        
    ##  [34] iterators_1.0.14       glue_1.6.2             gtable_0.3.1          
    ##  [37] gargle_1.2.1           ipred_0.9-13           zlibbioc_1.44.0       
    ##  [40] XVector_0.38.0         future.apply_1.10.0    BiocGenerics_0.44.0   
    ##  [43] scales_1.2.1           DBI_1.1.3              edgeR_3.40.1          
    ##  [46] Rcpp_1.0.10            viridisLite_0.4.1      xtable_1.8-4          
    ##  [49] proxy_0.4-27           bit_4.0.5              stats4_4.2.2          
    ##  [52] lava_1.7.0             prodlim_2019.11.13     htmlwidgets_1.6.1     
    ##  [55] httr_1.4.4             ellipsis_0.3.2         pkgconfig_2.0.3       
    ##  [58] XML_3.99-0.13          nnet_7.3-18            dbplyr_2.3.0          
    ##  [61] locfit_1.5-9.6         utf8_1.2.3             tidyselect_1.2.0      
    ##  [64] rlang_1.0.6            reshape2_1.4.4         AnnotationDbi_1.60.0  
    ##  [67] munsell_0.5.0          cellranger_1.1.0       tools_4.2.2           
    ##  [70] cachem_1.0.6           cli_3.6.0              generics_0.1.3        
    ##  [73] RSQLite_2.3.0          broom_1.0.1            evaluate_0.20         
    ##  [76] fastmap_1.1.0          yaml_2.3.7             ModelMetrics_1.2.2.2  
    ##  [79] knitr_1.42             bit64_4.0.5            fs_1.6.1              
    ##  [82] admisc_0.30            KEGGREST_1.38.0        future_1.31.0         
    ##  [85] xml2_1.3.3             compiler_4.2.2         rstudioapi_0.14       
    ##  [88] png_0.1-8              reprex_2.0.2           stringi_1.7.12        
    ##  [91] Matrix_1.5-1           vctrs_0.5.2            pillar_1.8.1          
    ##  [94] lifecycle_1.0.3        data.table_1.14.8      bitops_1.0-7          
    ##  [97] R6_2.5.1               IRanges_2.32.0         parallelly_1.34.0     
    ## [100] codetools_0.2-18       MASS_7.3-58.3          assertthat_0.2.1      
    ## [103] rprojroot_2.0.3        withr_2.5.0            S4Vectors_0.36.1      
    ## [106] GenomeInfoDbData_1.2.9 parallel_4.2.2         hms_1.1.2             
    ## [109] grid_4.2.2             rpart_4.1.19           timeDate_4021.107     
    ## [112] class_7.3-20           rmarkdown_2.20         googledrive_2.0.0     
    ## [115] pROC_1.18.0            Biobase_2.58.0         lubridate_1.9.0
