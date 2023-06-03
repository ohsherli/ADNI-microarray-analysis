Using ADNIMERGE diagnosis labels and comparing between labels
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
# BiocManager::install("sva")
library(sva)

library(stringr)
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
# import diagnosis info from adnimerge dataset
adnimerge_filePath <- "../data/ADNI_raw_files/Data___Database/ADNIMERGE.csv"
adnimerge <- read.csv(adnimerge_filePath)
# adnimerge %>% head
```

# Match ADNImerge subjects with list of subjects from gene expression data by study date

``` r
# split gene expression data by adni2 and adnigo
# adni_gene_expression_subjectInfo %>% head
adni_gene_expression_subjectInfo.split_phase <- adni_gene_expression_subjectInfo %>% 
  split(adni_gene_expression_subjectInfo$Phase)
# adni_gene_expression_subjectInfo.split_phase
```

## Get diagnosis for ADNI2 subjects

``` r
# for adni2, use viscode2 instead of viscode for Visit by adding additional column

## get viscode and viscode2 data from PTDEMOG file
PTDX <- read.csv("~/Desktop/NTU/EG_Lab/data/ADNI_raw_files/Diagnosis/DXSUM_PDXCONV_ADNIALL.csv")
# PTDX %>% head()
```

``` r
# extract patient RID 
?str_pad
PT_RID_info <- PTDX[,c("Phase", "RID", "SITEID", "VISCODE", "VISCODE2")]
# PT_RID_info %>% head

# add a column where RID is character instead of numeric
PT_RID_info <- PT_RID_info %>% 
  mutate(
    RID_chr = RID %>% as.character, .before = SITEID
  )

# pad elements of new RID_chr column with zeroes to make total of 4 digits
PT_RID_info$RID_chr <- str_pad(PT_RID_info$RID_chr, 4, pad = "0")
# PT_RID_info %>% head

PT_RID_info.splitPhase <- split(PT_RID_info, PT_RID_info$Phase)
PT_RID_info.splitPhase.ADNI2 <- PT_RID_info.splitPhase$ADNI2
```

``` r
# from ADNI2 gene expression data, add new VISCODE2 column based on RID_chr matched with subjectID and VISCODE
adni_gene_expression_subjectInfo.split_phase.ADNI2 <- adni_gene_expression_subjectInfo.split_phase$ADNI2 

RID_ADNI2_vector <- rep(0, nrow(adni_gene_expression_subjectInfo.split_phase.ADNI2))
RID_ADNI_list <- strsplit(adni_gene_expression_subjectInfo.split_phase.ADNI2$SubjectID, 
                          split = "_")

for (i in 1:length(RID_ADNI_list)){
  RID_ADNI_i <- RID_ADNI_list[[i]][3]
  RID_ADNI2_vector[i] <- RID_ADNI_i
}

adni_gene_expression_subjectInfo.split_phase.ADNI2 <- adni_gene_expression_subjectInfo.split_phase.ADNI2 %>%
  mutate(RID_chr = RID_ADNI2_vector, .after = SubjectID)

# adni_gene_expression_subjectInfo.split_phase.ADNI2 %>% head
```

``` r
# match new ADNI2 data with VISCODE2 in PT_RID_info

# initalize empty vector to store VISCODE2 values
ADNI2_geneExpression_VISCODE2_values <- rep(0, nrow(adni_gene_expression_subjectInfo.split_phase.ADNI2))
names(ADNI2_geneExpression_VISCODE2_values) <- rownames(adni_gene_expression_subjectInfo.split_phase.ADNI2)

for(i in 1:nrow(adni_gene_expression_subjectInfo.split_phase.ADNI2)){
  
  subject_RID_i <- adni_gene_expression_subjectInfo.split_phase.ADNI2$RID[i]
  viscode_i <- adni_gene_expression_subjectInfo.split_phase.ADNI2$Visit[i]
  
  # find index of RID in PT_RID_info
  ADNI2_geneExpression_PTDEMOG_RID_match_i <- which(
    PT_RID_info.splitPhase.ADNI2$RID_chr == subject_RID_i
  )
  # find index of viscode in PT_RID_info
  ADNI2_geneExpression_PTDEMOG_viscode_match_i <- which(
    PT_RID_info.splitPhase.ADNI2$VISCODE == viscode_i
  )
  # find intersection of RID and viscode to give row number
  ADNI2_RIDinfo_RID_viscode_match_i <- intersect(
    ADNI2_geneExpression_PTDEMOG_RID_match_i, 
    ADNI2_geneExpression_PTDEMOG_viscode_match_i
  )
  # select row and viscode2 value based on row number
  viscode2_i <- (PT_RID_info.splitPhase.ADNI2[ADNI2_RIDinfo_RID_viscode_match_i,])$VISCODE2
  
  # for rows with missing values, i.e. no VISCODE2 found, use NA
  if (length(viscode2_i)==0){
    viscode2_i = NA
  }
  
  # replace vector position i with viscode2 value
  ADNI2_geneExpression_VISCODE2_values[i] <- viscode2_i
  
}
```

``` r
# check results of running loop
# ADNI2_geneExpression_VISCODE2_values %>% head
ADNI2_geneExpression_VISCODE2_values %>% unique
```

    ## [1] "bl"  "m72" "m60" "m12" "m48" NA    "m06" "m84"

``` r
# get index of subjects with no viscode2 values
which(is.na(ADNI2_geneExpression_VISCODE2_values))
```

    ##  subject_73  subject_79 subject_124 subject_210 subject_361 subject_386 
    ##          47          52          76         135         223         243 
    ## subject_540 subject_636 subject_679 
    ##         335         385         407

``` r
 # subject_73  subject_79 subject_124 subject_210 subject_361 subject_386 subject_540 subject_636 subject_679 
 #         47          52          76         135         223         243         335         385         407 
```

``` r
# ADNI2_geneExpression_VISCODE2_values %>% head(18)
# adni_gene_expression_subjectInfo.split_phase.ADNI2 %>% head(18)
```

``` r
# match new ADNI2 data with VISCODE2 to ADNIMERGE for diagnosis labels

# split ADNIMERGE data by colprot (phase of study to minimise multiple labels across multiple study phases)
adnimerge_phaseSplit <- adnimerge %>% split(adnimerge$COLPROT)
adnimerge_adni2 <- adnimerge_phaseSplit$ADNI2

# adnimerge_adni2 %>% head
```

``` r
# initialize empty vector to store diagnosis labels for ADNI2 data
ADNI2_geneExpression_Dx_values <- rep(0, nrow(adni_gene_expression_subjectInfo.split_phase.ADNI2))
names(ADNI2_geneExpression_Dx_values) <- rownames(adni_gene_expression_subjectInfo.split_phase.ADNI2)

for(i in 1:nrow(adni_gene_expression_subjectInfo.split_phase.ADNI2)){
  
  # get subject RID to match with ADNIMERGE PTID
  subjectID_i <- adni_gene_expression_subjectInfo.split_phase.ADNI2$SubjectID[i]
  # get viscode2 for subject in gene expression data
  subject_viscode2_i <- ADNI2_geneExpression_VISCODE2_values[i]
  
  # for rows with no missing viscode2 values:
  if(is.na(subject_viscode2_i)){
    subject_ADNIMERGE_Dx_i <- NA
  } else {
    # find index of subjectID in ADNIMERGE data
    adnimerge_subjectID_match_i <- which(adnimerge_adni2$PTID == subjectID_i)
    # find index of viscode2 in ADNIMERGE data
    adnimerge_viscode2_match_i <- which(adnimerge_adni2$VISCODE == subject_viscode2_i)
    # find intersection of index of RID and index of viscode2 to give row number of diagnosis
    adnimerge_subjectID_viscode2_match_i <- intersect(adnimerge_subjectID_match_i, adnimerge_viscode2_match_i)
    # select row and diagnosis according to ADNIMERGE data
    subject_ADNIMERGE_Dx_i <- (adnimerge_adni2[adnimerge_subjectID_viscode2_match_i,])$DX
  }
  
  # store diagnosis value in vector
  ADNI2_geneExpression_Dx_values[i] <- subject_ADNIMERGE_Dx_i
}
```

``` r
# ADNI2_geneExpression_Dx_values %>% head
```

``` r
# fix the unknown diagnosis values by checking diagnosis before and after the visit

# get index of subjects with no viscode2 values
which(is.na(ADNI2_geneExpression_VISCODE2_values))
```

    ##  subject_73  subject_79 subject_124 subject_210 subject_361 subject_386 
    ##          47          52          76         135         223         243 
    ## subject_540 subject_636 subject_679 
    ##         335         385         407

``` r
 # subject_73  subject_79 subject_124 subject_210 subject_361 subject_386 subject_540 subject_636 subject_679 
 #         47          52          76         135         223         243         335         385         407
```

``` r
# subject_73
adni_gene_expression_subjectInfo.split_phase.ADNI2[47, ]
```

    ##            Phase Visit  SubjectID RID_chr X260.280 X260.230 RIN Affy.Plate
    ## subject_73 ADNI2   v04 072_S_4102    4102     2.06     1.64 6.7          8
    ##            YearofCollection
    ## subject_73             2011

``` r
# visit v04, subject ID 072_S_4102
subject_viscode_check <- which(PT_RID_info.splitPhase.ADNI2$RID_chr == "4102")
# [1]  163  351  860 2028
PT_RID_info.splitPhase.ADNI2[subject_viscode_check,]
```

    ##      Phase  RID RID_chr SITEID VISCODE VISCODE2
    ## 4506 ADNI2 4102    4102     34     v01       sc
    ## 4694 ADNI2 4102    4102     34     v03       bl
    ## 5203 ADNI2 4102    4102     34     v05      m06
    ## 6371 ADNI2 4102    4102     34     v11      m12

``` r
# viscode2 is between bl and m06
adnimerge_adni2[which(adnimerge_adni2$PTID == "072_S_4102"), c("RID", "PTID", "VISCODE", "DX_bl", "DX")]
```

    ##       RID       PTID VISCODE DX_bl       DX
    ## 8918 4102 072_S_4102     m12  EMCI Dementia
    ## 8919 4102 072_S_4102     m06  EMCI      MCI
    ## 8920 4102 072_S_4102     m03  EMCI         
    ## 8921 4102 072_S_4102      bl  EMCI      MCI

``` r
# at bl and m06, diagnosis is MCI, so use MCI as diagnosis
ADNI2_geneExpression_Dx_values[47] <- "MCI"
```

``` r
# subject_79
adni_gene_expression_subjectInfo.split_phase.ADNI2[52, ]
```

    ##            Phase Visit  SubjectID RID_chr X260.280 X260.230 RIN Affy.Plate
    ## subject_79 ADNI2   v04 006_S_4192    4192     2.08     0.28 5.9          5
    ##            YearofCollection
    ## subject_79             2011

``` r
# visit v04, subject ID 006_S_4192
subject_viscode_check <- which(PT_RID_info.splitPhase.ADNI2$RID_chr == "4192")
# [1]  355  652 1435 2240 3996
PT_RID_info.splitPhase.ADNI2[subject_viscode_check,]
```

    ##      Phase  RID RID_chr SITEID VISCODE VISCODE2
    ## 4698 ADNI2 4192    4192      4     v01       sc
    ## 4995 ADNI2 4192    4192      4     v03       bl
    ## 5778 ADNI2 4192    4192      4     v05      m06
    ## 6583 ADNI2 4192    4192      4     v11      m12
    ## 8339 ADNI2 4192    4192      4     v21      m24

``` r
# viscode2 is between bl and m06
adnimerge_adni2[which(adnimerge_adni2$PTID == "006_S_4192"), c("RID", "PTID", "VISCODE", "DX_bl", "DX")]
```

    ##       RID       PTID VISCODE DX_bl       DX
    ## 8568 4192 006_S_4192     m30    AD         
    ## 8569 4192 006_S_4192     m24    AD Dementia
    ## 8570 4192 006_S_4192     m18    AD         
    ## 8571 4192 006_S_4192     m12    AD Dementia
    ## 8572 4192 006_S_4192     m06    AD Dementia
    ## 8573 4192 006_S_4192     m03    AD         
    ## 8574 4192 006_S_4192      bl    AD Dementia

``` r
# at bl and m06, diagnosis is Dementia, so use Dementia as diagnosis
ADNI2_geneExpression_Dx_values[52] <- "Dementia"
```

``` r
# subject_124
adni_gene_expression_subjectInfo.split_phase.ADNI2[76, ]
```

    ##             Phase Visit  SubjectID RID_chr X260.280 X260.230 RIN Affy.Plate
    ## subject_124 ADNI2   v04 013_S_4268    4268     2.02     0.68 7.3          2
    ##             YearofCollection
    ## subject_124             2011

``` r
# visit v04, subject ID 013_S_4268
subject_viscode_check <- which(PT_RID_info.splitPhase.ADNI2$RID_chr == "4268")
# [1]  471  747 1619 2359 3844 4638
PT_RID_info.splitPhase.ADNI2[subject_viscode_check,]
```

    ##      Phase  RID RID_chr SITEID VISCODE VISCODE2
    ## 4814 ADNI2 4268    4268     10     v01       sc
    ## 5090 ADNI2 4268    4268     10     v03       bl
    ## 5962 ADNI2 4268    4268     10     v05      m06
    ## 6702 ADNI2 4268    4268     10     v11      m12
    ## 8187 ADNI2 4268    4268     10     v21      m24
    ## 8981 ADNI2 4268    4268     10     v31      m36

``` r
# viscode2 is between bl and m06
adnimerge_adni2[which(adnimerge_adni2$PTID == "013_S_4268"), c("RID", "PTID", "VISCODE", "DX_bl", "DX")]
```

    ##       RID       PTID VISCODE DX_bl  DX
    ## 8217 4268 013_S_4268     m36  EMCI  CN
    ## 8218 4268 013_S_4268     m30  EMCI    
    ## 8219 4268 013_S_4268     m24  EMCI  CN
    ## 8220 4268 013_S_4268     m18  EMCI    
    ## 8221 4268 013_S_4268     m12  EMCI  CN
    ## 8222 4268 013_S_4268     m06  EMCI MCI
    ## 8223 4268 013_S_4268     m03  EMCI    
    ## 8224 4268 013_S_4268      bl  EMCI MCI

``` r
# at bl and m06, diagnosis is MCI, so use MCI as diagnosis
# note reversion back to control from m12
ADNI2_geneExpression_Dx_values[76] <- "MCI"
```

``` r
# subject_210
adni_gene_expression_subjectInfo.split_phase.ADNI2[135, ]
```

    ##             Phase Visit  SubjectID RID_chr X260.280 X260.230 RIN Affy.Plate
    ## subject_210 ADNI2   v04 072_S_4103    4103     2.05      1.9 7.4          8
    ##             YearofCollection
    ## subject_210             2011

``` r
# visit v04, subject ID 072_S_4103
subject_viscode_check <- which(PT_RID_info.splitPhase.ADNI2$RID_chr == "4103")
# [1]  160  368 1479 2023 3703
PT_RID_info.splitPhase.ADNI2[subject_viscode_check,]
```

    ##      Phase  RID RID_chr SITEID VISCODE VISCODE2
    ## 4503 ADNI2 4103    4103     34     v01       sc
    ## 4711 ADNI2 4103    4103     34     v03       bl
    ## 5822 ADNI2 4103    4103     34     v05      m06
    ## 6366 ADNI2 4103    4103     34     v11      m12
    ## 8046 ADNI2 4103    4103     34     v21      m24

``` r
# viscode2 is between bl and m06
adnimerge_adni2[which(adnimerge_adni2$PTID == "072_S_4103"), c("RID", "PTID", "VISCODE", "DX_bl", "DX")]
```

    ##       RID       PTID VISCODE DX_bl DX
    ## 8911 4103 072_S_4103     m30    CN   
    ## 8912 4103 072_S_4103     m24    CN CN
    ## 8913 4103 072_S_4103     m18    CN   
    ## 8914 4103 072_S_4103     m12    CN CN
    ## 8915 4103 072_S_4103     m06    CN CN
    ## 8916 4103 072_S_4103     m03    CN   
    ## 8917 4103 072_S_4103      bl    CN CN

``` r
# at bl and m06, diagnosis is CN, so use CN as diagnosis
ADNI2_geneExpression_Dx_values[135] <- "CN"
```

``` r
# subject_361
adni_gene_expression_subjectInfo.split_phase.ADNI2[223, ]
```

    ##             Phase Visit  SubjectID RID_chr X260.280 X260.230 RIN Affy.Plate
    ## subject_361 ADNI2   v02 009_S_4564    4564     2.04     1.17 7.6          5
    ##             YearofCollection
    ## subject_361             2012

``` r
# visit v02, subject ID 009_S_4564
subject_viscode_check <- which(PT_RID_info.splitPhase.ADNI2$RID_chr == "4564")
# [1] 1023 4070
PT_RID_info.splitPhase.ADNI2[subject_viscode_check,]
```

    ##      Phase  RID RID_chr SITEID VISCODE VISCODE2
    ## 5366 ADNI2 4564    4564      6     v01       sc
    ## 8413 ADNI2 4564    4564      6     v03       bl

``` r
# viscode2 is between sc and bl
adnimerge_adni2[which(adnimerge_adni2$PTID == "009_S_4564"), c("RID", "PTID", "VISCODE", "DX_bl", "DX")]
```

    ##       RID       PTID VISCODE DX_bl  DX
    ## 4732 4564 009_S_4564      bl  LMCI MCI

``` r
# only diagnosis at bl available, diagnosis is MCI
ADNI2_geneExpression_Dx_values[223] <- "MCI"
```

``` r
# subject_386
adni_gene_expression_subjectInfo.split_phase.ADNI2[243, ]
```

    ##             Phase Visit  SubjectID RID_chr X260.280 X260.230 RIN Affy.Plate
    ## subject_386 ADNI2   v04 022_S_4266    4266     2.06     1.71 7.4          4
    ##             YearofCollection
    ## subject_386             2012

``` r
# visit v04, subject ID 022_S_4266
subject_viscode_check <- which(PT_RID_info.splitPhase.ADNI2$RID_chr == "4266")
# [1]  376  704 1589 2719 4048 5530
PT_RID_info.splitPhase.ADNI2[subject_viscode_check,]
```

    ##      Phase  RID RID_chr SITEID VISCODE VISCODE2
    ## 4719 ADNI2 4266    4266     16     v01       sc
    ## 5047 ADNI2 4266    4266     16     v03       bl
    ## 5932 ADNI2 4266    4266     16     v05      m06
    ## 7062 ADNI2 4266    4266     16     v11      m12
    ## 8391 ADNI2 4266    4266     16     v21      m24
    ## 9873 ADNI2 4266    4266     16     v41      m48

``` r
# viscode2 is between bl and m06
adnimerge_adni2[which(adnimerge_adni2$PTID == "022_S_4266"), c("RID", "PTID", "VISCODE", "DX_bl", "DX")]
```

    ##       RID       PTID VISCODE DX_bl DX
    ## 8225 4266 022_S_4266     m48    CN CN
    ## 8226 4266 022_S_4266     m24    CN CN
    ## 8227 4266 022_S_4266     m18    CN   
    ## 8228 4266 022_S_4266     m12    CN CN
    ## 8229 4266 022_S_4266     m06    CN CN
    ## 8230 4266 022_S_4266     m03    CN   
    ## 8231 4266 022_S_4266      bl    CN CN

``` r
# at bl and m06, diagnosis is CN, so use CN as diagnosis
ADNI2_geneExpression_Dx_values[243] <- "CN"
```

``` r
# subject_540
adni_gene_expression_subjectInfo.split_phase.ADNI2[335, ]
```

    ##             Phase Visit  SubjectID RID_chr X260.280 X260.230 RIN Affy.Plate
    ## subject_540 ADNI2   v04 114_S_4379    4379     2.05     1.38 5.7          5
    ##             YearofCollection
    ## subject_540             2012

``` r
# visit v04, subject ID 114_S_4379
subject_viscode_check <- which(PT_RID_info.splitPhase.ADNI2$RID_chr == "4379")
# [1]  604 1272 2462 3067
PT_RID_info.splitPhase.ADNI2[subject_viscode_check,]
```

    ##      Phase  RID RID_chr SITEID VISCODE VISCODE2
    ## 4947 ADNI2 4379    4379     42     v01       sc
    ## 5615 ADNI2 4379    4379     42     v03       bl
    ## 6805 ADNI2 4379    4379     42     v05      m06
    ## 7410 ADNI2 4379    4379     42     v11      m12

``` r
# viscode2 is between bl and m06
adnimerge_adni2[which(adnimerge_adni2$PTID == "114_S_4379"), c("RID", "PTID", "VISCODE", "DX_bl", "DX")]
```

    ##       RID       PTID VISCODE DX_bl       DX
    ## 6681 4379 114_S_4379     m18    AD         
    ## 6682 4379 114_S_4379     m12    AD Dementia
    ## 6683 4379 114_S_4379     m06    AD Dementia
    ## 6688 4379 114_S_4379     m03    AD         
    ## 6689 4379 114_S_4379      bl    AD Dementia

``` r
# at bl and m06, diagnosis is Dementia, so use Dementia as diagnosis
ADNI2_geneExpression_Dx_values[335] <- "Dementia"
```

``` r
# subject_636
adni_gene_expression_subjectInfo.split_phase.ADNI2[385, ]
```

    ##             Phase Visit  SubjectID RID_chr X260.280 X260.230 RIN Affy.Plate
    ## subject_636 ADNI2   v04 094_S_4282    4282     1.94     1.48 6.5          2
    ##             YearofCollection
    ## subject_636             2011

``` r
# visit v04, subject ID 094_S_4282
subject_viscode_check <- which(PT_RID_info.splitPhase.ADNI2$RID_chr == "4282")
# [1]  452  735 1291 3881
PT_RID_info.splitPhase.ADNI2[subject_viscode_check,]
```

    ##      Phase  RID RID_chr SITEID VISCODE VISCODE2
    ## 4795 ADNI2 4282    4282     37     v01       sc
    ## 5078 ADNI2 4282    4282     37     v03       bl
    ## 5634 ADNI2 4282    4282     37     v05      m06
    ## 8224 ADNI2 4282    4282     37     v21      m24

``` r
# viscode2 is between bl and m06
adnimerge_adni2[which(adnimerge_adni2$PTID == "094_S_4282"), c("RID", "PTID", "VISCODE", "DX_bl", "DX")]
```

    ##       RID       PTID VISCODE DX_bl       DX
    ## 8150 4282 094_S_4282     m30    AD         
    ## 8151 4282 094_S_4282     m24    AD Dementia
    ## 8152 4282 094_S_4282     m18    AD         
    ## 8153 4282 094_S_4282     m06    AD Dementia
    ## 8154 4282 094_S_4282     m03    AD         
    ## 8155 4282 094_S_4282      bl    AD Dementia

``` r
# at bl and m06, diagnosis is Dementia, so use Dementia as diagnosis
ADNI2_geneExpression_Dx_values[385] <- "Dementia"
```

``` r
# subject_679
adni_gene_expression_subjectInfo.split_phase.ADNI2[407, ]
```

    ##             Phase Visit  SubjectID RID_chr X260.280 X260.230 RIN Affy.Plate
    ## subject_679 ADNI2   v04 037_S_4432    4432     1.87     1.29 7.1          6
    ##             YearofCollection
    ## subject_679             2012

``` r
# visit v04, subject ID 037_S_4432
subject_viscode_check <- which(PT_RID_info.splitPhase.ADNI2$RID_chr == "4432")
# [1]  862 1185 1991 2785 4274 4854 5352
PT_RID_info.splitPhase.ADNI2[subject_viscode_check,]
```

    ##      Phase  RID RID_chr SITEID VISCODE VISCODE2
    ## 5205 ADNI2 4432    4432     26     v01       sc
    ## 5528 ADNI2 4432    4432     26     v03       bl
    ## 6334 ADNI2 4432    4432     26     v05      m06
    ## 7128 ADNI2 4432    4432     26     v11      m12
    ## 8617 ADNI2 4432    4432     26     v21      m24
    ## 9197 ADNI2 4432    4432     26     v31      m36
    ## 9695 ADNI2 4432    4432     26     v41      m48

``` r
# viscode2 is between bl and m06
adnimerge_adni2[which(adnimerge_adni2$PTID == "037_S_4432"), c("RID", "PTID", "VISCODE", "DX_bl", "DX")]
```

    ##        RID       PTID VISCODE DX_bl       DX
    ## 5748  4432 037_S_4432     m36  LMCI Dementia
    ## 5754  4432 037_S_4432     m24  LMCI Dementia
    ## 5755  4432 037_S_4432     m18  LMCI         
    ## 5756  4432 037_S_4432     m12  LMCI Dementia
    ## 5763  4432 037_S_4432     m06  LMCI Dementia
    ## 5764  4432 037_S_4432     m03  LMCI         
    ## 5765  4432 037_S_4432      bl  LMCI      MCI
    ## 12378 4432 037_S_4432     m48  LMCI Dementia

``` r
# at bl, diagnosis is MCI, but at m06, diagnosis is Dementia
# remove subject from analysis because diagnosis is unclear, tag as NA for now to prepare for removal
ADNI2_geneExpression_Dx_values[407] <- NA
```

## Get diagnosis information for ADNIGO subjects

``` r
adni_gene_expression_subjectInfo.split_phase$ADNIGO %>% head(10)
```

    ##             Phase Visit  SubjectID X260.280 X260.230 RIN Affy.Plate
    ## subject_1  ADNIGO   m48 116_S_1249     2.05     0.55 7.7          7
    ## subject_4  ADNIGO   m48 116_S_1232     2.03     1.52 6.8          7
    ## subject_8  ADNIGO    bl 003_S_2374     1.99     2.07 7.2          8
    ## subject_11 ADNIGO    bl 031_S_2018     1.97     1.71 7.6          1
    ## subject_18 ADNIGO   m48 128_S_0200     2.02     1.82 7.6          2
    ## subject_22 ADNIGO   m48 029_S_1218     2.01     1.88 6.9          5
    ## subject_25 ADNIGO    bl 021_S_2142     2.03     1.97 7.1          4
    ## subject_31 ADNIGO    bl 007_S_2394     2.03     1.01 7.3          9
    ## subject_32 ADNIGO   m36 128_S_1407     1.95     1.32 8.1          3
    ## subject_33 ADNIGO    bl 109_S_2200     2.06     1.71 7.4          5
    ##            YearofCollection
    ## subject_1              2011
    ## subject_4              2011
    ## subject_8              2011
    ## subject_11             2010
    ## subject_18             2010
    ## subject_22             2011
    ## subject_25             2010
    ## subject_31             2011
    ## subject_32             2010
    ## subject_33             2010

``` r
# since adnigo subjects already have viscode2 format for visit, can use visit column directly for matching with adnimerge dataset
# can match subjectID in gene expression data with PTID in adnimerge data
# adnimerge_phaseSplit$ADNIGO %>% head(10)
```

``` r
adni_gene_expression_subjectInfo.split_phase.ADNIGO <- adni_gene_expression_subjectInfo.split_phase$ADNIGO 
# adni_gene_expression_subjectInfo.split_phase.ADNIGO %>% head

adnimerge_ADNIGO <- adnimerge_phaseSplit$ADNIGO
```

``` r
# initialize empty vector to store diagnosis values
ADNIGO_geneExpression_Dx_values <- rep(0, nrow(adni_gene_expression_subjectInfo.split_phase.ADNIGO))
names(ADNIGO_geneExpression_Dx_values) <- rownames(adni_gene_expression_subjectInfo.split_phase.ADNIGO)

for (i in 1:length(ADNIGO_geneExpression_Dx_values)){
  
  subjectID_i <- adni_gene_expression_subjectInfo.split_phase.ADNIGO$SubjectID[i]
  visit_i <- adni_gene_expression_subjectInfo.split_phase.ADNIGO$Visit[i]
  
  # find index of subjectID in ADNIMERGE data
  adnimerge_subjectID_match_i <- which(adnimerge_ADNIGO$PTID == subjectID_i)
  # find index of visit in ADNIMERGE data
  adnimerge_visit_match_i <- which(adnimerge_ADNIGO$VISCODE == visit_i)
  # find intersection of index of subjectID and index of visit to give row number of diagnosis
  adnimerge_subjectID_visit_match_i <- intersect(adnimerge_subjectID_match_i, adnimerge_visit_match_i)
  # select row and diagnosis according to ADNIMERGE data
  subject_ADNIMERGE_Dx_i <- (adnimerge_ADNIGO[adnimerge_subjectID_visit_match_i,])$DX
  
  if(length(subject_ADNIMERGE_Dx_i) == 0){
    subject_ADNIMERGE_Dx_i = NA
  }
  
  # store diagnosis value into vector
  ADNIGO_geneExpression_Dx_values[i] = subject_ADNIMERGE_Dx_i
  
}
```

``` r
ADNIGO_geneExpression_Dx_values %>% unique
```

    ## [1] "CN"       "MCI"      "Dementia" ""

``` r
# [1] "CN"       "MCI"      "Dementia" ""        

which(ADNIGO_geneExpression_Dx_values == "")
```

    ## subject_60 
    ##         22

``` r
# subject_60 
#         22 

ADNIGO_geneExpression_Dx_values[22]
```

    ## subject_60 
    ##         ""

``` r
# check for diagnosis label for missing values
adni_gene_expression_subjectInfo.split_phase.ADNIGO[22,]
```

    ##             Phase Visit  SubjectID X260.280 X260.230 RIN Affy.Plate
    ## subject_60 ADNIGO   m03 073_S_2182     2.06     1.71 7.6          9
    ##            YearofCollection
    ## subject_60             2011

``` r
# subjectID 073_S_2182, visit m03
subject_adnimerge_match <- which(adnimerge_ADNIGO$PTID == "073_S_2182")
# [1] 752 753 754
adnimerge_ADNIGO[subject_adnimerge_match, c("PTID", "VISCODE", "DX_bl", "DX")]
```

    ##             PTID VISCODE DX_bl  DX
    ## 11824 073_S_2182      bl  EMCI MCI
    ## 11825 073_S_2182     m03  EMCI    
    ## 11826 073_S_2182     m06  EMCI MCI

``` r
# diagnosis at bl and m06 is MCI, so assign diagnosis MCI
ADNIGO_geneExpression_Dx_values[22] <- "MCI"
```

``` r
# create vector of new visit codes (for ADNI2 and ADNIGO collated)
ADNIGO_geneExpression_VISCODE2_values <- adni_gene_expression_subjectInfo.split_phase.ADNIGO$Visit
names(ADNIGO_geneExpression_VISCODE2_values) <- rownames(adni_gene_expression_subjectInfo.split_phase.ADNIGO)
viscode_adni2_adnigo_geneExpression <- c(ADNI2_geneExpression_VISCODE2_values, ADNIGO_geneExpression_VISCODE2_values)
viscode_adni2_adnigo_geneExpression <- viscode_adni2_adnigo_geneExpression[
  match(rownames(adni_gene_expression_values), 
        names(viscode_adni2_adnigo_geneExpression))]
viscode_adni2_adnigo_geneExpression %>% table
```

    ## .
    ##  bl m03 m06 m12 m36 m48 m60 m72 m84 
    ## 468   1   1   8   6 134 102  13   2

``` r
#  bl m03 m06 m12 m36 m48 m60 m72 m84 
# 468   1   1   8   6 134 102  13   2 
viscode_adni2_adnigo_geneExpression %>% unique
```

    ##  [1] "m48" "bl"  "m72" "m60" "m12" "m36" "m03" NA    "m06" "m84"

``` r
# create vector of new diagnosis codes (for ADNI2 and ADNIGO collated
Dx_adni2_adnigo_geneExpression <- c(ADNI2_geneExpression_Dx_values, ADNIGO_geneExpression_Dx_values)
Dx_adni2_adnigo_geneExpression <- Dx_adni2_adnigo_geneExpression[
  match(rownames(adni_gene_expression_values), 
        names(Dx_adni2_adnigo_geneExpression))
]
Dx_adni2_adnigo_geneExpression %>% table
```

    ## .
    ##       CN Dementia      MCI 
    ##      246      116      381

``` r
     #  CN Dementia      MCI 
     # 246      116      381 
Dx_adni2_adnigo_geneExpression %>% unique
```

    ## [1] "CN"       "Dementia" "MCI"      NA

``` r
# add viscode2 and diagnosis values back to the subjectInfo dataframe
adni_gene_expression_subjectInfo_adnimerge_Dx <- adni_gene_expression_subjectInfo %>%
  mutate(VISCODE2 = viscode_adni2_adnigo_geneExpression,  
         ADNIMERGE_Dx = Dx_adni2_adnigo_geneExpression, .after = SubjectID)

# adni_gene_expression_subjectInfo_adnimerge_Dx %>% head
```

## Add information about initial diagnosis, APOE genotype, age and sex

``` r
# adding information about age (baased on age data provided by Prof Nho Kwangsik)
ADNI_geneExpression_age <- read.csv("../data/ADNI_raw_files/ADNI_mRNA_Age.csv")
ADNI_geneExpression_age %>% head
```

    ##          IID Age
    ## 1 002_S_0413  82
    ## 2 002_S_0685  93
    ## 3 002_S_0729  69
    ## 4 002_S_1155  62
    ## 5 002_S_1261  76
    ## 6 002_S_1268  87

``` r
ADNI_geneExpression_age_ordered <- ADNI_geneExpression_age[match(adni_gene_expression_subjectInfo_adnimerge_Dx$SubjectID, ADNI_geneExpression_age$IID),]

adni_gene_expression_subjectInfo_adnimerge_Dx <- adni_gene_expression_subjectInfo_adnimerge_Dx %>%
  mutate(Age = ADNI_geneExpression_age_ordered$Age, .after=ADNIMERGE_Dx)

# adni_gene_expression_subjectInfo_adnimerge_Dx %>% head
```

``` r
# number of APOE4 alleles, age and sex from ADNIMERGE data sheet
adnimerge_selected_cols <- adnimerge[,c("COLPROT", "VISCODE", "PTID", "DX_bl", "PTGENDER", "APOE4", "PTRACCAT", "MMSE")]
# adnimerge_selected_cols %>% head
```

``` r
# initialize empty vectors
MMSE_vector <- rep(0, nrow(adni_gene_expression_subjectInfo_adnimerge_Dx))
DX_bl_vector <- rep(0, nrow(adni_gene_expression_subjectInfo_adnimerge_Dx))
PTGENDER_vector <- rep(0, nrow(adni_gene_expression_subjectInfo_adnimerge_Dx))
APOE4_vector <- rep(0, nrow(adni_gene_expression_subjectInfo_adnimerge_Dx))
PTRACCAT_vector <- rep(0, nrow(adni_gene_expression_subjectInfo_adnimerge_Dx))

for (i in 1:nrow(adni_gene_expression_subjectInfo_adnimerge_Dx)){
  
  # get information to match with adnimerge data
  phase_i <- adni_gene_expression_subjectInfo_adnimerge_Dx[i,]$Phase
  subjectID_i <- adni_gene_expression_subjectInfo_adnimerge_Dx[i,]$SubjectID
  viscode_i <- adni_gene_expression_subjectInfo_adnimerge_Dx[i,]$VISCODE2
  
  # match with adnimerge data
  phase_match_i <- which(adnimerge_selected_cols$COLPROT == phase_i)
  subjectID_match_i <- which(adnimerge_selected_cols$PTID == subjectID_i)
  viscode_match_i <- which(adnimerge_selected_cols$VISCODE == viscode_i)
  
  # for mmse, need to match phase, subjectID, and viscode
  row_intersect_phase_subjectID_viscode_i <- phase_match_i %>% 
    intersect(subjectID_match_i) %>% intersect(viscode_match_i)
  MMSE_i <- adnimerge_selected_cols[row_intersect_phase_subjectID_viscode_i,]$MMSE
  if (length(MMSE_i) == 0){
    MMSE_i = NA
  }
  MMSE_vector[i] <- MMSE_i
  
  # for the rest of the columns, only need a match for subjectID, use unique() to get only 1 value for each subject
  DX_bl_vector[i] <- adnimerge_selected_cols[subjectID_match_i,]$DX_bl %>% unique
  PTGENDER_vector[i] <- adnimerge_selected_cols[subjectID_match_i,]$PTGENDER %>% unique
  APOE4_vector[i] <- adnimerge_selected_cols[subjectID_match_i,]$APOE4 %>% unique
  PTRACCAT_vector[i] <- adnimerge_selected_cols[subjectID_match_i,]$PTRACCAT %>% unique
}
```

``` r
which(is.na(MMSE_vector))
```

    ##  [1]  60  73  79 124 210 361 386 487 540 636 679 683

``` r
# [1]  60  73  79 124 210 361 386 487 540 636 679 683

adni_gene_expression_subjectInfo_adnimerge_Dx[60,]
```

    ##             Phase Visit  SubjectID VISCODE2 ADNIMERGE_Dx Age X260.280 X260.230
    ## subject_60 ADNIGO   m03 073_S_2182      m03          MCI  66     2.06     1.71
    ##            RIN Affy.Plate YearofCollection
    ## subject_60 7.6          9             2011

``` r
# subject ID 073_S_2182, visit m03
# since MMSE at bl and m06 is 30, assume that MMSE at m03 is also 30
MMSE_vector[60] <- 30

adni_gene_expression_subjectInfo_adnimerge_Dx[73,]
```

    ##            Phase Visit  SubjectID VISCODE2 ADNIMERGE_Dx Age X260.280 X260.230
    ## subject_73 ADNI2   v04 072_S_4102     <NA>          MCI  67     2.06     1.64
    ##            RIN Affy.Plate YearofCollection
    ## subject_73 6.7          8             2011

``` r
# viscode2 is between bl and m06
# subject ID is 072_S_4102
# MMSE changes between bl and m06 from 30 to 25, so cannot assume MMSE at visit
# leave MMSE as NA

adni_gene_expression_subjectInfo_adnimerge_Dx[79,]
```

    ##            Phase Visit  SubjectID VISCODE2 ADNIMERGE_Dx Age X260.280 X260.230
    ## subject_79 ADNI2   v04 006_S_4192     <NA>     Dementia  82     2.08     0.28
    ##            RIN Affy.Plate YearofCollection
    ## subject_79 5.9          5             2011

``` r
# viscode2 is between bl and m06
# subject ID is 006_S_4192
# MMSE changes between bl and m06 from 19 to 22, so cannot assume MMSE at visit
# leave MMSE as NA

adni_gene_expression_subjectInfo_adnimerge_Dx[124,]
```

    ##             Phase Visit  SubjectID VISCODE2 ADNIMERGE_Dx Age X260.280 X260.230
    ## subject_124 ADNI2   v04 013_S_4268     <NA>          MCI  63     2.02     0.68
    ##             RIN Affy.Plate YearofCollection
    ## subject_124 7.3          2             2011

``` r
# viscode2 is between bl and m06
# subject ID is 013_S_4268
# MMSE changes between bl and m05 from 29 to 30

adni_gene_expression_subjectInfo_adnimerge_Dx[210,]
```

    ##             Phase Visit  SubjectID VISCODE2 ADNIMERGE_Dx Age X260.280 X260.230
    ## subject_210 ADNI2   v04 072_S_4103     <NA>           CN  71     2.05      1.9
    ##             RIN Affy.Plate YearofCollection
    ## subject_210 7.4          8             2011

``` r
# subject ID is 072_S_4103
# MMSE changes between bl and m05 from 29 to 28

adni_gene_expression_subjectInfo_adnimerge_Dx[361,]
```

    ##             Phase Visit  SubjectID VISCODE2 ADNIMERGE_Dx Age X260.280 X260.230
    ## subject_361 ADNI2   v02 009_S_4564     <NA>          MCI  65     2.04     1.17
    ##             RIN Affy.Plate YearofCollection
    ## subject_361 7.6          5             2012

``` r
# viscode2 is between sc and bl
# subject ID is 009_S_4564
# only MMSE available is at bl, MMSE = 27
# leave MMSE as NA

adni_gene_expression_subjectInfo_adnimerge_Dx[386,]
```

    ##             Phase Visit  SubjectID VISCODE2 ADNIMERGE_Dx Age X260.280 X260.230
    ## subject_386 ADNI2   v04 022_S_4266     <NA>           CN  71     2.06     1.71
    ##             RIN Affy.Plate YearofCollection
    ## subject_386 7.4          4             2012

``` r
# viscode2 is between bl and m06
# subject ID is 022_S_4266
# MMSE at bl and m06 = 30
MMSE_vector[386] <- 30

adni_gene_expression_subjectInfo_adnimerge_Dx[487,]
```

    ##              Phase Visit  SubjectID VISCODE2 ADNIMERGE_Dx Age X260.280 X260.230
    ## subject_487 ADNIGO   m48 128_S_0135      m48     Dementia  87     2.05     2.03
    ##             RIN Affy.Plate YearofCollection
    ## subject_487 6.8          2             2010

``` r
# Subject ID is 128_S_0135
# visit is m48
# no MMSE score available for visits before and after
# leave MMSE as NA

adni_gene_expression_subjectInfo_adnimerge_Dx[540,]
```

    ##             Phase Visit  SubjectID VISCODE2 ADNIMERGE_Dx Age X260.280 X260.230
    ## subject_540 ADNI2   v04 114_S_4379     <NA>     Dementia  88     2.05     1.38
    ##             RIN Affy.Plate YearofCollection
    ## subject_540 5.7          5             2012

``` r
# subject ID is 114_S_4379
# MMSE score at bl is 25, at m06 is 28
# leave MMSE as NA

adni_gene_expression_subjectInfo_adnimerge_Dx[636,]
```

    ##             Phase Visit  SubjectID VISCODE2 ADNIMERGE_Dx Age X260.280 X260.230
    ## subject_636 ADNI2   v04 094_S_4282     <NA>     Dementia  89     1.94     1.48
    ##             RIN Affy.Plate YearofCollection
    ## subject_636 6.5          2             2011

``` r
# subject ID is 094_S_4282
# MMSE score at bl is 23, at m06 is 21
# leave MMSE as NA

adni_gene_expression_subjectInfo_adnimerge_Dx[679,]
```

    ##             Phase Visit  SubjectID VISCODE2 ADNIMERGE_Dx Age X260.280 X260.230
    ## subject_679 ADNI2   v04 037_S_4432     <NA>         <NA>  62     1.87     1.29
    ##             RIN Affy.Plate YearofCollection
    ## subject_679 7.1          6             2012

``` r
# subject ID is 037_S_4432
# MMSE score at bl is 28, at m08 is 23
# leave MMSE as NA

adni_gene_expression_subjectInfo_adnimerge_Dx[683,]
```

    ##              Phase Visit  SubjectID VISCODE2 ADNIMERGE_Dx Age X260.280 X260.230
    ## subject_683 ADNIGO   m48 003_S_0907      m48           CN  92     2.04     0.74
    ##             RIN Affy.Plate YearofCollection
    ## subject_683 7.4          5             2010

``` r
# subject ID is 003_S_0907
# MMSE score consistently around 29 to 30
# assume MMSE is 30
MMSE_vector[683] <- 30
```

``` r
# adni_gene_expression_subjectInfo_adnimerge_Dx %>% head

adni_gene_expression_subjectInfo_allInfo <- adni_gene_expression_subjectInfo_adnimerge_Dx %>% mutate(
  DX_bl = DX_bl_vector, 
  PTGENDER = PTGENDER_vector, 
  APOE4 = APOE4_vector, 
  Race = PTRACCAT_vector, 
  MMSE = MMSE_vector, 
  .after = Age
)

# adni_gene_expression_subjectInfo_allInfo %>% head
```

``` r
write.csv(adni_gene_expression_subjectInfo_allInfo, file = "20230215_adni_gene_expression_subjectInfo_Dx.csv")
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
    ##  [1] plotly_4.10.1        Boruta_8.0.0         e1071_1.7-13        
    ##  [4] Metrics_0.1.4        tree_1.0-42          caret_6.0-93        
    ##  [7] lattice_0.20-45      randomForest_4.7-1.1 sva_3.46.0          
    ## [10] BiocParallel_1.32.5  genefilter_1.80.0    mgcv_1.8-41         
    ## [13] nlme_3.1-160         corrplot_0.92        polycor_0.8-1       
    ## [16] pheatmap_1.0.12      RColorBrewer_1.1-3   readxl_1.4.1        
    ## [19] forcats_0.5.2        stringr_1.5.0        dplyr_1.1.0         
    ## [22] purrr_1.0.1          readr_2.1.3          tidyr_1.3.0         
    ## [25] tibble_3.1.8         ggplot2_3.4.1        tidyverse_1.3.2     
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] backports_1.4.1        plyr_1.8.8             lazyeval_0.2.2        
    ##   [4] splines_4.2.2          listenv_0.9.0          GenomeInfoDb_1.34.9   
    ##   [7] digest_0.6.31          foreach_1.5.2          htmltools_0.5.4       
    ##  [10] fansi_1.0.4            magrittr_2.0.3         memoise_2.0.1         
    ##  [13] googlesheets4_1.0.1    tzdb_0.3.0             limma_3.54.1          
    ##  [16] recipes_1.0.3          globals_0.16.2         Biostrings_2.66.0     
    ##  [19] annotate_1.76.0        modelr_0.1.10          gower_1.0.0           
    ##  [22] matrixStats_0.63.0     hardhat_1.2.0          timechange_0.1.1      
    ##  [25] colorspace_2.1-0       blob_1.2.3             rvest_1.0.3           
    ##  [28] haven_2.5.1            xfun_0.37              crayon_1.5.2          
    ##  [31] RCurl_1.98-1.10        jsonlite_1.8.4         survival_3.5-5        
    ##  [34] iterators_1.0.14       glue_1.6.2             gtable_0.3.1          
    ##  [37] gargle_1.2.1           ipred_0.9-13           zlibbioc_1.44.0       
    ##  [40] XVector_0.38.0         future.apply_1.10.0    BiocGenerics_0.44.0   
    ##  [43] scales_1.2.1           DBI_1.1.3              edgeR_3.40.1          
    ##  [46] Rcpp_1.0.10            viridisLite_0.4.1      xtable_1.8-4          
    ##  [49] bit_4.0.5              proxy_0.4-27           stats4_4.2.2          
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
    ## [103] withr_2.5.0            S4Vectors_0.36.1       GenomeInfoDbData_1.2.9
    ## [106] parallel_4.2.2         hms_1.1.2              grid_4.2.2            
    ## [109] rpart_4.1.19           timeDate_4021.107      class_7.3-20          
    ## [112] rmarkdown_2.20         googledrive_2.0.0      pROC_1.18.0           
    ## [115] Biobase_2.58.0         lubridate_1.9.0
