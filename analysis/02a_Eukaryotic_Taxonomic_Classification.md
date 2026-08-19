---
title: "Eukaryotic Taxonomic Assignment" 
author: "Augustus Pendleton"
date: "16 April, 2025"
knitr: 
  opts_chunk: 
    fig.path: ../figures/02a_Eukaryotic_Taxonomic_Classificaiton/
editor_options: 
  chunk_output_type: console
---



# Purpose

The purpose of this document is to take the feature count (otu) table we generated with dada2 and assign taxonomy to our ASVs.

# Load packages, source files, and R objects



::: {.cell layout-align="center"}

```{.r .cell-code}
pacman::p_load(tidyverse, dada2, phyloseq, install = FALSE)

# Load in the colors and shapes
source("code/R/plotting_aesthetics.R")

source("code/R/notify_me.R")

# Load in our feature count table
load("data/01a_Eukaryotic_Denoising/seqtab_euk.RData")

# And load our meta_data
load("data/01a_Eukaryotic_Denoising/track_df_euk.RData")
```
:::


# Assign taxonomy using Dada/Silva

Here, we will use the silva database version 138!

::: {.cell layout-align="center"}

```{.r .cell-code}
# This takes awhile
taxa_wholeread <- assignTaxonomy(seqtab_euk, "code/mcnichol_split_scripts/pr2_database/pr2_version_5.1.0_SSU_dada2.fasta.gz",
                       multithread=30,
                       tryRC = TRUE,
                       taxLevels = c("Domain","Supergroup","Division","Subdivision", "Class","Order","Family","Genus","Species"))

seqtab_copy <- seqtab_euk

seqs <- colnames(seqtab_copy)

short_seqs <- substr(seqs, 1, 281)

colnames(seqtab_copy) <- short_seqs

taxa_forwardread <- assignTaxonomy(seqtab_copy, "code/mcnichol_split_scripts/pr2_database/pr2_version_5.1.0_SSU_dada2.fasta.gz",
                       multithread=30,
                       tryRC = TRUE,
                       taxLevels = c("Domain","Supergroup","Division","Subdivision", "Class","Order","Family","Genus","Species"))
```
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
identical(taxa_wholeread, taxa_forwardread)
```

::: {.cell-output .cell-output-stdout}
```
[1] FALSE
```
:::

```{.r .cell-code}
nrow(taxa_wholeread)
```

::: {.cell-output .cell-output-stdout}
```
[1] 1622
```
:::

```{.r .cell-code}
taxa_wholeread %>%
  data.frame() %>%
  summarize(across(everything(), \(x)sum(is.na(x))))
```

::: {.cell-output .cell-output-stdout}
```
  Domain Supergroup Division Subdivision Class Order Family Genus Species
1      2         17       18          37    71   126    179   299     376
```
:::

```{.r .cell-code}
taxa_forwardread %>%
  data.frame() %>%
  summarize(across(everything(), \(x)sum(is.na(x))))
```

::: {.cell-output .cell-output-stdout}
```
  Domain Supergroup Division Subdivision Class Order Family Genus Species
1      5         25       32          74   117   197    240   378     481
```
:::

```{.r .cell-code}
map(c(1:9), \(x)sum(taxa_wholeread[,x] == taxa_forwardread[,x], na.rm = TRUE))
```

::: {.cell-output .cell-output-stdout}
```
[[1]]
[1] 1616

[[2]]
[1] 1583

[[3]]
[1] 1576

[[4]]
[1] 1531

[[5]]
[1] 1485

[[6]]
[1] 1391

[[7]]
[1] 1333

[[8]]
[1] 1163

[[9]]
[1] 1044
```
:::

```{.r .cell-code}
taxa_wholeread %>%
  count(Domain)
```

::: {.cell-output .cell-output-error}
```
Error in UseMethod("count"): no applicable method for 'count' applied to an object of class "c('matrix', 'array', 'character')"
```
:::

```{.r .cell-code}
taxa_forwardread %>%
  count(Domain)
```

::: {.cell-output .cell-output-error}
```
Error in UseMethod("count"): no applicable method for 'count' applied to an object of class "c('matrix', 'array', 'character')"
```
:::
:::


Okay, using the whole read seems to be much more effective. Tbh I don't know if I understand the classifiers, cause it really seems like I would need to split this read, but something else to discuss. 


# Taxonomy Table 


::: {.cell layout-align="center"}

```{.r .cell-code}
##### Prepare tax table 
# Add the ASV sequences from the rownames to a column 
new_tax_tab <- taxa_wholeread %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASVseqs") 
head(new_tax_tab)
```

::: {.cell-output .cell-output-stdout}
```
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            ASVseqs
1 TTCCAGCTCCAATAGCGTATATTAAAGTTGTTGTGGTTAAAAAGCTCGTAGTTGGATCTCGGCGGACAGAAGTCGGTCCGAATGCTTTCGGTACTGGCTTGTTTGTCCGTGTGTTTTGCCAGAGGCAGGCGGGTGCTCTTCGCCGAGTGTCCGTTTGTGCTGGCAGGTTTACTTTGAAAAAATTAGAGTGCTCAAAGCAAGCTCTTACGCTTGAATATTCGTGCATGGAATAATAGAATAGGATGTCGTTTCTATTTTGTTGGTTTTCGGAAATCGACATATACTGCGAAAGCATTTGCCAAGAATGTTTTCATTAATCAAGAACGAAAGTTAGAGGTTCGAAGGCGATTAGATACCGCCCTAGTTCTAACCATAAACGATGCCAGCTAGCGATCTGCAGGCGTTTTATATTACGACCTTGCAGGGAGCTTCCGGGAAACCAAAGCTTGTGGGTTCCAGGGGAAGTATGGTTGCAAAGCTG
2 TTCCAGCTCCAATAGCGTATATTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTTGGATTTCTGGTAGGGAGTCAATGCCCGCTTGATGCGTGCGCAGAGGTTCCTTGCCATCCTTCTGTTAACCCTCGGACCTTCACTGGTCCGGTGGGAGCCCAGATATTTTACCTTGAGAAAATTAGAGTGTTTCAGGCAGGCTAGTGCCGGGATACATTAGCATGGAATAATAGAATAGGACTTTTTCCTCTTCTGTTGGTCTGAGGGGATGAAGTAATGATTAATAGGCTAACTTATGCGAAAGCATTTGCCAAGGATGTTTTCATTAATCAAGAACGAAAGTTAGGGGATCAAAGACGATCAGATACCGTCCTAGTCTTAACCATAAACTATGCCGACTCGGGATCGGGGGCAAGAATTATCTGCCTTCGGCACCGTATGAGAAATCAAAGTCTTTGGGTTCTGGGGGGAGTATGGTCGCAAGGCTG
3 TTCCAGCTCTAATAGCGTATATTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTCGGATGTCGGGCTCGGGCAAGCTGTCGGCCTTTGGTCGGACGGCAGGCTCGGGTCTTTCTGCCTGAGGAACCCGCTTGCTTTAACCGGCTGGCGGTGGACGCAGGTCGTTTACTTTGAAAAAATTAGAGTGTTCAAAGCAGGCTAGCGCTTGAATACATTAGCATGGAATAATGGAATAGGACTTTGGTGCTATTTTGTTGGTTTATGGGACCGAAGTAATGATTAACAACTTCTGCGAAAGCATTCGGCAAGGATGTTTTCATTGATCAAGAACGAAAGTTAGGGGATCGAAGACGATCAGATACCGTCGTAGTCTTAACCATAAACTATGCCGACTAGGGATCAGTGGATGACAATTTGCGTCTCCATTGGCACCTTGTGAGAAATCAAAGTTTTTGGGTTCCGGGGGGAGTATGGTCGCAAGGCTG
4 TTCCAGCTCTAATAGCATATATTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTCGGATTTTGGGCGCGTGTGAACTGTCAGCCTCTGGTTGGACGGTGCGCGCGTGCCTTTCTGCTTGGGAGCCCCATGGTGGTTAATTCTGTCGTGGTGGGTCACAGGCCGTTTACTTTGAAAAAATTAGAGTGTTCAAAGCAGGCTGTTGCTTGAATAGATTAGCATGGAATAATGGAATAGGACCTTGGTGCTATTTCGTTGGTTATGGGACTGGGGTAATGATTAATAACTTGTGCGAAAGCATTCGGCAAGGATGTTTTCATTGATCAAGAACGAAAGTTAGGGGATCGAAGACGATCAGATACCGTCGTAGTCTTAACCATAAACTATGCCAACTAGGGATCATTGAGTGTCAATTTACGACCTCATTGGCACCTTGTGAGAAATCAAAGTTTTTGGGTTCCGGGGGGAGTATGGTCGCAAGGCTG
5 TTCCAGCTCCAATAGCGTATATTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTTGGATCTCGGGCGCCGGTTTGTGGTACACTTTACTGTTGTTACTGCATGTCCGGCCCTAATCAATCGGTTGGTCGTTTATGCCCTTCACTGGGTGTTTTCGGTAGCCGGTACGTTTACTTTGAGAAAATTAGAGTGCTTAAAGCAGGCCTTACGCCTGAATAATGTTGCATGGAATAATAGAATAGGACCTCGGTTTTATTTTGTTGGTTTTGTAGAGCATTGAGGTAACAACTGCGAAAGCATTTGCCAAGAATGTTTTCATTAATCAAGAACGAAAGTTGGAGGTTCGAAGACGATTAGATACCGTCCTAGTTCCAACCGTAAACGATGCCAACTTACAATCCGTCGCACTTAATTAAATGAGTCGGCGGGCGGTATCCGGGAAACCAAAGTCTTTGGGTTCCGGGGGAAGTATGGTTGCAAAGCTG
6 TTCCAGCTCCAATAGCGTATATTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTTGGATTTCTGACAGGGTGACCTCGCCCGCTTATTGCGTGCGCGAAGCCGCTCTGTCATCCTTCTTTTAACCCTTTAGACCTTAATCAGTTTGCTGGGGGCTAAGATATTTTACCTTGAGAAAATTAGAGTGTTTCAGGCAGGCTAGCGCCGGAATACATTAGCATGGAATAATAGAATAGGACTTTTACCTCTTCTGTTGGTCTGAGGGGTTCAAGTAATGATTAATAGCTAACTTATGCGAAAGCATTTGCCAAGGATGTTTTCATTAATCAAGAACGAAAGTTAGGGGATCAAAGACGATCAGATACCGTCCTAGTCTTAACCATAAACTATGCCGACTCGGGATCGGGGGCGAGAGTTATCCGCCTTCGGCACCGTATGAGAAATCAAAGTCTTTGGGTTCTGGGGGGAGTATGGTCGCAAGGCTG
     Domain Supergroup     Division   Subdivision         Class           Order
1 Eukaryota     Obazoa Opisthokonta       Metazoa    Arthropoda       Crustacea
2 Eukaryota       TSAR    Alveolata    Ciliophora  Spirotrichea  Choreotrichida
3 Eukaryota  Cryptista  Cryptophyta Cryptophyta_X Cryptophyceae Cryptomonadales
4 Eukaryota  Cryptista  Cryptophyta Cryptophyta_X Cryptophyceae Cryptomonadales
5 Eukaryota     Obazoa Opisthokonta       Metazoa      Rotifera      Rotifera_X
6 Eukaryota       TSAR    Alveolata    Ciliophora  Spirotrichea  Choreotrichida
             Family           Genus              Species
1       Maxillopoda            <NA>                 <NA>
2    Strobilidiidae Rimostrombidium  Rimostrombidium_sp.
3 Cryptomonadales_X    Plagioselmis                 <NA>
4 Cryptomonadales_X     Cryptomonas  Cryptomonas_curvata
5       Rotifera_XX      Conochilus Conochilus_unicornis
6    Strobilidiidae Rimostrombidium  Rimostrombidium_sp.
```
:::

```{.r .cell-code}
# intution check 
stopifnot(new_tax_tab$ASVseqs == colnames(seqtab_euk))

asv_headers <- vector(dim(seqtab_euk)[2], mode = "character")

# loop through vector and fill it in with ASV names 

for (i in 1:dim(seqtab_euk)[2]) {
  asv_headers[i] <- paste(">EUK", i, sep = "_")
}

# intitution check
head(asv_headers, 20)
```

::: {.cell-output .cell-output-stdout}
```
 [1] ">EUK_1"  ">EUK_2"  ">EUK_3"  ">EUK_4"  ">EUK_5"  ">EUK_6"  ">EUK_7" 
 [8] ">EUK_8"  ">EUK_9"  ">EUK_10" ">EUK_11" ">EUK_12" ">EUK_13" ">EUK_14"
[15] ">EUK_15" ">EUK_16" ">EUK_17" ">EUK_18" ">EUK_19" ">EUK_20"
```
:::

```{.r .cell-code}
##### Rename ASVs in table then write out our ASV fasta file! 
#View(seqtab_euk)
asv_tab <- t(seqtab_euk)
#View(asv_tab)

## Rename our asvs! 
row.names(asv_tab) <- sub(">", "", asv_headers)

# Now let's add the ASV names 
rownames(new_tax_tab) <- rownames(asv_tab)
#View(new_tax_tab)

### Final prep of tax table. Add new column with ASV names 
asv_tax <- 
  new_tax_tab %>%
  # add rownames from count table for phyloseq handoff
  mutate(ASV = rownames(asv_tab)) %>%
  # Resort the columns with select
  dplyr::select(Domain:Species,ASV, ASVseqs)

# Intution check
stopifnot(asv_tax$ASV == rownames(asv_tax), rownames(asv_tax) == rownames(asv_tab))
```
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
# Write the table 
write.table(asv_tax, "data/02a_Eukaryotic_Taxonomic_Classification/ASV_final_taxonomy.tsv", sep = "\t", quote = FALSE, col.names = NA)
```
:::



## Handoff to phyloseq 

::: {.cell layout-align="center"}

```{.r .cell-code}
raw_euk_physeq <- phyloseq(otu_table(asv_tab, taxa_are_rows = TRUE),
                         tax_table(as.matrix(asv_tax)))

raw_euk_physeq
```

::: {.cell-output .cell-output-stdout}
```
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 1622 taxa and 76 samples ]
tax_table()   Taxonomy Table:    [ 1622 taxa by 11 taxonomic ranks ]
```
:::

```{.r .cell-code}
tax_table(raw_euk_physeq) %>% head()
```

::: {.cell-output .cell-output-stdout}
```
Taxonomy Table:     [6 taxa by 11 taxonomic ranks]:
      Domain      Supergroup  Division       Subdivision     Class          
EUK_1 "Eukaryota" "Obazoa"    "Opisthokonta" "Metazoa"       "Arthropoda"   
EUK_2 "Eukaryota" "TSAR"      "Alveolata"    "Ciliophora"    "Spirotrichea" 
EUK_3 "Eukaryota" "Cryptista" "Cryptophyta"  "Cryptophyta_X" "Cryptophyceae"
EUK_4 "Eukaryota" "Cryptista" "Cryptophyta"  "Cryptophyta_X" "Cryptophyceae"
EUK_5 "Eukaryota" "Obazoa"    "Opisthokonta" "Metazoa"       "Rotifera"     
EUK_6 "Eukaryota" "TSAR"      "Alveolata"    "Ciliophora"    "Spirotrichea" 
      Order             Family              Genus            
EUK_1 "Crustacea"       "Maxillopoda"       NA               
EUK_2 "Choreotrichida"  "Strobilidiidae"    "Rimostrombidium"
EUK_3 "Cryptomonadales" "Cryptomonadales_X" "Plagioselmis"   
EUK_4 "Cryptomonadales" "Cryptomonadales_X" "Cryptomonas"    
EUK_5 "Rotifera_X"      "Rotifera_XX"       "Conochilus"     
EUK_6 "Choreotrichida"  "Strobilidiidae"    "Rimostrombidium"
      Species                ASV    
EUK_1 NA                     "EUK_1"
EUK_2 "Rimostrombidium_sp."  "EUK_2"
EUK_3 NA                     "EUK_3"
EUK_4 "Cryptomonas_curvata"  "EUK_4"
EUK_5 "Conochilus_unicornis" "EUK_5"
EUK_6 "Rimostrombidium_sp."  "EUK_6"
      ASVseqs                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
EUK_1 "TTCCAGCTCCAATAGCGTATATTAAAGTTGTTGTGGTTAAAAAGCTCGTAGTTGGATCTCGGCGGACAGAAGTCGGTCCGAATGCTTTCGGTACTGGCTTGTTTGTCCGTGTGTTTTGCCAGAGGCAGGCGGGTGCTCTTCGCCGAGTGTCCGTTTGTGCTGGCAGGTTTACTTTGAAAAAATTAGAGTGCTCAAAGCAAGCTCTTACGCTTGAATATTCGTGCATGGAATAATAGAATAGGATGTCGTTTCTATTTTGTTGGTTTTCGGAAATCGACATATACTGCGAAAGCATTTGCCAAGAATGTTTTCATTAATCAAGAACGAAAGTTAGAGGTTCGAAGGCGATTAGATACCGCCCTAGTTCTAACCATAAACGATGCCAGCTAGCGATCTGCAGGCGTTTTATATTACGACCTTGCAGGGAGCTTCCGGGAAACCAAAGCTTGTGGGTTCCAGGGGAAGTATGGTTGCAAAGCTG"
EUK_2 "TTCCAGCTCCAATAGCGTATATTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTTGGATTTCTGGTAGGGAGTCAATGCCCGCTTGATGCGTGCGCAGAGGTTCCTTGCCATCCTTCTGTTAACCCTCGGACCTTCACTGGTCCGGTGGGAGCCCAGATATTTTACCTTGAGAAAATTAGAGTGTTTCAGGCAGGCTAGTGCCGGGATACATTAGCATGGAATAATAGAATAGGACTTTTTCCTCTTCTGTTGGTCTGAGGGGATGAAGTAATGATTAATAGGCTAACTTATGCGAAAGCATTTGCCAAGGATGTTTTCATTAATCAAGAACGAAAGTTAGGGGATCAAAGACGATCAGATACCGTCCTAGTCTTAACCATAAACTATGCCGACTCGGGATCGGGGGCAAGAATTATCTGCCTTCGGCACCGTATGAGAAATCAAAGTCTTTGGGTTCTGGGGGGAGTATGGTCGCAAGGCTG"
EUK_3 "TTCCAGCTCTAATAGCGTATATTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTCGGATGTCGGGCTCGGGCAAGCTGTCGGCCTTTGGTCGGACGGCAGGCTCGGGTCTTTCTGCCTGAGGAACCCGCTTGCTTTAACCGGCTGGCGGTGGACGCAGGTCGTTTACTTTGAAAAAATTAGAGTGTTCAAAGCAGGCTAGCGCTTGAATACATTAGCATGGAATAATGGAATAGGACTTTGGTGCTATTTTGTTGGTTTATGGGACCGAAGTAATGATTAACAACTTCTGCGAAAGCATTCGGCAAGGATGTTTTCATTGATCAAGAACGAAAGTTAGGGGATCGAAGACGATCAGATACCGTCGTAGTCTTAACCATAAACTATGCCGACTAGGGATCAGTGGATGACAATTTGCGTCTCCATTGGCACCTTGTGAGAAATCAAAGTTTTTGGGTTCCGGGGGGAGTATGGTCGCAAGGCTG"
EUK_4 "TTCCAGCTCTAATAGCATATATTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTCGGATTTTGGGCGCGTGTGAACTGTCAGCCTCTGGTTGGACGGTGCGCGCGTGCCTTTCTGCTTGGGAGCCCCATGGTGGTTAATTCTGTCGTGGTGGGTCACAGGCCGTTTACTTTGAAAAAATTAGAGTGTTCAAAGCAGGCTGTTGCTTGAATAGATTAGCATGGAATAATGGAATAGGACCTTGGTGCTATTTCGTTGGTTATGGGACTGGGGTAATGATTAATAACTTGTGCGAAAGCATTCGGCAAGGATGTTTTCATTGATCAAGAACGAAAGTTAGGGGATCGAAGACGATCAGATACCGTCGTAGTCTTAACCATAAACTATGCCAACTAGGGATCATTGAGTGTCAATTTACGACCTCATTGGCACCTTGTGAGAAATCAAAGTTTTTGGGTTCCGGGGGGAGTATGGTCGCAAGGCTG"
EUK_5 "TTCCAGCTCCAATAGCGTATATTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTTGGATCTCGGGCGCCGGTTTGTGGTACACTTTACTGTTGTTACTGCATGTCCGGCCCTAATCAATCGGTTGGTCGTTTATGCCCTTCACTGGGTGTTTTCGGTAGCCGGTACGTTTACTTTGAGAAAATTAGAGTGCTTAAAGCAGGCCTTACGCCTGAATAATGTTGCATGGAATAATAGAATAGGACCTCGGTTTTATTTTGTTGGTTTTGTAGAGCATTGAGGTAACAACTGCGAAAGCATTTGCCAAGAATGTTTTCATTAATCAAGAACGAAAGTTGGAGGTTCGAAGACGATTAGATACCGTCCTAGTTCCAACCGTAAACGATGCCAACTTACAATCCGTCGCACTTAATTAAATGAGTCGGCGGGCGGTATCCGGGAAACCAAAGTCTTTGGGTTCCGGGGGAAGTATGGTTGCAAAGCTG"
EUK_6 "TTCCAGCTCCAATAGCGTATATTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTTGGATTTCTGACAGGGTGACCTCGCCCGCTTATTGCGTGCGCGAAGCCGCTCTGTCATCCTTCTTTTAACCCTTTAGACCTTAATCAGTTTGCTGGGGGCTAAGATATTTTACCTTGAGAAAATTAGAGTGTTTCAGGCAGGCTAGCGCCGGAATACATTAGCATGGAATAATAGAATAGGACTTTTACCTCTTCTGTTGGTCTGAGGGGTTCAAGTAATGATTAATAGCTAACTTATGCGAAAGCATTTGCCAAGGATGTTTTCATTAATCAAGAACGAAAGTTAGGGGATCAAAGACGATCAGATACCGTCCTAGTCTTAACCATAAACTATGCCGACTCGGGATCGGGGGCGAGAGTTATCCGCCTTCGGCACCGTATGAGAAATCAAAGTCTTTGGGTTCTGGGGGGAGTATGGTCGCAAGGCTG"
```
:::

```{.r .cell-code}
otu_table(raw_euk_physeq) %>% head()
```

::: {.cell-output .cell-output-stdout}
```
OTU Table:          [6 taxa and 76 samples]
                     taxa are rows
      KH_D001 KH_D002 KH_D003 KH_D004 KH_D005 KH_D006 KH_D007 KH_D008 KH_D009
EUK_1      63     148     268      33     552       0    6962      23    1952
EUK_2       0       0       0       0       0       0       0       0       0
EUK_3      10       0      19       0      60      51      57      90      63
EUK_4      31       0       3       0     185     458      65     206     192
EUK_5       0       0       0       0       0       0       0       0       0
EUK_6      14       0       2       0      58      81      46      39      15
      KH_D010 KH_D011 KH_D012 KH_D013 KH_D014 KH_D015 KH_D016 KH_D017 KH_D018
EUK_1    1627       0       0       0       0     122      82     161      40
EUK_2       0      13      22     277       0       0       0       0       2
EUK_3     258     226     278    1262       0      32       9     118     215
EUK_4     387     184       0      66       0      11       8      17       0
EUK_5       0       0       0       0       0       0       0       0       0
EUK_6      86      73      67      23       0       0       9      10       0
      KH_D019 KH_D020 KH_D021 KH_D022 KH_D023 KH_D025 KH_D026 KH_D027 KH_D028
EUK_1      23      48     117       0       0       5     460     183    1025
EUK_2       0       0       0     121      43       0      27     195     374
EUK_3      19      39     542     288     195      70     145     178      80
EUK_4      19     192      26     352      64     114     113      70     154
EUK_5       0       0       0       0      72       0     200     364     668
EUK_6      12      12     238      76      23       0      25     160      76
      KH_D029 KH_D030 KH_D031 KH_D032 KH_D033 KH_D034 KH_D035 KH_D036 KH_D037
EUK_1       0      57       0       0       0       3      10       0    2059
EUK_2      15     647       0       0     428    1348     525       6     365
EUK_3       0     167       0      14     296      56     207       9     228
EUK_4      15     312       0       8     109      22     168      21     202
EUK_5       0     212       0       0      10       0      50       0    1328
EUK_6       0      49       0       0      95      83     101       0     104
      KH_D038 KH_D039 KH_D040 KH_D041 KH_D042 KH_D043 KH_D044 KH_D045 KH_D046
EUK_1     442       0    4997     212      59    2697     488     882      23
EUK_2       6     866      51     118      38       0     261     317     189
EUK_3       0      38     258     173      43     248     125      75     106
EUK_4      35      38      81      77      20      25     101      37      40
EUK_5      21       0     150     350      20      26     413     139      47
EUK_6       7       6      41      23      12      24     118      42      98
      KH_D047 KH_D051 KH_D052 KH_D053 KH_D054 KH_D055 KH_D056 KH_D057 KH_D058
EUK_1       0    2943       0     933     149       0     200       0     413
EUK_2     307     216       0     123     153       0      93       0     132
EUK_3      75     324      15     383     266       7     141      10      94
EUK_4     108     482     106      98     170      73      56      14      42
EUK_5     111     102       0      28      10       0      90       0      14
EUK_6       0      65       0      46      27       0      16       0     282
      KH_D059 KH_D060 KH_D061 KH_D062 KH_D063 KH_D064 KH_D065 KH_D066 KH_D067
EUK_1       0    2181       0     154    3160       0       0       0       0
EUK_2     293      55       0     358     156     496     287       0     306
EUK_3     203     141       0     226      90     150     169       0     146
EUK_4     339      30       0      51      22      91     142       0       0
EUK_5       0     351       0       0       0     104      65       0       0
EUK_6      66      21       0     665     291      87       7       0     137
      KH_D069 KH_D070 KH_D071 KH_D074 KH_D075 KH_D076 KH_D077 KH_D078 KH_D079
EUK_1      31     292      44       3     734    7061     333     400       0
EUK_2      65     317       0      18     108       4     108     518       0
EUK_3       0      89       0      51     137       0      44     359       0
EUK_4       0      34       0      15      88       0      44     199       0
EUK_5       0      76       0      21      15     172      77       0       0
EUK_6       5      53       0      11      44       0      37     153       0
      KH_D080 KH_D082 KH_D083 ZYMO_MOCK
EUK_1       0       4     674         0
EUK_2       0       0       4         0
EUK_3       0       0      11         0
EUK_4       0       0       0         0
EUK_5       0       0       0         0
EUK_6       0       0      42         0
```
:::

```{.r .cell-code}
save(raw_euk_physeq, file = "data/02a_Eukaryotic_Taxonomic_Classification/raw_euk_physeq.RData")
```
:::


# Session Information 

::: {.cell layout-align="center"}

```{.r .cell-code}
# Reproducibility
devtools::session_info()
```

::: {.cell-output .cell-output-stdout}
```
─ Session info ───────────────────────────────────────────────────────────────
 setting  value
 version  R version 4.3.3 (2024-02-29)
 os       Rocky Linux 9.4 (Blue Onyx)
 system   x86_64, linux-gnu
 ui       X11
 language (EN)
 collate  en_US.UTF-8
 ctype    en_US.UTF-8
 tz       America/New_York
 date     2025-04-16
 pandoc   3.1.1 @ /usr/lib/rstudio-server/bin/quarto/bin/tools/ (via rmarkdown)
 quarto   1.3.450 @ /usr/lib/rstudio-server/bin/quarto/bin/quarto

─ Packages ───────────────────────────────────────────────────────────────────
 ! package              * version    date (UTC) lib source
   abind                  1.4-8      2024-09-12 [1] RSPM (R 4.3.3)
   ade4                   1.7-23     2025-02-14 [1] RSPM (R 4.3.3)
   ape                    5.8-1      2024-12-16 [1] RSPM (R 4.3.3)
   Biobase                2.62.0     2023-10-24 [1] Bioconductor
   BiocGenerics           0.48.1     2023-11-01 [1] Bioconductor
   BiocParallel           1.36.0     2023-10-24 [1] Bioconductor
   biomformat             1.30.0     2023-10-24 [1] Bioconductor
   Biostrings             2.70.3     2024-03-13 [1] Bioconductor 3.18 (R 4.3.3)
   bitops                 1.0-9      2024-10-03 [1] RSPM (R 4.3.3)
   cachem                 1.1.0      2024-05-16 [1] RSPM (R 4.3.3)
   cli                    3.6.4      2025-02-13 [1] RSPM (R 4.3.3)
 P cluster                2.1.6      2023-12-01 [?] CRAN (R 4.3.3)
 P codetools              0.2-19     2023-02-01 [?] CRAN (R 4.3.3)
   colorspace             2.1-1      2024-07-26 [1] RSPM (R 4.3.3)
   crayon                 1.5.3      2024-06-20 [1] RSPM (R 4.3.3)
   dada2                * 1.30.0     2023-10-24 [1] Bioconductor
   data.table             1.17.0     2025-02-22 [1] RSPM (R 4.3.3)
   DelayedArray           0.28.0     2023-10-24 [1] Bioconductor
   deldir                 2.0-4      2024-02-28 [1] RSPM (R 4.3.3)
 P devtools               2.4.5      2022-10-11 [?] CRAN (R 4.3.2)
   digest                 0.6.37     2024-08-19 [1] RSPM (R 4.3.3)
   dplyr                * 1.1.4      2023-11-17 [1] RSPM (R 4.3.3)
 P ellipsis               0.3.2      2021-04-29 [?] CRAN (R 4.3.2)
   evaluate               1.0.3      2025-01-10 [1] RSPM (R 4.3.3)
   fastmap                1.2.0      2024-05-15 [1] RSPM (R 4.3.3)
   forcats              * 1.0.0      2023-01-29 [1] RSPM (R 4.3.3)
   foreach                1.5.2      2022-02-02 [1] RSPM (R 4.3.3)
   fs                     1.6.5      2024-10-30 [1] RSPM (R 4.3.3)
   generics               0.1.3      2022-07-05 [1] RSPM (R 4.3.3)
   GenomeInfoDb           1.38.8     2024-03-15 [1] Bioconductor 3.18 (R 4.3.3)
   GenomeInfoDbData       1.2.11     2025-04-11 [1] Bioconductor
   GenomicAlignments      1.38.2     2024-01-16 [1] Bioconductor 3.18 (R 4.3.3)
   GenomicRanges          1.54.1     2023-10-29 [1] Bioconductor
   ggplot2              * 3.5.2      2025-04-09 [1] RSPM (R 4.3.3)
   glue                   1.8.0      2024-09-30 [1] RSPM (R 4.3.3)
   gtable                 0.3.6      2024-10-25 [1] RSPM (R 4.3.3)
   hms                    1.1.3      2023-03-21 [1] RSPM (R 4.3.3)
   htmltools              0.5.8.1    2024-04-04 [1] RSPM (R 4.3.3)
 P htmlwidgets            1.6.4      2023-12-06 [?] CRAN (R 4.3.2)
 P httpuv                 1.6.15     2024-03-26 [?] CRAN (R 4.3.2)
   hwriter                1.3.2.1    2022-04-08 [1] RSPM (R 4.3.3)
   igraph                 2.1.4      2025-01-23 [1] RSPM (R 4.3.3)
   interp                 1.1-6      2024-01-26 [1] RSPM (R 4.3.3)
   IRanges                2.36.0     2023-10-24 [1] Bioconductor
   iterators              1.0.14     2022-02-05 [1] RSPM (R 4.3.3)
   jpeg                   0.1-11     2025-03-21 [1] RSPM (R 4.3.3)
   jsonlite               2.0.0      2025-03-27 [1] RSPM (R 4.3.3)
   knitr                  1.50       2025-03-16 [1] RSPM (R 4.3.3)
 P later                  1.4.2      2025-04-08 [?] RSPM (R 4.3.3)
 P lattice                0.22-5     2023-10-24 [?] CRAN (R 4.3.3)
   latticeExtra           0.6-30     2022-07-04 [1] RSPM (R 4.3.3)
   lifecycle              1.0.4      2023-11-07 [1] RSPM (R 4.3.3)
   lubridate            * 1.9.4      2024-12-08 [1] RSPM (R 4.3.3)
   magrittr               2.0.3      2022-03-30 [1] RSPM (R 4.3.3)
 P MASS                   7.3-60.0.1 2024-01-13 [?] CRAN (R 4.3.3)
 P Matrix                 1.6-5      2024-01-11 [?] CRAN (R 4.3.3)
   MatrixGenerics         1.14.0     2023-10-24 [1] Bioconductor
   matrixStats            1.5.0      2025-01-07 [1] RSPM (R 4.3.3)
   memoise                2.0.1      2021-11-26 [1] RSPM (R 4.3.3)
 P mgcv                   1.9-1      2023-12-21 [?] CRAN (R 4.3.3)
   mime                   0.13       2025-03-17 [1] RSPM (R 4.3.3)
 P miniUI                 0.1.1.1    2018-05-18 [?] CRAN (R 4.3.2)
   multtest               2.58.0     2023-10-24 [1] Bioconductor
   munsell                0.5.1      2024-04-01 [1] RSPM (R 4.3.3)
 P nlme                   3.1-164    2023-11-27 [?] CRAN (R 4.3.3)
 P pacman                 0.5.1      2019-03-11 [?] CRAN (R 4.3.2)
   permute                0.9-7      2022-01-27 [1] RSPM (R 4.3.3)
   phyloseq             * 1.46.0     2023-10-24 [1] Bioconductor
   pillar                 1.10.2     2025-04-05 [1] RSPM (R 4.3.3)
 P pkgbuild               1.4.7      2025-03-24 [?] RSPM (R 4.3.3)
   pkgconfig              2.0.3      2019-09-22 [1] RSPM (R 4.3.3)
 P pkgload                1.4.0      2024-06-28 [?] CRAN (R 4.3.2)
   plyr                   1.8.9      2023-10-02 [1] RSPM (R 4.3.3)
   png                    0.1-8      2022-11-29 [1] RSPM (R 4.3.3)
 P profvis                0.4.0      2024-09-20 [?] CRAN (R 4.3.2)
 P promises               1.3.2      2024-11-28 [?] CRAN (R 4.3.2)
   purrr                * 1.0.4      2025-02-05 [1] RSPM (R 4.3.3)
   R6                     2.6.1      2025-02-15 [1] RSPM (R 4.3.3)
   RColorBrewer           1.1-3      2022-04-03 [1] RSPM (R 4.3.3)
   Rcpp                 * 1.0.14     2025-01-12 [1] RSPM (R 4.3.3)
   RcppParallel           5.1.10     2025-01-24 [1] RSPM (R 4.3.3)
   RCurl                  1.98-1.17  2025-03-22 [1] RSPM (R 4.3.3)
   readr                * 2.1.5      2024-01-10 [1] RSPM (R 4.3.3)
 P remotes                2.5.0      2024-03-17 [?] CRAN (R 4.3.2)
   renv                   1.0.3      2023-09-19 [1] CRAN (R 4.3.3)
   reshape2               1.4.4      2020-04-09 [1] RSPM (R 4.3.3)
   rhdf5                  2.46.1     2023-11-29 [1] Bioconductor 3.18 (R 4.3.3)
   rhdf5filters           1.14.1     2023-11-06 [1] Bioconductor
   Rhdf5lib               1.24.2     2024-02-07 [1] Bioconductor 3.18 (R 4.3.3)
   rlang                  1.1.5      2025-01-17 [1] RSPM (R 4.3.3)
   rmarkdown              2.29       2024-11-04 [1] RSPM (R 4.3.3)
   Rsamtools              2.18.0     2023-10-24 [1] Bioconductor
   rstudioapi             0.17.1     2024-10-22 [1] RSPM (R 4.3.3)
   S4Arrays               1.2.1      2024-03-04 [1] Bioconductor 3.18 (R 4.3.3)
   S4Vectors              0.40.2     2023-11-23 [1] Bioconductor 3.18 (R 4.3.3)
   scales                 1.3.0      2023-11-28 [1] RSPM (R 4.3.3)
 P sessioninfo            1.2.3      2025-02-05 [?] RSPM (R 4.3.3)
 P shiny                  1.10.0     2024-12-14 [?] CRAN (R 4.3.2)
   ShortRead              1.60.0     2023-10-24 [1] Bioconductor
   SparseArray            1.2.4      2024-02-11 [1] Bioconductor 3.18 (R 4.3.3)
   stringi                1.8.7      2025-03-27 [1] RSPM (R 4.3.3)
   stringr              * 1.5.1      2023-11-14 [1] RSPM (R 4.3.3)
   SummarizedExperiment   1.32.0     2023-10-24 [1] Bioconductor
 P survival               3.5-8      2024-02-14 [?] CRAN (R 4.3.3)
   tibble               * 3.2.1      2023-03-20 [1] RSPM (R 4.3.3)
   tidyr                * 1.3.1      2024-01-24 [1] RSPM (R 4.3.3)
   tidyselect             1.2.1      2024-03-11 [1] RSPM (R 4.3.3)
   tidyverse            * 2.0.0      2023-02-22 [1] RSPM (R 4.3.3)
   timechange             0.3.0      2024-01-18 [1] RSPM (R 4.3.3)
   tzdb                   0.5.0      2025-03-15 [1] RSPM (R 4.3.3)
 P urlchecker             1.0.1      2021-11-30 [?] CRAN (R 4.3.2)
 P usethis                3.1.0      2024-11-26 [?] CRAN (R 4.3.2)
   vctrs                  0.6.5      2023-12-01 [1] RSPM (R 4.3.3)
   vegan                  2.6-10     2025-01-29 [1] RSPM (R 4.3.3)
   withr                  3.0.2      2024-10-28 [1] RSPM (R 4.3.3)
   xfun                   0.52       2025-04-02 [1] RSPM (R 4.3.3)
 P xtable                 1.8-4      2019-04-21 [?] CRAN (R 4.3.2)
   XVector                0.42.0     2023-10-24 [1] Bioconductor
   yaml                   2.3.10     2024-07-26 [1] RSPM (R 4.3.3)
   zlibbioc               1.48.2     2024-03-13 [1] Bioconductor 3.18 (R 4.3.3)

 [1] /local/workdir/arp277/AAM_ERI_2024/renv/library/R-4.3/x86_64-pc-linux-gnu
 [2] /home/arp277/.cache/R/renv/sandbox/R-4.3/x86_64-pc-linux-gnu/fd835031

 * ── Packages attached to the search path.
 P ── Loaded and on-disk path mismatch.

──────────────────────────────────────────────────────────────────────────────
```
:::
:::