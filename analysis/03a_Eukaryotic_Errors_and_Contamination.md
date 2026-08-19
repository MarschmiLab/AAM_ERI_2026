---
title: "Assess Errors and Contamination" 
author: "Augustus Pendleton"
date: "16 April, 2025"
knitr: 
  opts_chunk: 
    fig.path: ../figures/03a_Eukaryotic_Errors_and_Contamination/
editor_options: 
  chunk_output_type: console
---




# Goals of our Pre-Processing File

We will work here in this document to remove any unncessary ASVs or samples, such as: 

1. Mitochondria or plastids (unlikely)
2. Contamination within our blanks
3. Domain Bacteria
4. Splitting metazoans away

## Load Packages and Functions


::: {.cell layout-align="center"}

```{.r .cell-code}
# Load Packages
pacman::p_load(phyloseq, Biostrings, dada2, microViz, tidyverse, install = FALSE)

source("code/R/plotting_aesthetics.R")
```
:::


## Load the data!


::: {.cell layout-align="center"}

```{.r .cell-code}
# Load in the raw_phyloseq data from DADA2_workflow.Rmd
load("data/02a_Eukaryotic_Taxonomic_Classification/raw_euk_physeq.RData")

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
load("data/01a_Eukaryotic_Denoising/track_df_euk.RData")
```
:::


# Add metadata to physeq


::: {.cell layout-align="center"}

```{.r .cell-code}
sample_data(raw_euk_physeq) <- track_df_euk
```
:::



# Looking at bacteria


::: {.cell layout-align="center"}

```{.r .cell-code}
tax_df <- 
raw_euk_physeq %>%
  tax_table() %>%
  data.frame()

tax_df %>%
  count(Domain)
```

::: {.cell-output .cell-output-stdout}
```
          Domain    n
1       Bacteria    1
2      Eukaryota 1607
3 Eukaryota:nucl   12
4           <NA>    2
```
:::

```{.r .cell-code}
tax_df %>%
  filter(Domain == "Eukaryota:nucl")
```

::: {.cell-output .cell-output-stdout}
```
                 Domain     Supergroup         Division        Subdivision
EUK_32   Eukaryota:nucl Cryptista:nucl Cryptophyta:nucl Cryptophyta_X:nucl
EUK_188  Eukaryota:nucl Cryptista:nucl Cryptophyta:nucl Cryptophyta_X:nucl
EUK_316  Eukaryota:nucl Cryptista:nucl Cryptophyta:nucl Cryptophyta_X:nucl
EUK_488  Eukaryota:nucl Cryptista:nucl Cryptophyta:nucl Cryptophyta_X:nucl
EUK_615  Eukaryota:nucl Cryptista:nucl Cryptophyta:nucl Cryptophyta_X:nucl
EUK_709  Eukaryota:nucl Cryptista:nucl Cryptophyta:nucl Cryptophyta_X:nucl
EUK_834  Eukaryota:nucl Cryptista:nucl Cryptophyta:nucl Cryptophyta_X:nucl
EUK_901  Eukaryota:nucl Cryptista:nucl Cryptophyta:nucl Cryptophyta_X:nucl
EUK_1238 Eukaryota:nucl Cryptista:nucl Cryptophyta:nucl Cryptophyta_X:nucl
EUK_1331 Eukaryota:nucl Cryptista:nucl Cryptophyta:nucl Cryptophyta_X:nucl
EUK_1343 Eukaryota:nucl Cryptista:nucl Cryptophyta:nucl Cryptophyta_X:nucl
EUK_1531 Eukaryota:nucl Cryptista:nucl Cryptophyta:nucl Cryptophyta_X:nucl
                      Class                Order                 Family
EUK_32   Cryptophyceae:nucl Cryptomonadales:nucl Cryptomonadales_X:nucl
EUK_188  Cryptophyceae:nucl Cryptomonadales:nucl Cryptomonadales_X:nucl
EUK_316  Cryptophyceae:nucl Cryptomonadales:nucl Cryptomonadales_X:nucl
EUK_488  Cryptophyceae:nucl Pyrenomonadales:nucl  Pyrenomonadaceae:nucl
EUK_615  Cryptophyceae:nucl Cryptomonadales:nucl Cryptomonadales_X:nucl
EUK_709  Cryptophyceae:nucl Cryptomonadales:nucl Cryptomonadales_X:nucl
EUK_834  Cryptophyceae:nucl Cryptomonadales:nucl Cryptomonadales_X:nucl
EUK_901  Cryptophyceae:nucl Cryptomonadales:nucl   Hemiselmidaceae:nucl
EUK_1238 Cryptophyceae:nucl Cryptomonadales:nucl Cryptomonadales_X:nucl
EUK_1331 Cryptophyceae:nucl Cryptomonadales:nucl Cryptomonadales_X:nucl
EUK_1343 Cryptophyceae:nucl Cryptomonadales:nucl Cryptomonadales_X:nucl
EUK_1531 Cryptophyceae:nucl Cryptomonadales:nucl Cryptomonadales_X:nucl
                     Genus                           Species      ASV
EUK_32   Plagioselmis:nucl       Plagioselmis_lacustris:nucl   EUK_32
EUK_188  Plagioselmis:nucl       Plagioselmis_lacustris:nucl  EUK_188
EUK_316  Plagioselmis:nucl       Plagioselmis_lacustris:nucl  EUK_316
EUK_488    Rhodomonas:nucl               Rhodomonas_sp.:nucl  EUK_488
EUK_615   Cryptomonas:nucl    Cryptomonas_pyrenoidifera:nucl  EUK_615
EUK_709  Plagioselmis:nucl       Plagioselmis_lacustris:nucl  EUK_709
EUK_834   Cryptomonas:nucl                              <NA>  EUK_834
EUK_901    Chroomonas:nucl          Chroomonas_coerulea:nucl  EUK_901
EUK_1238  Cryptomonas:nucl          Cryptomonas_curvata:nucl EUK_1238
EUK_1331 Plagioselmis:nucl       Plagioselmis_lacustris:nucl EUK_1331
EUK_1343  Cryptomonas:nucl Cryptomonas_tetrapyrenoidosa:nucl EUK_1343
EUK_1531  Cryptomonas:nucl    Cryptomonas_pyrenoidifera:nucl EUK_1531
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   ASVseqs
EUK_32   TTCCAGCTCCAATAGCGTATACTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTCGAATGAAGGTTTCTTGAACTGCGCGAACTTCTAGGTGTAACCTAGACTAGTTTGTGCACAAGAACCTCCCTTTTGTAGTTGATTAAATCCAATTGGATCTAATTGTCTACATGCATTACTGTGAACAAAATAGAGTGTTCATGGCAGGTCAATTTAGGCTTGAATACCTCAGCATGGAATAATAGAATAGGACTTCTGGTTCTATTTTGTTGGTTTTTGGAGCCGAAGCAACTGCGAAAGCATCTGCCAAGGATGTTTTCATTGATCAAGAACGAAAGTTAGGGGATCGAAGACGATTAGATACCGTCGTAGTCCTAACCATAAACTATGCCGACTGGGGATTGGTAGATGTACCAAAAACGACTCTATCAGCACCCTATGAGAAATCAAAGTCTTTGGGTTCTGGGGGGAGTATGGTCGCAAGGCTG
EUK_188  TTCCAGCTCCAATAGCGTATACTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTCGAATGAAGGTTTCTTGAATTGCGCGAACTTCTAAGTGTAACTTAGACTAGTTTGCGCACAAGAACCTCCCTTATGTGTCTTAATAGATCCAATTTATTTTTGGATCTAATCAGTCACATGCATTACTGTGAACAAAATAGAGTGTTCATGGCAGGTCAATTTAGGCTTGAATACCTCAGCATGGAATAATAGAATAGGACTTCCGGTTCTATTTTGTTGGTTTTTGGACAACTGCGAAAGCATCTGCCAAGGATGTTTTCATTGATCAAGAACGAAAGTTAGGGGATCGAAGACGATTAGATACCGTCGTAGTCCTAACCATAAACTATGCCGACTGGGGATTGGTAGATGTACCAAAAACGACTCTATCAGCACCCTATGAGAAATCAAAGTCTTTGGGTTCTGGGGGGAGTATGGTCGCAAGGCTG
EUK_316  TTCCAGCTCCAATAGCGTATACTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTCGAATGAAGGTTTTTTTGTTATGCGGACTCTTAGGTGTTACCTAGGTGGTCTGCGTTCAATAACCTCCCTGGTGTGATACAAATTTTTTCGACTGAAAAAATTTGTTCACATGCATTACTGTGAACAAAATAGAGTGTTCATGGCAGGCCAAATTAGGCTTGAATACCTCAGCATGGAATAATAGAATAGGACTTCCGGTTCTATTTTGTTGGTTTTTGGAGCCGAAGTCAACTGCGAAAGCATCTGCCAAGGATGTTTTCATTGATCAAGAACGAAAGTTAGGGGATCGAAGACGATTAGATACCGTCGTAGTCCTAACCATAAACTATGCCGACTGGGGATTGGTAGATGTACCAAAAACGACTCTATCAGCACCCTATGAGAAATCAAAGTCTTTGGGTTCTGGGGGGAGTATGGTCGCAAGGCTG
EUK_488  TTCCAGCTCCAATAGCGTATACTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTCGAATTCGGACCTGTGGGTGCGGATTGATCTGTAAAGGTCCTTTGGCACACCAGGTCTAATTTTGTCGTTGCAAAGCCCCGCCTCCTTGGGCGTTGTGGTCTGCAAAGCGGCATCTTTTACTGTGAATAAAATAGAGTGTTCATGGCAGGCCTTTGAAAGGCTTAAATACCTCAGCATGGAATAATAGAATAGGACTTCCGGTTCTATTTTGTTGGTTTTTGGAGCCGACAACTGCGAAAGCATCTGCCAAGGATGTTTTCATTGATCAAGAACGAAAGTTAAGGGATCGAAGACGATTAGATACCGTCGTAGTCTTAACCATAAACTATGCCGACTGGGGATTGGTAGATGTACCAAAAACGACTCTATCAGCACCCTATGAGAAATCAAAGTCTATGGGTTCTGGGGGGAGTATGGTCGCAAGGCTG
EUK_615  TTCCAGCTCCAATAGCGTATACTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTCGAATTTTAGTTTTAGATTGTTGGCTCGTTATGAATAATAGCGTAGGCTAGAACAATATAAAACTTTGTTTGCTAATTTATAAAACGTCTCTTTAATGAGTCTTTGTAAATTATGCATTTTTTACTGTGAACAAAATAGAGTGTTCATGGCAGGCCTTAAAGGCTTATATACCTCAGCATGGAATAATAGAGTAGGACTTCTGGTTCTATTTTATTGGTTATTGGAGCCCAACTGCGAAAGCATCTGCCAAGGATGTTTTCATTGATCAAGAACGAAAGTTAGGGGATCGAAGACGATTAGATACCGTCGTAGTCTTAACCATAAACTATGCCGACTGGGGATTGGTAGATGTATCAAATACGACTCTATCAGTATCCTATGAGAAATCAAAGTCTCTGGGTTCTGGGGGGAGTATGGTCGCAAGGCTG
EUK_709  TTCCAGCTCCAATAGCGTATACTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTCGAATGAAGGTTTTTTTTGTTATGCGGACTCTTGGGTGTTACCTAGGTAGTCTGCGTTCAATAACCTCCCTGGTGTGATACAAATTTCTTCAACCGAAAAAATTTGTTCACATGCATTACTGTGAACAAAATAGAGTGTTCATGGCAGGCCAAATTAGGCTTGAATACCTCAGCATGGAATAATAGAATAGGACTTCCGGTTCTATTTTGTTGGTTTTTGGAGCCGAAGCAACTGCGAAAGCATCTGCCAAGGATGTTTTCATTGATCAAGAACGAAAGTTAGGGGATCGAAGACGATTAGATACCGTCGTAGTCCTAACCATAAACTATGCCGACTGGGGATTGGTAGATGTACCAAAAACGACTCTATCAGCACCCTATGAGAAATCAAAGTCTTTGGGTTCTGGGGGGAGTATGGTCGCAAGGCTG
EUK_834  TTCCAGCTCCAATAGCGTATACTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTCGAATATTAGTTTCCAATCGCTGCCTTGTTCTTTATTGTAACAGAGGCTTGGTGTGTTGGAATCTTTTTTTTGTCCAATGGAAGACATTCATTTTGAATGAAACCCATTGTGCATATATTACTGTGAATAAAATAGAGTGTTCATGGCAGGCCTTTATGGCTTATATACCTTAGCATGGAATAATAGAATAGGACTGTCGGTTCTATTTTATTGGTTATTGGAGCCGAACAACTGCGAAAGCATCTGCCAAGGATGTTTTCATTGATCAAGAACGAAAGTTAGGGGATCGAAGACGATTAGATACCGTCGTAGTCTTAACCATAAACTATGCCGACTGGGGATTGGTAGATGTATCAAATACAACTCTATCAGTATCCTATGAGAAATCAAAGTCTCTGGGTTCTGGGGGGAGTATGGTCGCAAGGCTG
EUK_901  TTCCAGCTCCAATAGTGTATATTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTCGAATTTTAAGATATCAGCGCCGATCTGTCTTTTAGGAGCCAGGCATTGGCGCAGATATCTTTCCTTGCCGACCACTGGCTAGTTTTTTATTAAATTAGCTTGGGGTATCGGTAATTTTTTACTGTGAATAAATTAGAGTGTTCATGGCAGGCCGAAAGGGCTTAGATATCTCAGCATGGAATAATAGAATAGGACTTCAGGTTCTATTTTGTTGGTTTTTTGGAGCCGCAACTGCGAAAGCATCTGCCAAGGATGTTTTCATTGATCAAGAACGAAAGTTAGGGGATCGAAGACGATTAGATACCGTCGTAGTCTTAACCATAAACTATGCCGACTGGGGATTGGTAGATGTATCAACGACAACTCTATCAGCACTCTATGAGAAATCAAAGTCTCTGGGTTCTGGGGGGAGTATGGTCGCAAGGCTG
EUK_1238 TTCCAGCTCCAATAGCGTATACTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTCGAATTTCAGTTTTAAACTGCTATCTTGTTATAATTTATAATAAAGGTTAGGTGAGTTTAGAACTTATTTTGCTAATTTGAAAATTACCTCTTAATTGAGTATTTTAAAATTATGCATTTTTTACTGTGAACAAAATAGAGTGTTCATGGCAGGCCTTTAAGGCTTATATACCTCAGCATGGAATAATAGAGTAGGACTTCCGGTTCTATTTTATTGGTTATTGGAGCCCAACTGCGAAAGCATCTGCCAAGGATGTTTTCATTGATCAAGAACGAAAGTTAGGGGATCGAAGACGATTAGATACCGTCGTAGTCTTAACCATAAACTATGCCGACTGGGGATTGGTAGATGTATCAAATACGACTCTATCAGTATCCTATGAGAAATCAAAGTCTCTGGGTTCTGGGGGGAGTATGGTCGCAAGGCTG
EUK_1331 AACGGAGGGGGCACTCGTACTCCGTTGTCAGAGGTGAAATTCTAGGATTGACGGAAGACGAACAACTGCGAAAGCATCTGCCAAGGATGTTTTCATTGATCAAGAACGAAAGTTAGGGGATCGAAGACGATTAGATACCGTCGTAGTCCTAACCATAAACTATGCCGACTGGGGATTGGTAGATGTACCAAAAACGACTCTATCAGCACCCTATGAGAAATCAAAGTCTTTGGGTTCTGGGGGGAGTATGGTCGCAAGGCTGAAACTTAAATAAATTGACGCAACTGCGAAAGCATCTGCCAAGGATGTTTTCATTGATCAAGAACGAAAGTTAGGGGATCGAAGACGATTAGATACCGTCGTAGTCCTAACCATAAACTATGCCGACTGGGGATTGGTAGATGTACCAAAAACGACTCTATCAGCACCCTATGAGAAATCAAAGTCTTTGGGTTCTGGGGGGAGTATGGTCGCAAGGCTG
EUK_1343 TTCCAGCTCCAATAGCGTATACTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTCGAATTTCAGTTCTAGATTGTTTCCTTGTTATGTATTATAATGTAGGATAAAACAGTTTAGAACTTTTTTTGTAAATTTATAAACATATTCTTCATTGAATTTTTATAAATATTGCAACTTTTACTGTGAACAAAATAGAGTGTTCATGGCAGGCCGTTATGGCTTATATACCTCAGCATGGAATAATAGAATAGGATTTCTGGTTCTATTTTATTGGTTATTGGAGCCCAACTGCGAAAGCATCTGCCAAGGATGTTTTCATTGATCAAGAACGAAAGTTAGGGGATCGAAGACGATTAGATACCGTCGTAGTCTTAACCATAAACTATGCCGACTGGGGATTGGTAGATGAAACAAATACGACTCTATCAGTATCCTATGAGAAATCAAAGTCTCTGGGTTCTGGGGGGAGTATGGTCGCAAGGCTG
EUK_1531 TTCCAGCTCCAATAGCGTATACTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTCGAATTTTAGTTTAAAACTGTTGCCTTGTTATGAATAGTAACATAGGATAAAACAATTTTGAACTTTGGTTGCTAATTTACGAAATGTTTCATCACTGAGCTTTTGTAAATTATGCATTTTTTACTGTGAACAAAATAGAGTGTTCATGGCAGGCCTTTAAGGCTTATATACCTCAGCATGGAATAATAGAGTAGGACTTCTGGTTCTATTTTATTGGTTATTGGAGCCCAACTGCGAAAGCATCTGCCAAGGATGTTTTCATTGATCAAGAACGAAAGTTAGGGGATCGAAGACGATTAGATACCGTCGTAGTCTTAACCATAAACTATGCCGACTGGGGATTGGTAGATGTATCAAATACGACTCTATCAGTATCCTATGAGAAATCAAAGTCTCTGGGTTCTGGGGGGAGTATGGTCGCAAGGCTG
```
:::

```{.r .cell-code}
tax_df %>%
  filter(Species == "Plagioselmis_lacustris")
```

::: {.cell-output .cell-output-stdout}
```
 [1] Domain      Supergroup  Division    Subdivision Class       Order      
 [7] Family      Genus       Species     ASV         ASVseqs    
<0 rows> (or 0-length row.names)
```
:::

```{.r .cell-code}
tax_df %>%
  filter(Species == "Cryptomonas_curvata")
```

::: {.cell-output .cell-output-stdout}
```
            Domain Supergroup    Division   Subdivision         Class
EUK_4    Eukaryota  Cryptista Cryptophyta Cryptophyta_X Cryptophyceae
EUK_1399 Eukaryota  Cryptista Cryptophyta Cryptophyta_X Cryptophyceae
EUK_1571 Eukaryota  Cryptista Cryptophyta Cryptophyta_X Cryptophyceae
                   Order            Family       Genus             Species
EUK_4    Cryptomonadales Cryptomonadales_X Cryptomonas Cryptomonas_curvata
EUK_1399 Cryptomonadales Cryptomonadales_X Cryptomonas Cryptomonas_curvata
EUK_1571 Cryptomonadales Cryptomonadales_X Cryptomonas Cryptomonas_curvata
              ASV
EUK_4       EUK_4
EUK_1399 EUK_1399
EUK_1571 EUK_1571
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   ASVseqs
EUK_4    TTCCAGCTCTAATAGCATATATTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTCGGATTTTGGGCGCGTGTGAACTGTCAGCCTCTGGTTGGACGGTGCGCGCGTGCCTTTCTGCTTGGGAGCCCCATGGTGGTTAATTCTGTCGTGGTGGGTCACAGGCCGTTTACTTTGAAAAAATTAGAGTGTTCAAAGCAGGCTGTTGCTTGAATAGATTAGCATGGAATAATGGAATAGGACCTTGGTGCTATTTCGTTGGTTATGGGACTGGGGTAATGATTAATAACTTGTGCGAAAGCATTCGGCAAGGATGTTTTCATTGATCAAGAACGAAAGTTAGGGGATCGAAGACGATCAGATACCGTCGTAGTCTTAACCATAAACTATGCCAACTAGGGATCATTGAGTGTCAATTTACGACCTCATTGGCACCTTGTGAGAAATCAAAGTTTTTGGGTTCCGGGGGGAGTATGGTCGCAAGGCTG
EUK_1399 TTCCAGCTCTAATAGCATATATTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTCGGATTTTGGGCGCGTGTGAACTGTCAGCCTCTGGTTGGACGGTGCGCGCGTGCCTTTCTGCTTGGGAGCCCCATGGTGGTTAATTCTGTCGTGGTGGGTCACAGGCCGTTTACTTTGAAAAACGATCAGATACCGTCGTAGTCTTAACCATAAACTATGCCAACTAGGGATCATTGAGTGTCAATTTACGACCTCATTGGCACCTTGTGAGAAATCAAAGTTTTTGGGTCTGCTTGGGAGCCCCATGGTGGTTAATTCTGTCGTGGTGGGTCACAGGCCGTTTACTTTGAAAAACGATCAGATACCGTCGTAGTCTTAACCATAAACTATGCCAACTAGGGATCATTGAGTGTCAATTTACGACCTCATTGGCACCTTGTGAGAAATCAAAGTTTTTGGGTTCCGGGGGGAGTATGGTCGCAAGGCTG
EUK_1571 TTCCAGCTCTAATAGCATATATTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTCGGATTTTGGGCGCGTGTGAACTGTCAGCCTCTGGTTGGACGGTGCACGCGTGCCGTTCATATTTCGTTGTCAGAGGTGAAATTCTTGGATTTACGAACGATGAACTTGTGCGAAAGCATTCGGCAAGGATGTTTTCATTGATCAAGAACGAAAGTTAGGGGATCGAAGACGATCAGATACCGTCGTAGTCTTAACCATAAACTATGCCAACTAGGGATCATTGAGTGTACTTGTGCGAAAGCATTCGGCAAGGATGTTTTCATTGATCAAGAACGAAAGTTAGGGGATCGAAGACGATCAGATACCGTCGTAGTCTTAACCATAAACTATGCCAACTAGGGATCATTGAGTGTCAATTTACGACCTCATTGGCACCTTGTGAGAAATCAAAGTTTTTGGGTTCCGGGGGGAGTATGGTCGCAAGGCTG
```
:::
:::


We'll definitley remove the one bacteria, and probably the two which have an NA domain, but I'm unsure what to do about the nucleomorphs. I'm assuming they'll make the tree look quite crazy (potentially interfering with downstream analyses). Are we also "counting" these species 2X then? Doesn't seem so for the Plagioformis, but it does for the Cryptomonas species. 

# Remove unnecessary taxa

We will need to remove the ASVs that are mitochondria or chloroplast

## Remove Control samples 

Next, we will need to remove ASVs that seem to be contamination based on our controls. But guess what? None of our control samples survived with enough reads from either cutadapt or filterandTrim, so those samples are null. 


::: {.cell layout-align="center"}

```{.r .cell-code}
control_samples <- c("KH_D024", "KH_D048", "KH_D049", "KH_D050", "KH_D084", "KH_D085", "KH_D085")

sample_names(raw_euk_physeq) %in% control_samples
```

::: {.cell-output .cell-output-stdout}
```
 [1] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
[13] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
[25] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
[37] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
[49] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
[61] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
[73] FALSE FALSE FALSE FALSE
```
:::
:::


# Evaluating Accuracy

I compare which ASVs I observe in our mock sequencing compared to which ASVs should be present


::: {.cell layout-align="center"}

```{.r .cell-code}
# Check the mock commmunity 
mock_physeq <- raw_euk_physeq %>%
  subset_samples(.,(DNA_ID=="ZYMO_MOCK")) %>%
  prune_taxa(taxa_sums(.) > 0, .)

tax_table(mock_physeq)[,1:9] # Look at what taxa are there
```

::: {.cell-output .cell-output-stdout}
```
Taxonomy Table:     [2 taxa by 9 taxonomic ranks]:
       Domain      Supergroup Division       Subdivision Class          
EUK_20 "Eukaryota" "Obazoa"   "Opisthokonta" "Fungi"     "Ascomycota"   
EUK_49 "Eukaryota" "Obazoa"   "Opisthokonta" "Fungi"     "Basidiomycota"
       Order              Family              Genus         
EUK_20 "Saccharomycotina" "Saccharomycetales" "Candida"     
EUK_49 "Agaricomycotina"  "Tremellomycetes"   "Cryptococcus"
       Species              
EUK_20 "Candida_wounanorum" 
EUK_49 "Cryptococcus_gattii"
```
:::

```{.r .cell-code}
otu_matrix <- otu_table(mock_physeq) %>% as.matrix() # Will use this matrix later

cat("DADA2 inferred", ntaxa(mock_physeq), "ASVs present in the Mock Community.")
```

::: {.cell-output .cell-output-stdout}
```
DADA2 inferred 2 ASVs present in the Mock Community.
```
:::

```{.r .cell-code}
#Who are they in the mock community? 

mock_seqs <- tax_table(mock_physeq) %>%
  as.data.frame() %>%
  select(ASV, ASVseqs)
```
:::


We're looking for Saccharomyces cerevisiae and Cryptococcus neoformans, which our classification at least doesn't quite match. However when I blast both those sequences, I get those pieces as 100% hits, so I feel good about our error rate (zero, supposedly). Let's remove the mock, non-Eukaryota domains, and save our new physeq.


::: {.cell layout-align="center"}

```{.r .cell-code}
preprocessed_physeq <- raw_euk_physeq %>%
  subset_taxa(Domain == "Eukaryota") %>%
  subset_samples(.,(DNA_ID!="ZYMO_MOCK")) %>%
  prune_taxa(taxa_sums(.) > 0, .) 

preprocessed_meta_physeq <- 
  preprocessed_physeq %>%
  subset_taxa(Subdivision == "Metazoa") %>%
  prune_samples(sample_sums(.)>0, .)

preprocessed_meta_physeq
```

::: {.cell-output .cell-output-stdout}
```
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 115 taxa and 70 samples ]
sample_data() Sample Data:       [ 70 samples by 6 sample variables ]
tax_table()   Taxonomy Table:    [ 115 taxa by 11 taxonomic ranks ]
```
:::

```{.r .cell-code}
preprocessed_euk_physeq <- 
  preprocessed_physeq %>%
  subset_taxa(Subdivision != "Metazoa") %>%
  prune_samples(sample_sums(.)>0, .)

preprocessed_euk_physeq
```

::: {.cell-output .cell-output-stdout}
```
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 1455 taxa and 73 samples ]
sample_data() Sample Data:       [ 73 samples by 6 sample variables ]
tax_table()   Taxonomy Table:    [ 1455 taxa by 11 taxonomic ranks ]
```
:::
:::



# Save preprocessed phyloseq!

::: {.cell layout-align="center"}

```{.r .cell-code}
save(preprocessed_euk_physeq, file = "data/03a_Eukaryotic_Errors_and_Contamination/preprocessed_euk_physeq.RData")

save(preprocessed_meta_physeq, file = "data/03a_Eukaryotic_Errors_and_Contamination/preprocessed_meta_physeq.RData")
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
   BiocGenerics         * 0.48.1     2023-11-01 [1] Bioconductor
   BiocParallel           1.36.0     2023-10-24 [1] Bioconductor
   biomformat             1.30.0     2023-10-24 [1] Bioconductor
   Biostrings           * 2.70.3     2024-03-13 [1] Bioconductor 3.18 (R 4.3.3)
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
   GenomeInfoDb         * 1.38.8     2024-03-15 [1] Bioconductor 3.18 (R 4.3.3)
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
   IRanges              * 2.36.0     2023-10-24 [1] Bioconductor
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
   microViz             * 0.12.7     2025-04-16 [1] https://david-barnett.r-universe.dev (R 4.3.3)
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
   S4Vectors            * 0.40.2     2023-11-23 [1] Bioconductor 3.18 (R 4.3.3)
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
   XVector              * 0.42.0     2023-10-24 [1] Bioconductor
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
