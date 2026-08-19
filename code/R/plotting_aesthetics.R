# This file indicates that colors, shapes, and themes across my documents
# let's keep the visualizations consistent! :) 
library(ggplot2)
triglav <- NatParksPalettes::natparks.pals("Triglav",n=6)
c(triglav)

NatParksPalettes::natparks.pals("Charmonix",n=6) %>%
  c()

erie_theme <- theme_classic() + 
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size=15),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 15),
        axis.line = element_line(linewidth = .6),
        strip.background = element_rect(colour = NA, fill = 'transparent'),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        panel.grid.major = element_blank(), #remove major gridlines
        panel.grid.minor = element_blank(), #remove minor gridlines
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent', color = NA) #transparent legend panel
  )

theme_set(erie_theme)

ggpreview <- function(...) {
  fname <- tempfile(fileext = ".png")
  ggsave(filename = fname, ...)
  system2("open", fname)
  invisible(NULL)
}

basin_colors <- 
  c(E = "#B9741F",
    C = "#990006",
    W = "#625D0A")

month_colors <- 
  c(May = "#F2C621",
    August = "#386EC2",
    September = "#E23B0E")

group_colors <- 
  c("1" = "#F2C621",
    "2" = "#B6AA0D",
    "3" = "#008FF8",
    "4" = "#E23B0E",
    "5" = "#196689")

# Set the phylum colors
phylum_colors <- c(
  p__Actinomycetota = "red1", 
  p__Pseudomonadota = "royalblue",
  p__Bacteroidota = "gold", 
  p__Cyanobacteriota = "limegreen",
  p__Verrucomicrobiota = "darkorchid1",
  p__Planctomycetota = "darkorange", 
  p__Chloroflexota="black", 
  p__Acidobacteriota = "navy", 
  p__Bdellovibrionota_473306 = "darkslategray2",
  p__Myxococcota_A_473307 = "#B5D6AA",
  p__Gemmatimonadota = "greenyellow",
  p__Armatimonadota = "deeppink1",
  p__Desulfobacterota_G_459543="magenta",
  c__Alphaproteobacteria = "plum2", 
  p__Gammaproteobacteria = "greenyellow",
  Other = "grey",
  Rare = "grey")

phylum_rjects <- c(c__Betaproteobacteria = "plum1", 
                   
                   
                   p__Crenarchaeota = "yellow",
                   
                   p__Deltaproteobacteria = "grey", 
                   
                   p__Firmicutes = "#3E9B96",
                   
                   p__Patescibacteria = "yellow",
                   
                   p__Nitrospirota = "palevioletred1",
                   
                
                   
                   p__Thermoplasmatota = "green")
# Taken from:
# https://stackoverflow.com/questions/9563711/r-color-palettes-for-many-data-classes
c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)


genus_colors <- c(
  "Rare" = "grey90",
  "OM27 clade" = "gold1",
  "betI-A" = "green4",
  "betI-B" = "green1",
  "betIV-A" = "palegreen2",
  "Candidatus Nitrotoga" = "#FB9A99",
  "acIV-C" = "blue1",
  "acIV-A" = "steelblue4",
  "acI-A" = "dodgerblue2",
  "acI-B" = "skyblue2",
  "hgcI clade" = "darkturquoise",
  "acI-C" = "turquoise1",
  "Cyanobium PCC-6307" = "olivedrab1",
  "bacIII-A" = "#FF7F00",
  "bacIII-B" = "#FDBF6F",
  "Fluviicola" = "khaki2",
  "bacI-A" = "yellow4",
  "bacII-A" = "yellow3",
  "Lacunisphaera" = "#CAB2D6",
  "LD19" = "#6A3D9A",
  "DEV114" = "orchid1",
  "SH3-11" = "deeppink1",
  "Chthoniobacter" = "maroon",
  "CL500-3" = "tan4",
  "Nitrospira" = "black",
  "alfV-A" = "#E31A1C"
)



class_colors <- c("(Rare)" = "grey80",
                  "Actinomycetes" = "blue1",
                  "Alphaproteobacteria" = "red",
                  "Cyanobacteriia" = "green4",
                  "Acidimicrobiia_401430" = "steelblue2",
                  "Kapaibacteriia" = "#FBBF6F",
                  "Gammaproteobacteria" = "deeppink1",
                  "Verrucomicrobiae" = "#6A3D9A",
                  "Anaerolineae" = "black",
                  "Oligoflexia" = "tan4",
                  "Bacteroidia" = "#FF7F00",
                  "Holophagae" = "olivedrab1",
                  "Planctomycetia" = "maroon",
                  "Limnocylindria" = "grey30",
                  "Gemmatimonadetes" = "gold1"
)

cyano_family_colors <- c("(Rare/Unknown)" = "grey80",
                         "Pseudanabaenaceae" = "blue1",
                         "Cyanobiaceae" = "red",
                         "Microcystaceae" = "green4",
                         "Nostocaceae" = "steelblue2",
                         "Microcoleaceae" = "#FF7F00",
                         "Phormidesmiceae" = "#643D9A",
                         "Microcystaceae_B" = "olivedrab1",
                         "Gloeocapsaceae" = "deeppink1",
                         "Synechococcaeae" = "maroon",
                         "Cyanobacteriaeae" = "gold1",
                         "Microcystaceae_A_24133" = "olivedrab3",
                         "Vampirovibrionaceae" = "black")
              
cyano_genus_colors <- 
  c(
    #Cyanobiaceae
    "Cyanobium_A_31525" = "#065563",
    "GFB01" = "#088096",
    "Synechococcus_D" = "#0eabc7",
    "Vulcanococcus" = "#11d3f5",
    "Unassigned Cyanobiaceae" = "#aaeefa",
    #Synechococcaceae
    "Synechococcus" = "#0e5bc7",
    #"Microcystaceae_A"
    "Microcystis" = "green4",
    "Gloeothece" = "green2",
    #Microcystaceae_B"
    "Snowella" = "olivedrab4",
    "Synechocystis" = "olivedrab3",
    #Pseudanabaenaceae
    "Pseudanabaena" = "#FF7F00",
    #Microcoleceae
    "Limnoraphis" = "maroon",
    "Planktothrix" = "red",
    #Nostocaceae
    "Dolichospermum" = "gold1",
    #Gloeocapsaceae
    "Gloeocapsa" = "#643D9A",
    #Phormidesmiaceae
    "Nodosilinea" = "orchid1",
    "(Rare)" = "grey90",
    "Unassigned" = "grey70",
    "(Rare/Unknown)" = "grey80")

cy_group_colors <- 
  c("Deep Desert" = "black",
  "Peak Bloom" = "green4",
  "Eastern Basin" = "#088096",
  "Aug. Central Basin" = "maroon",
  "Sep. Central Basin" = "red",
  "Bloom Adjacent 1" = "gold1",
  "Bloom Adjacent 2" = "olivedrab2")

class_colors_for_ancom <- c("Rare" = "grey80",
                            "Alphaproteobacteria" = "green4",
                            "Gammaproteobacteria"= "gold1",
                            "Bacteroidia" = "#FF7F00",
                            "Kapabacteria"= "#FDBF6F",
                            "Actinobacteria" = "blue4",
                            "Verrucomicrobiae" = "orchid1",
                            "Chlamydiae"= "deeppink1",
                            "Cyanobacteriia"     = "olivedrab1",
                            "Chlorobia" = "blue1",
                            "Vicinamibacteria" = "turquoise",
                            "Bacilli"=  "steelblue2",
                            "Chloroflexia"        = "black",
                            "SL56 marine group"   = "grey20",
                            "JG30-KF-CM66"= "grey40",
                            "Chloroflexia"= "grey60",
                            "Phycisphaerae" = "tan4",
                            "Planctomycetes"= "tan3",
                            "OM190"= "tan",
                            "Hydrogenedentia" = "#6A3D9A",
                            "Bdellovibrionia"= "yellow4",
                            "Oligoflexia"=  "steelblue4",
                            "Armatimonadia"= "red",
                            "Gemmatimonadetes"= "cornflowerblue",
                            "Fibrobacteria"= "maroon",
                            "Polyangia"  = "brown",
                            "Nanoarchaeia" = "purple"
)
