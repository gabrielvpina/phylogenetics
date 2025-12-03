library(dplyr)
library(ggplot2)
library(viridis)
library(ggtree)
library(seqinr)


# Checar fasta


# checar matriz de alinhamento
algn <- read.alignment(file = "mafft_result/mafft_res.fasta", format = "fasta")

# matriz de distancia
distance_matrix_identity <- dist.alignment(algn, matrix = "identity", gap = TRUE)
print(distance_matrix_identity)

mtx <- as.matrix(distance_matrix_identity)

heatmap(mtx, Colv = NA, Rowv = NA, scale="column")

# plotar árvore filogenétic
arvre <- read.tree("iqtree_result/sem_trimar/mafft_res.fasta.treefile") 


ggtree(arvre) + 
  geom_treescale() +
  geom_tiplab()



# ===============================================================
# GUIA PARA PROCESSO COMPLETO
# LINK: https://rpubs.com/mvillalobos/L01_Phylogeny
# ===============================================================

