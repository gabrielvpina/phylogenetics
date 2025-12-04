library(dplyr)
library(ggplot2)
library(viridis)
library(ggtree)
library(seqinr)
library(hrbrthemes)
library(phangorn)
library(ape)
library(knitr)

#=====================================================================
# Checar propriedades fasta
#=====================================================================
fastaLen <- read.table("sequencias.tsv", sep="\t", header=T)
# tamanho
mediaSeqs <- mean(fastaLen$Tamanho)
sdSeqs <- sd(fastaLen$Tamanho)
# conteudo gc
mediaGC <- mean(fastaLen$Conteudo_GC)
sdGC <- sd(fastaLen$Conteudo_GC)

## plot de tamanho das sequencias
seqLen <- ggplot(fastaLen, aes(x=Nome_da_Sequencia, y=Tamanho)) +
  geom_bar(stat = "identity", fill="steelblue", width=0.8) +
  xlab("Species sequence") +
  ylab("Tamanho da sequência (nt)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
print(seqLen)

seqGC <- ggplot(fastaLen, aes(x=Nome_da_Sequencia, y=Conteudo_GC)) +
  #geom_line(stat="identity",color="steelblue",size=1,group=TRUE) +
  geom_bar(stat = "identity", fill="steelblue", width=0.8) +
  xlab("Species sequence") +
  ylab("Conteúdo GC (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
print(seqGC)

# checar outliers
boxplot(fastaLen$Tamanho)

#=====================================================================
# Checar matriz de alinhamento
#=====================================================================
my_align <- read.alignment(file = "mafft_result/mafft_res.fasta", format = "fasta")

# matriz de distancia
distance_matrix_identity <- dist.alignment(algn, matrix = "identity", gap = TRUE)
library(knitr)
#print(distance_matrix_identity)

mtx <- as.matrix(distance_matrix_identity)

# Matriz simples de identidade (quanto mais escuro, menos similar)
heatmap(mtx, Colv = NA, Rowv = NA, scale="column")


# ===============================================================
# ANALISE FILOGENÉTICA
# GUIA PARA PROCESSO COMPLETO
# LINK: https://rpubs.com/mvillalobos/L01_Phylogeny
# ===============================================================


# Modelos de Substituição de Nucleotídeos
# Compute the p distance matrix for the aligned sequences
library(gt)
D_JK69 = dist.dna(as.DNAbin(my_align), model = "JC69")

# Convert the distance matrix to a matrix and round the values
D_JK69 = round(as.matrix(D_JK69),2)
D_JK69_df <- as.data.frame(D_JK69)

# Print the JC69 as an HTML table with a 20% width
#kable(D_JK69, format = "html", table.attr = "style='width:20%;'")
heatmap(D_JK69)


#=====================================================================
# Máxima parcimônia
#=====================================================================

suppressWarnings({
  suppressMessages({
    options(verbose = FALSE)
    
    # Distance Calculation
    # Compute the Hamming distance matrix for the given aligned sequences
    # This will serve as a measure of pairwise sequence dissimilarity
    D_hamming = dist.hamming(my_align)

    # Base Tree Construction
    # Compute the Neighbor Joining (NJ) tree as a base tree
    # NJ is a distance-based method for constructing phylogenetic trees
    nj_tree = nj(D_hamming)

    # Correct Negative Branch Lengths
    # Negative branch lengths are not biologically meaningful,
    # hence set any negative branch lengths to 0
    nj_tree$edge.length[which(nj_tree$edge.length < 0)] = 0

    # Rooting the Tree
    # Root the NJ tree at its midpoint to ensure ultrametricity,
    # and convert any multifurcating nodes to bifurcating nodes
    nj_tree = midpoint(multi2di(nj_tree))

    # Parsimony Optimization
    # Optimize the parsimony of the NJ tree
    # The optim.parsimony function refines the tree topology to minimize the parsimony score
    par_tree = optim.parsimony(nj_tree, my_align)

    
    #Edge estimation
    par_tree=acctran(par_tree, my_align)
    par_tree=midpoint(par_tree)
    
    # Calculating and Printing Parsimony
    # Calculate and print the parsimony score of the optimized tree
    # Lower parsimony scores indicate more parsimonious (simpler) trees
    # print(paste("Parsimony Score for the optimized tree:", parsimony(par_tree, my_align)))
  })
})

tree_MP <- plot_tree(par_tree, "Parsimony tree",180)
tree_MP
