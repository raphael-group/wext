# Load required modules
library(BiRewire)
library(igraph)
permute <- function (inputFile, outputDirectory, numPermutations, start_index=1, exact=TRUE, header=FALSE){
	# Read the graph and convert it to bipartite
	edgelist <- read.table(inputFile, heade=header)
	G <- graph.data.frame(edgelist)
	V(G)$type <- V(G)$name %in% edgelist[,1]
	
	# Extract the bipartite adjacency matrix
	M <- get.incidence(G)
	
	# Generate a permutation per 
	for (i in 0:(numPermutations-1) ){
		P <- birewire.rewire.bipartite(M, exact=exact)
		H <- graph.incidence(P)
		write.graph(H,
            file=sprintf("%s/permuted-edgelist-%d.txt", outputDirectory, i + start_index), format="ncol")
	}
}

# Parse command line args and run permute
args <- commandArgs(trailingOnly = TRUE)
permute(args[1], args[2], as.integer(args[3]), as.integer(args[4]))