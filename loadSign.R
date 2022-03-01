#' loadSign
#'
#' Load the output of SiGN-BN (HC+BS)
#'
#' @param fileName the result of SiGN-BN
#'
#' @return list of edges, nodes, strength, and bn (bnlearn)
#' @examples loadSign("result.txt")
#' @export
#'
loadSign <- function(fileName){
    returnList <- list()
    rawStr <- readChar(fileName, file.info(fileName)$size)
    edges <- read.csv(text=unlist(strsplit(rawStr, "\\[Edges]\n"))[2], sep="\t", header=F)
    nodes <- read.csv(text=unlist(strsplit(unlist(strsplit(rawStr, "\\[Edges]\n"))[1], "\\[Nodes]\n"))[2], sep="\t", header=F)
    
    changeName <- list()
    for (i in seq_len(dim(nodes)[1])){
        changeName[[as.character(nodes[i,]$V3)]] <- nodes[i,]$V1
    }
    
    edges$V1 <- sapply(edges$V1, function(x) changeName[[as.character(x)]])
    edges$V2 <- sapply(edges$V2, function(x) changeName[[as.character(x)]])
    
    signStr <- edges[,1:3]
    colnames(signStr) <- c("from","to","strength")
    attr(signStr, "nodes") <- unique(c(signStr$from, signStr$to))
    signStr = structure(signStr, method = "bootstrap", threshold = 0, class = c("bn.strength", class(signStr)))
    signBn <- averaged.network(signStr, threshold=bnlearn::inclusion.threshold(signStr))
    
    returnList[["str"]] <- signStr
    returnList[["edges"]] <- edges
    returnList[["nodes"]] <- nodes
    returnList[["av"]] <- signBn
    return(returnList)
}


#' loadIngor
#'
#' Load the output of ingor
#'
#' @param fileName the result of ingor
#'
#' @return list of strength, and bn (bnlearn)
#' @examples loadIngor("result.txt")
#' @export
#'
loadIngor <- function(fileName){
    returnList <- list()
    
    df <- read.table(fileName, header=1)
    ingStr <- df[,c(1,2,6)]
    colnames(ingStr) <- c("from","to","strength")
    attr(ingStr, "nodes") <- unique(c(ingStr$from, ingStr$to))
    ingStr = structure(ingStr, method = "bootstrap", threshold = 0, class = c("bn.strength", class(ingStr)))
    ingBn <- averaged.network(ingStr, threshold=bnlearn::inclusion.threshold(ingStr))
    
    returnList[["str"]] <- ingStr
    returnList[["av"]] <- ingBn
    return(returnList)
}