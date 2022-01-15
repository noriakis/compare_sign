# compare_sign

Comparing the Bayesian network inferred by SiGN-BN (*B*-spline nonparametric regression) HC+BS and bnlearn (Gaussian).

## Download
https://ytlab.jp/sign/signbn/download.html

## Command
```
./signbn.1.8.3 -N 500 input.txt -o bs/result
```

## R
### Load SiGN results and perform boot.strength() in bnlearn
```R
## import SiGN result
source("loadSign.R")
signBn <- loadSign("result") # "result" is an output file of SiGN

## same input as SiGN
load("bwmod.rda")
inp <- bwmod$MEs
inp$ME0 <- NULL
bs <- boot.strength(inp, R=500, algorithm = "hc")
ave <- averaged.network(bs)

## compare
graphviz.compare(ave, signBn$av)
```
<img src="https://github.com/noriakis/compare_sign/blob/main/images/compare.png?raw=true" width="800px">

```R
tpEdges <- as.ggplot(as.grob(~graphviz.plot(as.bn(intersection(bnlearn::as.igraph(ave), bnlearn::as.igraph(signBn$av))))))
tpEdges
```
<img src="https://github.com/noriakis/compare_sign/blob/main/images/tpEdges.png?raw=true" width="800px">

## Extract TP edges between BNs, and annotate the modules by clusterProfiler

```R
## Extract TP edges
tp <- as.bn(intersection(bnlearn::as.igraph(ave), bnlearn::as.igraph(signBn$av)))$arcs
MEs <- unique(c(tp[,1], tp[,2]))

## ME number
MEnum <- unlist(strsplit(MEs,"ME"))
MEnum <- MEnum[MEnum!=""]

## Check EA
pathNum <- c() # No enrichment
pathName <- list()
for (i in MEnum){
    ensg <- names(bwmod$colors)
    candEnsg <- ensg[bwmod$colors==i]
    candEntre <- AnnotationDbi::select(org.Hs.eg.db, candEnsg, c("ENTREZID"), c("ENSEMBL"))$ENTREZID
    candEntre <- candEntre[!is.na(candEntre)]
    res <- ReactomePA::enrichPathway(candEntre, pAdjustMethod = "bonferroni")
    num <- dim(subset(res@result, p.adjust<0.05))[1]
    print(num)
    if (num==0){
        pathNum <- c(pathNum, i)
        pathName[[paste0("ME",i)]] <- paste0("ME",i)
    } else {
        pathName[[paste0("ME",i)]] <- paste0("ME", i, " (",res@result$Description[1],")")
    }
}
```
## Make ggraph, and make a wordcloud for the non-enriched modules
```R
## Extract igraph
ig <- intersection(bnlearn::as.igraph(ave), bnlearn::as.igraph(signBn$av))
ig <- delete.vertices(ig, !names(V(ig)) %in% MEs)

vPathName <- sapply(names(V(ig)), function(x) pathName[[x]])
vPathName[sapply(vPathName, function(x) is.null(x))] <- NA
V(ig)$path <- unlist(vPathName)

## Draw in ggraph
library(ggraph)
gg <- ggraph(ig, layout="sugiyama") +
    geom_edge_diagonal(
        angle_calc = "along",
        label_dodge=unit(3,'mm'),
        arrow=arrow(length=unit(4, 'mm')),
        end_cap=circle(5, 'mm')) +
    geom_node_point()+
    geom_node_text(aes_(label=~stringr::str_wrap(path, width = 25)),
                   check_overlap=TRUE, repel=TRUE, size = 4,
                   color = "black",
                   bg.color = "white", segment.color="white",
                   bg.r = .15)+
    theme_graph()

## Generate wordcloud for non-enriched modules
mywc <- list()
for (i in pathNum){
    ensg <- names(bwmod$colors)
    candEnsg <- ensg[bwmod$colors==i]
    candEntre <- AnnotationDbi::select(org.Hs.eg.db, candEnsg, c("ENTREZID"), c("ENSEMBL"))$ENTREZID
    candEntre <- candEntre[!is.na(candEntre)]
    mywc[[as.character(i)]] <- wcGeneSummary(candEntre, max.words=200, random.order=FALSE, excludeFreq=15000,
                  colors=palettetown::pokepal(sample(1:151,1)), shape="circle", rot.per=0.4)$wc +
        ggtitle(paste0("ME",i))
}

wrapped <- patchwork::wrap_plots(mywc)
gg / wrapped + plot_layout(height=c(3,7))
```
<img src="https://github.com/noriakis/compare_sign/blob/main/images/tpWithWc.png?raw=true" width="800px">


## References
- Scutari, M. (2010). Learning Bayesian Networks with the bnlearn R Package. Journal of Statistical Software, 35(3), 1–22. https://doi.org/10.18637/jss.v035.i03
- Tamada, Y., Shimamura, T., Yamaguchi, R., Imoto, S., Nagasaki, M., and Miyano, S. (2011). SiGN: Large-scale gene network estimation environment for high performance computing, Genome Informatics, 25 (1), 40-52.
- Chen,X. et al. (2019) 5-methylcytosine promotes pathogenesis of bladder cancer through stabilizing mRNAs. Nat. Cell Biol., 21, 978–990.
