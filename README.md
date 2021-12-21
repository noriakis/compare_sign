# compare_sign

## Download
https://ytlab.jp/sign/signbn/download.html

## Command
```
./signbn.1.8.3 -N 500 input.txt -o bs/result
```

## R
```R
## import SiGN result
source("loadSign.R")
signBn <- loadSign("result")

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

## References
- Scutari, M. (2010). Learning Bayesian Networks with the bnlearn R Package. Journal of Statistical Software, 35(3), 1–22. https://doi.org/10.18637/jss.v035.i03
- Tamada, Y., Shimamura, T., Yamaguchi, R., Imoto, S., Nagasaki, M., and Miyano, S. (2011). SiGN: Large-scale gene network estimation environment for high performance computing, Genome Informatics, 25 (1), 40-52.
