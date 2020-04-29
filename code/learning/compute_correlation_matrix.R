# function to generate the feature correlation matrix


library(Hmisc)
library(RcmdrMisc)

calc_corr_matrix <- function(data_corr){
  r <- rcorr(as.matrix(data_corr), type="spearman")

  adjusted <- p.adjust(r$P, method = "holm")
  r$P <- adjusted

  flattenCorrMatrix <- function(cormat, pmat) {
    ut <- upper.tri(cormat)
    data.frame(
      row = rownames(cormat)[row(cormat)[ut]],
      column = rownames(cormat)[col(cormat)[ut]],
      cor  =(cormat)[ut],
      p = pmat[ut]
    )
  }

  new_r <- flattenCorrMatrix(r$r, r$P) %>% 
    filter(cor==1) %>% 
    filter(p<0.01) %>% 
    write_csv("data/process/sig_flat_corr_matrix.csv")
}


