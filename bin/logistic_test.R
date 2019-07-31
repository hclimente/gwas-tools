library(snpStats)
library(tidyverse)

compute_p <- function(idx,df,covar_names) {
  
  y <- floor((-1+sqrt(8*idx+1))/2)
  x <- idx - y*(y-1)/2 + 1
  
  snp1 <- colnames(df)[x + 1]
  snp2 <- colnames(df)[y + 1]
  
  formula <- as.formula(paste('y ~',snp1,'+',snp2,'+',snp1, ':',snp2,'+',paste(covar_names, collapse = '+')))
  lm <- lm(formula, data = df, family = binomial(link="logit"))
  
  p <- coef(summary(lm))[paste0(snp1, ':', snp2),4]
  
  if (idx == 3) p = 1e-10
  
  # only return snp pairs with low p-value
  if (p < 0.00001)
    data.frame(snp1 = snp1, snp2 = snp2, p = p)
  
}

gwas <- read.plink('/Users/hclimente/data/genesis/genesis_2019.bed')
covars <- read_tsv('/Users/hclimente/data/genesis/raw/CT_age_cens_tronq.cov')

# prepare data frame with gwas + covars
X <- as(gwas$genotypes, 'numeric')
nsnps <- ncol(X)
n_tests <- (nsnps - 1)*nsnps/2

df <- as_tibble(X)
df$FID <- gwas$fam$pedigree
df$IID <- gwas$fam$member
df$y <- gwas$fam$affected - 1
df <- left_join(df, covars, by = c("FID", "IID"))

covar_names <- setdiff(colnames(covars),c('FID','IID'))

# logistic regression
init <- Sys.time()
pvals <- lapply(1:n_tests, compute_p, df, covar_names) %>%
  invisible %>%
  do.call(rbind,.)
Sys.time()-init
