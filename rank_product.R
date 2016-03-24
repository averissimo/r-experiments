#
rm(list = ls())
#
source('rank_product/pvalue_exact.r')
source('rank_product/rank_product_1.R')
#
set.seed(1985)

original_vector <- 1:1000000

order1 <- sample(original_vector)
order2 <- sample(original_vector)
order3 <- sample(original_vector)

rank_prod <- as.double(order1) * as.double(order2) * as.double(order3)

res <- rankprodbounds(rank_prod, length(original_vector), 3, 'geometric')

# print(res < 0.05)
sum_under_0.05 <- sum(res<0.05)
print(paste0(sum_under_0.05, ' in ', length(original_vector), ' (', 100 * sum_under_0.05 / length(original_vector),'%)'))

# source("http://bioconductor.org/biocLite.R")
# biocLite('qvalue')

library('qvalue')

resq <- qvalue(res)

summary(resq)

