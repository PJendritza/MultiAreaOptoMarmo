# Fig8_b_stats
# R code for pairwise comparisons of proportions with multiple comp. corr. (Benjamini Hochberg)

nhits  <- c(156.0520, 108.0000, 109.0000, 66.0000, 89.0000, 53.0000, 44.0000)
ntrls <- c(467, 115, 123, 118, 109, 119, 116)

pairwise.prop.test(nhits , ntrls, p.adjust.method = "BH")


# run here: https://rdrr.io/snippets/
#
# Result:
# 
# Pairwise comparisons using Pairwise comparison of proportions
# 
# data:  nhits out of ntrls
# 
#   1       2       3       4       5       6
# 2 < 2e-16 -       -       -       -       -
# 3 < 2e-16 0.250   -       -       -       -
# 4 1.9e-05 1.7e-10 5.7e-08 -       -       -
# 5 < 2e-16 0.013   0.221   9.5e-05 -       -
# 6 0.041   4.0e-15 2.6e-12 0.129   3.6e-08 -
# 7 0.419   < 2e-16 4.0e-15 0.013   1.7e-10 0.389
# 
# P value adjustment method: BH