#PERMANOVA for community-level multivariate comparisons
# Load libraries
library(microbiome)
library(ggplot2)
library(dplyr)

# Probiotics intervention example data 
pseq <- ps_final # Rename the example data

# Pick relative abundances (compositional) and sample metadata 
pseq.rel <- microbiome::transform(pseq, "compositional")
otu <- abundances(pseq.rel)
meta <- meta(pseq.rel)


p <- plot_landscape(pseq.rel, method = "NMDS", distance = "bray", col = "type", size = 3)
print(p)


#PERMANOVA significance test for group-level differences
#Now let us evaluate whether the group (probiotics vs. placebo) has a significant effect on overall gut microbiota composition. Perform PERMANOVA:

# samples x species as input
library(vegan)
permanova <- adonis(t(otu) ~ type,
                    data = meta, permutations=99, method = "bray")

# P-value
print(as.data.frame(permanova$aov.tab)["type", "Pr(>F)"])

#Checking the homogeneity condition
#Check that variance homogeneity assumptions hold (to ensure the reliability of the results):
# Note the assumption of similar multivariate spread among the groups
# ie. analogous to variance homogeneity
# Here the groups have signif. different spreads and
# permanova result may be potentially explained by that.
dist <- vegdist(t(otu))
anova(betadisper(dist, meta$type))

#Investigate the top factors
#Show coefficients for the top taxa separating the groups

coef <- coefficients(permanova)["type1",]
top.coef <- coef[rev(order(abs(coef)))[1:20]]
par(mar = c(3, 14, 2, 1))
barplot(sort(top.coef), horiz = T, las = 1, main = "Top taxa")






