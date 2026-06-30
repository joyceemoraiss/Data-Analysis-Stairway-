library(ape)
library(phytools)
library(picante)
library(phylocom)
library(vegan)
library(PCPS)
library (SYNCSA)
library (devtools)
library (betapart)
library (reshape2)

#Load community data, enviromental, species pool list and the original tree
data <- read.xlsx(here::here("")) %>% 
envir<-read.table(here::here(""))
spp<-read.table(file.choose(""))
tree<-read.tree(file.choose(""))

#Math tree with species pool and prune tree
comparison<-name.check(tree,spp,data.names=spp[,1])
comparison #check if all species in your pool appears in original tree
comparison$tree_not_data# Check the species present in the original tree and absent from the species list.
pruned.tree<-drop.tip(tree,comparison$tree_not_data) #remove the species present in the original tree and absent from the species list.

# Ensure the tree is ultrametric, rooted and dicotomic
is.ultrametric(tree) 
is.rooted(tree) 
is.binary(tree)

#Phylogenetic distance 
dist_phylo<-cophenetic(tree)
dist_phylo

## Faith's Phylogenetic Diversity (PD) and  Standardized Effect Size of PD (SES.PD)
pd <- pd(dist_phylo, tree = tree, TRUE)
pd_values<-ses.pd(data,tree)

## Net Relatedness Index (NRI) using MPD
nri_comm<-ses.mpd(data,dist_phylo,null.model="independentswap",runs=499)
nri.valor<-nri_comm$mpd.obs.z*-1
# Nearest Taxon Index (NTI) using MNTD
nti_comm<-ses.mntd(data,dist_phylo,null.model="independentswap",runs=499)
nti.valor<-nti_comm$mntd.obs.z*-1

#organize values
ses.pd_strata<-data_ses.pd<-read.table("")
ses.pd_habitat<-data_ses.pd<-read.table("")

nti_strata<-read.table("")
nti_habitat<-read.table("") 


nri_strata<-read.table("")
nri_habitat<-read.table("")

# T-tests (for vertical strata: canopy vs understory)
#PD
residuals_pdstrata <- lm(ses.pd ~ strata, data = ses.pd_strata)
residuals_model <- residuals(residuals_pdstrata)
shapiro.test(residuals_model)
t.test(ses.pd ~ strata, data = ses.pd_strata, var.equal = TRUE)
#NRI
residuals_nristrata <- lm(nri ~ strata, data = nri_strata)
residuals_modelnri <- residuals(residuals_nritrata)
shapiro.test(residuals_modelnri)
t.test(nri ~ strata, data = nri_strata, var.equal = TRUE)
#NTI
residuals_ntistrata <- lm(nti ~ strata, data = nti_strata)
residuals_ntistrata
residuals_modelnti <- residuals(residuals_ntitrata)
shapiro.test(residuals_modelnti)
t.test(nti ~ strata, data = nti_strata, var.equal = TRUE)

#Two-way ANOVA Habitat Type
#PD
anova_pd <- aov(ses.pd ~ habitat, data = ses.pd_habitat) 
shapiro.test(residuals(anova_pd))
bartlett.test(ses.pd ~ habitat, data = ses.pd_habitat)
anova(anova_pd)
#NRI
anova_nri<- aov(nri ~ habitat, data = nri_habitat) 
shapiro.test(residuals(anova_nri))
bartlett.test(nri ~ habitat, data = nri_habitat)
anova(anova_nri)
#NTI
anova_nti<- aov(nti ~ habitat, data = nti_habitat) 
shapiro.test(residuals(anova_nti))
bartlett.test(nti ~ habitat, data = nti_habitat)
anova(anova_nti)

#PCPS

# organizing the data 
org <- organize.pcps(comm = data, phylodist = dist.phylo, envir = envir)
envir <- org$environmental
comm <- org$community
dist.phylo <- org$phylodist

## extracting the P matrix - community matrix based on phylogenetic information
P_comm <- SYNCSA::matrix.p(as.matrix(data), as.matrix(dist.phylo))$matrix.P
P_comm

# performing the PCPS
pcps_comm <- pcps(data, dist.phylo, method = "bray", squareroot = T)
pcps_comm$values %>% 
filter(Relative_eig > 0.05) # os 4 primeiros eixos que possuem mais de 5% de explicação dos 

# Autovectors (PCPS) 
pcps <- pcps_comm$vectors
pcps

# correlation between PCPS axis and matrix.P
pcps_correl <- pcps_comm$correlations
pcps_correl

# score vectors for sites and species for axis with more than 5% of explanation
scores_pcps <- scores.pcps(pcps_comm, choices = c(1, 2)) # axis 1 and 2
scores_pcps3 <- scores.pcps(pcps_comm, choices = c(3, 4)) # axis 3 and 4

plot(pcps_comm, choices = c(1, 2))

## LME - PCPS (Debastiani & Duarte 2014): testing the effect of Strata (with strata as fixed variable and habitat type as random)
# PCPS 1 
res1 <- pcps.sig(comm = data, phylodist = dist.phylo, FUN = FUN.LME.marginal, 
                 formula = pcps.1 ~ strata , envir = envir, 
                 random = ~ 1|habitat, choices = 1, runs = 999)

summary(res1$model)

# PCPS 2
res2 <- pcps.sig(comm = data, phylodist = dist.phylo, FUN = FUN.LME.marginal, 
                 formula = pcps.2 ~ strata , envir = envir, 
                 random = ~ 1|habitat, choices = 2, runs = 999)

summary(res2$model)

# PCPS 3
res3 <- pcps.sig(comm = data, phylodist = dist.phylo, FUN = FUN.LME.marginal, 
                 formula = pcps.3 ~ strata , envir = envir, 
                 random = ~ 1|habitat, choices = 3, runs = 999)

summary(res3$model)

# PCPS 4
res4 <- pcps.sig(comm = data, phylodist = dist.phylo, FUN = FUN.LME.marginal, 
                 formula = pcps.4 ~ strata , envir = envir, 
                 random = ~ 1|habitat, choices = 4, runs = 999)

summary(res4$model)


