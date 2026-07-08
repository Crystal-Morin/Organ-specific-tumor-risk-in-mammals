################################################################################
##            ORGAN-SPECIFIC TUMOR RISK STATISTICAL ANALYSES                  ##
################################################################################

rm(list=ls())
setwd("C:/Users/cmorin09/Desktop/Organ-specific tumor risk")

library(ape)
library(geiger)
library(car)
library(phytools)
library(nlme)
library(emmeans)
library(brms)
library(readr)
library(dplyr)
library(stringr)

Prev <- read_delim("Final_Prevalence_Data.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)

phy <- read.nexus("mammalian_tree.tre")
phy <- drop.tip(phy, phy$tip.label[!phy$tip.label%in%unique(Prev$Species)])
phy <- force.ultrametric(phy, method = c("extend"))
phy_matrix <- ape::vcv.phylo(phy, corr = TRUE)


################################################################################
### PHYLOGENETIC SIGNAL AND ORDER COMPARISON

phylosig_ordercomp <- function(Prevdata, phylotree){
  dataPrev <- Prevdata[!is.na(Prevdata$Vertebrate),]; row.names(dataPrev) <- dataPrev$Species
  phyPrev <- treedata(phylotree,dataPrev, warnings = FALSE)$phy
  
  ### Phylogenetic signal
  Prev <- Prevdata$Prevalence ; names(Prev) <- Prevdata$Species
  Prev <- Prev[!is.na(Prev)]
  print(paste("Number of species =",length(Prev)))
  phylo <- treedata(phylotree,Prev, warnings = F)$phy
  print("Phylogenetic signal")
  print(phylosig(phylo, Prev, method = "lambda", test = T))
  
  ### Linear regression exploring order differences in tumor prevalence
  mx <- lm(Prevalence*100 ~ Order, Prevdata[Prevdata$Order%in%names(table(Prevdata$Order))[which(table(Prevdata$Order)>3)],])
  newdata = expand.grid(Order=names(table(Prevdata$Order))[which(table(Prevdata$Order)>3)])
  print("Order average prevalence")
  print(cbind(newdata,
              as.data.frame(emmeans(mx,~Order))))
  print("Pairwise order differences in prevalence")
  print(as.data.frame(pairs(emmeans(mx,~Order))))
}

# TO REPETE FOR EVERY STUDIED ORGAN SYSTEMS
phylosig_ordercomp(subset(Prev, Group == "All Organs" & Species != "Phascolarctos_cinereus"),phy)
#phylosig_ordercomp(subset(Prev, Group == "All Organs"),phy)


################################################################################
### BAYESIAN REGRESSIONS (DIET AND PETO'S PARADOX)

levels_diet <- function(data){
  cols <- c("Animal", 'Invertebrate', "Vertebrate", "Herptile", "Mammal")
  for (col in cols) {
    data[[col]] <- factor(data[[col]], levels=c("Rarely or Never","Primary or secondary food item"))
  }
  data
}

Prev <- levels_diet(Prev)

bmrs_stat <- function(data){
  priors <- default_prior(Occurence | trials(N_Death) ~ log(Life_Expectancy) + log(Body_Mass) + Necropsy_Effort
                          + (1|gr(Species, cov = phy_matrix)), 
                          data = data, 
                          data2 = list(phy_matrix = phy_matrix),
                          family = binomial(link = "logit"))
  mod <- brm(formula = Occurence | trials(N_Death) ~ log(Life_Expectancy) + log(Body_Mass) + Necropsy_Effort
             + (1|gr(Species, cov = phy_matrix)), 
             data = data,
             family = beta_binomial(),
             data2 = list(phy_matrix = phy_matrix),
             prior = priors,
             iter = 4000, chains = 4, cores = 4,
             control = list(adapt_delta = 0.99))
  
  print(summary(mod))
  
  priors <- default_prior(Occurence | trials(N_Death) ~ Animal + log(Life_Expectancy) + log(Body_Mass) + Necropsy_Effort
                          + (1|gr(Species, cov = phy_matrix)), 
                          data = data, 
                          data2 = list(phy_matrix = phy_matrix),
                          family = binomial(link = "logit"))
  mod <- brm(formula = Occurence | trials(N_Death) ~ Animal + log(Life_Expectancy) + log(Body_Mass) + Necropsy_Effort
             + (1|gr(Species, cov = phy_matrix)), 
             data = data,
             family = beta_binomial(),
             data2 = list(phy_matrix = phy_matrix),
             prior = priors,
             iter = 4000, chains = 4, cores = 4,
             control = list(adapt_delta = 0.99))
  
  print(summary(mod))
  
  priors <- default_prior(Occurence | trials(N_Death) ~ Vertebrate+ log(Life_Expectancy) + log(Body_Mass) + Necropsy_Effort
                          + (1|gr(Species, cov = phy_matrix)), 
                          data = data, 
                          data2 = list(phy_matrix = phy_matrix),
                          family = binomial(link = "logit"))
  mod <- brm(formula = Occurence | trials(N_Death) ~ Vertebrate + log(Life_Expectancy) + log(Body_Mass) + Necropsy_Effort
             + (1|gr(Species, cov = phy_matrix)), 
             data = data,
             family = beta_binomial(),
             data2 = list(phy_matrix = phy_matrix),
             prior = priors,
             iter = 4000, chains = 4, cores = 4,
             control = list(adapt_delta = 0.99))
  
  print(summary(mod))
  
  priors <- default_prior(Occurence | trials(N_Death) ~ Invertebrate + log(Life_Expectancy) + log(Body_Mass) + Necropsy_Effort
                          + (1|gr(Species, cov = phy_matrix)), 
                          data = data, 
                          data2 = list(phy_matrix = phy_matrix),
                          family = binomial(link = "logit"))
  mod <- brm(formula = Occurence | trials(N_Death) ~ Invertebrate + log(Life_Expectancy) + log(Body_Mass) + Necropsy_Effort
             + (1|gr(Species, cov = phy_matrix)), 
             data = data,
             family = beta_binomial(),
             data2 = list(phy_matrix = phy_matrix),
             prior = priors,
             iter = 4000, chains = 4, cores = 4,
             control = list(adapt_delta = 0.99))
  
  print(summary(mod))
  
  priors <- default_prior(Occurence | trials(N_Death) ~ Herptile + log(Life_Expectancy) + log(Body_Mass) + Necropsy_Effort
                          + (1|gr(Species, cov = phy_matrix)), 
                          data = data, 
                          data2 = list(phy_matrix = phy_matrix),
                          family = binomial(link = "logit"))
  mod <- brm(formula = Occurence | trials(N_Death) ~ Herptile + log(Life_Expectancy) + log(Body_Mass) + Necropsy_Effort
             + (1|gr(Species, cov = phy_matrix)), 
             data = data,
             family = beta_binomial(),
             data2 = list(phy_matrix = phy_matrix),
             prior = priors,
             iter = 4000, chains = 4, cores = 4,
             control = list(adapt_delta = 0.99))
  
  print(summary(mod))
  
  priors <- default_prior(Occurence | trials(N_Death) ~ Mammal + log(Life_Expectancy) + log(Body_Mass) + Necropsy_Effort
                          + (1|gr(Species, cov = phy_matrix)), 
                          data = data, 
                          data2 = list(phy_matrix = phy_matrix),
                          family = binomial(link = "logit"))
  mod <- brm(formula = Occurence | trials(N_Death) ~ Mammal + log(Life_Expectancy) + log(Body_Mass) + Necropsy_Effort
             + (1|gr(Species, cov = phy_matrix)), 
             data = data,
             family = beta_binomial(),
             data2 = list(phy_matrix = phy_matrix),
             prior = priors,
             iter = 4000, chains = 4, cores = 4,
             control = list(adapt_delta = 0.99))
  
  print(summary(mod))
  
}


# TO REPETE FOR EVERY STUDIED ORGAN SYSTEMS
bmrs_stat(subset(Prev, Group == "All Organs" & Species != "Phascolarctos_cinereus"))
#bmrs_stat(subset(Prev, Group == "All Organs"))


################################################################################
### XI TEST FOR TUMOR PROPORTION ANALYZE

Prevalence_Group <- subset(Prev, Group!="All Organs" & Species!="Phascolarctos_cinereus")

sp <- unique(Prevalence_Group$Species)
occ <- c()
for (k in 1:length(sp)){
  occ <- c(occ, sum(subset(Prevalence_Group, Species==sp[k])$Occurence))
}
Occs <- data.frame(sp,occ)
colnames(Occs)<-c("Species","Occurence_tot")
Prevalence_Group <- merge(Prevalence_Group,Occs,by="Species",all.x=TRUE)
rm(Occs,k,occ,sp)

#Tumor cases threshold
Prevalence_Group <- subset(Prevalence_Group, Occurence_tot>=10)


Group_Abundance <- data.frame(row.names = unique(Prevalence_Group$Species))
Groups <- unique(Prevalence_Group$Group)
for(k in 1:length(Groups)){
  dgroup <- subset(Prevalence_Group, Group==Groups[k])
  Group_Abundance <- cbind(Group_Abundance, dgroup$Occurence)
}
colnames(Group_Abundance) <- Groups
rm(dgroup,k,Groups)

Group_Abundance$"Other Organs" <- Group_Abundance$Bone + Group_Abundance$`Circulatory system` +
  Group_Abundance$`Nervous system` + Group_Abundance$Skin + Group_Abundance$`Fat tissues` +
  Group_Abundance$`Exocrine system` + Group_Abundance$`Diaphragm`

Group_Abundance <- subset(Group_Abundance, select = c("Lymphatic system",'Digestive system',
                                                      'Uro-genital system','Endocrine system',
                                                      "Respiratory system","Other Organs"))


chisq.test(Group_Abundance)


################################################################################
### CORRELATION WITH VINCZE ET AL (2021) 


#FOR THIS PART, PLEASE DOWNLOAD VINCZE ET AL'S CMR ESTIMATION DATA 
#FROM https://github.com/OrsolyaVincze/VinczeEtal2021Nature/blob/main/data.csv
CMR_data <- read_csv("CMR_Data.csv")

Prevalence_General <- subset(Prev, Group=="All Organs" & Species != "Phascolarctos_cinereus")

Comp_data <- merge(subset(Prevalence_General, select=c("Species","Prevalence")),
                   subset(CMR_data, select=c("Species","CMR")),
                   by="Species", all.x=TRUE)
Comp_data <- na.omit(Comp_data)

row.names(Comp_data) <- Comp_data$Species
Comp_data <- Comp_data[,-1]

result = cor.test(Comp_data$Prevalence, Comp_data$CMR, method = "spearman")
print(result)
