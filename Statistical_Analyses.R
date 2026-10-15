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

Prev <- read.csv("Final_Prevalence_Data.csv",sep=";")

phy <- read.nexus("mammalian_tree.tre")
phy <- drop.tip(phy, phy$tip.label[!phy$tip.label%in%unique(Prev$Species)])
phy <- force.ultrametric(phy, method = c("extend"))
phy_matrix <- ape::vcv.phylo(phy, corr = TRUE)


################################################################################
### CORRELATION WITH VINCZE ET AL (2021) 


#FOR THIS PART, PLEASE DOWNLOAD VINCZE ET AL'S CMR ESTIMATION DATA 
#FROM https://github.com/OrsolyaVincze/VinczeEtal2021Nature/blob/main/data.csv
CMR_data <- read_csv("Vincze_2022_CMR.csv")

Prevalence_General <- subset(Prev, Group=="All Organs")

Comp_data <- merge(subset(Prevalence_General, select=c("Species","Prevalence")),
                   subset(CMR_data, select=c("Species","CMR")),
                   by="Species", all.x=TRUE)
Comp_data <- na.omit(Comp_data)

row.names(Comp_data) <- Comp_data$Species
Comp_data <- Comp_data[,-1]

result = cor.test(Comp_data$Prevalence, Comp_data$CMR, method = "spearman")
print(result)


################################################################################
### PHYLOGENETIC SIGNAL AND ORDER COMPARISON

rm(CMR_data,Comp_data,Prevalence_General,result)

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

#Prev2 <- subset(Prev, Species!="Phascolarctos_cinereus") #add to run without Koala

# TO REPETE FOR EVERY STUDIED ORGAN SYSTEMS
phylosig_ordercomp(subset(Prev, Group == "All Organs"),phy)


################################################################################
### BAYESIAN REGRESSIONS (DIET AND PETO'S PARADOX)

Prev$Animal<- factor(Prev$Animal, levels= c("Rarely or never", "Primarily food item"))
Prev$Vertebrate<- factor(Prev$Vertebrate, levels= c("Rarely or never", "Primarily food item"))
Prev$Invertebrate<- factor(Prev$Invertebrate, levels= c("Rarely or never", "Primarily food item"))
Prev$Herptile<- factor(Prev$Herptile, levels= c("Rarely or never", "Primarily food item"))
Prev$Mammal<- factor(Prev$Mammal, levels= c("Rarely or never", "Primarily food item"))

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

#Prev2 <- subset(Prev, Species!="Phascolarctos_cinereus") #add to run without Koala

# TO REPETE FOR EVERY STUDIED ORGAN SYSTEMS
bmrs_stat(subset(Prev, Group == "All Organs"))


################################################################################
### XI TEST FOR TUMOR PROPORTION ANALYZE

Prevalence_Group <- subset(Prev, Group!="All Organs")

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
### ICM AND CMR PGLS REGRESSIONS

# Summary function GLS
sumGLS <- function(model){ 
  su <- as.data.frame(summary(model)$tTable)
  su[,1:3] <- round(su[,1:3],2); su[,4] <- round(su[,4],4)
  su[,1] <- paste(su[,1],'_(', su[,2],')', sep='')
  su <- su[,c(1,3,4)]
  x <- nrow(su)
  #su[x+1,1] <- paste("ShaP=",round(as.numeric(unlist(shapiro.test(resid(model)))[[2]][1]),4), sep="")
  su[x+1,2] <- paste("AIC=", round(AIC(model),2), sep="")
  su[x+1,3] <- paste("n=", length(resid(model)), sep='')
  su[x+1,1] <-paste('lambda=',round(as.numeric(model[[1]])[1], 2), sep='') # put lambda value in bottom right corner
  row.names(su)[nrow(su)] <- "ModelStats"
  return(su)
}

Statis_ICM <- function(PrevGroupsdata,phylotree){
  Prev <- PrevGroupsdata
  row.names(Prev) <- Prev$Species
  
  ### Association between cancer risk and diet animal content
  # Models of non-zero prevalence in function of life expectancy, body mass and diet 
  dataPrev <- Prev[!is.na(Prev$ICM),]; row.names(dataPrev) <- dataPrev$Species
  phyPrev <- treedata(phylotree,dataPrev, warnings = FALSE)$phy
  dataPrev$occurrence <- ifelse(dataPrev$ICM==0,0,1)
  dataPrev <- dataPrev[as.character(phyPrev$tip.label), ] # organise data to follow same order as phylogeny
  # Drop tips on tree:
  phy_nonzeroPrev <- drop.tip(phyPrev, which(!dataPrev$occurrence==1))
  # logistic data:
  data_nonzeroPrev <- dataPrev[dataPrev$Species%in%phy_nonzeroPrev$tip.label, ] # organise data to follow same order as phylogeny
  data_nonzeroPrev$ICM <- as.numeric(data_nonzeroPrev$ICM)
  data_nonzeroPrev$Life_Expectancy <- as.numeric(data_nonzeroPrev$Life_Expectancy)
  data_nonzeroPrev$N_Death <- as.numeric(data_nonzeroPrev$N_Death)
  data_nonzeroPrev$Body_Mass <- as.numeric(data_nonzeroPrev$Body_Mass)
  
  #models
  DietMod0 <- gls(logit(ICM) ~  log(Body_Mass)+log(Life_Expectancy), data_nonzeroPrev,
                  cor=corPagel(0,phy_nonzeroPrev, form = ~Species, fixed=F), weights = ~1/log(N_Death))
  DietMod1 <- gls(logit(ICM) ~  log(Body_Mass)+log(Life_Expectancy)+as.factor(Animal), data_nonzeroPrev,
                  cor=corPagel(0,phy_nonzeroPrev, form = ~Species, fixed=F), weights = ~1/log(N_Death))
  DietMod2 <- gls(logit(ICM) ~  log(Body_Mass)+log(Life_Expectancy)+as.factor(Vertebrate), data_nonzeroPrev,
                  cor=corPagel(0,phy_nonzeroPrev, form = ~Species, fixed=F), weights = ~1/log(N_Death))
  DietMod3 <- gls(logit(ICM) ~  log(Body_Mass)+log(Life_Expectancy)+as.factor(Invertebrate), data_nonzeroPrev,
                  cor=corPagel(0,phy_nonzeroPrev, form = ~Species, fixed=F), weights = ~1/log(N_Death))
  DietMod4 <- gls(logit(ICM) ~  log(Body_Mass)+log(Life_Expectancy)+as.factor(Mammal), data_nonzeroPrev,
                  cor=corPagel(0,phy_nonzeroPrev, form = ~Species, fixed=T), weights = ~1/log(N_Death))
  DietMod5 <- gls(logit(ICM) ~  log(Body_Mass)+log(Life_Expectancy)+as.factor(Herptile), data_nonzeroPrev,
                  cor=corPagel(0,phy_nonzeroPrev, form = ~Species, fixed=F), weights = ~1/log(N_Death))
  
  X <- as.data.frame(array(,c(1,3))); names(X) <- names(sumGLS(DietMod0))
  print("Summary of the GLS models assessing for effect of diet animal content")
  print(rbind(X, sumGLS(DietMod0), X, sumGLS(DietMod1),X, sumGLS(DietMod2), X, 
              sumGLS(DietMod3),X, sumGLS(DietMod5), X, sumGLS(DietMod4)))
}

Statis_CMR <- function(PrevGroupsdata,phylotree){
  Prev <- PrevGroupsdata
  row.names(Prev) <- Prev$Species
  
  ### Association between cancer risk and diet animal content
  # Models of non-zero prevalence in function of life expectancy, body mass and diet 
  dataPrev <- Prev[!is.na(Prev$CMR),]; row.names(dataPrev) <- dataPrev$Species
  phyPrev <- treedata(phylotree,dataPrev, warnings = FALSE)$phy
  dataPrev$occurrence <- ifelse(dataPrev$CMR==0,0,1)
  dataPrev <- dataPrev[as.character(phyPrev$tip.label), ] # organise data to follow same order as phylogeny
  # Drop tips on tree:
  phy_nonzeroPrev <- drop.tip(phyPrev, which(!dataPrev$occurrence==1))
  # logistic data:
  data_nonzeroPrev <- dataPrev[dataPrev$Species%in%phy_nonzeroPrev$tip.label, ] # organise data to follow same order as phylogeny
  data_nonzeroPrev$CMR <- as.numeric(data_nonzeroPrev$CMR)
  data_nonzeroPrev$Life_Expectancy <- as.numeric(data_nonzeroPrev$Life_Expectancy)
  data_nonzeroPrev$N_Death <- as.numeric(data_nonzeroPrev$N_Death)
  data_nonzeroPrev$Body_Mass <- as.numeric(data_nonzeroPrev$Body_Mass)
  
  #models
  DietMod0 <- gls(logit(CMR) ~  log(Body_Mass)+log(Life_Expectancy), data_nonzeroPrev,
                  cor=corPagel(0,phy_nonzeroPrev, form = ~Species, fixed=F), weights = ~1/log(N_Death))
  DietMod1 <- gls(logit(CMR) ~  log(Body_Mass)+log(Life_Expectancy)+as.factor(Animal), data_nonzeroPrev,
                  cor=corPagel(0,phy_nonzeroPrev, form = ~Species, fixed=F), weights = ~1/log(N_Death))
  DietMod2 <- gls(logit(CMR) ~  log(Body_Mass)+log(Life_Expectancy)+as.factor(Vertebrate), data_nonzeroPrev,
                  cor=corPagel(0,phy_nonzeroPrev, form = ~Species, fixed=F), weights = ~1/log(N_Death))
  DietMod3 <- gls(logit(CMR) ~  log(Body_Mass)+log(Life_Expectancy)+as.factor(Invertebrate), data_nonzeroPrev,
                  cor=corPagel(0,phy_nonzeroPrev, form = ~Species, fixed=F), weights = ~1/log(N_Death))
  DietMod4 <- gls(logit(CMR) ~  log(Body_Mass)+log(Life_Expectancy)+as.factor(Mammal), data_nonzeroPrev,
                  cor=corPagel(0,phy_nonzeroPrev, form = ~Species, fixed=T), weights = ~1/log(N_Death))
  DietMod5 <- gls(logit(CMR) ~  log(Body_Mass)+log(Life_Expectancy)+as.factor(Herptile), data_nonzeroPrev,
                  cor=corPagel(0,phy_nonzeroPrev, form = ~Species, fixed=F), weights = ~1/log(N_Death))
  
  X <- as.data.frame(array(,c(1,3))); names(X) <- names(sumGLS(DietMod0))
  print("Summary of the GLS models assessing for effect of diet animal content")
  print(rbind(X, sumGLS(DietMod0), X, sumGLS(DietMod1),X, sumGLS(DietMod2), X, 
              sumGLS(DietMod3),X, sumGLS(DietMod5), X, sumGLS(DietMod4)))
}

#Prev2 <- subset(Prev, Species!="Phascolarctos_cinereus") #add to run without Koala

# TO REPETE FOR EVERY STUDIED ORGAN SYSTEMS
Statis_ICM(subset(Prev, Group=="All Organs"),phy)
Statis_CMR(subset(Prev, Group=="All Organs"),phy)
