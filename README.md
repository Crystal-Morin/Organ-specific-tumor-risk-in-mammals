### This repository provides data and code to reproduce results presented in the article:
# "ORGAN-SPECIFIC PATTERNS OF TUMOR SUSCEPTIBILITY ACROSS MAMMALS"
### by Crystal Morin, David Sayag, Morgane Tidière, Karin Lemberger, Jean-Francois Lemaitre, Samuel Pavard, Fernando Colchero, Baptiste Mulot, Orsolya Vincze, Mathieu Giraudeau.

Files included are:
mammalian_tree.tre - Consensus phylogeny tree

Statistical_Analyses.R - R code to reproduce results presented in the article.

Final_Prevalence_Data.csv - Database containing species-specific necropsy effort, diet, body mass, adult life expectancy and organ-specific tumor risk. 
Included variables are as follows:
  
      Species     	      - Latin scientific name of the taxa.
	  Order       	      - Phylogenetic order of the species.
      Life_Expectancy	  - Average number of days lived after sexual maturity was reached, i.e. remaining life expectancy.
	  Body_Mass   	      - Average body mass of adults (kg).
	  RDI        	      - "Relevant death information", number of adult individuals per species that have been necropsied after 2013 (year at which the Medical Module started to be used).
	  N_Death	       	  - Number of individuals per species that died after 2013.
      Necropsy_Effort     - Estimation of the necrospy effort, = RDI/N_death.
	  Animal	          - Animal prey regular presence (Primary or secondary food item) or absence (Rarely or Never) in the diet. Following: Kissling et al. Ecol. Evol. 4, 2913–2930 (2014).
	  Vertebrate	      - Vertebrate prey regular presence (Primary or secondary food item) or absence (Rarely or Never) in the diet. Following: Kissling et al. Ecol. Evol. 4, 2913–2930 (2014).
	  Invertebrate	      - Invertebrate prey regular presence (Primary or secondary food item) or absence (Rarely or Never) in the diet. Following: Kissling et al. Ecol. Evol. 4, 2913–2930 (2014).
	  Herptile	          - Reptile and Amphibian prey regular presence (Primary or secondary food item) or absence (Rarely or Never) in the diet. Following: Kissling et al. Ecol. Evol. 4, 2913–2930 (2014).
	  Mammal              - Mammalian prey regular presence (Primary or secondary food item) or absence (Rarely or Never) in the diet. Following: Kissling et al. Ecol. Evol. 4, 2913–2930 (2014).
      Group      	      - Organ system, as described in the article. "All Organs" refers to all systems merged together.
      Occurence	          - Total number of individuals per species that were diagnosed with at least one case of tumor in the corresponding organ system. For "All Organs", total number of individuals per species that were diagnoses with at least one tumor case in any organ system.
      Proportion 	      - Per species, proportion of tumor cases that were diagnosed in the corresponding organ system, = corresponding organ system occurrence / sum of occurrences across all organ systems except "All Organs. Undefined for "All Organs". The variable takes values between 0 and 1.
	  Prevalence 	      - Per species, proportion of tumor cases that were diagnosed in the corresponding organ system, = corresponding organ system occurrence / N_death. The variable takes values between 0 and 1.


### Please note that in order to reproduce all statistical analyses, Cancer Mortality Risk data from Vincze et al. Nature (2021) is needed. This data can be dowloaded from :
[https://github.com/OrsolyaVincze/VinczeEtal2021Nature/blob/main/data.csv]
