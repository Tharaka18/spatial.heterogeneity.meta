####################################################################################################################################################################################################################################

# Article title: Crop and landscape heterogeneity increase biodiversity in agricultural landscapes: A global review and meta-analysis

# Journal: Ecology Letters

# Authors: Tharaka S. Priyadarshana (ORCID: 0000-0003-3962-5465), Emily A. Martin (0000-0001-5785-9105), Cl?lia Sirami (0000-0003-1741-3082), Ben A. Woodcock (0000-0003-0300-9951), Eben Goodale (0000-0003-3403-2847), Carlos Mart?nez-N??ez (0000-0001-7814-4985), Myung-Bok Lee (0000-0003-2680-5707), Emilio Pagani-N??ez (0000-0001-8839-4005), Chlo? A. Raderschall (0000-0003-2005-1705), Llu?s Brotons (0000-0002-4826-4457), Anushka Rege (0000-0002-8383-0258), Annie Ouin (0000-0001-7046-2719), Teja Tscharntke (0000-0002-4482-3178), Eleanor M. Slade (0000-0002-6108-1196)

# Contact: Tharaka S. Priyadarshana; tharakas001@e.ntu.edu.sg; tharakas.priyadarshana@gmail.com	

# All these data files are also publicly available in the Dryad Digital Repository, at https://doi.org/10.5061/dryad.dbrv15f7j (Priyadarshana et al. 2024). 

# Citation for data: Priyadarshana, T.S., Martin, E.A., Sirami, C., Woodcock, B.A., Goodale, E., Mart?nez-N??ez, C., et al. (2024). Crop and landscape heterogeneity increase biodiversity in agricultural landscapes: A global review and meta-analysis [dataset]. Dryad Digital Repository, https://doi.org/10.5061/dryad.dbrv15f7j. 

####################################################################################################################################################################################################################################

# There are eight datafiles containing the following metadata. To replicate our results, please run each data file with the corresponding R scripts

# Specific datafiles

## For invertebrate, vertebrate, and animal analyses: "dat.csv"

## For sensitivity analyses: "dat.csv"; "datCropArea.csv", "datSemiNaturalArea.csv", and "datOtherArea.csv". 

## For functional group analyses: "datPollinators.csv", "datPredators.csv", and "datPests.csv" 
 
## For climatic region and cropping system analyses: "dat.csv" and "datCroppingSystems.csv", 

####################################################################################################################################################################################################################################

# The spreadsheets contain following metadata that require analysis using the specified R scripts provided (see below). 

## Author: Citation for each study

## Study_ID: Identifier for each study

## Effect_Size_ID: Identifier for effect size

## Study_Country: Study country

## Study_Region: Study region

## Biome: Tropical/sub-tropical vs temperate farmlands

## Study_System: Dominant crops

## Verts_Inverts_Plants: Taxa categorisation into invertebrates, vertebrates, and plants

## Verts_Inverts_Plants_Without_Pests: Taxa categorisation into invertebrates, vertebrates, and plants, but excluding pest invertebrates and vertebrates   

## Animals_Plants: Taxa categorisation into animals (i.e.  invertebrates and vertebrates without pests) and plants

## Functional_Groups: Functional groups

## Orders: Taxonomic orders, including birds

## Taxa_Category_2: More details about the study taxa

## Biodiversity_Category: Biodiversity metrics (i.e. abundance, richness and Shannon diversity)

## Heterogeneity_Measure: Spatial heterogeneity components

## Heterogeneity_Mesure_Fine: Computation of the spatial heterogeneity components

## Heterogeneity_Type: Spatial heterogeneity component (i.e. crop compositional heterogeneity, crop configurational heterogeneity, landscape compositional heterogeneity, and landscape configurational heterogeneity)
Heterogeneity	Spatial heterogeneity types (i.e. compositional heterogeneity = crop and landscape compositional heterogeneity; configurational heterogeneity = crop and landscape configurational heterogeneity )

## Heterogeneity_Type_L_C: Land-use type (i.e. crop heterogeneity = crop compositional and configurational heterogeneity; landscape heterogeneity = landscape compositional and configurational heterogeneity)

## Pearson_Correlation: Effect size (i.e. Pearson correlation coefficient)

## Sample_Size: Sample size

## Scale: Spatial scales that the heterogeneity components are defined at

## Scale_Cat: Spatial scale (A_LocalScale [i.e. < 0.5 km radius area], LandscapeScale_1 [i.e. ? 0.5 km, but < 1 km], and LandscapeScale_2 [i.e. ? 1 km radius area])

## Seminatural_Cover_Mean: Mean semi-natural area (%)

## Crop_Cover_Mean: Mean cultivated area (%)

## Other_Cover_Mean: Mean other land-use area (%)

## Annual_Perennial: Annual vs perennial crops

## Article_Tittle: Article tittle 

## Notes: Additional notes

####################################################################################################################################################################################################################################

# R scripts

## 1 Publication_Bias (a R script testing for publication bias, outliers, and influential studies; this script should be run with ?dat? data file).

## 2 Effects_Of_Semi-Natural_Area (a R script testing for the influence of mean semi-natural area on the estimated average relationships between spatial heterogeneity and biodiversity; this script should be run with ?datSemiNaturalArea? data file)

## 3 Effects_Of_Crop_Area (a R script testing for the influence of mean cultivated area on the estimated average relationships between spatial heterogeneity and biodiversity; this script should be run with ?datCropArea? data file)

## 4 Effects_Of_Other_Land_Use_Area (a R script testing for the influence of mean other land-use area on the estimated average relationships between spatial heterogeneity and biodiversity; this script should be run with ?datOtherArea? data file)

## 5 Invertebrate_Abundance (a R script testing for the effect of spatial heterogeneity on invertebrate abundance; this script should be run with ?dat? data file).
 
## 6 Invertebrate_Richness (a R script testing for the effect of spatial heterogeneity on invertebrate richness; this script should be run with ?dat? data file).

## 7 Invertebrate_Diversity (a R script testing for the effect of spatial heterogeneity on invertebrate Shannon diversity; this script should be run with ?dat? data file).

## 8 Vertebrate_Abundance (a R script testing for the effect of spatial heterogeneity on vertebrate abundance; this script should be run with ?dat? data file).

## 9 Vertebrate_Richness (a R script testing for the effect of spatial heterogeneity on vertebrate richness; this script should be run with ?dat? data file).

## 10 Vertebrate_Diversity (a R script testing for the effect of spatial heterogeneity on vertebrate Shannon diversity; this script should be run with ?dat? data file).

## 11 Plant_Abundance (a R script testing for the effect of spatial heterogeneity on plant abundance; this script should be run with ?dat? data file).

## 12 Plant_Richness (a R script testing for the effect of spatial heterogeneity on plant richness; this script should be run with ?dat? data file).

## 13 Plant_Diversity (a R script testing for the effect of spatial heterogeneity on plant Shannon diversity; this script should be run with ?dat? data file).

## 14 Pollinator_Abundance (a R script testing for the effect of spatial heterogeneity on pollinator abundance; this script should be run with ?datPollinators? data file).

## 15 Pollinator_Richness (a R script testing for the effect of spatial heterogeneity on pollinator richness; this script should be run with ?datPollinators? data file).

## 16 Pollinator_Diversity (a R script testing for the effect of spatial heterogeneity on pollinator Shannon diversity; this script should be run with ?datPollinators? data file).

## 17 Predator_Abundance (a R script testing for the effects of spatial heterogeneity on predator abundance; this script should be run with ?datPredators? data file).

## 18 Predator_Richness (a R script testing for the effect of spatial heterogeneity on predator richness; this script should be run with ?datPredators? data file).

## 19 Predator_Diversity (a R script testing for the effect of spatial heterogeneity on predator Shannon diversity; this script should be run with ?datPredators? data file).

## 20 Pest_Abundance (a R script testing for the effect of spatial heterogeneity on pest abundance; this script should be run with ?datPests? data file).

## 21 Pest_Richness (a R script testing for the effect of spatial heterogeneity on pest richness; this script should be run with ?datPests? data file).

## 22 Araneae_Abundance (a R script testing for the effect of spatial heterogeneity on Araneae abundance; this script should be run with ?dat? data file).

## 23 Araneae_Richness (a R script testing for the effect of spatial heterogeneity on Araneae richness; this script should be run with ?dat? data file).

## 24 Araneae_Diversity (a R script testing for the effect of spatial heterogeneity on Araneae Shannon diversity; this script should be run with ?dat? data file).

## 25 Bird_Abundance (a R script testing for the effect of spatial heterogeneity on bird abundance; this script should be run with ?dat? data file).

## 26 Bird_Richness (a R script testing for the effect of spatial heterogeneity on bird richness; this script should be run with ?dat? data file).

## 27 Bird_Diversity (a R script testing for the effect of spatial heterogeneity on bird Shannon diversity; this script should be run with ?dat? data file).

## 28 Coleoptera_Abundance (a R script testing for the effect of spatial heterogeneity on Coleoptera abundance; this script should be run with ?dat? data file).

## 29 Coleoptera_Richness (a R script testing for the effect of spatial heterogeneity on Coleoptera richness; this script should be run with ?dat? data file).

## 30 Coleoptera_Diversity (a R script testing for the effect of spatial heterogeneity on Coleoptera Shannon diversity; this script should be run with ?dat? data file).

## 31 Diptera_Abundance (a R script testing for the effect of spatial heterogeneity on Diptera abundance; this script should be run with ?dat? data file).

## 32 Diptera_Richness (a R script testing for the effect of spatial heterogeneity on Diptera richness; this script should be run with ?dat? data file).

## 33 Diptera_Diversity (a R script testing for the effect of spatial heterogeneity on Diptera Shannon diversity; this script should be run with ?dat? data file).

## 34 Hymenoptera_Abundance (a R script testing for the effect of spatial heterogeneity on Hymenoptera abundance; this script should be run with ?dat? data file).

## 35 Hymenoptera_Richness (a R script testing for the effect of spatial heterogeneity on Hymenoptera richness; this script should be run with ?dat? data file).

## 36 Hymenoptera_Diversity (a R script testing for the effect of spatial heterogeneity on Hymenoptera Shannon diversity; this script should be run with ?dat? data file).

## 37 Lepidoptera_Abundance (a R script testing for the effect of spatial heterogeneity on Lepidoptera abundance; this script should be run with ?dat? data file).

## 38 Lepidoptera_Richness (a R script testing for the effect of spatial heterogeneity on Lepidoptera richness; this script should be run with ?dat? data file).

## 39 Lepidoptera_Diversity (a R script testing for the effect of spatial heterogeneity on Lepidoptera Shannon diversity; this script should be run with ?dat? data file).

## 40 Animal_Abundance (a R script testing for the effect of spatial heterogeneity on animal [i.e. invertebrate and vertebrate without pests] abundance; this script should be run with ?dat? data file).

## 41 Animal_Richness (a R script testing for the effect of spatial heterogeneity on animal [i.e. invertebrate and vertebrate without pests] richness; this script should be run with ?dat? data file).

## 42 Animal_Diversity (a R script testing for the effects of spatial heterogeneity metrics on animal [i.e. invertebrate and vertebrate without pests] Shannon diversity; this script should be run with ?dat? data file).

## 43 Tropical_vs_Temperate_Farmlands (a R script testing for the variabilities of the effect of spatial heterogeneity on animals [i.e. invertebrate and vertebrate without pests] in tropical an temperate farmlands; this script should be run with ?dat? data file).

## 44 Annual_vs_Perennial_Crops (a R script testing for the variabilities of the effect of spatial heterogeneity on animals [i.e. invertebrate and vertebrate without pests] in annual and perennial farmlands; this script should be run with ?datCroppingSystems? data file).

## 45 Invertebrate_Scale (a R script testing for the variabilities of the effect of spatial heterogeneity on invertebrate in local-level and landscape-level scales; this script should be run with ?dat? data file). 

## 46 Vertebrate_Scale (a R script testing for the variabilities of the effect of spatial heterogeneity on vertebrate in local-level and landscape-level scales; this script should be run with ?dat? data file).

## 47 Animal_Scale (a R script testing for the variabilities of the effect of spatial heterogeneity on animals [i.e. invertebrate and vertebrate without pests] in local-level and landscape-level scales; this script should be run with ?dat? data file). 

## 48 Pollinator_Scale (a R script testing for the variabilities of the effect of spatial heterogeneity on pollinators in local-level and landscape-level scales; this script should be run with ?datPollinators? data file).

## 49 Predator_Scale (a R script testing for the variabilities of the effect of spatial heterogeneity on predators in local-level and landscape-level scales; this script should be run with ?datPredators? data file).

####################################################################################################################################################################################################################################