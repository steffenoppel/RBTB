# RBTB
Analysis of Red-billed Tropicbird tracking data in collaboration with St Helena Marine Section and Uni Barcelona. Basis for manuscript submitted to Marine Ecology Progress Series in 2018:

Foraging ecology of tropicbirds breeding in two contrasting marine environments

Ngoné Diop, Laura Zango, Annalea Beard, Cheikh Tidiane Ba, Papa Ibnou Ndiaye, Leeann Henry, Elizabeth Clingham, Steffen Oppel, Jacob González-Solís


Two different R files contain code for the following steps:

Seabird_ENVIRON_DATA_Preparation.R: This file reads in data downloaded from Movebank and Copernicus to create a regular grid over both study areas and associate all background and bird locations with the 13 different environmental variables in space and time. Contains code for spatial GIS operations, overlays, distance claculations, terrain calculations, and general data manipulation to return a flat data.frame that can be used as input for modelling [which is done in the second script]

RBTB_environmental_models_MEPS_manuscript.r: This file contains the statistical models to compare foraging and commuting locations, and the RandomForest model to compare background and foraging locations, assess variable importance and spatial autocorrelation, predictive ability, and create plots of variable importance and variable relationships.

