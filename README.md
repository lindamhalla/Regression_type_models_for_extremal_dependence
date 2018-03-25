# Regression_type_models_for_extremal_dependence

In this repository, new VGAM families for various angular density families are implemented along with examples for their use requiring the R package VGAM 1.04.

Bivariate_Dirichlet_model.R
 - dirichlet2d.vgam : VGAM family function for the bivariate Dirichlet angular density.
Bivariate_Husler-Reiss_model.R
 - hr2d.vgam : VGAM family function for the bivariate Hüsler—Reiss angular density.
Bivariate_logistic_model.R
 - log2d.vgam : VGAM family function for the bivariate logistic angular density.
Trivariate_Pairwise_Beta_model.R
 - pBeta3d.vgam : VGAM family function for the trivariate pairwise beta angular density.

The R packages mev and BMAmevt are needed to simulate from the angular densities. In the examples of use, the VGAMs are fitted, for the sake of illustration, in a well specified framework where the responses come from the angular density.
