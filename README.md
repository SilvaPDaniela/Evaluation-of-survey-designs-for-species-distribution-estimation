# Evaluation-of-survey-designs-for-species-distribution-estimation

R Code corresponding to the manuscript: Adapting the sampling design of research surveys to improve the biomass estimation of non-target species - the case study of Raja clavata

The code corresponding to the modelling process is not included, which implementation was based on Krainski et al. (2019).

Script_1_Selection_procedure: Having as input a predicted thinned grid, we start defining the constraints, the minimum distance and the measures for locations selection. The points.selection.function() is a loop function for select locations guaranting designs constraints using the measures previous defined. The result is a set of alternative designs, each containing the coordinates of the selected locations for each measure.

Script_2_Assessment_of_selection_measures: Having as input a database of predicted and observed abundance values for a selected benchmark survey, compute the assessment criteria for each alternative design and species through a loop function.








References:

E. T. Krainski et al. (2019). Advanced Spatial Modeling with Stochastic Partial Differential Equations Using R and INLA. Chapman & Hall/CRC Press. Boca Raton, FL.

Lindgren, F., Rue, H., and Lindström, J. (2011). An explicit link between gaussian fields and gaussian markov random fields: the stochastic partial differential equation approach. Statistical Methodology, 73(4):423–498.

Zuur, A. F., Ieno, E. N., and Saveliev, A. A. (2017). Beginner’s Guide to Spatial, Temporal, and Spatial-Temporal Ecological Data Analysis with R-INLA, volume I: Using GLM and GLMM, chapter 20: GAM with correlation and zero-inflation in R-INLA using owl data. Highland Statistics Ltd.
