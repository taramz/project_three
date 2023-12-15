# Project Three
 
Empirical investigations often span two distinct populations—the study population and the target population. The alignment of these populations is challenging, especially when constructing predictive models intended for application in the target population. Discrepancies between them can compromise model efficacy, emphasizing the need for transportability analysis.

Datasets from the Framingham Study and NHANES are leveraged to explore cardiovascular disease (CVD) prediction across diverse populations. The Framingham Study, initiated in 1948, is a longitudinal cohort investigating CVD epidemiology and risk factors. NHANES evaluates the health and nutritional status of the U.S. population, integrating self-reported and comprehensive examination data.

Population disparities between Framingham and NHANES underscore the imprudence of extrapolating model performance across datasets. Transportability analysis becomes imperative to understand model generalizability. The Brier score, a metric for probabilistic predictions, is employed, but adjustments are made to evaluate transportability. Multiple imputation addresses NHANES' missing data, and gender-specific average Brier scores are computed.

The ADEMP framework guides simulation considerations, encompassing aim, data generating mechanism, estimand, methods, and performance measures. Simulation involves diverse parameter variations for gender-specific data generation, aiming to replicate NHANES demographics.

The evaluation of simulated datasets using bias and mean squared error reveals nuanced differences between genders. Optimal covariance values for minimizing bias and MSE differ for men and women. The Brier score comparison and summary statistics show congruence with NHANES for women but slight discrepancies for men, potentially attributed to data variability.

The data generation process incorporates missing variables, potentially contributing to disparities in Brier scores. Computational constraints limit simulations, raising the possibility of unstable estimates. These considerations highlight areas for refining the data generation process and underscore the importance of meticulous transportability analysis in model development and application across diverse populations.

## Directory:
- project3.R: the code outlining all the data pre-processing performed, as well as the initial EDA that did not make it into the final analysis.
- simulation_code.R: the code where all the simulations were performed to determine the best parameters.
- sim_results_men_1.rds: the simulation results for men.
- sim_results_women_1.rds: the simulation results for women.
- diabetes_mod.rds: the model used to determine the probability of an individual having diabetes.
- Project_Three.Rmd: the .Rmd file including all the code and text from the analysis.
- Project_Three.pdf: the exported .pdf file.
- references.bib: the references used in this analysis, in the BibTex format.

## Acknowledgements:

This analysis has been performed with help from Dr. Jon Steingrimsson from the School of Public Health at Brown University.
