# Integrated modelling of labile and glycated hemoglobin with glucose for enhanced diabetes detection and short-term monitoring

These files correspond to the code developed for the project titled "Integrated Modelling of Labile and Glycated Hemoglobin with Glucose." Please note that we are currently in the process of finalizing details, and this file will be updated with more information shortly.

The files are divided into 4 folders, each used in an aspect of the project:

• Kinetic Variable Study: This folder contains the code used to select and compute values for simulations related to different patients. 
  1) The 'Obtain variables' subfolder includes three files:

    - initial_common_var_simulations.m: This file is used to calculate the best 'common kinetic variables' for all patients involved in this study.
    - kinetic_variables_study_commentated.m: This file calls the function 'kin_var_sim' to generate simulations with varying numbers of common and patient-specific kinetic variables.
    - kin_var_sim.m: This file computes our values.

  2) Oat Sensitivity: This subfolder contains files that compute and plot errors for each group of kinetic variables.

  3) The Akaike Information Criterion: This subfolder contains files that compute the Akaike information criterion value, determining the number of variables ultimately used.

• Sensitivity: This folder contains files that compute the global sensitivity of our function using the Sobol index method.

• Stability: This folder contains files that compute an approximation suggesting the stability of our equations through Floquet theory.

• R Dimensionality Reduction and Regression: This folder contains files used for U-maps and logistic regression.

For more information, please contact the first and corresponding authors at:

Gabriel F. Calvo: gabriel.fernandez@uclm.es
David G. Aragonés: david.aragones@uclm.es
José Antonio Romero: joseantonio.romero@uclm.es
