import sklearn as sklearn
from sklearn import datasets, linear_model
import pandas as pd
import numpy as np

covariates = pd.read_csv(r"C:\Users\mahmo\Downloads\MMIT\RPII\CYP2D6\Cov_Correction\TAD_Age_CYP2D6Short2_All.txt", sep='\t', header = None)
metabolic_ratio = pd.read_csv(r"C:\Users\mahmo\Downloads\MMIT\RPII\CYP2D6\Cov_Correction\LnMR_CYP2D6_All.txt", sep='\t', header = None)


# First initialize an instance of the linear model class
correction_model = sklearn.linear_model.LinearRegression(fit_intercept=True)

# Fit the model
correction_model.fit(covariates, metabolic_ratio)

# Now calculate the outcomes
predicted_effects = correction_model.predict(covariates)

# And subtract those from the expression matrix
corrected_metabolic_ratio = metabolic_ratio - predicted_effects

np.savetxt(r"C:\Users\mahmo\Downloads\MMIT\RPII\CYP2D6\Cov_Correction\Corrected_LnMR_CYP2D6_ALL_Age_TAD_CYP2D6Short2.txt", corrected_metabolic_ratio, delimiter='\t')   # X is an array
