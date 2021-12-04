
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from statsmodels.formula.api import ols
from statsmodels.formula.api import mixedlm
from statsmodels.tools.eval_measures import rmse
import statsmodels
import statsmodels.api as sm
import numpy as np
from patsy.contrasts import Sum
from scipy.stats.distributions import chi2

dataset = pd.read_csv('data_p3 - Copy.csv')

             
#Using maximum likelihood estimation to estimate the betas --- will use mixed linear effects modeling.
mod= mixedlm('Y ~ 1+ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19', data=dataset, groups=dataset["patient"])
model=mod.fit(reml=False)
results=model.summary()
print(results)

#Backward Elimination
mod1= mixedlm('Y ~ 1+ X1+X2+X4+X5+X8+X9+X11+X12+X13+X14+X17+X19', data=dataset, groups=dataset["group"])
model1=mod1.fit(reml=False)
results1=model1.summary()
print(results1)


mod2= mixedlm('Y ~ 1+X2+X8+X9+X12+X13+X14', data=dataset, groups=dataset["group"])
model2=mod2.fit(reml=False)
results2=model2.summary()
print(results2)


mod3= mixedlm('Y ~ 1+X2+X9+X13', data=dataset, groups=dataset["group"])
model3=mod3.fit(reml=False)
results3=model3.summary()
print(results3)
