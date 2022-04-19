import os
import numpy as np
import pandas as pd

DAYS_IN_MONTH = 365/12

# current path
CURRENT_PATH               =    os.path.dirname(__file__)

# mouse clock paths
MEER_CLOCK_PATH            =    os.path.join(CURRENT_PATH, 'Clocks/Mouse/MeerClock.tsv')
PETKOVICH_CLOCK_PATH       =    os.path.join(CURRENT_PATH, 'Clocks/Mouse/PetkovichClock.tsv')
THOMPSON_CLOCK_PATH        =    os.path.join(CURRENT_PATH, 'Clocks/Mouse/ThompsonClock-all.tsv')
WANG_CLOCK_PATH            =    os.path.join(CURRENT_PATH, 'Clocks/Mouse/WangLiver.tsv')
      
# human clock paths      
HORVATH_CLOCK_PATH         =    os.path.join(CURRENT_PATH, 'Clocks/Human/Horvath2013.tsv')
HANNUM_CLOCK_PATH          =    os.path.join(CURRENT_PATH, 'Clocks/Human/Hannum2013.tsv')
PHENOAGE_CLOCK_PATH        =    os.path.join(CURRENT_PATH, 'Clocks/Human/PhenoAge.tsv')
DUNEDIN_CLOCK_PATH         =    os.path.join(CURRENT_PATH, 'Clocks/Human/Dunedin.tsv')
PEDBE_CLOCK_PATH           =    os.path.join(CURRENT_PATH, 'Clocks/Human/PedBE.tsv')
ZHANG_CLOCK_PATH           =    os.path.join(CURRENT_PATH, 'Clocks/Human/Zhang2019.tsv')

# mouse clock constants   
MEER_CLOCK_INTERCEPT       =    234.64
PETKOVICH_CLOCK_INTERCEPT  =    0.
THOMPSON_CLOCK_INTERCEPT   =    30.3172
WANG_CLOCK_INTERCEPT       =    5.827926399

# human clock constants
HORVATH_CLOCK_INTERCEPT    =    0.695507258
HANNUM_CLOCK_INTERCEPT     =    0.
PHENOAGE_CLOCK_INTERCEPT   =    60.664
DUNEDIN_CLOCK_INTERCEPT    =    -0.0692980470852994
PEDBE_CLOCK_INTERCEPT      =    0.
ZHANG_CLOCK_INTERCEPT      =    65.792950

# CLOCK BUILDING BLOCKS ###########################################################

# A simple linear clock
# (1) selects relevant CpG sites
# (2) imputes missing values
# (3) takes the vector dot product of methylation levels with model slopes
# (4) adds model intercept for predicted age
class LinearClock():
    def __init__(self, sites, slopes, intercept):
        self.sites = sites # should be a numpy array of strings
        self.slopes = slopes # should be a numpy array of floats
        self.intercept = intercept # should be a float
        
    def _impute(self, betas, method):
        return betas.fillna(method(betas))
    
    # Accepts a dataframe with a very specific form
    # Indices contain mouse chromosome names
    # Values are in [0, 1] or are nan
    def get_age(self, 
                betas, 
                impute_method=lambda x: x.mean(), 
                scale_betas=lambda x: x, 
                scale_age=lambda x: x):
    
        betas = betas.reindex(self.sites) # Isolate relevant sites
        assert betas.shape[0] == len(self.sites) # Assert that we have all the necessary sites
        betas = scale_betas(betas) # Scale from [0, 1] domain to another domain
        betas = self._impute(betas, impute_method) # Impute nan values
        age = (betas * self.slopes.reshape(-1,1)).sum(axis = 0) + self.intercept # Perform linear computation
        return pd.DataFrame(scale_age(age), columns=['PredictedAge']) # Scale to appropriate unit
    
    def get_num_nan(self, betas):
        return pd.DataFrame(betas.reindex(self.sites).isna().sum(), columns=['Number of NaNs'])

class LinearClockFromPath(LinearClock):
    def __init__(self, path, intercept):
        _df = pd.read_csv(path, sep='\t', index_col=0)
        super().__init__(_df.index.values, _df.iloc[:,0].astype(float).values, intercept)

# A common clock pattern which uses the "trafo" function defined in Horvath2013
class TrafoClock(LinearClock):
    def __init__(self, sites, slopes, intercept):
        super().__init__(sites, slopes, intercept)
    
    def _trafo(self, arr, adult_age=20):
        return (arr > 0) * (arr * (adult_age + 1) + adult_age) + (arr <= 0) * (np.exp(arr + np.log(adult_age + 1)) - 1)

# MOUSE CLOCKS ####################################################################

# Input has to be in [0, 100]
class MeerClock(LinearClockFromPath):
    def __init__(self):
        super().__init__(MEER_CLOCK_PATH, MEER_CLOCK_INTERCEPT)
        
    # Adjust age output to months
    def get_age(self, betas):
        return super().get_age(betas,
                               scale_betas = lambda x: 100 * x,
                               scale_age = lambda x: x / DAYS_IN_MONTH)

# Input has to be in [0,1]
class ThompsonClock(LinearClockFromPath):
    def __init__(self):
        super().__init__(THOMPSON_CLOCK_PATH, THOMPSON_CLOCK_INTERCEPT)

# Input has to be in [0, 100]
# Trained on liver tissue; performs best on liver tissue
class WangClock(LinearClockFromPath):
    def __init__(self):
        super().__init__(WANG_CLOCK_PATH, WANG_CLOCK_INTERCEPT)
        
    # Adjust by DAYS_IN_MONTH and take 2 to the power of this new value
    def get_age(self, betas):
        return 2**(super().get_age(betas, 
                                   scale_betas = lambda x: 100 * x,
                                   scale_age = lambda x: x / DAYS_IN_MONTH))

# Petkovich is non-linear
class PetkovichClock(LinearClockFromPath):
    def __init__(self):
        super().__init__(PETKOVICH_CLOCK_PATH, PETKOVICH_CLOCK_INTERCEPT)
        self.a = 0.1666
        self.b = 0.4185 
        self.c = -1.712
        
    def _petkovich_nonlinearity(self, val): return ((val - self.c) / self.a) ** (1 / self.b)

    def get_age(self, betas):
        return super().get_age(betas,
                               scale_age = lambda x: self._petkovich_nonlinearity(x) / DAYS_IN_MONTH)

# HUMAN CLOCKS ####################################################################

class HorvathClock(LinearClockFromPath, TrafoClock):
    def __init__(self):
        super().__init__(HORVATH_CLOCK_PATH, HORVATH_CLOCK_INTERCEPT)
    
    def get_age(self, betas):
        return super().get_age(betas, scale_age=self._trafo)
    

class HannumClock(LinearClockFromPath):
    def __init__(self):
        super().__init__(HANNUM_CLOCK_PATH, HANNUM_CLOCK_INTERCEPT)
        
class PhenoAgeClock(LinearClockFromPath):
    def __init__(self):
        super().__init__(PHENOAGE_CLOCK_PATH, PHENOAGE_CLOCK_INTERCEPT)
        
class DunedinClock(LinearClockFromPath):
    def __init__(self):
        super().__init__(DUNEDIN_CLOCK_PATH, DUNEDIN_CLOCK_INTERCEPT)
        
class PedBEClock(LinearClockFromPath, TrafoClock):
    def __init__(self):
        super().__init__(PEDBE_CLOCK_PATH, PEDBE_CLOCK_INTERCEPT)
        
    def get_age(self, betas):
        return super().get_age(betas, scale_age=self._trafo)
        
class ZhangClock(LinearClockFromPath):
    def __init__(self):
        super().__init__(ZHANG_CLOCK_PATH, ZHANG_CLOCK_INTERCEPT)
