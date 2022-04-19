# Library of Existing Epigenetic Clocks

A simple python library with a collection of a handful of existing epigenetic clocks which can be easily imported and applied to properly formatted input.


## Usage
To import this library, download this directory and use:

```
import sys
sys.path.append('/path/to/dir')
from clocks_library import PetkovichClock 

# import desired clock from following available clocks:
#
# MOUSE CLOCKS
# MeerClock
# PetkovichClock
# ThompsonClock
# WangClock
#
# HUMAN CLOCKS
# HorvathClock
# HannumClock
# PhenoAgeClock
# DunedinClock
# PedBEClock
# ZhangClock
```

The `get_age()` method for each class accepts a dataframe whose indices are valid CpG identifiers and whose values are real numbers between 0 and 1 inclusive (or `nan`). 

To predict age using a given clock, do:

```
# import a dataframe containing methylation data with CpG identifiers as the index
# and values between 0 and 1 inclusive (or nan)
df = pd.read_csv('input_data.csv', index_col=0)
pvc = PetkovichClock()
predicted_ages = pvc.get_age(df)
```

## TODOS:
- Potential additions:
  - Stubbs clock
  - Kerepesi clock
  - AltumAge clock
  - PCAge clock
