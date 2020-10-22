# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 15:55:46 2020

@author: Benjamin
"""

import pandas as pd

# Test dataframe for TL inversion

b = np.arange(0,10)
A = np.array([b*0.13]*3).transpose()
A = np.r_[b,b]

df = pd.DataFrame(b,columns=['obs0'])
df.append(b,columns=['obs1'])