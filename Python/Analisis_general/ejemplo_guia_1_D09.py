# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 16:57:25 2021

@author: yolib

Ejemplo del análisis de una guía en concreto para ver si comparte los alelos 
más comunes y el efecto de la normalización. 
"""


#%%

import os
import pandas as pd
import numpy as np
import seaborn as sns 
import matplotlib.pyplot as plt
#%% 

os.chdir('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos')
guides_data = pd.read_csv("guides_sequence_2.csv")

os.chdir('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos\\results')
guides_files = os.listdir()

data = pd.read_csv(guides_files[216])
data_sorted = data.sort_values('RL384-00020_H18', ascending = False)
#data_10_sorted = data_sorted.iloc[:10, :len(data.columns)-1]
data_5_sorted = data_sorted.iloc[:5, :]

freq_table = pd.melt(data_5_sorted, id_vars="Alleles", var_name="Sample", value_name="Frequency")

#%%
import seaborn as sns
sns.set_theme(style="whitegrid")

g = sns.catplot(
    data=freq_table, kind="bar",
    x="Alleles", y="Frequency", hue="Sample",
    ci="sd", palette="dark", alpha=.6, height=6
)

plt.title("Frequency of the first 4 alleles (guide 1_D09)", fontsize=15)

#%% Ahora, normalizados con respecto al PC

df_no_variant = freq_table[freq_table["Alleles"] == "no variant"]
df_no_variant["Freq_corte"] = list(100 - freq_table[freq_table["Alleles"] == "no variant"]["Frequency"])


freq_table_normalized = pd.merge(freq_table, df_no_variant.loc[:, ["Freq_corte", "Sample"]], on='Sample')
freq_table_normalized["Freq_normalized"] = 100*(freq_table_normalized["Frequency"]/freq_table_normalized["Freq_corte"])
freq_table_normalized = freq_table_normalized.drop(freq_table_normalized[freq_table_normalized['Alleles']=="no variant"].index)

sns.set_theme(style="whitegrid")

g = sns.catplot(
    data=freq_table_normalized, kind="bar",
    x="Alleles", y="Freq_normalized", hue="Sample",
    ci="sd", palette="dark", alpha=.6, height=6
)

plt.ylabel("Normalized frequency", size = 12)
plt.title("Frequency of the firts 4 alleles normalized with respect to the cut-off percentage \n (guide 1_D09)",
          fontsize = 15)