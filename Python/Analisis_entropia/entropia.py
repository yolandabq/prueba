# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 17:16:23 2021

@author: yolib

Script para analizar, mediante el parámetro de Entropía de Shannon, la diversidad del perfil mutacional que genera cada guía.

"""

#%%

import os
import pandas as pd
import numpy as np
import seaborn as sns 
import matplotlib.pyplot as plt
import glob
from collections import Counter
import math
os.chdir('C:\\Users\\yolib\\OneDrive\\Documentos\\TFM\\Documento TFM\\Codigo\\Python\\Funciones')
import utils_functions
import entropy_functions as enp

#%% Cargo las secuencias guía

os.chdir('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos')
guides_data_no_replic = pd.read_csv("guides_sequence_2_no_replic.csv")
guides_replic_data = pd.read_csv("guides_sequence_2.csv")

frames = [guides_replic_data, guides_data_no_replic]
guides_data = pd.concat(frames)


#%% Calculo la entropía de cada guía

#entropy_list = []
entropy_dict = {}
porc_corte_dict = {}
most_common_dict = {}

for carpeta in ["results", "results_no_replic"]:
    
    os.chdir('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos\\' + carpeta)
    guides_files = os.listdir()
    guides_id = list(set([int(x.split("_")[0]) for x in guides_files]))

    for g_id in guides_id:
        
        #print(g_id)
        #print(glob.glob(str(g_id) + "_*"))
        files = glob.glob(str(g_id) + "_*")
        porc_deletion = 0
        
        if len(files) >= 2:
            
            df_list_results = []
            df_list_insertion = []

            for file in files:  # si hay más de un archivo de la misma guía: 
            
                df_list_results.append(pd.read_csv(file))
                
            df_results = utils_functions.merge_df_results_alleles(df_list_results)
            df_results = df_results.fillna(0)
            
        else: 
            
            df_results = pd.read_csv(files[0])
            
        if "no variant" in list(df_results["Alleles"]):
            porc_corte = (100-np.mean(list(df_results[df_results["Alleles"] == "no variant"].iloc[0, :-1])))
            #porc_corte_list.append(porc_corte)
            entropy_dict[str(g_id)] = enp.entropy(df_results)
            #entropy_dict[str(g_id)] = entropy_25(df_results)
            
            df_results = df_results[df_results["Alleles"] != "no variant"] 
        
        else: 
            porc_corte = 100
            entropy_dict[str(g_id)] = enp.entropy(df_results)
            #entropy_dict[str(g_id)] = entropy_25(df_results)
            
        norm = 1/porc_corte
        
        df_results = df_results.sort_values(by = list(df_results.columns)[:-1], ascending = False)
        most_common = np.mean(df_results.iloc[0, :-1])
        most_common_dict[str(g_id)] = most_common * norm
        porc_corte_dict[g_id] = porc_corte
        
        
#%%

print("La entropía media es: " + str(np.mean(list(entropy_dict.values()))) + " +- " + str(np.std(list(entropy_dict.values()))))
# La entropía media es: 5.684740579965953 +- 1.0715673630409448

values = [] #in same order as traversing keys
keys = [] #also needed to preserve order
for key in entropy_dict.keys():
  keys.append(key)
  values.append(entropy_dict[key])

fig, ax = plt.subplots(figsize=(10, 8))

sns.histplot(values, stat = "density", ax = ax)
plt.xlabel("Entropy value", size = 12)
plt.title("Distribution of entropy values", fontsize = 15)
#plt.title("Distribution of entropy values \n (using the 25 most common outcomes)", fontsize = 15)

list_of_key = list(entropy_dict.keys())
list_of_value = list(entropy_dict.values())

"""
position = list_of_value.index(0.15272878019499037)
print(list_of_key[position])
files = glob.glob(str(list_of_key[position]) + "_*")
result_df = pd.read_csv('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos\\results_no_replic\\' + files[0])
"""
#%%
fig, ax = plt.subplots(figsize=(10, 8))
ax.scatter(y = list(entropy_dict.values()),x = list(most_common_dict.values()))

#%% Correlación de Pearson entre la entropía y la eficiencia de corte

from numpy.random import randn
from numpy.random import seed
from scipy.stats import pearsonr
import scipy.stats 
from matplotlib.font_manager import FontProperties

y = list(entropy_dict.values())
x = list(porc_corte_dict.values())
sns.set_style("white")
slope, intercept, r, p, stderr = scipy.stats.linregress(x, y)
line = f'Regression line: y={intercept:.2f}+{slope:.2f}x, r={r:.2f}'

fontP = FontProperties()
fontP.set_size('medium')
# these are matplotlib.patch.Patch properties
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

fig, ax = plt.subplots(figsize=(12, 8))
ax.plot(x, y, linewidth=0, marker='o', label='Data points', 
        markersize=5, c='cornflowerblue')
ax.plot(x, intercept + slope * np.array(x), 
        label=line, c = "black", linestyle=":")
ax.set_xlabel("Efficiency", size = 13)
ax.set_ylabel("Entropy", size = 13)
plt.title("Correlation between the entropy of the mutational profile \n and the cutting efficiency", size = 15)
ax.legend(facecolor='white', bbox_to_anchor=(1.05, 1), loc='upper left', prop=fontP)
text = str('R = ') + str(round(r,3))
ax.text(1.2, 0.85, text, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', bbox=props)

plt.show()

#%% Correlación de Pearson entre la entropía y la freq del alelo más común 

from numpy.random import randn
from numpy.random import seed
from scipy.stats import pearsonr
import scipy.stats 
from matplotlib.font_manager import FontProperties

y = list(entropy_dict.values())
x = list(most_common_dict.values())
sns.set_style("white")
slope, intercept, r, p, stderr = scipy.stats.linregress(x, y)
line = f'Regression line: y={intercept:.2f}+{slope:.2f}x, r={r:.2f}'

fontP = FontProperties()
fontP.set_size('medium')
# these are matplotlib.patch.Patch properties
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

fig, ax = plt.subplots(figsize=(12, 8))
ax.plot(x, y, linewidth=0, marker='o', label='Data points', 
        markersize=5, c='cornflowerblue')
ax.plot(x, intercept + slope * np.array(x), 
        label=line, c = "black", linestyle=":")
ax.set_xlabel("Frequency of the most common outcome", size = 13)
ax.set_ylabel("Entropy", size = 13)
plt.title("Correlation between the entropy of the mutational profile \n and the frequency of the most frequent outcome", size = 15)
ax.legend(facecolor='white', bbox_to_anchor=(1.05, 1), loc='upper left', prop=fontP)
text = str('R = ') + str(round(r,3))
ax.text(1.2, 0.85, text, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', bbox=props)

plt.show()



#%% Puedo separar las guías según su valor de entropía. 

list_of_key = list(entropy_dict.keys())
list_of_value = list(entropy_dict.values())

index_menor_4 = [list_of_value.index(x) for x in list_of_value if x <= 4]
key_menor_4 = [list_of_key[x] for x in index_menor_4]

index_entre_4_6 = [list_of_value.index(x) for x in list_of_value if x > 4 and x <= 6]
key_entre_4_6 = [list_of_key[x] for x in index_entre_4_6]

index_mayor_6 = [list_of_value.index(x) for x in list_of_value if x > 6]
key_mayor_6 = [list_of_key[x] for x in index_mayor_6]

#%%

#seq_menor_4 = [list(guides_data[guides_data["Guide_ID"] == int(guide)]["Sequence"])[0][25:34].upper() for guide in key_menor_4]
#seq_entre_4_6 = [list(guides_data[guides_data["Guide_ID"] == int(guide)]["Sequence"])[0][25:34].upper() for guide in key_entre_4_6]
#seq_mayor_6 = [list(guides_data[guides_data["Guide_ID"] == int(guide)]["Sequence"])[0][25:34].upper() for guide in key_mayor_6]



seq_menor_4 = [list(guides_data[guides_data["Guide_ID"] == int(guide)]["Sequence"])[0].upper() for guide in key_menor_4]
seq_entre_4_6 = [list(guides_data[guides_data["Guide_ID"] == int(guide)]["Sequence"])[0].upper() for guide in key_entre_4_6]
seq_mayor_6 = [list(guides_data[guides_data["Guide_ID"] == int(guide)]["Sequence"])[0].upper() for guide in key_mayor_6]

#all_guides = seq_menor_4 + seq_entre_4_6 + seq_mayor_6





