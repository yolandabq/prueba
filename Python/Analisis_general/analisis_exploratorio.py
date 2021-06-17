# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 12:28:47 2021

@author: yolib

Análisis exploratorio. En este script se realizan varios análisis para explorar características generales de los datos. 

* Distribución de la frecuencia del alelo más común
* Eficiencia de corte 
* Distribución de la frecuencia de las deleciones
* Distribución del número de alelos obtenidos tras utilizar cada guía 
* % GC de las guías 
* Suma acumulativa de los alelos 


"""

#%%
import os
import pandas as pd
import numpy as np
import seaborn as sns 
import matplotlib.pyplot as plt
import glob
os.chdir('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos\\python_scripts')
import utils_functions

#%%
os.chdir('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos')
guides_data = pd.read_csv("guides_sequence_2.csv")

#%% Frecuencia del alelo más común 

most_common_freq = {}
most_common_allele = {}
frequ_1_1I_most_common = {}

for carpeta in ["results"]:
    
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
            df_results = df_results[df_results["Alleles"] != "no variant"]
        
        else: 
            porc_corte = 100
            
        norm = 1/porc_corte
        
        df_results["Average"] = df_results.loc[:, list(df_results.columns)[:-1]].mean(axis=1)*norm
        
        df_results = df_results.sort_values(by=list(df_results.columns)[:-1], ascending = False)

        most_common_freq[g_id] = df_results.iloc[0, -1]
        most_common_allele[g_id] = df_results.iloc[0, -2]
        
        if df_results.iloc[0, -2] == "1:1I": 
            frequ_1_1I_most_common[g_id] = df_results.iloc[0, -1]

#%% # Frecuencia del alelo más común y la distribución de 1:1I cuando es el alelo más común

most_common_freq_values = most_common_freq.values()

fig, ax = plt.subplots(figsize=(10, 8))
sns.set_style("whitegrid")
sns.histplot(most_common_freq_values, stat = "probability", ax = ax, legend = False)
plt.title("Distribution of the frequency of the most common allele", fontsize = 15)
plt.xlabel("Frequency of the most common allele", size = 12)
plt.show()

most_common_freq_values = frequ_1_1I_most_common.values()

fig, ax = plt.subplots(figsize=(10, 8))
sns.set_style("whitegrid")
sns.histplot(most_common_freq_values, stat = "probability", ax = ax, legend = False)
plt.title("Distribution of the frequency of 1:1I when it is the most common allele", fontsize = 15)
plt.xlabel("Frequency of 1:1I as the most common allele", size = 12)

#%% Porcentaje de corte (cutting efficiency)

porc_corte_dict = {}
porc_deletion_dict = {}

for carpeta in ["results"]:
    
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
        
        else: 
            porc_corte = 100
            
        norm = 1/porc_corte
        alleles_indel = [x for x in list(df_results["Alleles"]) if "no variant" not in x and "SNV" not in x and "," not in x]
            
        alleles_deletion = [x for x in alleles_indel if x[-1] == "D"]
        deletions_index = list(df_results.loc[df_results['Alleles'].isin(alleles_deletion)].index)
        
        for row in deletions_index:
            #if np.mean(list(df_results.loc[row, list(df_results.columns)[:-1]])) >= 0.5:
            porc_deletion = porc_deletion + norm*np.mean(list(df_results.loc[row, list(df_results.columns)[:-1]]))
       
        porc_deletion_dict[g_id] = porc_deletion
        porc_corte_dict[g_id] = porc_corte
        
        
#%% Representamos la distribución del cutting efficiency

porc_corte_values = list(porc_corte_dict.values())

fig, ax = plt.subplots(figsize=(8, 6))
sns.set_style("white")
sns.histplot(porc_corte_values, binwidth=5, stat = "probability", ax = ax, legend = False)
plt.title("Cutting efficiency distribution", fontsize = 15)
plt.xlabel("Cutting efficiency", size = 12)
text = "Mean = " + str(round(np.mean(porc_corte_values),2))
plt.text(80,0.075,text,fontsize = 13)

#%% Representamos la distribución de la frecuencia de deleciones

porc_deletion_values = list(porc_deletion_dict.values())

fig, ax = plt.subplots(figsize=(10, 8))
sns.set_style("whitegrid")
sns.boxplot(list(porc_deletion_values), ax = ax,  orient="v")
plt.title("Deletion frequency distribution", fontsize = 15)
plt.xlabel("Deletion frequency", size = 12)

np.mean(porc_deletion_values)

#%%

porc_deletion_values = list(porc_deletion_dict.values())

fig, ax = plt.subplots(figsize=(10, 8))
sns.set_style("whitegrid")
#plt.style.use('classic')
sns.histplot(porc_deletion_values,  stat = "probability", ax = ax, legend = False, color= "#1f77b4")
plt.title("Deletion frequency distribution", fontsize = 15)
plt.xlabel("Deletion frequency", size = 12)

arr_porc_deletion_values = np.array(porc_deletion_values)


print("Guías con menos de 0.6: " + str(len(arr_porc_deletion_values[arr_porc_deletion_values < 0.6])*100/len(arr_porc_deletion_values)))
print("Guías con menos de 0.7: " + str(len(arr_porc_deletion_values[arr_porc_deletion_values < 0.7])*100/len(arr_porc_deletion_values)))
print("Guías con menos de 0.9: " + str(len(arr_porc_deletion_values[arr_porc_deletion_values < 0.9])*100/len(arr_porc_deletion_values)))

print("Guías con más de 0.7: " + str(len(arr_porc_deletion_values[arr_porc_deletion_values > 0.7])*100/len(arr_porc_deletion_values)))


#calculate interquartile range 
q3, q1 = np.percentile(arr_porc_deletion_values, [75 ,25])
iqr = q3 - q1

#display interquartile range 
iqr




#%% Número de alelos obtenidos tras usar cada guía (filtramos alelos cuya freq normalizada > 1/1000)

n_alleles_dict = {}

for carpeta in ["results"]:
    
    os.chdir('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos\\' + carpeta)
    guides_files = os.listdir()
    guides_id = list(set([int(x.split("_")[0]) for x in guides_files]))

    for g_id in guides_id:
        
        #print(g_id)
        #print(glob.glob(str(g_id) + "_*"))
        files = glob.glob(str(g_id) + "_*")
       
        
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
        
        else: 
            porc_corte = 100
            
        norm = 1/porc_corte    
        df_results["Average"] = df_results.loc[:, list(df_results.columns)[:-1]].mean(axis=1)*norm
        df_results = df_results[df_results["Average"] >= 1/1000]
        
        alleles = df_results["Alleles"]
        
        if "no variant" in alleles:
            alleles.remove("no variant")
        
        n_alleles_dict[g_id] = len(alleles)
        
#%%

n_alleles_values = list(n_alleles_dict.values())
#print(np.mean(n_alleles_values))

fig, ax = plt.subplots(figsize=(8, 6))
sns.set_style("white")
sns.histplot(n_alleles_values,  stat = "probability", ax = ax, legend = False)
plt.title("Distribution of the number of outcomes > 1/1000", fontsize = 15)
plt.xlabel("Number of outcomes", size = 12)
text = "Mean = " + str(round(np.mean(n_alleles_values),2))
plt.text(185,0.15,text,fontsize = 13)


#%% Exploramos también el porcentaje de GC de las guías para ver si se relaciona con el % de corte. 

guides_id = list(set(list(guides_data["Guide_ID"])))
GC = {}
strand = {}

for k in guides_id: 
    sequence = (list(guides_data[guides_data["Guide_ID"] == k]["Sequence"])[0])[27:33].upper()
    GC[k] = (sequence.count("G") + sequence.count("C"))/len(sequence)
    strand[k] = (list(guides_data[guides_data["Guide_ID"] == k]["Strand"])[0])
    

df_data = pd.DataFrame(index = guides_id, 
                                 columns = ["GC", "Porc_corte"])
            
for k in guides_id: 
    df_data.loc[k, "GC"] = float(GC[k])  
    df_data.loc[k, "Porc_corte"] = float(porc_corte_dict[k]  )
    
plt.scatter(x = list(GC.values()), y = list(porc_corte_dict.values()))

df_data = pd.DataFrame(index = guides_id, 
                                 columns = ["Strand", "Porc_corte"])
            
for k in guides_id: 
    df_data.loc[k, "Strand"] = (strand[k])  
    df_data.loc[k, "Porc_corte"] = float(porc_corte_dict[k]  )

#sns.swarmplot(data = df_data, x = "Strand", y = "Porc_corte")
#%%

data = pd.DataFrame(columns = ["Cut Efficiency", "Nº Alleles", "Deletion frequency", "GC"], index = list(n_alleles_dict.keys()))

for k in n_alleles_dict.keys(): 
    data.loc[k, :] = [porc_corte_dict[k], n_alleles_dict[k], porc_deletion_dict[k], GC[k]]

#%%
fig, ax = plt.subplots(figsize=(10, 8))
sns.scatterplot(y = "Nº Alleles", x = "Cut Efficiency", data = data, ax = ax)

#%% Correlación de Pearson
# https://machinelearningmastery.com/how-to-use-correlation-to-understand-the-relationship-between-variables/

from numpy.random import randn
from numpy.random import seed
from scipy.stats import pearsonr
import scipy.stats 
from matplotlib.font_manager import FontProperties

y = list(data["Nº Alleles"])
x = list(data["Cut Efficiency"])

slope, intercept, r, p, stderr = scipy.stats.linregress(x, y)
line = f'Regression line: y={intercept:.2f}+{slope:.2f}x, r={r:.2f}'

fontP = FontProperties()
fontP.set_size('medium')
# these are matplotlib.patch.Patch properties
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

fig, ax = plt.subplots(figsize=(15, 10))
ax.plot(x, y, linewidth=0, marker='o', label='Data points', 
        markersize=5, c='darkmagenta')
ax.plot(x, intercept + slope * np.array(x), 
        label=line, c = "black", linestyle=":")
ax.set_xlabel("Cutting efficiency", size = 10)
ax.set_ylabel("Nº Alleles", size = 10)
plt.title("Correlation between cutting efficiency and Nº alleles", size = 15)
ax.legend(facecolor='white', bbox_to_anchor=(1.05, 1), loc='upper left', prop=fontP)
text = str('R = ') + str(round(r,3))
ax.text(1.2, 0.9, text, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', bbox=props)

plt.show()

#%% Suma acumulativa de los alelos 

import operator
import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots(figsize=(12, 4))

n_al = 31
row = 0


for carpeta in ["results"]:
    
    os.chdir('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos\\' + carpeta)
    guides_files = os.listdir()
    guides_id = list(set([int(x.split("_")[0]) for x in guides_files]))

    n_samples = len(guides_id)
    freq = np.zeros((n_samples,n_al ))
    
    for g_id in guides_id[0:n_samples]:
        
        #print(g_id)
        #print(glob.glob(str(g_id) + "_*"))
        files = glob.glob(str(g_id) + "_*")
       
        
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
            df_results = df_results[df_results["Alleles"] != "no variant"]
            #porc_corte_list.append(porc_corte)
        
        else: 
            porc_corte = 100
            
        norm = 1/porc_corte    
        df_results["Average"] = df_results.loc[:, list(df_results.columns)[:-1]].mean(axis=1)*norm
        df_results = df_results.sort_values(by = "Average", ascending = False)
        
        y = [0]
        suma = 0
        
        
        for allele in list(df_results.index)[0:n_al-1]:  
            
            if df_results.loc[allele, "Alleles"] != "no variant": 
                suma = suma + df_results.loc[allele, "Average"]
                y.append(suma)
        
        freq[row] = y 
        row = row + 1
        plt.plot(np.arange(0, len(y)), y, color = 'lightgrey')
    
    plt.plot(np.arange(0, len(y)), np.mean(freq, axis=0), color = 'royalblue') # la línea azul es la media
    plt.ylim((0,1))
    #plt.title("Cumulative sum of allele frequencies", fontsize = 15)
    
    plt.ylabel("Cumulative frequency of alleles", size = 12 )
    plt.xlabel("Number of alleles", size = 12 )
    
    
    style = dict(size=10, color='gray')
    ax.text(n_al-1.5, 0.8, str(round(np.mean(freq, axis=0)[-1],2)), transform=ax.transData)
    ax.text(24.5, 0.775, str(round(np.mean(freq, axis=0)[25],2)), transform=ax.transData)
    ax.text(19.5, 0.75, str(round(np.mean(freq, axis=0)[20],2)), transform=ax.transData)
    
    #ax.text('30', 0.8, "hs", **style)
    
    plt.show()
