# -*- coding: utf-8 -*-
"""
Created on Sat Mar 13 18:10:09 2021

@author: yolib


Script para comprobar si el proceso de reparación es aleatorio o es replicable. 
Para ello, usamos las guías con las que se han realizado 3 o más réplicas realizamos 2 comprobaciones: 
    1. Que la divergencia de KL sea menor entre réplicas que entre muestras que han usado diferentes gRNAs.
    2. Que el porcentaje de mutaciones frameshift (que no son múltipo de 3) coincide entre las réplicas. 
    
Voy a normalizar con respecto al % de corte (que será 100 - "No Variant"), ya que el % de corte puede variar 
por factores experimentales. Divido la frecuencia de cada alelo entre ese valor. 


"""

#%%

import os
import pandas as pd
import numpy as np
os.chdir('C:\\Users\\yolib\\OneDrive\\Documentos\\TFM\\Documento TFM\\Codigo\\Python\\Funciones')
import replicate_functions as rpl

#%% Obtengo los alelos de las guías con más de 3 réplicas

os.chdir('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos\\results')

guides_files = os.listdir()

lista_dict_min_freq = [] # diccionario de cada una de las muestras
lista_samples_min_freq = []
lista_guide_index_freq = []

for guide in guides_files: 
    data = pd.read_csv(guide)
    n_samples = data.shape[1] - 1
    #print(n_samples)
    #print(list(data.columns))
    #str(guide.split('.')[0] + )
    for sample in list(data.columns)[0:n_samples]:
        #print(guide.split('.')[0] + '_' + sample)
        sample_name_min_freq = str(guide.split('.')[0] + '_' + sample)
        sample_dict_min_freq = {}
        for k,line in data.iterrows():
            sample_dict_min_freq[line['Alleles']] = line[sample]
          
        lista_dict_min_freq.append(sample_dict_min_freq)    
        lista_samples_min_freq.append(sample_name_min_freq)
        lista_guide_index_freq.append(sample_name_min_freq.split('_')[0])

#%% Cálculo de la divergencia de KL entre todas las muestras.
# Como ya lo he calculado y tengo guardados los valores, la casilla está comentada, ya que este paso es lento. 
"""
KL_array_min_freq = np.zeros((len(lista_dict_min_freq), len(lista_dict_min_freq)))

#%%

for sample_p1 in np.arange(len(lista_dict_min_freq)):
    #print(sample_i)
   
    for sample_p2 in np.arange(len(lista_dict_min_freq)):
        KL_array_min_freq[sample_p1, sample_p2] = rpl.symmetricKL(lista_dict_min_freq[sample_p1], lista_dict_min_freq[sample_p2], 0.001)
        
    
np.savetxt('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos\\KL_normalizado.csv', KL_array_min_freq)
"""

#%% Representación del Heatmap. 
# Ordeno las muestras de forma que las réplicas estén una detrás de otra

os.chdir('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos\\')

KL_array_min_freq = np.loadtxt('KL_normalizado.csv') # archivo con la matriz de los valores de divergencia de KL. 

import seaborn as sns; sns.set_theme()
import matplotlib.pyplot as plt


#cmap = sns.diverging_palette(0, 230, 90, 60, as_cmap=True)

i = 0; j = 100
samples_names = [p.split('_')[0] for p in lista_samples_min_freq[i:j]]

KL_100_min_freq = KL_array_min_freq[i:j, i:j]
sns.set(font_scale=1)
ax = sns.heatmap(KL_100_min_freq, linewidths=0, # si no quiero que haya entre los cuadrados, pongo 0 
                 cmap = 'bone',
                 vmin=0, 
                 vmax=4, # vmax entre 2 y 4 para que se aprecie bien la diferencia de colores (si no se solapan)
                 #annot=True,
                 #annot_kws={'size':8},
                 #xticklabels=5,
                 #xticklabels=True,
                 #yticklabels=5
                 #yticklabels = samples_names
                 )


labels = [samples_names[index] for index in np.arange(0, j, 3)]

ax.set_xticks(np.arange(0, j, 3))
ax.set_xticklabels(labels, rotation=90,fontsize = 8)

ax.set_yticks(np.arange(0, j, 3))
ax.set_yticklabels(labels,rotation=0,fontsize = 8)

plt.title("Heatmap of the values of KL", fontsize=15)

plt.show()

#%% 
# Voy a guardar, por una parte, los valores pertenecientes a las réplicas, y, por otra 
# parte, los que corresponden a valores de KL entre guías random.
# Voy a hacer un dataframe que me indique que es cada valor. 0 significa que es él mismo, 
# 1 significa réplica y 2 significa random.
# No uso el nombre de la guía, si no el número de antes (del nombre del archivo)


import pandas as pd
#data = pd.DataFrame(columns=lista_samples_min_freq, index = lista_samples_min_freq)
data = pd.DataFrame(columns=lista_samples_min_freq, index = lista_samples_min_freq)


for row in np.arange(len(data.columns)): 
#for row in np.arange(10):
    #for col in np.arange(len(data.columns)):
    for col in np.arange(row):    
        if data.columns[row] == data.index[col]: 
            data.loc[[data.index[row]], data.columns[col]] = 0
            
            
        elif data.columns[row].split('_')[0] == data.index[col].split('_')[0]: 
            data.loc[[data.index[row]], data.columns[col]] = 1
            
        else: 
            data.loc[[data.index[row]], data.columns[col]] = 2
    
data = data.to_numpy()

cond_rep = (data == 1)
cond_no_rep = (data == 2)

replic = np.extract(cond_rep, KL_array_min_freq)
no_replic = np.extract(cond_no_rep, KL_array_min_freq)

np.savetxt('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos\\repl_norm.csv', replic, delimiter=',')
np.savetxt('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos\\no_replicas_norm.csv', no_replic, delimiter=',')


#%%  Represento los datos

#colors = ['#78C850', '#6890F0']
          
sns.set_style("whitegrid") 
boxplot = sns.boxplot(data= [replic, no_replic], palette = "Blues", 
                      showfliers = False, # si no quiero los outliers
                      width=0.7)
boxplot.axes.set_title("Boxplot of KL values", fontsize=15)
boxplot.set_xlabel("Samples", fontsize=12)
boxplot.set_ylabel("KL divergence", fontsize=12)

plt.xticks([0,1], ["Replicates", "Not Replicates"])

plt.show()

#%%

fig, axs = plt.subplots(2)
fig.suptitle('Histograma de los valores de KL')
#plt.ylabel('Frecuencia')

axs[0].hist(x=replic, bins = np.arange(0, 1, 0.01)) #, bins=intervalos, color='#F2AB6D', rwidth=0.85)

axs[1].hist(x=no_replic,bins = np.arange(0, max(no_replic), 0.5)) #, bins=intervalos, color='#F2AB6D', rwidth=0.85)

#plt.title('Histograma de KL')
plt.xlabel('KL')

axs[0].set(ylabel="Frecuencia en réplicas")
axs[1].set(ylabel="Frecuencia en no réplicas")

#plt.xticks(intervalos)

#%%

import statsmodels.api as sm
import pylab

sm.qqplot(replic, line='q')
pylab.title("QQ plot Replicas")
pylab.show()


#%% Como no tienen distribución normal, no debería un t test. Voy a probar a hacer un test no paramétrico
# Mann-Whitney U test
# Aunque como la diferencia es tan grande, no importa si hago t test o no paramétrico. 

from numpy.random import seed
import random
from scipy.stats import mannwhitneyu
# seed the random number generator
seed(1)
# generate two independent samples

random_replic = random.sample(list(replic), 1000)
random_no_replic = random.sample(list(no_replic), 1000)

# compare samples
stat, p = mannwhitneyu(random_replic, random_no_replic)
print('Statistics=%.3f, p=%.3f' % (stat, p))
# interpret
alpha = 0.05
if p > alpha:
	print('Same distribution (fail to reject H0)')
else:
	print('Different distribution (reject H0)')

#%%

porc_in_frame = pd.DataFrame(columns = lista_samples_min_freq) # nombres de las muestras

for sample_i in np.arange(len(lista_dict_min_freq)):
#for sample_i in np.arange(2):  
    if "no variant" in lista_dict_min_freq[sample_i].keys() and lista_dict_min_freq[sample_i]["no variant"] < 95: 
        # Solo cojo los experimento con más de un 5% de corte. 
        #print(lista_samples_min_freq[sample_i])
        porc_in_frame[porc_in_frame.columns[sample_i]] = [rpl.in_frame(lista_dict_min_freq[sample_i])]
        
porc_in_frame = porc_in_frame.dropna(axis = 1)
#%%

replica_1 = []
replica_2 = []

for col_i in np.arange(len(porc_in_frame.columns)): 
    for col_j in np.arange(col_i + 1,len(porc_in_frame.columns)):
        if porc_in_frame.columns[col_i].split('_')[0] == porc_in_frame.columns[col_j].split('_')[0]:
            replica_1.append(float(porc_in_frame[porc_in_frame.columns[col_i]]))
            replica_2.append(float(porc_in_frame[porc_in_frame.columns[col_j]]))
        else: 
            break
            
#%%

dif = np.array(replica_1) - np.array(replica_2)
sns.histplot(dif, bins = np.arange(-20,20,1))
plt.title("Distribution of the difference between the percentage \n of in frame mutations between replicates ")
plt.xlabel("Difference in the percentage of in frame mutations \n between replicates", size = 10)

#%%
#plt.hist(porc_in_frame.loc[0])

plt.scatter(x=replica_1, y=replica_2, marker='o', s=10, c='darkmagenta')
plt.xlabel("% In frame mutations (replica 1)", size = 10)
plt.ylabel("% In frame mutations (replica 2)", size = 10)
plt.title("In frame mutations reproducibility", size = 15)
#plt.xlim([0, 100])
#plt.ylim([0, 100])

#%% Correlación de Pearson
# https://machinelearningmastery.com/how-to-use-correlation-to-understand-the-relationship-between-variables/

from numpy.random import randn
from numpy.random import seed
from scipy.stats import pearsonr
import scipy.stats
from matplotlib.font_manager import FontProperties

slope, intercept, r, p, stderr = scipy.stats.linregress(replica_1, replica_2)
line = f'Regression line: y={intercept:.2f}+{slope:.2f}x, r={r:.2f}'

fontP = FontProperties()
fontP.set_size('small')
# these are matplotlib.patch.Patch properties
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

fig, ax = plt.subplots(figsize=(10, 8))
ax.plot(replica_1, replica_2, linewidth=0, marker='o', label='Data points', 
        markersize=5, c='cornflowerblue')
ax.plot(replica_1, intercept + slope * np.array(replica_1), 
        label=line, c = "black", linestyle=":")
ax.set_xlabel("% In frame mutations (replica 1)", size = 12)
ax.set_ylabel("% In frame mutations (replica 2)", size = 12)
plt.title("In frame mutations reproducibility", size = 15)
ax.legend(facecolor='white', bbox_to_anchor=(1.05, 1), loc='upper left', prop=fontP)
text = str('R = ') + str(round(r,3))
ax.text(1.2, 0.9, text, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', bbox=props)

plt.show()