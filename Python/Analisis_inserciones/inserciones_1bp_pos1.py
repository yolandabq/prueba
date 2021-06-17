# -*- coding: utf-8 -*-
"""
Created on Tue May  4 17:12:29 2021

@author: yolib

Script para analizar qué nucleótidos se suelen insertar tras la base X en la posición -1

"""

#%%
import os
import pandas as pd
import numpy as np
import seaborn as sns 
import matplotlib.pyplot as plt
import glob 
from collections import Counter
from Bio.Seq import Seq

os.chdir('C:\\Users\\yolib\\OneDrive\\Documentos\\TFM\\Documento TFM\\Codigo\\Python\\Funciones')

import utils_functions

#%% Cargo los datos

os.chdir('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos\\results_no_replic_ins')
guides_ins_no_replic_files = os.listdir()

#path_no_replic = 'C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos\\results_no_replic_ins\\' + guides_ins_no_replic_files[0]

os.chdir('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos\\results_ins')
guides_ins_files = os.listdir()

os.chdir('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos\\results')
guides_results_files = os.listdir()


os.chdir('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos')
guides_data_no_replic = pd.read_csv("guides_sequence_2_no_replic.csv")
guides_replic_data = pd.read_csv("guides_sequence_2.csv")

frames = [guides_replic_data, guides_data_no_replic]
guides_data = pd.concat(frames)


#%%  Guardamos el nucleótido insertado según su nucleótido adyacente en la posición -1
# Análisis cualitativo, no tenemos en cuenta las frecuencias. 

nt_list = {} # Guardo las secuencias dependiendo de la posición -1
nt_list["A"] = []
nt_list["C"] = []
nt_list["T"] = []
nt_list["G"] = []

#nt_list["no_comun"] = []


guides_id = list(set([int(x.split("_")[0]) for x in guides_results_files]))

for carpeta in ["results", "results_no_replic"]:
    
    os.chdir('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos\\' + carpeta)
    guides_files = os.listdir()
    guides_id = list(set([int(x.split("_")[0]) for x in guides_files]))

    for g_id in guides_id:
        
        #print(g_id)
        #print(glob.glob(str(g_id) + "_*"))
        files = glob.glob(str(g_id) + "_*")
        
        
        sequence = list(guides_data[guides_data["Guide_ID"] == g_id]["Sequence"])[0]
        nt_PC = sequence[29].upper()
        nt_PC_1 = sequence[30].upper()
        strand = list(guides_data[guides_data["Guide_ID"] == g_id]["Strand"])[0]
           
        
        if len(files) >= 2:
            
            df_list_results = []
            df_list_insertion = []

            for file in files:  # si hay más de un archivo de la misma guía: 
            
                df_list_results.append(pd.read_csv(file))
                df_list_insertion.append(pd.read_csv('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos\\' + carpeta + '_ins\\'+ file.split(".")[0] + "_ins.csv"))
            
            df_results = utils_functions.merge_df_results_alleles(df_list_results)
            df_insertion_data = pd.concat(df_list_insertion, ignore_index = True)
            suma_counts = dict(df_insertion_data.groupby(["sample_name"])["count"].sum())
            df_insertion_data['proportion'] = [df_insertion_data.loc[x, "count"]/suma_counts[df_insertion_data.loc[x, "sample_name"]] for x in df_insertion_data.index]

        else: 
            
            df_results = pd.read_csv(files[0])
            df_insertion_data = pd.read_csv('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos\\' + carpeta + '_ins\\'+ files[0].split(".")[0] + "_ins.csv")
            suma_counts = dict(df_insertion_data.groupby(["sample_name"])["count"].sum())
            df_insertion_data['proportion'] = [df_insertion_data.loc[x, "count"]/suma_counts[df_insertion_data.loc[x, "sample_name"]] for x in df_insertion_data.index]

        # me quito así el problema de tener que abrir varios archivos y tal
        # Ahora quiero, si hay -1:1I o 2:1I, ver si se puede cambiar el alineamiento y, si se puede, poner como nombre 1:1I. 
        
        if "-1:1I" in list(df_results["Alleles"]) and np.mean(list(df_results[df_results["Alleles"] == "-1:1I"].iloc[0, :-1])) > 0:
            
            nt = df_insertion_data[df_insertion_data["cigar_label"] == '-1:1I'] 
            
            if strand == '+': 
                
                if nt_PC in list(nt['seq']): 
                    
                    index = df_insertion_data[df_insertion_data["cigar_label"] == '-1:1I'][df_insertion_data["seq"] == nt_PC]["cigar_label"].index
                    df_insertion_data.loc[index, "cigar_label"] = "1:1I"
                    
            elif strand == "-": 
                
                if str(Seq(list(nt_PC)[0]).reverse_complement()) in list(nt['seq']): 
                    
                    index = df_insertion_data[df_insertion_data["cigar_label"] == '-1:1I'][df_insertion_data["seq"] == str(Seq(list(nt_PC)[0]).reverse_complement())]["cigar_label"].index
                    df_insertion_data.loc[index, "cigar_label"] = "1:1I"
                
        if "2:1I" in list(df_results["Alleles"]) and np.mean(list(df_results[df_results["Alleles"] == "2:1I"].iloc[0, :-1])) > 0:
            
            nt = df_insertion_data[df_insertion_data["cigar_label"] == '-1:1I']
            
            if strand == '+': 
                
                if nt_PC_1 in list(nt['seq']): 
                    
                    index = df_insertion_data[df_insertion_data["cigar_label"] == '-1:1I'][df_insertion_data["seq"] == nt_PC_1]["cigar_label"].index
                    df_insertion_data.loc[index, "cigar_label"] = "1:1I"
                    
            elif strand == "-": 
                
                if str(Seq(list(nt_PC)[0]).reverse_complement()) in list(nt['seq']): 
                    
                    index = df_insertion_data[df_insertion_data["cigar_label"] == '-1:1I'][df_insertion_data["seq"] == str(Seq(list(nt_PC)[0]).reverse_complement())]["cigar_label"].index
                    df_insertion_data.loc[index, "cigar_label"] = "1:1I"
        
        
        if "1:1I" in list(df_results["Alleles"]):
            
                samples_1_nt = []
                
                for sample in list(set(df_insertion_data["sample_name"])):
                    try: 
                        nt = df_insertion_data[df_insertion_data["cigar_label"] == '1:1I'][df_insertion_data["sample_name"] == sample].sort_values(by = "proportion", ascending = False)
                        samples_1_nt.append(nt.iloc[0, 1])
                        
                    except IndexError: 
                        samples_1_nt.append("NaN")
                        
        mut_common_1 = Counter(samples_1_nt)
        
        if len(mut_common_1) == 1: 
            
            if strand == '+': 
                nt_list[nt_PC].append(list(mut_common_1)[0])
                
            elif strand == '-': 
                nt_list[nt_PC].append(str(Seq(list(mut_common_1)[0]).reverse_complement()))
                    
        else: 
            
            nt_list[nt_PC].append("not common")
             
#%%

#colors = sns.color_palette("hls", 14)
colors = ['#ff9999','#66b3ff','#99ff99','#ffcc99', '#9999ff']

fig, ax1 = plt.subplots(figsize=(9, 6), subplot_kw=dict(aspect="equal"))

nt = "T" # A, C, G o T. Cambiar para obtener los otros

# Data to plot
labels = []
porc = []

for x, y in sorted(dict(Counter(nt_list[nt])).items()):
    labels.append(x)
    porc.append(y)
    
# Plot
patches, texts, autotexts = plt.pie(porc, colors = colors, startangle=90,
        labels=labels, autopct='%1.1f%%',  textprops={'fontsize': 9}, pctdistance=0.9)
        
"""
#draw circle
centre_circle = plt.Circle((0,0),0.70,fc='white')
fig = plt.gcf()
fig.gca().add_artist(centre_circle)
# Equal aspect ratio ensures that pie is drawn as a circle
ax1.axis('equal')  
plt.tight_layout()
"""

"""
ax.legend([str(x + ": " + str(round(dic_porc_mayor_alelo_comun[x], 2)) + "%") for x in dic_porc_mayor_alelo_comun.keys()],
          title="Alelo mayoritario",
          loc="center left",
          bbox_to_anchor=(1, 0, 0.5, 1))
"""
plt.title("Frequency of the inserted nucleotide when " + nt + " is at position -1")

plt.axis('equal')

plt.show()

#%% Análisis cuantitativo, ahora sí que guardo las frecuencias, no solo cual es el nt insertado más común según la base -1.

df_nt_menos_1 = pd.DataFrame(columns = ["nt_menos_1", "nt_inserted", "frequency"])

guides_id = list(set([int(x.split("_")[0]) for x in guides_results_files]))

for carpeta in ["results", "results_no_replic"]:
    
    os.chdir('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos\\' + carpeta)
    guides_files = os.listdir()
    guides_id = list(set([int(x.split("_")[0]) for x in guides_files]))

    for g_id in guides_id:
        
        #print(g_id)
        #print(glob.glob(str(g_id) + "_*"))
        files = glob.glob(str(g_id) + "_*")
        
        porc_i_ins = 0
        
        sequence = list(guides_data[guides_data["Guide_ID"] == g_id]["Sequence"])[0]
        nt_PC = sequence[29].upper()
        nt_PC_1 = sequence[30].upper()
        strand = list(guides_data[guides_data["Guide_ID"] == g_id]["Strand"])[0]
           
        
        if len(files) >= 2:
            
            df_list_results = []
            df_list_insertion = []

            for file in files:  # si hay más de un archivo de la misma guía: 
            
                df_list_results.append(pd.read_csv(file))
                df_list_insertion.append(pd.read_csv('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos\\' + carpeta + '_ins\\'+ file.split(".")[0] + "_ins.csv"))
            
            df_results = utils_functions.merge_df_results_alleles(df_list_results)
            df_insertion_data = pd.concat(df_list_insertion, ignore_index = True)
            suma_counts = dict(df_insertion_data.groupby(["sample_name"])["count"].sum())
            df_insertion_data['proportion'] = [df_insertion_data.loc[x, "count"]/suma_counts[df_insertion_data.loc[x, "sample_name"]] for x in df_insertion_data.index]

        else: 
            
            df_results = pd.read_csv(files[0])
            df_insertion_data = pd.read_csv('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos\\' + carpeta + '_ins\\'+ files[0].split(".")[0] + "_ins.csv")
            suma_counts = dict(df_insertion_data.groupby(["sample_name"])["count"].sum())
            df_insertion_data['proportion'] = [df_insertion_data.loc[x, "count"]/suma_counts[df_insertion_data.loc[x, "sample_name"]] for x in df_insertion_data.index]
            
            
        if "no variant" in list(df_results["Alleles"]):
            porc_corte = (100-np.mean(list(df_results[df_results["Alleles"] == "no variant"].iloc[0, :-1])))
            #porc_corte_list.append(porc_corte)
        
        else: 
            porc_corte = 100
            
        norm = 1/porc_corte  
        
        # me quito así el problema de tener que abrir varios archivos y tal
        # Ahora quiero, si hay -1:1I o 2:1I, ver si se puede cambiar el alineamiento y, si se puede, poner como nombre 1:1I. 
        
        if "-1:1I" in list(df_results["Alleles"]) and np.mean(list(df_results[df_results["Alleles"] == "-1:1I"].iloc[0, :-1])) > 0:
            
            nt = df_insertion_data[df_insertion_data["cigar_label"] == '-1:1I'] 
            
            if strand == '+': 
                
                if nt_PC in list(nt['seq']): 
                    
                    index = df_insertion_data[df_insertion_data["cigar_label"] == '-1:1I'][df_insertion_data["seq"] == nt_PC]["cigar_label"].index
                    df_insertion_data.loc[index, "cigar_label"] = "1:1I"
                    porc_i_ins = porc_i_ins + norm*np.mean(list(df_results[df_results["Alleles"] == "-1:1I"].iloc[0, :-1]))
                    
                    
            elif strand == "-": 
                
                if str(Seq(list(nt_PC)[0]).reverse_complement()) in list(nt['seq']): 
                    
                    index = df_insertion_data[df_insertion_data["cigar_label"] == '-1:1I'][df_insertion_data["seq"] == str(Seq(list(nt_PC)[0]).reverse_complement())]["cigar_label"].index
                    df_insertion_data.loc[index, "cigar_label"] = "1:1I"
                    porc_i_ins = porc_i_ins + norm*np.mean(list(df_results[df_results["Alleles"] == "-1:1I"].iloc[0, :-1]))
                    
                
                
        if "2:1I" in list(df_results["Alleles"]) and np.mean(list(df_results[df_results["Alleles"] == "2:1I"].iloc[0, :-1])) > 0:
            
            nt = df_insertion_data[df_insertion_data["cigar_label"] == '-1:1I']
            
            if strand == '+': 
                
                if nt_PC_1 in list(nt['seq']): 
                    
                    index = df_insertion_data[df_insertion_data["cigar_label"] == '2:1I'][df_insertion_data["seq"] == nt_PC_1]["cigar_label"].index
                    df_insertion_data.loc[index, "cigar_label"] = "1:1I"
                    porc_i_ins = porc_i_ins + norm*np.mean(list(df_results[df_results["Alleles"] == "2:1I"].iloc[0, :-1]))
                    
                    
                    
            elif strand == "-": 
                
                if str(Seq(list(nt_PC)[0]).reverse_complement()) in list(nt['seq']): 
                    
                    index = df_insertion_data[df_insertion_data["cigar_label"] == '2:1I'][df_insertion_data["seq"] == str(Seq(list(nt_PC)[0]).reverse_complement())]["cigar_label"].index
                    df_insertion_data.loc[index, "cigar_label"] = "1:1I"
                    porc_i_ins = porc_i_ins + norm*np.mean(list(df_results[df_results["Alleles"] == "2:1I"].iloc[0, :-1]))
                    
        if "1:1I" in list(df_results["Alleles"]): 
            
            porc_i_ins = porc_i_ins + norm*np.mean(list(df_results[df_results["Alleles"] == "1:1I"].iloc[0, :-1]))

            
        if "1:1I" in list(df_insertion_data["cigar_label"]):
            
                nt_df = df_insertion_data[df_insertion_data["cigar_label"] == '1:1I'].loc[:, ["seq", "sample_name", "proportion"]].groupby(["seq"]).mean()  
                #nt_df["frequency"] = nt_df["proportion"]*(porc_i_ins)
                
                for n in nt_df.index:
                    
                    if strand == "+":
                    
                        df_nt_menos_1 = df_nt_menos_1.append(pd.Series([nt_PC, n, nt_df.loc[n, "proportion"]], 
                                                                     index=df_nt_menos_1.columns),
                                                                     ignore_index=True)
                
                    else: 
                        df_nt_menos_1 = df_nt_menos_1.append(pd.Series([nt_PC, str(Seq(n).reverse_complement()), nt_df.loc[n, "proportion"]], 
                                                                     index=df_nt_menos_1.columns),
                                                                     ignore_index=True)
                
#%%

df_nt_menos_1 = df_nt_menos_1[df_nt_menos_1["nt_inserted"] != "N"]

plt.figure(figsize=(12, 10))

ax = sns.catplot(x = "nt_menos_1", y = "frequency", data = df_nt_menos_1, hue="nt_inserted", 
             dodge=True, kind="box")

plt.title("Frequency of the inserted nucleotide \n depending on position -1", fontsize = 15)

plt.xlabel("Nucleotide at position -1", size = 12)
plt.ylabel("1bp insertion frequency (%)", size = 12)

ax._legend.set_title("Nucleotide \n inserted")



"""
ax = sns.catplot(x = "nt_inserted", y = "frequency", data = df_nt_menos_1, hue="nt_menos_1", 
             dodge=True, size = 2 , kind = "strip")



ax = sns.boxplot(x = "nt_inserted", y = "frequency", data = df_nt_menos_1, hue="nt_menos_1", 
             dodge=True, whis=np.inf)
ax = sns.stripplot(x = "nt_inserted", y = "frequency", data = df_nt_menos_1, hue="nt_menos_1", 
             dodge=True, size = 2 )

"""

df_nt_menos_1.to_csv('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos\\analisis_ins\\ins_nt_menos1.csv', header=True, index=False)


#%% Vamos a ver cuántas guías tienen su alelo mayoritario 1:1I según su nt en la pos -1
# Tendria que normalizarlo con respecto al número de guías con ese nt en -1
# Es decir, del número de guías que tienen A en menos 1, cuántas de ellas tienen como alelo mayoritario 1:1I

nt_list_mayor = {} # Guardo las secuencias dependiendo de la posición X que elija 
nt_list_mayor["A"] = []
nt_list_mayor["C"] = []
nt_list_mayor["T"] = []
nt_list_mayor["G"] = []
#nt_list["no_comun"] = []

number_guides = {}
number_guides["A"] = 0
number_guides["C"] = 0
number_guides["G"] = 0
number_guides["T"] = 0


guides_id = list(set([int(x.split("_")[0]) for x in guides_results_files]))

for carpeta in ["results", "results_no_replic"]:
    
    os.chdir('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos\\' + carpeta)
    guides_files = os.listdir()
    guides_id = list(set([int(x.split("_")[0]) for x in guides_files]))

    for g_id in guides_id:
        
        #print(g_id)
        #print(glob.glob(str(g_id) + "_*"))
        files = glob.glob(str(g_id) + "_*")
        
        
        sequence = list(guides_data[guides_data["Guide_ID"] == g_id]["Sequence"])[0]
        nt_PC = sequence[29].upper()
        number_guides[nt_PC] = number_guides[nt_PC] + 1
            
        nt_PC_1 = sequence[30].upper()
        strand = list(guides_data[guides_data["Guide_ID"] == g_id]["Strand"])[0]
           
        
        if len(files) >= 2:
            
            df_list_results = []
            df_list_insertion = []

            for file in files:  # si hay más de un archivo de la misma guía: 
            
                df_list_results.append(pd.read_csv(file))
                df_list_insertion.append(pd.read_csv('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos\\' + carpeta + '_ins\\'+ file.split(".")[0] + "_ins.csv"))
            
            df_results = utils_functions.merge_df_results_alleles(df_list_results)
            df_insertion_data = pd.concat(df_list_insertion, ignore_index = True)
            suma_counts = dict(df_insertion_data.groupby(["sample_name"])["count"].sum())
            df_insertion_data['proportion'] = [df_insertion_data.loc[x, "count"]/suma_counts[df_insertion_data.loc[x, "sample_name"]] for x in df_insertion_data.index]

        else: 
            
            df_results = pd.read_csv(files[0])
            df_insertion_data = pd.read_csv('C:\\Users\\yolib\\Documents\\TFM\\Linf_T\\Datos\\' + carpeta + '_ins\\'+ files[0].split(".")[0] + "_ins.csv")
            suma_counts = dict(df_insertion_data.groupby(["sample_name"])["count"].sum())
            df_insertion_data['proportion'] = [df_insertion_data.loc[x, "count"]/suma_counts[df_insertion_data.loc[x, "sample_name"]] for x in df_insertion_data.index]

        # me quito así el problema de tener que abrir varios archivos y tal
        
        
        # Ordeno el dataframe de resultados y quito el "no variants"
        
        df_results = df_results.sort_values(by = list(df_results.columns)[0:-1], ascending = False)
        
        if "no variant" in list(df_results["Alleles"]): 
            df_results = df_results.drop(df_results[df_results['Alleles']=='no variant'].index)

        if df_results.iloc[0, -1] == "1:1I": 
            nt_list_mayor[nt_PC].append(g_id)
            
        elif df_results.iloc[0, -1] == "-1:1I":
                
            nt = df_insertion_data[df_insertion_data["cigar_label"] == '-1:1I']
            
            if strand == '+': 
                
                if nt_PC in list(nt['seq']): 
                    
                    nt_list_mayor[nt_PC].append(g_id)
                    
            elif strand == "-": 
                
                if str(Seq(list(nt_PC)[0]).reverse_complement()) in list(nt['seq']): 
                    
                    nt_list_mayor[nt_PC].append(g_id)
                    
        elif df_results.iloc[0, -1] == "2:1I":
            
            nt = df_insertion_data[df_insertion_data["cigar_label"] == '-1:1I']
            
            if strand == '+': 
                
                if nt_PC_1 in list(nt['seq']): 
                    
                    nt_list_mayor[nt_PC].append(g_id)
                    
            elif strand == "-": 
                
                if str(Seq(list(nt_PC)[0]).reverse_complement()) in list(nt['seq']): 
                    
                    nt_list_mayor[nt_PC].append(g_id)
        
        
#%%

plt.figure(figsize=(9, 6))

# Data to plot
labels = []
porc = []

for x, y in nt_list_mayor.items():
    labels.append(x)
    porc.append(len(y)/number_guides[x]) # normalizo entre el número de guías con X nt en la pos -1
    # Es decir, estoy calculando, por ejemplo, de las 657 guías que tienen una A en la pos -1, el 40% (273 guías)
    # tienen como alelo mayoritario 1:1I

df_nt_mayor = pd.DataFrame(columns=["Nucleotide -1", "Percent of guides"])
df_nt_mayor["Nucleotide -1"] = labels
df_nt_mayor["Percent of guides"] = porc

plot = sns.barplot(x='Nucleotide -1', y='Percent of guides', data=df_nt_mayor)
plt.title("Percent of guides with the most common allele 1:1I")


for p in plot.patches:
    plot.annotate(format(p.get_height(), '.1f'), 
                   (p.get_x() + p.get_width() / 2., p.get_height()), 
                   ha = 'center', va = 'center', 
                   xytext = (0, 5), 
                   textcoords = 'offset points')

plt.show()



