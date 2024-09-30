# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 11:45:43 2024

@author: Elisabetta
"""

import pandas as pd
from scipy import sparse #This is necessary to create the sparse matrix, which is a lighter matrix in which zero values have been removed.
import numpy as np
import os
import seaborn as sb #This is needed for the graphical rapresentation of the monte carlo results
import matplotlib.pyplot as plt #This is needed for the graphical rapresentation of the monte carlo results
from matplotlib.ticker import FuncFormatter #This is needed for the graphical rapresentation of the monte carlo results

#%% - CASES COMPARISON EXIOBASE SMALL - Sector "EU28-Agriculture-Forestry-Fishing"

folder = r"C:\Users\Elisabetta\Desktop\Risultati_analisi_incertezza_Aalborg\Exiobase_small_sector_1"

file_1 = "CASE_1_uniform_0.1_MC_simulations_EU28-Agriculture-Forestry-Fishing.csv"
file_2 = "CASE_2_uniform_0.2_MC_simulations_EU28-Agriculture-Forestry-Fishing.csv"
file_3 = "CASE_3_uniform_0.3_MC_simulations_EU28-Agriculture-Forestry-Fishing.csv"
file_4 = "CASE_4_log-normal_1.1_MC_simulations_EU28-Agriculture-Forestry-Fishing.csv"
file_5 = "CASE_5_log-normal_1.2_MC_simulations_EU28-Agriculture-Forestry-Fishing.csv"
file_6 = "CASE_6_log-normal_1.3_MC_simulations_EU28-Agriculture-Forestry-Fishing.csv"

comp_1 = pd.read_csv(os.path.join(folder, file_1))
comp_2 = pd.read_csv(os.path.join(folder, file_2))
comp_3 = pd.read_csv(os.path.join(folder, file_3))
comp_4 = pd.read_csv(os.path.join(folder, file_4))
comp_5 = pd.read_csv(os.path.join(folder, file_5))
comp_6 = pd.read_csv(os.path.join(folder, file_6))

cases_comparision_exio_small =pd.DataFrame({
    'Unif, 10%': comp_1['kg CO2eq'],
    'Unif, 20%': comp_2['kg CO2eq'],
    'Unif, 30%': comp_3['kg CO2eq'],
    'log-normal, 1.1': comp_4['kg CO2eq'],
    'log-normal, 1.2': comp_5['kg CO2eq'],
    'log-normal, 1.3': comp_6['kg CO2eq']
})

def scientific_format(x, pos):
    return f'{x:.2e}'
formatter = FuncFormatter(scientific_format)

plt.figure(figsize=(12, 8))  # Adjust the figure size
sb.boxplot(data=cases_comparision_exio_small, palette="Set2")

# Add labels and title
plt.xlabel('Scenarios')
plt.ylabel('kg CO2eq')
plt.title('EU28-Agriculture-Forestry-Fishing (exiobase_small)')
plt.gca().yaxis.set_major_formatter(formatter)

#Save the plot
filename = "Monte_Carlo_Comparison_between_cases.png"  # You can change the extension to .pdf, .jpg, etc.

# Create the folder if it doesn't exist
os.makedirs(folder, exist_ok=True)

# Save the plot
plt.savefig(os.path.join(folder, filename))



#%% - CASES COMPARISON EXIOBASE SMALL - Sector "EU28-Basic_iron_and_steel_and_of_ferro-alloys_and_first_products_thereof"

folder = r"C:\Users\Elisabetta\Desktop\Risultati_analisi_incertezza_Aalborg\Exiobase_small_sector_2"

file_1 = "CASE_1_uniform_0.1_MC_simulations_EU28-Basic_iron_and_steel_and_of_ferro-alloys_and_first_products_thereof.csv"
file_2 = "CASE_2_uniform_0.2_MC_simulations_EU28-Basic_iron_and_steel_and_of_ferro-alloys_and_first_products_thereof.csv"
file_3 = "CASE_3_uniform_0.3_MC_simulations_EU28-Basic_iron_and_steel_and_of_ferro-alloys_and_first_products_thereof.csv"
file_4 = "CASE_4_log-normal_1.1_MC_simulations_EU28-Basic_iron_and_steel_and_of_ferro-alloys_and_first_products_thereof.csv"
file_5 = "CASE_5_log-normal_1.2_MC_simulations_EU28-Basic_iron_and_steel_and_of_ferro-alloys_and_first_products_thereof.csv"
file_6 = "CASE_6_log-normal_1.3_MC_simulations_EU28-Basic_iron_and_steel_and_of_ferro-alloys_and_first_products_thereof.csv"

comp_1 = pd.read_csv(os.path.join(folder, file_1))
comp_2 = pd.read_csv(os.path.join(folder, file_2))
comp_3 = pd.read_csv(os.path.join(folder, file_3))
comp_4 = pd.read_csv(os.path.join(folder, file_4))
comp_5 = pd.read_csv(os.path.join(folder, file_5))
comp_6 = pd.read_csv(os.path.join(folder, file_6))

cases_comparision_exio_small =pd.DataFrame({
    'Unif, 10%': comp_1['kg CO2eq'],
    'Unif, 20%': comp_2['kg CO2eq'],
    'Unif, 30%': comp_3['kg CO2eq'],
    'log-normal, 1.1': comp_4['kg CO2eq'],
    'log-normal, 1.2': comp_5['kg CO2eq'],
    'log-normal, 1.3': comp_6['kg CO2eq']
})

def scientific_format(x, pos):
    return f'{x:.2e}'
formatter = FuncFormatter(scientific_format)

plt.figure(figsize=(12, 8))  # Adjust the figure size
sb.boxplot(data=cases_comparision_exio_small, palette="Set2")

# Add labels and title
plt.xlabel('Scenarios')
plt.ylabel('kg CO2eq')
plt.title('EU28-Basic_iron_and_steel_and_of_ferro-alloys_and_first_products_thereof (exiobase_small)')
plt.gca().yaxis.set_major_formatter(formatter)

#Save the plot
filename = "Monte_Carlo_Comparison_between_cases.png"  # You can change the extension to .pdf, .jpg, etc.

# Create the folder if it doesn't exist
os.makedirs(folder, exist_ok=True)

# Save the plot
plt.savefig(os.path.join(folder, filename))


#%% - CASES COMPARISON EXIOBASE SMALL - Sector "EU28-Biodiesels"

folder = r"C:\Users\Elisabetta\Desktop\Risultati_analisi_incertezza_Aalborg\Exiobase_small_sector_3"

file_1 = "CASE_1_uniform_0.1_MC_simulations_EU28-Biodiesels.csv"
file_2 = "CASE_2_uniform_0.2_MC_simulations_EU28-Biodiesels.csv"
file_3 = "CASE_3_uniform_0.3_MC_simulations_EU28-Biodiesels.csv"
file_4 = "CASE_4_log-normal_1.1_MC_simulations_EU28-Biodiesels.csv"
file_5 = "CASE_5_log-normal_1.2_MC_simulations_EU28-Biodiesels.csv"
file_6 = "CASE_6_log-normal_1.3_MC_simulations_EU28-Biodiesels.csv"

comp_1 = pd.read_csv(os.path.join(folder, file_1))
comp_2 = pd.read_csv(os.path.join(folder, file_2))
comp_3 = pd.read_csv(os.path.join(folder, file_3))
comp_4 = pd.read_csv(os.path.join(folder, file_4))
comp_5 = pd.read_csv(os.path.join(folder, file_5))
comp_6 = pd.read_csv(os.path.join(folder, file_6))

cases_comparision_exio_small =pd.DataFrame({
    'Unif, 10%': comp_1['kg CO2eq'],
    'Unif, 20%': comp_2['kg CO2eq'],
    'Unif, 30%': comp_3['kg CO2eq'],
    'log-normal, 1.1': comp_4['kg CO2eq'],
    'log-normal, 1.2': comp_5['kg CO2eq'],
    'log-normal, 1.3': comp_6['kg CO2eq']
})

def scientific_format(x, pos):
    return f'{x:.2e}'
formatter = FuncFormatter(scientific_format)

plt.figure(figsize=(12, 8))  # Adjust the figure size
sb.boxplot(data=cases_comparision_exio_small, palette="Set2")

# Add labels and title
plt.xlabel('Scenarios')
plt.ylabel('kg CO2eq')
plt.title('EU28-Biodiesels (exiobase_small)')
plt.gca().yaxis.set_major_formatter(formatter)

#Save the plot
filename = "Monte_Carlo_Comparison_between_cases.png"  # You can change the extension to .pdf, .jpg, etc.

# Create the folder if it doesn't exist
os.makedirs(folder, exist_ok=True)

# Save the plot
plt.savefig(os.path.join(folder, filename))


#%% - CASES COMPARISON EXIOBASE SMALL - Sector "RoW-Services"

folder = r"C:\Users\Elisabetta\Desktop\Risultati_analisi_incertezza_Aalborg\Exiobase_small_sector_4"

file_1 = "CASE_1_uniform_0.1_MC_simulations_RoW-Services.csv"
file_2 = "CASE_2_uniform_0.2_MC_simulations_RoW-Services.csv"
file_3 = "CASE_3_uniform_0.3_MC_simulations_RoW-Services.csv"
file_4 = "CASE_4_log-normal_1.1_MC_simulations_RoW-Services.csv"
file_5 = "CASE_5_log-normal_1.2_MC_simulations_RoW-Services.csv"
file_6 = "CASE_6_log-normal_1.3_MC_simulations_RoW-Services.csv"

comp_1 = pd.read_csv(os.path.join(folder, file_1))
comp_2 = pd.read_csv(os.path.join(folder, file_2))
comp_3 = pd.read_csv(os.path.join(folder, file_3))
comp_4 = pd.read_csv(os.path.join(folder, file_4))
comp_5 = pd.read_csv(os.path.join(folder, file_5))
comp_6 = pd.read_csv(os.path.join(folder, file_6))

cases_comparision_exio_small =pd.DataFrame({
    'Unif, 10%': comp_1['kg CO2eq'],
    'Unif, 20%': comp_2['kg CO2eq'],
    'Unif, 30%': comp_3['kg CO2eq'],
    'log-normal, 1.1': comp_4['kg CO2eq'],
    'log-normal, 1.2': comp_5['kg CO2eq'],
    'log-normal, 1.3': comp_6['kg CO2eq']
})

def scientific_format(x, pos):
    return f'{x:.2e}'
formatter = FuncFormatter(scientific_format)

plt.figure(figsize=(12, 8))  # Adjust the figure size
sb.boxplot(data=cases_comparision_exio_small, palette="Set2")

# Add labels and title
plt.xlabel('Scenarios')
plt.ylabel('kg CO2eq')
plt.title('RoW-Services (exiobase_small)')
plt.gca().yaxis.set_major_formatter(formatter)

#Save the plot
filename = "Monte_Carlo_Comparison_between_cases.png"  # You can change the extension to .pdf, .jpg, etc.

# Create the folder if it doesn't exist
os.makedirs(folder, exist_ok=True)

# Save the plot
plt.savefig(os.path.join(folder, filename))



#%% - SECTORS COMPARISON EXIOBASE SMALL - CASE 1 (unifrom distribution, range 10%)

folder_1 = r"C:\Users\Elisabetta\Desktop\Risultati_analisi_incertezza_Aalborg\Exiobase_small_sector_1"
folder_2 = r"C:\Users\Elisabetta\Desktop\Risultati_analisi_incertezza_Aalborg\Exiobase_small_sector_2"
folder_3 = r"C:\Users\Elisabetta\Desktop\Risultati_analisi_incertezza_Aalborg\Exiobase_small_sector_3"
folder_4 = r"C:\Users\Elisabetta\Desktop\Risultati_analisi_incertezza_Aalborg\Exiobase_small_sector_4"

file_1 = "CASE_1_uniform_0.1_MC_simulations_EU28-Agriculture-Forestry-Fishing.csv"
file_2 = "CASE_1_uniform_0.1_MC_simulations_EU28-Basic_iron_and_steel_and_of_ferro-alloys_and_first_products_thereof.csv"
file_3 = "CASE_1_uniform_0.1_MC_simulations_EU28-Biodiesels.csv"
file_4 = "CASE_1_uniform_0.1_MC_simulations_RoW-Services.csv"


comp_1 = pd.read_csv(os.path.join(folder_1, file_1))
comp_2 = pd.read_csv(os.path.join(folder_2, file_2))
comp_3 = pd.read_csv(os.path.join(folder_3, file_3))
comp_4 = pd.read_csv(os.path.join(folder_4, file_4))


cases_comparision_exio_small =pd.DataFrame({
    'EU28-Agriculture\n-Forestry\n-Fishing': comp_1['kg CO2eq'],
    'EU28-Basic_iron\nand_steel_and_of\nferro-alloys_and\nfirst_products_thereof': comp_2['kg CO2eq'],
    'EU28-Biodiesels': comp_3['kg CO2eq'],
    'RoW-Services': comp_4['kg CO2eq'],
})

def scientific_format(x, pos):
    return f'{x:.2e}'
formatter = FuncFormatter(scientific_format)

plt.figure(figsize=(12, 8))  # Adjust the figure size
sb.boxplot(data=cases_comparision_exio_small, palette="Set2")

# Add labels and title
plt.xlabel('Scenarios')
plt.ylabel('kg CO2eq')
plt.title('Case 1 - uniform distribution, range 10% (exiobase_small)')
plt.gca().yaxis.set_major_formatter(formatter)

#Save the plot
filename = "Monte_Carlo_Comparison_between_sectors.png"  # You can change the extension to .pdf, .jpg, etc.

# Create the folder if it doesn't exist
os.makedirs(folder, exist_ok=True)

# Save the plot
plt.savefig(os.path.join(folder, filename))


#%% - SECTORS COMPARISON EXIOBASE SMALL - CASE 2 (unifrom distribution, range 20%)

folder_1 = r"C:\Users\Elisabetta\Desktop\Risultati_analisi_incertezza_Aalborg\Exiobase_small_sector_1"
folder_2 = r"C:\Users\Elisabetta\Desktop\Risultati_analisi_incertezza_Aalborg\Exiobase_small_sector_2"
folder_3 = r"C:\Users\Elisabetta\Desktop\Risultati_analisi_incertezza_Aalborg\Exiobase_small_sector_3"
folder_4 = r"C:\Users\Elisabetta\Desktop\Risultati_analisi_incertezza_Aalborg\Exiobase_small_sector_4"

file_1 = "CASE_2_uniform_0.2_MC_simulations_EU28-Agriculture-Forestry-Fishing.csv"
file_2 = "CASE_2_uniform_0.2_MC_simulations_EU28-Basic_iron_and_steel_and_of_ferro-alloys_and_first_products_thereof.csv"
file_3 = "CASE_2_uniform_0.2_MC_simulations_EU28-Biodiesels.csv"
file_4 = "CASE_2_uniform_0.2_MC_simulations_RoW-Services.csv"

comp_1 = pd.read_csv(os.path.join(folder_1, file_1))
comp_2 = pd.read_csv(os.path.join(folder_2, file_2))
comp_3 = pd.read_csv(os.path.join(folder_3, file_3))
comp_4 = pd.read_csv(os.path.join(folder_4, file_4))

cases_comparision_exio_small =pd.DataFrame({
    'EU28-Agriculture\n-Forestry\n-Fishing': comp_1['kg CO2eq'],
    'EU28-Basic_iron\nand_steel_and_of\nferro-alloys_and\nfirst_products_thereof': comp_2['kg CO2eq'],
    'EU28-Biodiesels': comp_3['kg CO2eq'],
    'RoW-Services': comp_4['kg CO2eq'],
})

def scientific_format(x, pos):
    return f'{x:.2e}'
formatter = FuncFormatter(scientific_format)

plt.figure(figsize=(12, 8))  # Adjust the figure size
sb.boxplot(data=cases_comparision_exio_small, palette="Set2")

# Add labels and title
plt.xlabel('Scenarios')
plt.ylabel('kg CO2eq')
plt.title('Case 2 - uniform distribution, range 20% (exiobase_small)')
plt.gca().yaxis.set_major_formatter(formatter)

#Save the plot
filename = "Monte_Carlo_Comparison_between_sectors.png"  # You can change the extension to .pdf, .jpg, etc.

# Create the folder if it doesn't exist
os.makedirs(folder, exist_ok=True)

# Save the plot
plt.savefig(os.path.join(folder, filename))


#%% - SECTORS COMPARISON EXIOBASE SMALL - CASE 3 (unifrom distribution, range 30%)

folder_1 = r"C:\Users\Elisabetta\Desktop\Risultati_analisi_incertezza_Aalborg\Exiobase_small_sector_1"
folder_2 = r"C:\Users\Elisabetta\Desktop\Risultati_analisi_incertezza_Aalborg\Exiobase_small_sector_2"
folder_3 = r"C:\Users\Elisabetta\Desktop\Risultati_analisi_incertezza_Aalborg\Exiobase_small_sector_3"
folder_4 = r"C:\Users\Elisabetta\Desktop\Risultati_analisi_incertezza_Aalborg\Exiobase_small_sector_4"

file_1 = "CASE_3_uniform_0.3_MC_simulations_EU28-Agriculture-Forestry-Fishing.csv"
file_2 = "CASE_3_uniform_0.3_MC_simulations_EU28-Basic_iron_and_steel_and_of_ferro-alloys_and_first_products_thereof.csv"
file_3 = "CASE_3_uniform_0.3_MC_simulations_EU28-Biodiesels.csv"
file_4 = "CASE_3_uniform_0.3_MC_simulations_RoW-Services.csv"

comp_1 = pd.read_csv(os.path.join(folder_1, file_1))
comp_2 = pd.read_csv(os.path.join(folder_2, file_2))
comp_3 = pd.read_csv(os.path.join(folder_3, file_3))
comp_4 = pd.read_csv(os.path.join(folder_4, file_4))

cases_comparision_exio_small =pd.DataFrame({
    'EU28-Agriculture\n-Forestry\n-Fishing': comp_1['kg CO2eq'],
    'EU28-Basic_iron\nand_steel_and_of\nferro-alloys_and\nfirst_products_thereof': comp_2['kg CO2eq'],
    'EU28-Biodiesels': comp_3['kg CO2eq'],
    'RoW-Services': comp_4['kg CO2eq'],
})

def scientific_format(x, pos):
    return f'{x:.2e}'
formatter = FuncFormatter(scientific_format)

plt.figure(figsize=(12, 8))  # Adjust the figure size
sb.boxplot(data=cases_comparision_exio_small, palette="Set2")

# Add labels and title
plt.xlabel('Scenarios')
plt.ylabel('kg CO2eq')
plt.title('Case 3 - uniform distribution, range 30% (exiobase_small)')
plt.gca().yaxis.set_major_formatter(formatter)

#Save the plot
filename = "Monte_Carlo_Comparison_between_sectors.png"  # You can change the extension to .pdf, .jpg, etc.

# Create the folder if it doesn't exist
os.makedirs(folder, exist_ok=True)

# Save the plot
plt.savefig(os.path.join(folder, filename))


#%% - SECTORS COMPARISON EXIOBASE SMALL - CASE 4 (log-normal distribution, GSD 1.1)

folder_1 = r"C:\Users\Elisabetta\Desktop\Risultati_analisi_incertezza_Aalborg\Exiobase_small_sector_1"
folder_2 = r"C:\Users\Elisabetta\Desktop\Risultati_analisi_incertezza_Aalborg\Exiobase_small_sector_2"
folder_3 = r"C:\Users\Elisabetta\Desktop\Risultati_analisi_incertezza_Aalborg\Exiobase_small_sector_3"
folder_4 = r"C:\Users\Elisabetta\Desktop\Risultati_analisi_incertezza_Aalborg\Exiobase_small_sector_4"

file_1 = "CASE_4_log-normal_1.1_MC_simulations_EU28-Agriculture-Forestry-Fishing.csv"
file_2 = "CASE_4_log-normal_1.1_MC_simulations_EU28-Basic_iron_and_steel_and_of_ferro-alloys_and_first_products_thereof.csv"
file_3 = "CASE_4_log-normal_1.1_MC_simulations_EU28-Biodiesels.csv"
file_4 = "CASE_4_log-normal_1.1_MC_simulations_RoW-Services.csv"

comp_1 = pd.read_csv(os.path.join(folder_1, file_1))
comp_2 = pd.read_csv(os.path.join(folder_2, file_2))
comp_3 = pd.read_csv(os.path.join(folder_3, file_3))
comp_4 = pd.read_csv(os.path.join(folder_4, file_4))

cases_comparision_exio_small =pd.DataFrame({
    'EU28-Agriculture\n-Forestry\n-Fishing': comp_1['kg CO2eq'],
    'EU28-Basic_iron\nand_steel_and_of\nferro-alloys_and\nfirst_products_thereof': comp_2['kg CO2eq'],
    'EU28-Biodiesels': comp_3['kg CO2eq'],
    'RoW-Services': comp_4['kg CO2eq'],
})

def scientific_format(x, pos):
    return f'{x:.2e}'
formatter = FuncFormatter(scientific_format)

plt.figure(figsize=(12, 8))  # Adjust the figure size
sb.boxplot(data=cases_comparision_exio_small, palette="Set2")

# Add labels and title
plt.xlabel('Scenarios')
plt.ylabel('kg CO2eq')
plt.title('Case 4 - log-normal distribution, GSD 1.1 (exiobase_small)')
plt.gca().yaxis.set_major_formatter(formatter)

#Save the plot
filename = "Monte_Carlo_Comparison_between_sectors.png"  # You can change the extension to .pdf, .jpg, etc.

# Create the folder if it doesn't exist
os.makedirs(folder, exist_ok=True)

# Save the plot
plt.savefig(os.path.join(folder, filename))


#%% - SECTORS COMPARISON EXIOBASE SMALL - CASE 5 (log-normal distribution, GSD 1.2)

folder_1 = r"C:\Users\Elisabetta\Desktop\Risultati_analisi_incertezza_Aalborg\Exiobase_small_sector_1"
folder_2 = r"C:\Users\Elisabetta\Desktop\Risultati_analisi_incertezza_Aalborg\Exiobase_small_sector_2"
folder_3 = r"C:\Users\Elisabetta\Desktop\Risultati_analisi_incertezza_Aalborg\Exiobase_small_sector_3"
folder_4 = r"C:\Users\Elisabetta\Desktop\Risultati_analisi_incertezza_Aalborg\Exiobase_small_sector_4"

file_1 = "CASE_5_log-normal_1.2_MC_simulations_EU28-Agriculture-Forestry-Fishing.csv"
file_2 = "CASE_5_log-normal_1.2_MC_simulations_EU28-Basic_iron_and_steel_and_of_ferro-alloys_and_first_products_thereof.csv"
file_3 = "CASE_5_log-normal_1.2_MC_simulations_EU28-Biodiesels.csv"
file_4 = "CASE_5_log-normal_1.2_MC_simulations_RoW-Services.csv"

comp_1 = pd.read_csv(os.path.join(folder_1, file_1))
comp_2 = pd.read_csv(os.path.join(folder_2, file_2))
comp_3 = pd.read_csv(os.path.join(folder_3, file_3))
comp_4 = pd.read_csv(os.path.join(folder_4, file_4))

cases_comparision_exio_small =pd.DataFrame({
    'EU28-Agriculture\n-Forestry\n-Fishing': comp_1['kg CO2eq'],
    'EU28-Basic_iron\nand_steel_and_of\nferro-alloys_and\nfirst_products_thereof': comp_2['kg CO2eq'],
    'EU28-Biodiesels': comp_3['kg CO2eq'],
    'RoW-Services': comp_4['kg CO2eq'],
})

def scientific_format(x, pos):
    return f'{x:.2e}'
formatter = FuncFormatter(scientific_format)

plt.figure(figsize=(12, 8))  # Adjust the figure size
sb.boxplot(data=cases_comparision_exio_small, palette="Set2")

# Add labels and title
plt.xlabel('Scenarios')
plt.ylabel('kg CO2eq')
plt.title('Case 5 - log-normal distribution, GSD 1.2 (exiobase_small)')
plt.gca().yaxis.set_major_formatter(formatter)

#Save the plot
filename = "Monte_Carlo_Comparison_between_sectors.png"  # You can change the extension to .pdf, .jpg, etc.

# Create the folder if it doesn't exist
os.makedirs(folder, exist_ok=True)

# Save the plot
plt.savefig(os.path.join(folder, filename))


#%% - SECTORS COMPARISON EXIOBASE SMALL - CASE 6 (log-normal distribution, GSD 1.3)

folder_1 = r"C:\Users\Elisabetta\Desktop\Risultati_analisi_incertezza_Aalborg\Exiobase_small_sector_1"
folder_2 = r"C:\Users\Elisabetta\Desktop\Risultati_analisi_incertezza_Aalborg\Exiobase_small_sector_2"
folder_3 = r"C:\Users\Elisabetta\Desktop\Risultati_analisi_incertezza_Aalborg\Exiobase_small_sector_3"
folder_4 = r"C:\Users\Elisabetta\Desktop\Risultati_analisi_incertezza_Aalborg\Exiobase_small_sector_4"

file_1 = "CASE_6_log-normal_1.3_MC_simulations_EU28-Agriculture-Forestry-Fishing.csv"
file_2 = "CASE_6_log-normal_1.3_MC_simulations_EU28-Basic_iron_and_steel_and_of_ferro-alloys_and_first_products_thereof.csv"
file_3 = "CASE_6_log-normal_1.3_MC_simulations_EU28-Biodiesels.csv"
file_4 = "CASE_6_log-normal_1.3_MC_simulations_RoW-Services.csv"

comp_1 = pd.read_csv(os.path.join(folder_1, file_1))
comp_2 = pd.read_csv(os.path.join(folder_2, file_2))
comp_3 = pd.read_csv(os.path.join(folder_3, file_3))
comp_4 = pd.read_csv(os.path.join(folder_4, file_4))

cases_comparision_exio_small =pd.DataFrame({
    'EU28-Agriculture\n-Forestry\n-Fishing': comp_1['kg CO2eq'],
    'EU28-Basic_iron\nand_steel_and_of\nferro-alloys_and\nfirst_products_thereof': comp_2['kg CO2eq'],
    'EU28-Biodiesels': comp_3['kg CO2eq'],
    'RoW-Services': comp_4['kg CO2eq'],
})

def scientific_format(x, pos):
    return f'{x:.2e}'
formatter = FuncFormatter(scientific_format)

plt.figure(figsize=(12, 8))  # Adjust the figure size
sb.boxplot(data=cases_comparision_exio_small, palette="Set2")

# Add labels and title
plt.xlabel('Scenarios')
plt.ylabel('kg CO2eq')
plt.title('Case 6 - log-normal distribution, GSD 1.3 (exiobase_small)')
plt.gca().yaxis.set_major_formatter(formatter)

#Save the plot
filename = "Monte_Carlo_Comparison_between_sectors.png"  # You can change the extension to .pdf, .jpg, etc.

# Create the folder if it doesn't exist
os.makedirs(folder, exist_ok=True)

# Save the plot
plt.savefig(os.path.join(folder, filename))




#%% - CASES COMPARISON EXIOBASE BIG - Sector "CH-Beverages"

folder = r"C:\Users\Elisabetta\Desktop\Risultati_analisi_incertezza_Aalborg\Exiobase_big_sector_1"

file_1 = "CASE_1_uniform_0.1_MC_simulations_CH-Beverages.csv"
file_2 = "CASE_2_uniform_0.2_MC_simulations_CH-Beverages.csv"
file_3 = "CASE_3_uniform_0.3_MC_simulations_CH-Beverages.csv"
file_4 = "CASE_4_log-normal_1.1_MC_simulations_CH-Beverages.csv"
file_5 = "CASE_5_log-normal_1.2_MC_simulations_CH-Beverages.csv"
file_6 = "CASE_6_log-normal_1.3_MC_simulations_CH-Beverages.csv"

comp_1 = pd.read_csv(os.path.join(folder, file_1))
comp_2 = pd.read_csv(os.path.join(folder, file_2))
comp_3 = pd.read_csv(os.path.join(folder, file_3))
comp_4 = pd.read_csv(os.path.join(folder, file_4))
comp_5 = pd.read_csv(os.path.join(folder, file_5))
comp_6 = pd.read_csv(os.path.join(folder, file_6))

cases_comparision_exio_small =pd.DataFrame({
    'Unif, 10%': comp_1['kg CO2eq'],
    'Unif, 20%': comp_2['kg CO2eq'],
    'Unif, 30%': comp_3['kg CO2eq'],
    'log-normal, 1.1': comp_4['kg CO2eq'],
    'log-normal, 1.2': comp_5['kg CO2eq'],
    'log-normal, 1.3': comp_6['kg CO2eq']
})

def scientific_format(x, pos):
    return f'{x:.2e}'
formatter = FuncFormatter(scientific_format)

plt.figure(figsize=(12, 8))  # Adjust the figure size
sb.boxplot(data=cases_comparision_exio_small, palette="Set2")

# Add labels and title
plt.xlabel('Scenarios')
plt.ylabel('kg CO2eq')
plt.title('CH-Beverages (exiobase_big)')
plt.gca().yaxis.set_major_formatter(formatter)

#Save the plot
filename = "Monte_Carlo_Comparison_between_cases.png"  # You can change the extension to .pdf, .jpg, etc.

# Create the folder if it doesn't exist
os.makedirs(folder, exist_ok=True)

# Save the plot
plt.savefig(os.path.join(folder, filename))



#%% - CASES COMPARISON EXIOBASE SMALL - Sector "SE-SE-Basic iron and steel and of ferro-alloys and first products thereof"

folder = r"C:\Users\Elisabetta\Desktop\Risultati_analisi_incertezza_Aalborg\Exiobase_big_sector_2"

file_1 = "CASE_1_uniform_0.1_MC_simulations_SE-Basic iron and steel and of ferro-alloys and first products thereof.csv"
file_2 = "CASE_2_uniform_0.2_MC_simulations_SE-Basic iron and steel and of ferro-alloys and first products thereof.csv"
file_3 = "CASE_3_uniform_0.3_MC_simulations_SE-Basic iron and steel and of ferro-alloys and first products thereof.csv"
file_4 = "CASE_4_log-normal_1.1_MC_simulations_SE-Basic iron and steel and of ferro-alloys and first products thereof.csv"
file_5 = "CASE_5_log-normal_1.2_MC_simulations_SE-Basic_iron_and_steel_and_of_ferro-alloys_and_first_products_thereof.csv"
file_6 = "CASE_6_log-normal_1.3_MC_simulations_SE-Basic_iron_and_steel_and_of_ferro-alloys_and_first_products_thereof.csv"

comp_1 = pd.read_csv(os.path.join(folder, file_1))
comp_2 = pd.read_csv(os.path.join(folder, file_2))
comp_3 = pd.read_csv(os.path.join(folder, file_3))
comp_4 = pd.read_csv(os.path.join(folder, file_4))
comp_5 = pd.read_csv(os.path.join(folder, file_5))
comp_6 = pd.read_csv(os.path.join(folder, file_6))

cases_comparision_exio_small =pd.DataFrame({
    'Unif, 10%': comp_1['kg CO2eq'],
    'Unif, 20%': comp_2['kg CO2eq'],
    'Unif, 30%': comp_3['kg CO2eq'],
    'log-normal, 1.1': comp_4['kg CO2eq'],
    'log-normal, 1.2': comp_5['kg CO2eq'],
    'log-normal, 1.3': comp_6['kg CO2eq']
})

def scientific_format(x, pos):
    return f'{x:.2e}'
formatter = FuncFormatter(scientific_format)

plt.figure(figsize=(12, 8))  # Adjust the figure size
sb.boxplot(data=cases_comparision_exio_small, palette="Set2")

# Add labels and title
plt.xlabel('Scenarios')
plt.ylabel('kg CO2eq')
plt.title('SE-Basic iron and steel and of ferro-alloys and first products thereof (exiobase_big)')
plt.gca().yaxis.set_major_formatter(formatter)

#Save the plot
filename = "Monte_Carlo_Comparison_between_cases.png"  # You can change the extension to .pdf, .jpg, etc.

# Create the folder if it doesn't exist
os.makedirs(folder, exist_ok=True)

# Save the plot
plt.savefig(os.path.join(folder, filename))


#%% - CASES COMPARISON EXIOBASE SMALL - Sector "DE-Biodiesels"

folder = r"C:\Users\Elisabetta\Desktop\Risultati_analisi_incertezza_Aalborg\Exiobase_big_sector_3"

file_1 = "CASE_1_uniform_0.1_MC_simulations_DE-Biodiesels.csv"
file_2 = "CASE_2_uniform_0.2_MC_simulations_DE-Biodiesels.csv"
file_3 = "CASE_3_uniform_0.3_MC_simulations_DE-Biodiesels.csv"
file_4 = "CASE_4_log-normal_1.1_MC_simulations_DE-Biodiesels.csv"
file_5 = "CASE_5_log-normal_1.2_MC_simulations_DE-Biodiesels.csv"
file_6 = "CASE_6_log-normal_1.3_MC_simulations_DE-Biodiesels.csv"

comp_1 = pd.read_csv(os.path.join(folder, file_1))
comp_2 = pd.read_csv(os.path.join(folder, file_2))
comp_3 = pd.read_csv(os.path.join(folder, file_3))
comp_4 = pd.read_csv(os.path.join(folder, file_4))
comp_5 = pd.read_csv(os.path.join(folder, file_5))
comp_6 = pd.read_csv(os.path.join(folder, file_6))

cases_comparision_exio_small =pd.DataFrame({
    'Unif, 10%': comp_1['kg CO2eq'],
    'Unif, 20%': comp_2['kg CO2eq'],
    'Unif, 30%': comp_3['kg CO2eq'],
    'log-normal, 1.1': comp_4['kg CO2eq'],
    'log-normal, 1.2': comp_5['kg CO2eq'],
    'log-normal, 1.3': comp_6['kg CO2eq']
})

def scientific_format(x, pos):
    return f'{x:.2e}'
formatter = FuncFormatter(scientific_format)

plt.figure(figsize=(12, 8))  # Adjust the figure size
sb.boxplot(data=cases_comparision_exio_small, palette="Set2")

# Add labels and title
plt.xlabel('Scenarios')
plt.ylabel('kg CO2eq')
plt.title('DE-Biodiesels (exiobase_big)')
plt.gca().yaxis.set_major_formatter(formatter)

#Save the plot
filename = "Monte_Carlo_Comparison_between_cases.png"  # You can change the extension to .pdf, .jpg, etc.

# Create the folder if it doesn't exist
os.makedirs(folder, exist_ok=True)

# Save the plot
plt.savefig(os.path.join(folder, filename))


#%% - CASES COMPARISON EXIOBASE SMALL - Sector "CN-Railway transportation services"

folder = r"C:\Users\Elisabetta\Desktop\Risultati_analisi_incertezza_Aalborg\Exiobase_big_sector_4"

file_1 = "CASE_1_uniform_0.1_MC_simulations_CN-Railway transportation services.csv"
file_2 = "CASE_2_uniform_0.2_MC_simulations_CN-Railway transportation services.csv"
file_3 = "CASE_3_uniform_0.3_MC_simulations_CN-Railway transportation services.csv"
file_4 = "CASE_4_log-normal_1.1_MC_simulations_CN-Railway transportation services.csv"
file_5 = "CASE_5_log-normal_1.2_MC_simulations_CN-Railway_transportation_services.csv"
file_6 = "CASE_6_log-normal_1.3_MC_simulations_CN-Railway_transportation_services.csv"

comp_1 = pd.read_csv(os.path.join(folder, file_1))
comp_2 = pd.read_csv(os.path.join(folder, file_2))
comp_3 = pd.read_csv(os.path.join(folder, file_3))
comp_4 = pd.read_csv(os.path.join(folder, file_4))
comp_5 = pd.read_csv(os.path.join(folder, file_5))
comp_6 = pd.read_csv(os.path.join(folder, file_6))

cases_comparision_exio_small =pd.DataFrame({
    'Unif, 10%': comp_1['kg CO2eq'],
    'Unif, 20%': comp_2['kg CO2eq'],
    'Unif, 30%': comp_3['kg CO2eq'],
    'log-normal, 1.1': comp_4['kg CO2eq'],
    'log-normal, 1.2': comp_5['kg CO2eq'],
    'log-normal, 1.3': comp_6['kg CO2eq']
})

def scientific_format(x, pos):
    return f'{x:.2e}'
formatter = FuncFormatter(scientific_format)

plt.figure(figsize=(12, 8))  # Adjust the figure size
sb.boxplot(data=cases_comparision_exio_small, palette="Set2")

# Add labels and title
plt.xlabel('Scenarios')
plt.ylabel('kg CO2eq')
plt.title('CN-Railway transportation services (exiobase_big)')
plt.gca().yaxis.set_major_formatter(formatter)

#Save the plot
filename = "Monte_Carlo_Comparison_between_cases.png"  # You can change the extension to .pdf, .jpg, etc.

# Create the folder if it doesn't exist
os.makedirs(folder, exist_ok=True)

# Save the plot
plt.savefig(os.path.join(folder, filename))


