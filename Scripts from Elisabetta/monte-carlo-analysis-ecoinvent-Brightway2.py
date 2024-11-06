# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 13:30:41 2024

@author: Elisabetta

SIMULATION WITH ECOINVENT IN BRIGHTWAY 2

  - S1: Test without the foreground system
  - S2: Test with the foreground system (TSF) and the pedigree matrix for the foreground

"""

#%% importing packages

import brightway2 as bw
import os
import pandas as pd
import numpy as np
from lci_to_bw2 import *
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter 

#%% S1: Test without the foreground system

electricity = bw.Database("ecoinvent-3.9-consequential").get("a435daf86e858232e75f32771d785f83")
functional_unit_noTSF = {electricity: 1000}
mymethod = ('EF v3.1 no LT',  'climate change no LT',  'global warming potential (GWP100) no LT')
lca_noTSF = bw.LCA(functional_unit_noTSF, mymethod)  # LCA calculations with method
lca_noTSF.lci()
lca_noTSF.lcia()

print(lca.score)

mc_noTSF = bw.MonteCarloLCA({electricity: 1000}, mymethod)  # Monte Carlo class
mc_results_noTSF = [next(mc) for x in range(500)]

# Look at the MC results
plt.hist(mc_results_noTSF, density=True)  # From matplotlib package. Use bins = int(500/15) to increase number of bars
plt.ylabel("Probability")
plt.xlabel('lca.score')


#%% S2: Test with the foreground system (TSF) and the pedigree matrix for the foreground

os.getcwd()
os.chdir(r"C:\Users\Elisabetta\OneDrive - Alma Mater Studiorum Università di Bologna\Documents\UNIVERSITA\Progetto di dottorato\Periodo_all_estero_2024\UA_Brightway")

bw.projects.report()
bw.projects.set_current('UA_IOT_mini')
bw.databases

mydb = pd.read_csv('TSF_inventory_matrix_ecoinvent.csv', header=[0], delimiter=";")
#mydb.head()

# Create a dict that can be written as database
bw2_db = lci_to_bw2(mydb) # a function from the lci_to_bw2 module
bw2_db

t_db = bw.Database('TSF') # it works because the database name in the excel file is the same
# shut down all other notebooks using the same project
t_db.write(bw2_db)

#Choose the LCIA method:
#list(bw.methods)

#Set the funztional unit and run the LCA:
functional_unit = {t_db.get('Sewage sludge treatment'): 1000}
mymethod = ('EF v3.1 no LT',  'climate change no LT',  'global warming potential (GWP100) no LT')
lca = bw.LCA(functional_unit, mymethod)  # LCA calculations with method
lca.lci()
lca.lcia()

print(lca.score)

mc = bw.MonteCarloLCA({t_db.get('Sewage sludge treatment'): 1000}, mymethod)  # Monte Carlo class
mc_results = [next(mc) for x in range(500)]

# Look at the MC results
plt.hist(mc_results, density=True)  # From matplotlib package. Use bins = int(500/15) to increase number of bars
plt.ylabel("Probability")
plt.xlabel('lca.score')

#%% Comparison between the uncertainty range of the results with and without TSF

def descriptive_stats(data):
    return {
        'count': np.count_nonzero(data),
        'mean': np.mean(data),
        'std_dev': np.std(data),
        'min': np.min(data),
        'max': np.max(data),
        '25th_percentile': np.percentile(data, 25),
        '50th_percentile': np.percentile(data, 50),  # Median
        '75th_percentile': np.percentile(data, 75),
    }

# Call the function
stats = descriptive_stats(mc_results)
stats_noTSF = descriptive_stats(mc_results_noTSF)

# Display the results
for stat_name, value in stats.items():
    print(f"{stat_name}: {value}")

for stat_name, value in stats_noTSF.items():
    print(f"{stat_name}: {value}")


comparison_unc_with_and_without_TSF = pd.DataFrame({
    'noTSF': mc_results_noTSF,
    'with_TSF': mc_results
})

def scientific_format(x, pos):
    return f'{x:.2e}'
formatter = FuncFormatter(scientific_format)

plt.figure(figsize=(12, 8))  # Adjust the figure size

# Creating the box plot
plt.boxplot([comparison_unc_with_and_without_TSF['noTSF'], comparison_unc_with_and_without_TSF['with_TSF']],
            labels=['noTSF', 'with_TSF'], 
            patch_artist=True,
            boxprops=dict(facecolor='lightblue', color='black'),
            medianprops=dict(color='red'),
            whiskerprops=dict(color='black'),
            capprops=dict(color='black'),
            flierprops=dict(marker='o', color='black', alpha=0.5))

# Add labels and title
plt.xlabel('Scenarios')
plt.ylabel('kg CO2eq')
plt.title('Uncertainty Range Comparison')

# Setting the y-axis to use scientific notation
plt.gca().yaxis.set_major_formatter(FuncFormatter(scientific_format))

# Show the plot
plt.tight_layout()
plt.show()