# -*- coding: utf-8 -*-
"""
Created on Fri May  2 10:26:34 2025

@author: Elisabetta
"""

#%% --- Importing necessary packages ---
import bw2data as bd
import bw2calc as bc
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import ecoinvent_interface


bd.projects.set_current('UA_TSF_2024')
bd.databases

#bi.import_ecoinvent_release(
#    version='3.11',
#    system_model='consequential',
#    username ='XXX', # use your own
#    password='XXX' # use your own
#)

#bd.databases

#activity_name = "market group for heat, district or industrial, natural gas"
#for activity in bd.Database("ecoinvent-3.11-consequential").search(activity_name):  
#    print(activity)
#    print(activity['code'])

#mycode = '41094fad451a2d8d99cf69e998e81f86'
#myact = bd.Database("ecoinvent-3.11-consequential").get(mycode)

#print(myact['name'])


#%% --- Setting the working directory ---

os.getcwd()
os.chdir(r"C:\Users\Elisabetta\Desktop\RISULTATI AALBORG\20250509_Risultati_II_fase") 

from lci_to_bw2 import *

#%% --- Configuration and Setup ---

# Set project path and current project
PROJECT_PATH = r"C:\Users\Elisabetta\Desktop\RISULTATI AALBORG\20250509_Risultati_II_fase" 
PROJECT_NAME = 'UA_TSF_2024'
os.chdir(PROJECT_PATH)
bd.projects.set_current(PROJECT_NAME)


#%% --- Utility Functions ---

def descriptive_stats(data):
    """Return descriptive statistics of a given dataset."""
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

def scientific_format(x, pos):
    """Formatter function to display values in scientific notation."""
    return f'{x:.2e}'


def monte_carlo_simulation(functional_unit, method, iterations=500):
    """Run Monte Carlo simulation and return results."""
    mc = bc.LCA(functional_unit, method, use_distributions=True) # ATTENTION: you should add your uncertainty and then run this
    mc.lci()
    mc.lcia()
    return [mc.score for _ in zip(range(500), mc)]


def plot_monte_carlo_results(results, xlabel='lca.score'):
    """Plot histogram for Monte Carlo results."""
    plt.hist(results, density=True, bins=int(len(results)/15), color='lightblue', edgecolor='black')
    plt.ylabel("Probability")
    plt.xlabel(xlabel)
    plt.show()


def run_scenario_with_TSF():
    """Run LCA with the foreground system (TSF) using the provided foreground database."""
    
    # Load and register foreground system data
    mydb = pd.read_csv('TSF_inventory_matrix_ecoinvent_v311.csv', header=0, delimiter=";")
    bw2_db = lci_to_bw2(mydb)
    
    t_db = bd.Database('TSF')
    t_db.write(bw2_db)
    
    # Define functional unit and LCIA method
    functional_unit = {t_db.get('Sewage sludge treatment'): 100000000}
    method = ('ecoinvent-3.11','EF v3.1 no LT', 'climate change no LT', 'global warming potential (GWP100) no LT')
    
    # Run LCA
    lca = bc.LCA(functional_unit, method)
    lca.lci()
    lca.lcia()
    print("S2 LCA score (with TSF):", lca.score)
        
    # Monte Carlo simulation
    mc_results = monte_carlo_simulation(functional_unit, method)
    
    # Plotting the results
    # plot_monte_carlo_results(mc_results)
    
    return mc_results

#%%
if __name__ == "__main__":
    
    # S2: Test with the foreground system (TSF) and the pedigree matrix for the foreground
    mc_results_with_TSF = run_scenario_with_TSF()
    
    # Calculate descriptive statistics
    stats_with_TSF = descriptive_stats(mc_results_with_TSF)
    
    print("\nDescriptive statistics for scenario with TSF:")
    for stat_name, value in stats_with_TSF.items():
        print(f"{stat_name}: {value}")


#%% --- Export the results ---

output = pd.DataFrame(mc_results_with_TSF, columns=["kg CO2 eq"])
output.to_csv('20250804_mc_results_with_TSF.csv', index=False)
