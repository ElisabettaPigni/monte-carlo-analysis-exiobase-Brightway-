# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 13:30:41 2024

@author: Elisabetta

SIMULATION WITH ECOINVENT IN BRIGHTWAY 2
- S1: Test without the foreground system
- S2: Test with the foreground system (TSF) and the pedigree matrix for the foreground
"""

# Importing necessary packages
import bw2data as bd
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from lci_to_bw2 import *


# --- Configuration and Setup ---

# Set project path and current project
PROJECT_PATH = r"C:\Users\Elisabetta\OneDrive - Alma Mater Studiorum Università di Bologna\Documents\UNIVERSITA\Progetto di dottorato\Periodo_all_estero_2024\UA_Brightway"
PROJECT_NAME = 'UA_IOT_mini'
os.chdir(PROJECT_PATH)
bd.projects.set_current(PROJECT_NAME)


#  --- Utility Functions ---

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
    mc = bd.LCA(functional_unit, method, use_distributions=True) # ATTENTION: you should add your uncertainty and then run this
    return [next(mc) for _ in range(iterations)]

def plot_monte_carlo_results(results, xlabel='lca.score'):
    """Plot histogram for Monte Carlo results."""
    plt.hist(results, density=True, bins=int(len(results)/15), color='lightblue', edgecolor='black')
    plt.ylabel("Probability")
    plt.xlabel(xlabel)
    plt.show()


def run_scenario_noTSF():
    """Run LCA without the foreground system (TSF)."""
    electricity = bd.Database("ecoinvent-3.9-consequential").get("a435daf86e858232e75f32771d785f83")
    functional_unit_noTSF = {electricity: 1000}
    mymethod = ('EF v3.1 no LT', 'climate change no LT', 'global warming potential (GWP100) no LT')
    
    lca_noTSF = bd.LCA(functional_unit_noTSF, mymethod)
    lca_noTSF.lci()
    lca_noTSF.lcia()
    
    print("S1 LCA score (no TSF):", lca_noTSF.score)
    
    # Monte Carlo simulation
    mc_results_noTSF = monte_carlo_simulation(functional_unit_noTSF, mymethod)
    
    # Plotting the results
    plot_monte_carlo_results(mc_results_noTSF)
    
    return mc_results_noTSF


def run_scenario_with_TSF():
    """Run LCA with the foreground system (TSF) using the provided foreground database."""
    
    # Load and register foreground system data
    mydb = pd.read_csv('TSF_inventory_matrix_ecoinvent.csv', header=0, delimiter=";")
    bw2_db = lci_to_bw2(mydb)
    
    t_db = bd.Database('TSF')
    t_db.write(bw2_db)
    
    # Define functional unit and LCIA method
    functional_unit = {t_db.get('Sewage sludge treatment'): 1000}
    mymethod = ('EF v3.1 no LT', 'climate change no LT', 'global warming potential (GWP100) no LT')
    
    # Run LCA
    lca = bd.LCA(functional_unit, mymethod)
    lca.lci()
    lca.lcia()
    print("S2 LCA score (with TSF):", lca.score)
    
    # Monte Carlo simulation
    mc_results = monte_carlo_simulation(functional_unit, mymethod)
    
    # Plotting the results
    plot_monte_carlo_results(mc_results)
    
    return mc_results


if __name__ == "__main__":

    # S1: Test without the foreground system
    mc_results_noTSF = run_scenario_noTSF()
    
    # S2: Test with the foreground system (TSF) and the pedigree matrix for the foreground
    mc_results_with_TSF = run_scenario_with_TSF()
    
    # Comparison between the uncertainty range of the results with and without TSF
    
    # Calculate descriptive statistics
    stats_noTSF = descriptive_stats(mc_results_noTSF)
    stats_with_TSF = descriptive_stats(mc_results_with_TSF)
    
    print("\nDescriptive statistics for scenario without TSF:")
    for stat_name, value in stats_noTSF.items():
        print(f"{stat_name}: {value}")
    
    print("\nDescriptive statistics for scenario with TSF:")
    for stat_name, value in stats_with_TSF.items():
        print(f"{stat_name}: {value}")
    
    # Prepare data for comparison plot
    comparison_unc_with_and_without_TSF = pd.DataFrame({
        'noTSF': mc_results_noTSF,
        'with_TSF': mc_results_with_TSF
    })
    
    # Plot comparison
    plt.figure(figsize=(12, 8))
    plt.boxplot(
        [comparison_unc_with_and_without_TSF['noTSF'], comparison_unc_with_and_without_TSF['with_TSF']],
        labels=['noTSF', 'with_TSF'], 
        patch_artist=True,
        boxprops=dict(facecolor='lightblue', color='black'),
        medianprops=dict(color='red'),
        whiskerprops=dict(color='black'),
        capprops=dict(color='black'),
        flierprops=dict(marker='o', color='black', alpha=0.5)
    )
    plt.xlabel('Scenarios')
    plt.ylabel('kg CO2eq')
    plt.title('Uncertainty Range Comparison')
    plt.gca().yaxis.set_major_formatter(FuncFormatter(scientific_format))
    plt.tight_layout()
    plt.show()