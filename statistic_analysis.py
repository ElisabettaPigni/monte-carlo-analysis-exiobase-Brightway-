import seaborn as sb #This is needed for the graphical rapresentation of the monte carlo results
import matplotlib.pyplot as plt #This is needed for the graphical rapresentation of the monte carlo results
from matplotlib.ticker import FuncFormatter #This is needed for the graphical rapresentation of the monte carlo results
import pandas as pd
import numpy as np
import os

class StatisticAnalysis:
    # Custom formatter function
    def scientific_format(self, x, pos):
        return f'{x:.2e}'


    def simu_plot(self):
        # Create a formatter object
        formatter = FuncFormatter(self.scientific_format)
        
        # Load the Monte Carlo results
        mc_results = pd.read_csv("/Users/bp45th/Documents/github/monte-carlo-analysis-exiobase-Brightway-/output/CASE_1_MC_simulations.csv", header=0)

        # Rename the column for clarity
        mc_results.columns = ["kg CO2eq"]

        # Save summary statistics to a CSV file
        summary_stats = mc_results.describe()
        summary_stats.to_csv("CASE_5_summary_statistics.csv")

        # Create and save the boxplot
        plt.figure(figsize=(10, 6))
        sb.boxplot(data=mc_results, y="kg CO2eq")
        plt.title('Climate change (kg CO2eq)')
        plt.gca().yaxis.set_major_formatter(formatter)
        plt.savefig("CASE_5_boxplot.png")
        plt.close()

        # Create and save the violin plot
        plt.figure(figsize=(10, 6))
        sb.violinplot(data=mc_results, y="kg CO2eq")
        plt.title('Climate change (kg CO2eq)')
        plt.gca().yaxis.set_major_formatter(formatter)
        plt.savefig("CASE_5_violinplot.png")
        plt.close()

        # Create and save the scatter plot
        plt.figure(figsize=(10, 6))
        sb.scatterplot(data=mc_results, y="kg CO2eq", x=mc_results.index)
        plt.title('Climate change (kg CO2eq)')
        plt.xlabel('MC Simulations')
        plt.gca().yaxis.set_major_formatter(formatter)
        plt.savefig("CASE_5_scatterplot.png")
        plt.close()

        print("Summary statistics and plots saved successfully.")


    def matrix_plot(self, A, A_IO, B):
        matrix_IO = A_IO[:8,:8]
        matrix_p = A[:8,:8]
        I_1 = np.eye(8)
        leontief_inverse = np.linalg.inv(I_1 - matrix_IO)
        inverse_matrix = np.linalg.inv(matrix_p)
        matrix_B = B[:8,:8] 

        fig, axs = plt.subplots(2, 3, figsize=(15, 10))

        sb.heatmap(matrix_IO, ax=axs[0, 0], annot=True, cbar=False)
        axs[0, 0].set_title('Input-Output Matrix (A_IO)')
        sb.heatmap(leontief_inverse, ax=axs[0, 2], annot=True, cbar=False)
        axs[0, 2].set_title('Leontief Inverse ((I-A_IO)^-1)')
        sb.heatmap(matrix_B, ax=axs[0, 1], annot=True, cbar=False)
        axs[0, 1].set_title('Direct stressor/impact coefficients (S)')

        sb.heatmap(matrix_p, ax=axs[1, 0], annot=True, cbar=False)
        axs[1, 0].set_title('Technological Matrix (A_p)')
        sb.heatmap(inverse_matrix, ax=axs[1, 2], annot=True, cbar=False)
        axs[1, 2].set_title('Inverse Matrix (A_p^-1)')
        sb.heatmap(matrix_B, ax=axs[1, 1], annot=True, cbar=False)
        axs[1, 1].set_title('Emission Factors (B)')

        plt.show()

