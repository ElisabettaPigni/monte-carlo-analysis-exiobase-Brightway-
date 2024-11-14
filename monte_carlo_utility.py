'''
 # @ Author: Ning An
 # @ Create Time: 2024-10-29 20:19:31
 # @ Contributor(s): Elisabetta Pigni
 # @ Modified by: Ning An
 # @ Modified time: 2024-11-06 10:01:48
 '''

import bw_processing as bwp
import bw2calc as bc
from scipy import sparse
import pandas as pd
import numpy as np
from random import sample
import os
from constants import *
import matplotlib.pyplot as plt
import seaborn as sb
from matplotlib.ticker import FuncFormatter
from constants import *
import textwrap
import re


class SimulationScript:
    def check(self, ):
        """
        Before run the simulation, check everything is ready.
        """
        pass

    def get_activities(self, A_file_path: str) -> list:
        """
        Form activities by combing <country_name> and <sector_name>.
        """
        A_raw = pd.read_csv(A_file_path, sep='\t', low_memory=False)
        countries = A_raw['region'].drop_duplicates().iloc[2:].tolist()
        sectors = list(A_raw.iloc[2:,1].drop_duplicates())
        activities = [ x + '-' + y for x in countries for y in sectors]
        return activities

    def choose_activities(self, activities: list, amount: int) -> list:
        """
        Choose activities randomly.
        """
        chosen_activities = []
        indices = sample(range(len(activities) - 1), amount)
        chosen_activities = [(activities[i], i) for i in indices]
        return chosen_activities
    
    def get_index(self, activities: list, activity_name: str) -> int:
        """
        Get corresponding index for an activity.
        """
        index = activities.index(activity_name)
        return index
    
    def build_bw_matrix_add(self, A_file_path, S_file_path):
        pass

    def build_bw_matrix(self, A_file_path, S_file_path):
        """
        Build matrix for brightway
        """
        activities = self.get_activities(A_file_path)

        A_raw = pd.read_csv(A_file_path, sep='\t', low_memory=False)
        A_IO = A_raw.iloc[2:,2:].astype('float').values
        I = np.identity(len(A_IO))
        A_ = I - A_IO
        A = -A_
        np.fill_diagonal(A, -A.diagonal())

        Asparse = sparse.coo_array(A)
        a_data = Asparse.data
        a_indices = np.array([tuple(coord) for coord in np.transpose(Asparse.nonzero())], dtype=bwp.INDICES_DTYPE)
        a_flip = np.array([False if i[0] == i[1] else True for i in a_indices ])

        S_raw = pd.read_csv(S_file_path, header=[0,1], index_col=[0], sep='\t', low_memory=False)
        GHG_rows = ["CO2 - combustion - air",
                    "CO2 - non combustion - Cement production - air",
                    "CO2 - non combustion - Lime production - air",
                    "CO2 - waste - biogenic - air", 
                    "CO2 - waste - fossil - air",
                    "CO2 - agriculture - peat decay - air", 
                    "CH4 - agriculture - air",
                    "CH4 - waste - air",
                    "CH4 - combustion - air",
                    "CH4 - non combustion - Extraction/production of (natural) gas - air",
                    "CH4 - non combustion - Extraction/production of crude oil - air",
                    "CH4 - non combustion - Mining of antracite - air",
                    "CH4 - non combustion - Mining of bituminous coal - air",
                    "CH4 - non combustion - Mining of coking coal - air",
                    "CH4 - non combustion - Mining of lignite (brown coal) - air",
                    "CH4 - non combustion - Mining of sub-bituminous coal - air",
                    "CH4 - non combustion - Oil refinery - air",
                    "N2O - combustion - air",
                    "N2O - agriculture - air",
                    "HFC - air",
                    "SF6 - air"]
        S = S_raw.loc[GHG_rows]

        B = S.to_numpy()
        Bsparse = sparse.coo_array(B)
        b_data = Bsparse.data
        b_indices_remap = [[i[0] + len(activities),i[1]] for i in np.transpose(Bsparse.nonzero())]
        b_indices = np.array([tuple(coord) for coord in b_indices_remap], dtype=bwp.INDICES_DTYPE)

        CFs = [1., 29.8, 273., 29.8, 29.8, 29.8, 29.8, 29.8, 29.8, 29.8, 29.8, 1., 1., 25200., 14600., 29.8, 1., 273., 29.8, 1., 1.]
        C = np.matrix(np.zeros((len(CFs), len(CFs))))
        C_diag = np.matrix(CFs)
        np.fill_diagonal(C, C_diag)

        Csparse = sparse.coo_array(C)
        c_data =  Csparse.data 
        c_indices_remap = [[i[1] + len(activities),i[1]+ len(activities)] for i in np.transpose(Csparse.nonzero())]
        c_indices = np.array([tuple(coord) for coord in c_indices_remap], dtype=bwp.INDICES_DTYPE) 
        return A, A_, A_IO, B, C, a_data, b_data, c_data, a_indices, b_indices, c_indices, a_flip
    
    def perform_baseline(self, index, a_data, b_data, c_data, a_indices, b_indices, c_indices, a_flip, A, A_, B, C, directory, k, t, myact):
        """
        Perform baseline simulation.
        """
        dp_static = bwp.create_datapackage()
        dp_static.add_persistent_vector(
            matrix='technosphere_matrix',
            indices_array=a_indices,
            data_array=a_data,
            flip_array=a_flip,
        )
        dp_static.add_persistent_vector(
            matrix='biosphere_matrix',
            indices_array=b_indices,
            data_array=b_data,
        )
        dp_static.add_persistent_vector(
            matrix='characterization_matrix',
            indices_array=c_indices,
            data_array=c_data,
        )

        lca = bc.LCA(
            demand={index: 1},
            data_objs=[dp_static],
        )
        lca.lci()
        lca.lcia()

        os.makedirs(directory, exist_ok=True)
        filename = os.path.join(directory, f"CASE_{k}_{t}_MC_simulations_{myact}.csv")

        print('ioscore',np.sum(C.dot(B.dot((np.linalg.inv(A_)).dot(f)))))
        print('LCA score: ', lca.score)  ## TODO: better to save the score in the result file too, so that we don't miss it.
        
        with open(filename, "w") as file:
            file.write("kg CO2eq\n") # Write the header
            file.write(f"{lca.score}")
            print(f"Baseline result saved to {filename}.")

    def add_uncertainty(self, t, u, a_data, b_data, c_data, a_indices, b_indices, c_indices, a_flip):
        """
        Add uncertainty to matrices.
        """
        if t == "uniform":
            results_a = []
            results_b = []
            
            min_val_a = a_data - (a_data * u)
            max_val_a = a_data + (a_data * u)

            for i in range(len(a_data)):
                if not a_flip[i]:
                    parameters_a = (0, a_data[i], np.NaN, np.NaN, np.NaN, np.NaN, False)
                else:
                    parameters_a = (4, np.NaN, np.NaN, np.NaN, min_val_a[i], max_val_a[i], False)
                results_a.append(parameters_a)
            
            a_uncertainty = np.array(results_a, 
                dtype=[('uncertainty_type', 'i4'), ('loc', 'f4'), ('scale', 'f4'), 
                    ('shape', 'f4'), ('minimum', 'f4'), ('maximum', 'f4'), ('negative', 'b')]
            )

            min_val_b = b_data - (b_data * u)
            max_val_b = b_data + (b_data * u)

            for i in range(len(b_data)):
                parameters_b = (4, np.NaN, np.NaN, np.NaN, min_val_b[i], max_val_b[i], False)
                results_b.append(parameters_b)
            
            b_uncertainty = np.array(results_b, 
                dtype=[('uncertainty_type', 'i4'), ('loc', 'f4'), ('scale', 'f4'), 
                    ('shape', 'f4'), ('minimum', 'f4'), ('maximum', 'f4'), ('negative', 'b')]
            )
        elif t == "log-normal":
            results_a = []
            results_b = []

            mu_a = np.log(a_data)
            sigma = np.log(u)
            
            for i in range(len(a_data)):
                if not a_flip[i]:
                    parameters_a = (0, a_data[i], np.NaN, np.NaN, np.NaN, np.NaN, False)
                else:
                    parameters_a = (2, mu_a[i], sigma, np.NaN, np.NaN, np.NaN, False)
                results_a.append(parameters_a)

            a_uncertainty = np.array(results_a, 
                dtype=[('uncertainty_type', 'i4'), ('loc', 'f4'), ('scale', 'f4'), 
                    ('shape', 'f4'), ('minimum', 'f4'), ('maximum', 'f4'), ('negative', 'b')]
            )

            mu_b = np.log(b_data)
            for i in range(len(b_data)):
                parameters_b = (2, mu_b[i], sigma, np.NaN, np.NaN, np.NaN, False)
                results_b.append(parameters_b)

            b_uncertainty = np.array(results_b, 
                dtype=[('uncertainty_type', 'i4'), ('loc', 'f4'), ('scale', 'f4'), 
                    ('shape', 'f4'), ('minimum', 'f4'), ('maximum', 'f4'), ('negative', 'b')]
            )

        dp_stochastic = bwp.create_datapackage()
        dp_stochastic.add_persistent_vector(
            matrix='technosphere_matrix',
            indices_array=a_indices,
            data_array=a_data,
            flip_array=a_flip,
            distributions_array=a_uncertainty,
        )
        dp_stochastic.add_persistent_vector(
            matrix='biosphere_matrix',
            indices_array=b_indices,
            data_array=b_data,
            distributions_array=b_uncertainty,
        )
        dp_stochastic.add_persistent_vector(
            matrix='characterization_matrix',
            indices_array=c_indices,
            data_array=c_data,
        )
        return dp_stochastic
    
    def perform_simu(self, index, dp_stochastic, directory, k, myact, t, u):
        """
        Perform Monte Carlo simulation and save the lca score.
        """
        lca = bc.LCA(
            demand={index: 1},
            data_objs=[dp_stochastic],
            use_distributions=True,
        )
        lca.lci()
        lca.lcia()

        print('LCA score: ', lca.score)

        os.makedirs(directory, exist_ok=True)
        filename = os.path.join(directory, f"CASE_{k}_{t}_{u}_MC_simulations_{myact}.csv")
        
        # TODO: these can be passed as parameters.
        batch_size = 50
        num_batches = 10

        with open(filename, "w") as file:
            file.write("kg CO2eq\n")
            for p in range(num_batches):
                batch_results = [lca.score for _ in zip(range(batch_size), lca)]
                df_batch = pd.DataFrame(batch_results, columns=["kg CO2eq"])
                df_batch.to_csv(file, header=False, index=False)
                print(f"Batch {p} saved to {filename}.")

        print(f"Results saved to {filename}.")

    def sort_file(self, file_list):
        """
        Sort files in a folder when file name has both string and number.
        """
        for file in file_list:
            match = re.search(r'\d+', file)
            return int(match.group()) if match else float('inf')

    def concate_files(self, folder_path):
        """
        Merge all columns from all files in a folder(exis=0).
        """
        data = pd.DataFrame()
        sorted_files = sorted(os.listdir(folder_path), key=self.sort_file)
        for file in sorted_files:
            if "CASE" in file:
                file_path = os.path.join(folder_path, file)
                df = pd.read_csv(file_path)
                old_column_name = "kg CO2eq"
                df.rename(columns={old_column_name: os.path.splitext(file)[0]}, inplace=True)
                data = pd.concat([data, df], axis=1)
            data.to_csv(f"{folder_path}/all_results.csv", index=False)

    def collect_data(self, folder_path, database_type):
        """
        Merge all columns from all files in a folder(axis=1), ready for plot drawing.
        Note: This function is used when you have a folder hierarchy.
        """
        data = pd.DataFrame()
        sorted_folders = sorted(os.listdir(folder_path), key=self.sort_file)
        for folder in sorted_folders:
            if database_type in folder:
                path = os.path.join(folder_path, folder)
                sorted_files = sorted(os.listdir(path), key=self.sort_file)
                for file in sorted_files:
                    if file.endswith(".csv"):
                        file_path = os.path.join(path, file)
                        print(f"Reading file: {file}")
                        df = pd.read_csv(file_path)
                        df["case"] = "_".join(file.split("_")[2:4])
                        df["sector"] = file.split("_")[-1].split(".")[0]
                        data = pd.concat([data, df], ignore_index=True)

        print(f"Check cases: {data['case'].unique()}")
        print(f"Check row numbers: {len(data)}")
        print(f"Check column numbers: {len(data.columns)}")
        return data

    def collect_data_direct(self, folder_path):
        """
        Merge all columns from all files in a folder(axis=1), ready for plot drawing.
        Note: This function is used when you have all files in one folder.
        """
        data = pd.DataFrame()
        sorted_folders = sorted(os.listdir(folder_path), key=self.sort_file)
        for file in sorted_folders:
            if file.endswith(".csv"):
                file_path = os.path.join(folder_path, file)
                df = pd.read_csv(file_path)
                df["case"] = "_".join(file.split("_")[2:4]) # TODO: maybe better to get it by regular expression. get the uncertainty type and number
                match = re.search(r'simulations_(.*)\.csv', file)
                df["sector"] = match.group(1) if match else print("Sector name not founded.")
                data = pd.concat([data, df], ignore_index=True)

        print(f"Check cases: {data['case'].unique()}")
        print(f"Check row numbers: {len(data)}")
        print(f"Check column numbers: {len(data.columns)}")
        return data

    def draw_plot(self, data, compare_type, database_name, sector_names, save_path):
        if not os.path.exists(save_path):
            os.makedirs(save_path)

        def scientific_format(x, pos):
            return f'{x:.2e}'
        formatter = FuncFormatter(scientific_format)

        if compare_type == "cases": # means one plot includes all uncertainty cases for one sector
            for sector in data["sector"].unique():
                filtered_data = data[data["sector"] == sector].copy()
                plt.figure(figsize=(12, 10))
                plt.xlabel("Scenarios")
                plt.ylabel("kg CO2eq")
                sb.boxplot(x=filtered_data["case"], y=filtered_data["kg CO2eq"], data=filtered_data, order=sector_names, hue="case", palette="Set2")
                plt_title = "_".join([sector, database_name])
                plt.title(plt_title)
                plt_name = f"MC_{plt_title}_{compare_type}.png"
                plt.gca().yaxis.set_major_formatter(formatter)
                plt.savefig(os.path.join(save_path, plt_name))
                plt.close() # always remember to free memory
        elif compare_type == "sectors": # means one plot includes all sectors
            for case in data["case"].unique():
                filtered_data = data[data["case"] == case].copy()
                plt.figure(figsize=(12, 10))
                plt.xlabel("Scenarios")
                plt.ylabel("kg CO2eq")
                sb.boxplot(x=filtered_data["sector"], y=filtered_data["kg CO2eq"], data=filtered_data, order=sector_names, hue="sector", palette="Set2")
                plt_title = "_".join([case, database_name])
                plt.title(plt_title)
                plt.xticks(
                    ticks=range(len(sector_names)), 
                    labels=["\n".join(textwrap.wrap(label, width=21)) for label in sector_names],
                )
                plt_name = f"MC_{plt_title}_{compare_type}.png"
                plt.gca().yaxis.set_major_formatter(formatter)
                plt.savefig(os.path.join(save_path, plt_name))
                plt.close() # always remember to free memory
