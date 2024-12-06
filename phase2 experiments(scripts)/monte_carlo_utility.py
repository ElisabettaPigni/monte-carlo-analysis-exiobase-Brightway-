import bw_processing as bwp
import bw2calc as bc
from scipy import sparse
import pandas as pd
import numpy as np
from random import sample
import os
import matplotlib.pyplot as plt
import seaborn as sb
from matplotlib.ticker import FuncFormatter
import textwrap
import re
import bw2data as bd
from typing import List, Tuple, Any


class SimulationScript:
    def __init__(self):
        """
        self.metadata: save technosphere data, including column index and gsd
        """
        self.metadata = []

    def file_preprocessing(self, file_name, delimiter: str, column_name: str, expacted_order: list):
        """
        Preprocess a file and return a DataFrame with the desired order.
        """
        df = pd.read_csv(file_name, delimiter=delimiter)
        df_sorted = df.set_index(column_name).reindex(expacted_order).reset_index()

        return df_sorted

    def get_index(self, activities: list, activity_name: str) -> int:
        """
        Get corresponding index for an activity.
        """
        index = activities.index(activity_name)
        return index
    
    def form_tech_matrix(self, raw_tech: pd.DataFrame):
        identity_matrix = np.identity(len(raw_tech))
        tech_matrix = - (identity_matrix - raw_tech)
        np.fill_diagonal(tech_matrix, -tech_matrix.diagonal())

        return tech_matrix

    def extend_matrix(self, original_matrix, extend_data: pd.DataFrame, names: list, is_technosphere=True):
        """
        Concatenage additional column and line to the matrix.
        * names: this is the list of activities or emissions.
        """
        if is_technosphere:
            row = np.zeros([original_matrix.shape[1]]).reshape(1, -1)
            column = np.zeros([original_matrix.shape[0]])
            for act, data in zip(extend_data.iloc[:, 0], extend_data.iloc[:, 1]):
                column[self.get_index(names, act)] = data
            column = np.nan_to_num(column, nan=0)
            column = np.insert(column, 0, 1)
            column = np.array([column]).T
            extended_matrix = np.concatenate((column, np.concatenate((row, original_matrix), axis=0)), axis=1)
        else:
            column = np.zeros([original_matrix.shape[0]]).reshape(1, -1).T
            extended_matrix = np.concatenate((column, original_matrix), axis=1)

        return extended_matrix
    
    def form_bio_matrix(self, bio_df, emissions) -> np.ndarray:
        """
        Get biosphere matrix data.
        """
        bio_matrix = bio_df.loc[emissions].to_numpy()

        return bio_matrix
    
    # ATTENTION: have to make sure the file has the same order as biosphere emissions.
    def form_cf_matrix(self, emission_file: str, method: tuple) -> pd.DataFrame:
        emission_code = pd.read_csv(emission_file, delimiter=",") 
        codes = emission_code.iloc[:, -1]

        bw_method = bd.Method(method)
        method_df = pd.DataFrame(bw_method.load(), columns=["database_code", "cf_number"])
        method_df[["database", "code"]] = method_df["database_code"].to_list()
        cf_selected = method_df[method_df["code"].isin(codes)][["code", "cf_number"]]
        cf_dict = cf_selected.set_index("code")["cf_number"].to_dict()
        missing_codes = list(set(codes.unique()) - set(cf_selected["code"]))
        
        cf_matrix = []
        if not missing_codes:
            for code in codes:
                cf_matrix.append(cf_dict.get(code))
        else:
            miss_dict = emission_code[["ecoinvent name", "brightway code"]].set_index("brightway code")["ecoinvent name"].to_dict()
            fixed_codes = []
            for code in missing_codes:
                name = miss_dict.get(code)
                if "Carbon dioxide" in name:
                    cf_dict[code] = 1.0
                    fixed_codes.append(code)
            if missing_codes != fixed_codes:
                print(f"CF data imcomplete, missing: {missing_codes}")
            else:
                for code in codes:
                    cf_matrix.append(cf_dict.get(code))

        cf_matrix = np.diagflat(cf_matrix)
            
        return cf_matrix

    def get_country_sector(self, activity):
        """
        Separate the country and sector.
        """
        country, sector = activity.split("-", 1)

        return country, sector

    def map_pedigree_uncertainty(self, country_file, sector_file, region_sector_file):
        """
        Build dictionaries to mapping specific uncertainty.
        """
        country_dfs = pd.read_csv(country_file, delimiter=";")
        sector_dfs = pd.read_csv(sector_file, delimiter=";")
        region_sector_dfs = pd.read_csv(region_sector_file, delimiter=";")

        country_region = country_dfs.set_index(country_dfs.columns[0])[country_dfs.columns[1]].to_dict()
        sector_seccat = sector_dfs.set_index(sector_dfs.columns[0])[sector_dfs.columns[1]].to_dict()

        return country_region, sector_seccat, region_sector_dfs

    def find_pedigree_uncertainty(self, activity, country_region, sector_seccat, region_sector_dfs):
        """
        Search for uncertainty for specific activity or biosphere flow.
        """
        country, sector = self.get_country_sector(activity)
        region_category =  country_region.get(country, None)
        sector_category = sector_seccat.get(sector, None)

        if region_category !=  None and sector_category != None:
            gsd = float(region_sector_dfs[(region_sector_dfs.iloc[:, 0] == region_category) & (region_sector_dfs.iloc[:, 1] == sector_category)]["GSD"].iloc[0])
        else:
            # print("No GSD found.")
            gsd = None

        return gsd
    
    def calc_specific_uncertainty(self, data, uncertainty):
        loc = np.log(data)
        scale = np.log(uncertainty)

        return loc, scale
    
    def get_uncertainty_negative(self, df, column_name):
        df
        [column_name]
    
    # TODO: design for pedigree and specific uncertainty
    def add_uncertainty(self, bw_data, bw_indices, bw_flip):
        bw_uncertainties = []
        if bw_flip is not None: # technosphere
            k = 0
            for data, indices, flip in zip(bw_data, bw_indices, bw_flip):
                uncertainty = list(self.metadata[indices[1]].items())[0][1][1]
                if uncertainty is not None:
                    if indices[1] == 0:
                        uncertainty = uncertainty[k]
                        k += 1
                        if uncertainty == 0:
                            parameters_a = (0, data, np.NaN, np.NaN, np.NaN, np.NaN, False)
                        else:
                            loc = np.log(data)
                            scale = np.log(uncertainty)
                            if not flip:
                                parameters_a = (0, data, np.NaN, np.NaN, np.NaN, np.NaN, False)
                            else:
                                parameters_a = (2, loc, scale, np.NaN, np.NaN, np.NaN, False)
                    else:
                        loc = np.log(data)
                        scale = np.log(uncertainty)
                        if not flip:
                            parameters_a = (0, data, np.NaN, np.NaN, np.NaN, np.NaN, False)
                        else:
                            parameters_a = (2, loc, scale, np.NaN, np.NaN, np.NaN, False)
                else:
                    parameters_a = (0, data, np.NaN, np.NaN, np.NaN, np.NaN, False)
                bw_uncertainties.append(parameters_a)
        else:
            k = 0
            for data, indices in zip(bw_data, bw_indices):
                uncertainty = list(self.metadata[indices[1]].items())[0][1][1]
                if uncertainty is not None:
                    if indices[1] == 0:
                        uncertainty = uncertainty[k]
                        k += 1
                        if uncertainty == 0:
                            parameters_a = (0, data, np.NaN, np.NaN, np.NaN, np.NaN, False)
                        else:
                            loc = np.log(data)
                            scale = np.log(uncertainty)
                            parameters_a = (2, loc, scale, np.NaN, np.NaN, np.NaN, False)
                    else:
                        loc = np.log(data)
                        scale = np.log(uncertainty)
                        parameters_a = (2, loc, scale, np.NaN, np.NaN, np.NaN, False)
                else:
                    parameters_a = (0, data, np.NaN, np.NaN, np.NaN, np.NaN, False)
                bw_uncertainties.append(parameters_a)

        return np.array(bw_uncertainties, dtype=bwp.UNCERTAINTY_DTYPE)

    
    def prepare_bw_matrix(self, tech_matrix, bio_matrix, cf_matrix, activities):
        """
        Transform matrices data to bw matrices data, ready for the datapackages.
        """
        tech_sparse = sparse.coo_array(tech_matrix)
        tech_coors = np.column_stack(tech_sparse.nonzero())
        
        max_coor = tech_coors[np.argmax(np.sum(tech_coors, axis=1))]
        tech_data = tech_sparse.data
        tech_indices = np.array([tuple(coor) for coor in tech_coors], dtype=bwp.INDICES_DTYPE)
        tech_flip = np.array([False if i[0] == i[1] else True for i in tech_indices])

        bio_sparse = sparse.coo_array(bio_matrix)
        bio_coors = np.column_stack(bio_sparse.nonzero())
        bio_data = bio_sparse.data
        bio_indices = np.array([tuple([coord[0] + max_coor[0] + 1, coord[1]]) for coord in bio_coors], dtype=bwp.INDICES_DTYPE)
        
        cf_sparse = sparse.coo_array(cf_matrix)
        cf_coors = np.column_stack(cf_sparse.nonzero())
        cf_data =  cf_sparse.data
        cf_indices = np.array([tuple([coord[0] + max_coor[0] + 1, coord[1] + max_coor[1] + 1]) for coord in cf_coors], dtype=bwp.INDICES_DTYPE)

        tech_poss = list(set(coord[1] for coord in tech_indices))
        for act, tech_pos in zip(activities, tech_poss):
            self.metadata.append({tech_pos: act})

        return [
            (tech_data, tech_indices, tech_flip),
            (bio_data, bio_indices),
            (cf_data, cf_indices)
        ]

    # TODO: uncertainty is None currently
    # Is it common people add ucnertainty to cf matrix?
    def prepare_datapackage(self, datapackage_data: List[Tuple[Any, ...]], uncertainty=None):
        tech_data, tech_indices, tech_flip = datapackage_data[0]
        bio_data, bio_indices = datapackage_data[1]
        cf_data, cf_indices = datapackage_data[2]
        tech_uncertainty, bio_uncertainty, cf_uncerainty = uncertainty[0], uncertainty[1], uncertainty[2]

        dp = bwp.create_datapackage()
        dp.add_persistent_vector(
            matrix='technosphere_matrix',
            indices_array=tech_indices,
            data_array=tech_data,
            flip_array=tech_flip,
            distributions_array=tech_uncertainty,
        )
        dp.add_persistent_vector(
            matrix='biosphere_matrix',
            indices_array=bio_indices,
            data_array=bio_data,
            distributions_array=bio_uncertainty,
        )
        dp.add_persistent_vector(
            matrix='characterization_matrix',
            indices_array=cf_indices,
            data_array=cf_data,
            distributions_array=cf_uncerainty,
        )

        return dp

    def perform_simulation(self, index, datapackage):
        lca = bc.LCA(
            demand={index: 1},
            data_objs=[datapackage],
        )
        lca.lci()
        lca.lcia()

        return lca.score

    def manual_lca(self, A, B, C, A_, index):
        f = np.zeros(len(A))
        f[index] = 1
        lca_score = np.sum(C.dot(B.dot((np.linalg.inv(A_)).dot(f))))

        return lca_score

    def sort_file(self, filename):
        """
        Sort files in a folder when file name has both string and number.
        """
        match = re.search(r'CASE_(\d+)_', filename)
        if match:
            return int(match.group(1))
        else:
            return float('inf')

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
                case_name = " ".join(file.split("_")[2:4])
                df["case"] = case_name if "static" not in case_name else "deterministic"
                match = re.search(r'simulations_(.*)\.csv', file)
                df["sector"] = match.group(1).replace("_", " ") if match else print("Sector name not founded.")
                data = pd.concat([data, df], ignore_index=True)

        print(f"Check cases: {data['case'].unique()}")
        print(f"Check row numbers: {len(data)}")
        print(f"Check column numbers: {len(data.columns)}")
        return data

    def draw_plot(self, data, compare_type, database_name, sector_names, save_path):
        font = {'size': 16}
        plt.rc('font', **font)
        if not os.path.exists(save_path):
            os.makedirs(save_path)

        def scientific_format(x, pos):
            return f'{x:.2e}'
        formatter = FuncFormatter(scientific_format)

        if compare_type == "cases": # means one plot includes all uncertainty cases for one sector
            for sector in data["sector"].unique():
                filtered_data = data[data["sector"] == sector].copy()
                plt.figure(figsize=(16, 14))
                plt.xlabel("Cases", labelpad=20)
                plt.ylabel("kg CO\u2082eq", labelpad=20)
                sb.boxplot(x=filtered_data["case"], y=filtered_data["kg CO2eq"], data=filtered_data, order=sector_names, hue="case", palette="Set2")
                plt_title = " ".join([sector, database_name])
                plt.title(plt_title, labelpad=20)
                plt_name = f"MC_{plt_title}_{compare_type}.png"
                plt.gca().yaxis.set_major_formatter(formatter)
                plt.savefig(os.path.join(save_path, plt_name))
                plt.close() # always remember to free memory
        elif compare_type == "sectors": # means one plot includes all sectors
            for case in data["case"].unique():
                filtered_data = data[data["case"] == case].copy()
                plt.figure(figsize=(16, 14))
                plt.xlabel("Activities", labelpad=20)
                plt.ylabel("kg CO\u2082eq", labelpad=20)
                sb.boxplot(x=filtered_data["sector"], y=filtered_data["kg CO2eq"], data=filtered_data, order=sector_names, hue="sector", palette="Set2")
                plt_title = " ".join([case, database_name])
                plt.title(plt_title)
                plt.xticks(
                    ticks=range(len(sector_names)), 
                    labels=["\n".join(textwrap.wrap(label, width=21)) for label in sector_names],
                )
                plt_name = f"MC_{plt_title}_{compare_type}.png"
                plt.gca().yaxis.set_major_formatter(formatter)
                plt.savefig(os.path.join(save_path, plt_name))
                plt.close() # always remember to free memory
