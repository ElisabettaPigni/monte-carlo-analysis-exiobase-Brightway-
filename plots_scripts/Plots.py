# -*- coding: utf-8 -*-
"""
Created on Wed Apr 16 10:37:41 2025

@author: Elisabetta
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import seaborn as sb
import textwrap
import re
import matplotlib.patches as mpatches

os.getcwd()
os.chdir(r"C:\Users\Elisabetta\Desktop\RISULTATI AALBORG\Plots_script_and_graphs_phase_I_and_II") 

#%% --- Phase 1 - Boxplot exiobase --- 

def draw_plot_exio(filename, compare_type="cases", database_name="", save_path="plots_phase1/"):
    # Load csv file
    df = pd.read_csv(filename)

    # Sector renaming map
    sector_rename = {
        "CN-Railway_transportation_services": "CN, Railway services",
        "DE-Biodiesels": "DE, Biodiesel",
        "CH-Beverages": "CH, Beverages",
        "SE-Basic_iron_and_steel_and_of_ferro-alloys_and_first_products_thereof": "SE, Basic iron",
        "DE-Paraffin_Waxes": "DE, Paraffin waxes",
        "RU-Food_waste_for_treatment__incineration": "RU, Food waste treatment",
        "NO-Other_services_(93)": "NO, Other services",
        "AT-Office_machinery_and_computers_(30)": "AT, Office machinery"
    }

    # Identify Monte Carlo columns
    nominal_columns = df.columns[:8]  # first 8 columns = nominal values
    mc_columns = df.columns[8:]      # remaining columns = MC simulations

    long_data = []

    # Add nominal values first
    for i, col in enumerate(nominal_columns):
        match = re.match(r"CASE_\d+_static_MC_simulations_(.+)", col)
        if match:
            sector = match.group(1)
            sector_clean = sector_rename.get(sector, sector.replace("_", " "))
            nominal_value = df[col].dropna().iloc[0] / 1000  # scale to tons
            long_data.append({
                "sector": sector_clean,
                "case": "nominal",
                "kg CO2eq": nominal_value
            })

    # Add MC values
    for col in mc_columns:
        match = re.match(r"CASE_\d+_(.+?)_MC_simulations_(.+)", col)
        if match:
            case_raw = match.group(1).replace("_", " ")
            sector = match.group(2)
            sector_clean = sector_rename.get(sector, sector.replace("_", " "))
            case = case_raw.strip()
            for value in df[col].dropna():
                long_data.append({
                    "sector": sector_clean,
                    "case": case,
                    "kg CO2eq": value / 1000  # scale to tons
                })

    # Convert to DataFrame
    plot_df = pd.DataFrame(long_data)

    # Create output directory
    os.makedirs(save_path, exist_ok=True)

    if compare_type == "cases":
        for sector in plot_df["sector"].unique():
            filtered_data = plot_df[plot_df["sector"] == sector]

            # Custom ordering: nominal, then uniform, then log-normal
            cases = list(filtered_data["case"].unique())
            uniform_cases = sorted([c for c in cases if "uniform" in c.lower()])
            lognormal_cases = sorted([c for c in cases if "log" in c.lower()])
            other_cases = sorted([c for c in cases if c.lower() not in [*uniform_cases, *lognormal_cases, "nominal"]])
            ordered_cases = ["nominal"] + uniform_cases + lognormal_cases + other_cases

            plt.figure(figsize=(16, 14))
            plt.xlabel("Cases", labelpad=20, fontsize=24)
            plt.ylabel("t CO₂eq / MEUR", labelpad=20, fontsize=24)
            sb.boxplot(
                x="case", y="kg CO2eq", data=filtered_data,
                order=ordered_cases,
                hue="case", palette="Set2"
            )
            plt.title(sector, fontsize=28)  # Removed database_name from title
            plt.xticks(rotation=45, ha='right', fontsize=20)  # Diagonal case names
            plt.yticks(fontsize=20)  # Regular (non-scientific) formatting
            plt_name = f"MC_{sector}_{database_name}_{compare_type}.png".replace("/", "_")
            plt.savefig(os.path.join(save_path, plt_name), bbox_inches='tight')
            plt.close()

    else:
        raise NotImplementedError("Only 'cases' comparison is currently implemented.")
        
# Create the plots for phase 1, different uncertainty cases for different exiobase sectors
draw_plot_exio(
    filename=r"C:\Users\Elisabetta\Desktop\RISULTATI AALBORG\Ning_all_experiments\20241206_EXIOBASE_phase1\all_results.csv",
    compare_type="cases",
    database_name="Phase_1_noTSF",
    save_path="20250807_exio_plots"
)        

#%% --- Phase 1 - boxplot exiobase aggregated ---

def draw_plot_exio_agg(filename, compare_type="cases", database_name="", save_path="plots_phase1/"):
    # Load csv file
    df = pd.read_csv(filename)

    # Sector renaming map
    sector_rename = {
        "RoW-Services": "RoW, Services",
        "EU28-Biodiesels": "EU28, Biodiesels",
        "EU28-Agriculture-Forestry-Fishing": "EU28, AFOLU",
        "EU28-Basic_iron_and_steel_and_of_ferro-alloys_and_first_products_thereof": "EU28, Basic iron",
        "EU28-Energy": "EU28, Energy",
        "RoW-Waste_management": "RoW, Waste management",
        "EU28-Services": "EU28, Services",
        "EU28-Industry": "EU28, Industry"
    }

    # Identify Monte Carlo columns
    nominal_columns = df.columns[:8]  # first 8 columns = nominal values
    mc_columns = df.columns[8:]      # remaining columns = MC simulations

    long_data = []

    # Add nominal values first
    for i, col in enumerate(nominal_columns):
        match = re.match(r"CASE_\d+_static_MC_simulations_(.+)", col)
        if match:
            sector = match.group(1)
            sector_clean = sector_rename.get(sector, sector.replace("_", " "))
            nominal_value = df[col].dropna().iloc[0] / 1000  # scale to tons
            long_data.append({
                "sector": sector_clean,
                "case": "nominal",
                "kg CO2eq": nominal_value
            })

    # Add MC values
    for col in mc_columns:
        match = re.match(r"CASE_\d+_(.+?)_MC_simulations_(.+)", col)
        if match:
            case_raw = match.group(1).replace("_", " ")
            sector = match.group(2)
            sector_clean = sector_rename.get(sector, sector.replace("_", " "))
            case = case_raw.strip()
            for value in df[col].dropna():
                long_data.append({
                    "sector": sector_clean,
                    "case": case,
                    "kg CO2eq": value / 1000  # scale to tons
                })

    # Convert to DataFrame
    plot_df = pd.DataFrame(long_data)

    # Create output directory
    os.makedirs(save_path, exist_ok=True)

    if compare_type == "cases":
        for sector in plot_df["sector"].unique():
            filtered_data = plot_df[plot_df["sector"] == sector]

            # Custom ordering: nominal, then uniform, then log-normal
            cases = list(filtered_data["case"].unique())
            uniform_cases = sorted([c for c in cases if "uniform" in c.lower()])
            lognormal_cases = sorted([c for c in cases if "log" in c.lower()])
            other_cases = sorted([c for c in cases if c.lower() not in [*uniform_cases, *lognormal_cases, "nominal"]])
            ordered_cases = ["nominal"] + uniform_cases + lognormal_cases + other_cases

            plt.figure(figsize=(16, 14))
            plt.xlabel("Cases", labelpad=20, fontsize=24)
            plt.ylabel("t CO₂eq / MEUR", labelpad=20, fontsize=24)
            sb.boxplot(
                x="case", y="kg CO2eq", data=filtered_data,
                order=ordered_cases,
                hue="case", palette="Set2"
            )
            plt.title(sector, fontsize=28)  # Removed database_name from title
            plt.xticks(rotation=45, ha='right', fontsize=20)  # Diagonal case names
            plt.yticks(fontsize=20)  # Regular (non-scientific) formatting
            plt_name = f"MC_{sector}_{database_name}_{compare_type}.png".replace("/", "_")
            plt.savefig(os.path.join(save_path, plt_name), bbox_inches='tight')
            plt.close()

    else:
        raise NotImplementedError("Only 'cases' comparison is currently implemented.")

draw_plot_exio_agg(
    filename=r"C:\Users\Elisabetta\Desktop\RISULTATI AALBORG\Ning_all_experiments\20241206_EXIOBASE_AGGREGATED_phase1\all_results.csv",
    compare_type="cases",
    database_name="Phase_1_noTSF",
    save_path="20250807_exio_agg_plots"
)  


#%% --- Sectors comparison ---

def combine_sector_plots_multi(
    file_sector_map,  # Dict: filename -> list of sector names
    database_name="",
    save_path="PIPPOLO/",
    compare_type="cases"
):
    import pandas as pd
    import seaborn as sb
    import matplotlib.pyplot as plt
    import os
    import re

    # Renaming map
    sector_rename = {
        "CN-Railway_transportation_services": "CN, Railway services",
        "DE-Biodiesels": "DE, Biodiesel",
        "CH-Beverages": "CH, Beverages",
        "SE-Basic_iron_and_steel_and_of_ferro-alloys_and_first_products_thereof": "SE, Basic iron",
        "DE-Paraffin_Waxes": "DE, Paraffin waxes",
        "RU-Food_waste_for_treatment__incineration": "RU, Food waste treatment",
        "NO-Other_services_(93)": "NO, Other services",
        "AT-Office_machinery_and_computers_(30)": "AT, Office machinery",
        "RoW-Services": "RoW, Services",
        "EU28-Biodiesels": "EU28, Biodiesels",
        "EU28-Agriculture-Forestry-Fishing": "EU28, AFOLU",
        "EU28-Basic_iron_and_steel_and_of_ferro-alloys_and_first_products_thereof": "EU28, Basic iron",
        "EU28-Energy": "EU28, Energy",
        "RoW-Waste_management": "RoW, Waste management",
        "EU28-Services": "EU28, Services",
        "EU28-Industry": "EU28, Industry"
    }

    long_data = []

    for filename, sector_names in file_sector_map.items():
        df = pd.read_csv(filename)
        nominal_columns = df.columns[:8]
        mc_columns = df.columns[8:]

        for i, col in enumerate(nominal_columns):
            match = re.match(r"CASE_\d+_static_MC_simulations_(.+)", col)
            if match:
                sector = match.group(1)
                sector_clean = sector_rename.get(sector, sector.replace("_", " "))
                if sector_clean in sector_names:
                    nominal_value = df[col].dropna().iloc[0] / 1000
                    long_data.append({
                        "sector": sector_clean,
                        "case": "nominal",
                        "kg CO2eq": nominal_value
                    })

        for col in mc_columns:
            match = re.match(r"CASE_\d+_(.+?)_MC_simulations_(.+)", col)
            if match:
                case_raw = match.group(1).replace("_", " ")
                sector = match.group(2)
                sector_clean = sector_rename.get(sector, sector.replace("_", " "))
                for sector_filter in sector_names:
                    if sector_clean == sector_filter:
                        for value in df[col].dropna():
                            long_data.append({
                                "sector": sector_clean,
                                "case": case_raw.strip(),
                                "kg CO2eq": value / 1000
                            })

    # Convert to DataFrame
    plot_df = pd.DataFrame(long_data)
    os.makedirs(save_path, exist_ok=True)

    # Maintain sector order as provided in file_sector_map
    all_sectors = []
    for _, sector_names in file_sector_map.items():
        for name in sector_names:
            if name not in all_sectors:
                all_sectors.append(name)

    num_sectors = len(all_sectors)

    if num_sectors == 0:
        raise ValueError("No data matched. Please check sector names and filenames.")

    fig, axs = plt.subplots(1, num_sectors, figsize=(10 * num_sectors, 8), sharey=True)
    if num_sectors == 1:
        axs = [axs]

    global_ymin = plot_df["kg CO2eq"].min()
    global_ymax = plot_df["kg CO2eq"].max()

    for ax, sector in zip(axs, all_sectors):
        filtered_data = plot_df[plot_df["sector"] == sector]

        cases = list(filtered_data["case"].unique())
        uniform_cases = sorted([c for c in cases if "uniform" in c.lower()])
        lognormal_cases = sorted([c for c in cases if "log" in c.lower()])
        other_cases = sorted([c for c in cases if c.lower() not in [*uniform_cases, *lognormal_cases, "nominal"]])
        ordered_cases = ["nominal"] + uniform_cases + lognormal_cases + other_cases

        sb.boxplot(
            x="case", y="kg CO2eq", data=filtered_data,
            order=ordered_cases,
            ax=ax, palette="Set2"
        )
        ax.set_title(sector, fontsize=20)
        ax.set_xlabel("Cases", fontsize=16)
        ax.set_ylabel("t CO₂eq / MEUR", fontsize=16)
        ax.tick_params(axis='x', rotation=45, labelsize=16)
        ax.tick_params(axis='y', labelsize=16)

    for ax in axs:
        ax.set_ylim(global_ymin, global_ymax)

    plt.tight_layout()
    out_filename = f"Combined_MC_{'_'.join([s.replace(', ', '_') for s in all_sectors])}_{database_name}_{compare_type}.png"
    plt.savefig(os.path.join(save_path, out_filename), bbox_inches="tight")
    plt.close()

    print(f"✅ Combined plot saved as: {os.path.join(save_path, out_filename)}")


# Compare the biodiesel sector between exio and exio aggregated
combine_sector_plots_multi(
    file_sector_map={
        "exiobase_all_results.csv": ["DE, Biodiesel"],
        "exiobase_agg_all_results.csv": ["EU28, Biodiesels"],
    },
    database_name="Phase_1_noTSF",
    save_path="20250807_combined_plots"
)

# Compare the industry sector between exio and exio aggregated
combine_sector_plots_multi(
    file_sector_map={
        "exiobase_all_results.csv": ["AT, Office machinery"],
        "exiobase_agg_all_results.csv": ["EU28, Industry"],
    },
    database_name="Phase_1_noTSF",
    save_path="20250807_combined_plots"
)

# Compare the AFOLU sector between exio and exio aggregated
combine_sector_plots_multi(
    file_sector_map={
        "exiobase_all_results.csv": ["CH, Beverages"],
        "exiobase_agg_all_results.csv": ["EU28, AFOLU"],
    },
    database_name="Phase_1_noTSF",
    save_path="20250807_combined_plots"
)

# Compare the service sector (CN vs RoW)
combine_sector_plots_multi(
    file_sector_map={
        "exiobase_all_results.csv": ["CN, Railway services"],
        "exiobase_agg_all_results.csv": ["RoW, Services"],
    },
    database_name="Phase_1_noTSF",
    save_path="20250807_combined_plots"
)

# Compare the basic iron sector
combine_sector_plots_multi(
    file_sector_map={
        "exiobase_all_results.csv": ["SE, Basic iron"],
        "exiobase_agg_all_results.csv": ["EU28, Basic iron"],
    },
    database_name="Phase_1_noTSF",
    save_path="20250807_combined_plots"
)

# Compare food waste management
combine_sector_plots_multi(
    file_sector_map={
        "exiobase_all_results.csv": ["RU, Food waste treatment"],
        "exiobase_agg_all_results.csv": ["RoW, Waste management"],
    },
    database_name="Phase_1_noTSF",
    save_path="20250807_combined_plots"
)

# Compare the energy sector
combine_sector_plots_multi(
    file_sector_map={
        "exiobase_all_results.csv": ["DE, Paraffin waxes"],
        "exiobase_agg_all_results.csv": ["EU28, Energy"],
    },
    database_name="Phase_1_noTSF",
    save_path="20250807_combined_plots"
)

# Compare services (NO vs EU28)
combine_sector_plots_multi(
    file_sector_map={
        "exiobase_all_results.csv": ["NO, Other services"],
        "exiobase_agg_all_results.csv": ["EU28, Services"],
    },
    database_name="Phase_1_noTSF",
    save_path="20250807_combined_plots"
)



#%% --- Phase 2 - Boxplot: exiobase, exiobase agg, ecoinvent with the foreground system ---


filepath_exio_agg_nominal = r"C:\Users\Elisabetta\Desktop\RISULTATI AALBORG\Ning_all_experiments\20250807_EXIOBASE AGGREGATED_phase2\CASE_1_static_MC_simulations_exrta_column.csv"
filepath_exio_agg = r"C:\Users\Elisabetta\Desktop\RISULTATI AALBORG\Ning_all_experiments\20250807_EXIOBASE AGGREGATED_phase2\CASE_2_exrta_column_MC_simulations_pedigree.csv"

filepath_exio_nominal = r"C:\Users\Elisabetta\Desktop\RISULTATI AALBORG\Ning_all_experiments\20250807_EXIOBASE_phase2\CASE_1_static_MC_simulations_extra_column.csv"
filepath_exio = r"C:\Users\Elisabetta\Desktop\RISULTATI AALBORG\Ning_all_experiments\20250807_EXIOBASE_phase2\CASE_2_extra_column_MC_simulations_pedigree.csv"

filepath_ecoinvent_nominal = r"C:\Users\Elisabetta\Desktop\RISULTATI AALBORG\Ecoinvent_script_and_mc_results\20250806_static_results_with_TSF.csv"
filepath_ecoinvent = r"C:\Users\Elisabetta\Desktop\RISULTATI AALBORG\Ecoinvent_script_and_mc_results\20250806_mc_results_with_TSF.csv"

exio_nominal = pd.read_csv(filepath_exio_nominal).rename(columns={"kg CO2eq": "exio nominal"})
exio = pd.read_csv(filepath_exio).rename(columns={"kg CO2eq": "exio"})

exio_agg_nominal = pd.read_csv(filepath_exio_agg_nominal).rename(columns={"kg CO2eq": "exio agg nominal"})
exio_agg = pd.read_csv(filepath_exio_agg).rename(columns={"kg CO2eq": "exio agg"})

ecoinvent_nominal = pd.read_csv(filepath_ecoinvent_nominal).rename(columns={"kg CO2 eq": "ecoinvent nominal"})
ecoinvent = pd.read_csv(filepath_ecoinvent).rename(columns={"kg CO2 eq": "ecoinvent"})

phase2_results = pd.concat([exio_agg_nominal, exio_nominal, ecoinvent_nominal, exio_agg, exio, ecoinvent], axis=1)

def draw_phase2_boxplot(dataframe, save_path="phase2/"):
    # Scale to Million tons
    df_Mtons = dataframe / 1e9

    # Separate nominal and MC columns
    nominal_columns = df_Mtons.columns[:3]
    mc_columns = df_Mtons.columns[3:]
    column_order = list(nominal_columns) + list(mc_columns)  # for consistent x-axis order

    # Prepare long-form DataFrame for seaborn
    long_data = []

    # Add nominal values (1 value per column)
    for col in nominal_columns:
        value = df_Mtons[col].dropna().iloc[0]
        long_data.append({
            "Database": col,
            "Value": value,
            "Type": "nominal"
        })

    # Add Monte Carlo simulation values
    for col in mc_columns:
        for val in df_Mtons[col].dropna():
            long_data.append({
                "Database": col,
                "Value": val,
                "Type": "Monte Carlo"
            })

    plot_df = pd.DataFrame(long_data)

    # Create plot
    plt.figure(figsize=(12, 8))

    # Boxplot for Monte Carlo values
    sb.boxplot(
        data=plot_df[plot_df["Type"] == "Monte Carlo"],
        x="Database", y="Value", hue="Database", palette="Set2",
        order=column_order, width=0.6, dodge=False
    )
    
    legend = plt.gca().get_legend()
    if legend:
        legend.remove()

    # Overlay nominal values (black diamonds)
    nominal_df = plot_df[plot_df["Type"] == "nominal"]
    sb.stripplot(
        data=nominal_df, x="Database", y="Value",
        size=10, color='black', marker='D',
        order=column_order
    )

    # Aesthetics
    plt.ylabel("Mt CO₂eq", fontsize=20, labelpad=20)
    plt.xticks(fontsize=16)
    new_labels = [label.replace(" ", "\n") for label in column_order]
    plt.gca().set_xticklabels(new_labels, fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel("")  # No x-axis label

    # Turn off scientific notation
    plt.gca().yaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x:,.2f}'))

    # Save plot
    os.makedirs(save_path, exist_ok=True)
    plt.tight_layout()
    plt.savefig(os.path.join(save_path, "phase2_boxplot_PROVA2.png"), bbox_inches='tight')
    plt.close()

# Example usage
draw_phase2_boxplot(phase2_results, save_path="20250807_phase2/")


#%% --- descriptive statistics of all scenarios and sectors ---

exio_NO_TSF = pd.read_csv("exiobase_all_results.csv")
exio_agg_NO_TSF = pd.read_csv("exiobase_agg_all_results.csv")

# Calculate the natural logarithm for the log-normal results
lognorm_1 = np.log(exio_NO_TSF.iloc[:, 32:56])
lognorm_2 = np.log(exio_agg_NO_TSF.iloc[:, 32:56])

unif_1 = exio_NO_TSF.iloc[:, 0:32]
unif_2 = exio_agg_NO_TSF.iloc[:, 0:32]

exiobase_without_TSF = pd.concat([unif_1, lognorm_1], axis=1)
exiobase_agg_without_TSF = pd.concat([unif_2, lognorm_2], axis=1)

def descriptive_stats(data, label):
    """Return descriptive statistics including IQR for a given dataset."""
    q1 = np.percentile(data, 25)
    q3 = np.percentile(data, 75)
    stats = {
        'Dataset': label,
        'count': np.count_nonzero(data),
        'mean': np.mean(data),
        'std_dev': np.std(data),
        'min': np.min(data),
        'max': np.max(data),
        '25th_percentile': q1,
        '50th_percentile': np.percentile(data, 50),  # Median
        '75th_percentile': q3,
        'IQR': q3 - q1
    }
    return stats


# Collect stats for each dataset
stats_list = [
    descriptive_stats(exio, "exio"),
    descriptive_stats(exio_agg, "exio_agg"),
    descriptive_stats(ecoinvent, "ecoinvent"),
    descriptive_stats(exiobase_without_TSF.iloc[:,8], "uniform_0.1_CN, Railway services"),
    descriptive_stats(exiobase_without_TSF.iloc[:,9], "uniform_0.1_DE, Biodiesel"),
    descriptive_stats(exiobase_without_TSF.iloc[:,10], "uniform_0.1_CH, Beverages"),
    descriptive_stats(exiobase_without_TSF.iloc[:,11], "uniform_0.1_SE, Basic iron"),
    descriptive_stats(exiobase_without_TSF.iloc[:,12], "uniform_0.1_DE, Paraffin waxes"),
    descriptive_stats(exiobase_without_TSF.iloc[:,13], "uniform_0.1_RU, Food waste treatment"),
    descriptive_stats(exiobase_without_TSF.iloc[:,14], "uniform_0.1_NO, Other services"),
    descriptive_stats(exiobase_without_TSF.iloc[:,15], "uniform_0.1_AT, Office machinery"),
    descriptive_stats(exiobase_without_TSF.iloc[:,16], "uniform_0.2_CN, Railway services"),
    descriptive_stats(exiobase_without_TSF.iloc[:,17], "uniform_0.2_DE, Biodiesel"),
    descriptive_stats(exiobase_without_TSF.iloc[:,18], "uniform_0.2_CH, Beverages"),
    descriptive_stats(exiobase_without_TSF.iloc[:,19], "uniform_0.2_SE, Basic iron"),
    descriptive_stats(exiobase_without_TSF.iloc[:,20], "uniform_0.2_DE, Paraffin waxes"),
    descriptive_stats(exiobase_without_TSF.iloc[:,21], "uniform_0.2_RU, Food waste treatment"),
    descriptive_stats(exiobase_without_TSF.iloc[:,22], "uniform_0.2_NO, Other services"),
    descriptive_stats(exiobase_without_TSF.iloc[:,23], "uniform_0.2_AT, Office machinery"),
    descriptive_stats(exiobase_without_TSF.iloc[:,24], "uniform_0.3_CN, Railway services"),
    descriptive_stats(exiobase_without_TSF.iloc[:,25], "uniform_0.3_DE, Biodiesel"),
    descriptive_stats(exiobase_without_TSF.iloc[:,26], "uniform_0.3_CH, Beverages"),
    descriptive_stats(exiobase_without_TSF.iloc[:,27], "uniform_0.3_SE, Basic iron"),
    descriptive_stats(exiobase_without_TSF.iloc[:,28], "uniform_0.3_DE, Paraffin waxes"),
    descriptive_stats(exiobase_without_TSF.iloc[:,29], "uniform_0.3_RU, Food waste treatment"),
    descriptive_stats(exiobase_without_TSF.iloc[:,30], "uniform_0.3_NO, Other services"),
    descriptive_stats(exiobase_without_TSF.iloc[:,31], "uniform_0.3_AT, Office machinery"),
    descriptive_stats(exiobase_without_TSF.iloc[:,32], "lognorm_1.106_CN, Railway services"),
    descriptive_stats(exiobase_without_TSF.iloc[:,33], "lognorm_1.106_DE, Biodiesel"),
    descriptive_stats(exiobase_without_TSF.iloc[:,34], "lognorm_1.106_CH, Beverages"),
    descriptive_stats(exiobase_without_TSF.iloc[:,35], "lognorm_1.106_SE, Basic iron"),
    descriptive_stats(exiobase_without_TSF.iloc[:,36], "lognorm_1.106_DE, Paraffin waxes"),
    descriptive_stats(exiobase_without_TSF.iloc[:,37], "lognorm_1.106_RU, Food waste treatment"),
    descriptive_stats(exiobase_without_TSF.iloc[:,38], "lognorm_1.106_NO, Other services"),
    descriptive_stats(exiobase_without_TSF.iloc[:,39], "lognorm_1.106_AT, Office machinery"),
    descriptive_stats(exiobase_without_TSF.iloc[:,40], "lognorm_1.225_CN, Railway services"),
    descriptive_stats(exiobase_without_TSF.iloc[:,41], "lognorm_1.225_DE, Biodiesel"),
    descriptive_stats(exiobase_without_TSF.iloc[:,42], "lognorm_1.225_CH, Beverages"),
    descriptive_stats(exiobase_without_TSF.iloc[:,43], "lognorm_1.225_SE, Basic iron"),
    descriptive_stats(exiobase_without_TSF.iloc[:,44], "lognorm_1.225_DE, Paraffin waxes"),
    descriptive_stats(exiobase_without_TSF.iloc[:,45], "lognorm_1.225_RU, Food waste treatment"),
    descriptive_stats(exiobase_without_TSF.iloc[:,46], "lognorm_1.225_NO, Other services"),
    descriptive_stats(exiobase_without_TSF.iloc[:,47], "lognorm_1.225_AT, Office machinery"),
    descriptive_stats(exiobase_without_TSF.iloc[:,48], "lognorm_1.363_CN, Railway services"),
    descriptive_stats(exiobase_without_TSF.iloc[:,49], "lognorm_1.363_DE, Biodiesel"),
    descriptive_stats(exiobase_without_TSF.iloc[:,50], "lognorm_1.363_CH, Beverages"),
    descriptive_stats(exiobase_without_TSF.iloc[:,51], "lognorm_1.363_SE, Basic iron"),
    descriptive_stats(exiobase_without_TSF.iloc[:,52], "lognorm_1.363_DE, Paraffin waxes"),
    descriptive_stats(exiobase_without_TSF.iloc[:,53], "lognorm_1.363_RU, Food waste treatment"),
    descriptive_stats(exiobase_without_TSF.iloc[:,54], "lognorm_1.363_NO, Other services"),
    descriptive_stats(exiobase_without_TSF.iloc[:,55], "lognorm_1.363_AT, Office machinery"),
    descriptive_stats(exiobase_agg_without_TSF.iloc[:,8], "uniform_0.1_RoW, Services"),
    descriptive_stats(exiobase_agg_without_TSF.iloc[:,9], "uniform_0.1_EU28, Biodiesels"),
    descriptive_stats(exiobase_agg_without_TSF.iloc[:,10], "uniform_0.1_EU28, AFOLU"),
    descriptive_stats(exiobase_agg_without_TSF.iloc[:,11], "uniform_0.1_EU28, Basic iron"),
    descriptive_stats(exiobase_agg_without_TSF.iloc[:,12], "uniform_0.1_EU28, Energy"),
    descriptive_stats(exiobase_agg_without_TSF.iloc[:,13], "uniform_0.1_RoW, Waste management"),
    descriptive_stats(exiobase_agg_without_TSF.iloc[:,14], "uniform_0.1_EU28, Services"),
    descriptive_stats(exiobase_agg_without_TSF.iloc[:,15], "uniform_0.1_EU28, Industry"),
    descriptive_stats(exiobase_agg_without_TSF.iloc[:,16], "uniform_0.2_RoW, Services"),
    descriptive_stats(exiobase_agg_without_TSF.iloc[:,17], "uniform_0.2_EU28, Biodiesels"),
    descriptive_stats(exiobase_agg_without_TSF.iloc[:,18], "uniform_0.2_EU28, AFOLU"),
    descriptive_stats(exiobase_agg_without_TSF.iloc[:,19], "uniform_0.2_EU28, Basic iron"),
    descriptive_stats(exiobase_agg_without_TSF.iloc[:,20], "uniform_0.2_EU28, Energy"),
    descriptive_stats(exiobase_agg_without_TSF.iloc[:,21], "uniform_0.2_RoW, Waste management"),
    descriptive_stats(exiobase_agg_without_TSF.iloc[:,22], "uniform_0.2_EU28, Services"),
    descriptive_stats(exiobase_agg_without_TSF.iloc[:,23], "uniform_0.2_EU28, Industry"),
    descriptive_stats(exiobase_agg_without_TSF.iloc[:,24], "uniform_0.3_RoW, Services"),
    descriptive_stats(exiobase_agg_without_TSF.iloc[:,25], "uniform_0.3_EU28, Biodiesels"),
    descriptive_stats(exiobase_agg_without_TSF.iloc[:,26], "uniform_0.3_EU28, AFOLU"),
    descriptive_stats(exiobase_agg_without_TSF.iloc[:,27], "uniform_0.3_EU28, Basic iron"),
    descriptive_stats(exiobase_agg_without_TSF.iloc[:,28], "uniform_0.3_EU28, Energy"),
    descriptive_stats(exiobase_agg_without_TSF.iloc[:,29], "uniform_0.3_RoW, Waste management"),
    descriptive_stats(exiobase_agg_without_TSF.iloc[:,30], "uniform_0.3_EU28, Services"),
    descriptive_stats(exiobase_agg_without_TSF.iloc[:,31], "uniform_0.3_EU28, Industry"),
    descriptive_stats(exiobase_agg_without_TSF.iloc[:,32], "lognorm_1.106_RoW, Services"),
    descriptive_stats(exiobase_agg_without_TSF.iloc[:,33], "lognorm_1.106_EU28, Biodiesels"),
    descriptive_stats(exiobase_agg_without_TSF.iloc[:,34], "lognorm_1.106_EU28, AFOLU"),
    descriptive_stats(exiobase_agg_without_TSF.iloc[:,35], "lognorm_1.106_EU28, Basic iron"),
    descriptive_stats(exiobase_agg_without_TSF.iloc[:,36], "lognorm_1.106_EU28, Energy"),
    descriptive_stats(exiobase_agg_without_TSF.iloc[:,37], "lognorm_1.106_RoW, Waste management"),
    descriptive_stats(exiobase_agg_without_TSF.iloc[:,38], "lognorm_1.106_EU28, Services"),
    descriptive_stats(exiobase_agg_without_TSF.iloc[:,39], "lognorm_1.106_EU28, Industry"),
    descriptive_stats(exiobase_agg_without_TSF.iloc[:,40], "lognorm_1.225_RoW, Services"),
    descriptive_stats(exiobase_agg_without_TSF.iloc[:,41], "lognorm_1.225_EU28, Biodiesels"),
    descriptive_stats(exiobase_agg_without_TSF.iloc[:,42], "lognorm_1.225_EU28, AFOLU"),
    descriptive_stats(exiobase_agg_without_TSF.iloc[:,43], "lognorm_1.225_EU28, Basic iron"),
    descriptive_stats(exiobase_agg_without_TSF.iloc[:,44], "lognorm_1.225_EU28, Energy"),
    descriptive_stats(exiobase_agg_without_TSF.iloc[:,45], "lognorm_1.225_RoW, Waste management"),
    descriptive_stats(exiobase_agg_without_TSF.iloc[:,46], "lognorm_1.225_EU28, Services"),
    descriptive_stats(exiobase_agg_without_TSF.iloc[:,47], "lognorm_1.225_EU28, Industry"),
    descriptive_stats(exiobase_agg_without_TSF.iloc[:,48], "lognorm_1.363_RoW, Services"),
    descriptive_stats(exiobase_agg_without_TSF.iloc[:,49], "lognorm_1.363_EU28, Biodiesels"),
    descriptive_stats(exiobase_agg_without_TSF.iloc[:,50], "lognorm_1.363_EU28, AFOLU"),
    descriptive_stats(exiobase_agg_without_TSF.iloc[:,51], "lognorm_1.363_EU28, Basic iron"),
    descriptive_stats(exiobase_agg_without_TSF.iloc[:,52], "lognorm_1.363_EU28, Energy"),
    descriptive_stats(exiobase_agg_without_TSF.iloc[:,53], "lognorm_1.363_RoW, Waste management"),
    descriptive_stats(exiobase_agg_without_TSF.iloc[:,54], "lognorm_1.363_EU28, Services"),
    descriptive_stats(exiobase_agg_without_TSF.iloc[:,55], "lognorm_1.363_EU28, Industry"),
]

# Convert to DataFrame
stats_df_complete = pd.DataFrame(stats_list)

