import os

# ---------- PARAMETERS FOR SIMULATION ---------- 
# EXIOBASE FILE PATH
EXIOBASE_INPUT = os.path.join(os.getcwd(), "IOT_2022_pxp")
EXIOBASE_OUTPUT = os.path.join(os.getcwd(), "big_output")
EXIOBASE_A_FILE = os.path.join(os.getcwd(), "IOT_2022_pxp", "A.txt")
EXIOBASE_S_FILE = os.path.join(os.getcwd(), "IOT_2022_pxp", "satellite", "S.txt")

# EXIOBASE AGGREGATED FILE PATH
EXIOBASE_AGGREGATED_INPUT = os.path.join(os.getcwd(), "exiobase_2022_small")
EXIOBASE_AGGREGATED_OUTPUT = os.path.join(os.getcwd(), "small_output")
EXIOBASE_AGGREGATED_A_FILE = os.path.join(os.getcwd(), "exiobase_2022_small", "A.txt")
EXIOBASE_AGGREGATED_S_FILE = os.path.join(os.getcwd(), "exiobase_2022_small", "satellite", "S.txt")

# SELECTED ACTIVITIES
SELECTED_EXIOBASE = [("CN-Railway transportation services", 6156), ("DE-Biodiesels", 1093), ("CH-Beverages", 7651), ("SE-Basic iron and steel and of ferro-alloys and first products thereof", 4903), ("DE-Paraffin Waxes", 1081), ("RU-Food waste for treatment: incineration", 7375), ("NO-Other services (93)", 8397), ("AT-Office machinery and computers (30)", 118)]
SELECTED_AGGREGATED = [("RoW-Services", 68), ("EU28-Biodiesels", 11), ("EU28-Agriculture-Forestry-Fishing", 0), ("EU28-Basic iron and steel and of ferro-alloys and first products thereof", 13), ("EU28-Energy", 1), ("RoW-Waste management", 71), ("EU28-Services", 30), ("EU28-Industry", 3)]

# UNCERTAINTY
DIST_TYPE = ["static", "uniform", "log-normal", "pedigree"] # Define the types of distribution
U_UNIFORM = [0.1, 0.2, 0.3] # Define the uncertainty for uniform distribution
U_LOG = [1.106, 1.225, 1.363] # Define the uncertainty for log-normal distribution

# COMBINED PARAMETERS FOR PHASE 1 EXPERIMENTS (WITHOUT EXTRA DATA)
# COMBINED_PARAMETERS = [(EXIOBASE_AGGREGATED_A_FILE, EXIOBASE_AGGREGATED_S_FILE, SELECTED_AGGREGATED, EXIOBASE_AGGREGATED_OUTPUT, "EXIOBASE AGGREGATED"), (EXIOBASE_A_FILE, EXIOBASE_S_FILE, SELECTED_EXIOBASE, EXIOBASE_OUTPUT, "EXIOBASE")]

# ---------- CONSTANTS FOR CASE STUDY ----------
METHOD = ('EF v3.1', 'climate change', 'global warming potential (GWP100)')
METHOD_IPCC = ('IPCC 2013', 'climate change', 'global warming potential (GWP100)')
GHG = ["CO2 - combustion - air",
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
        "SF6 - air"]

CFS = [1., 1., 1., 1., 1., 1., 27.0, 29.8, 29.8, 29.8, 29.8, 29.8, 29.8, 29.8, 29.8, 29.8, 29.8, 273., 273., 25200.]

COUNTRY_FILE = "../gsd_data/Grouping_reg.csv"
SECTOR_FILE = "../gsd_data/Grouping_sec.csv"
GSD_FILE = "../gsd_data/GSD_sec_reg.csv"
GSD_SMALL_FILE = "../gsd_data/GSD_background_pedigree_exiobase_small.csv"

# ----- EXTEND FILE FOR FIRST PEIDGREE EXPERIMENTS
# BIG_EXTEND = "../extend_data/20241121_foreground_system_big.csv"
# SMALL_EXTEND = "../extend_data/20241121_foreground_system_small.csv"

# ----- EXTEND FILE FOR SECODN PEDIGREE EXPERIMENTS
BIG_EXTEND = "../extend_data/20250806_foreground_system_big.csv"
SMALL_EXTEND = "../extend_data/20250806_foreground_system_small.csv"

# COMBINED PARAMETERS FOR PHASE 2 EXPERIMENTS (WITH EXTRA DATA)
# COMBINED_PARAMETERS = [(EXIOBASE_AGGREGATED_A_FILE, EXIOBASE_AGGREGATED_S_FILE, SELECTED_AGGREGATED, EXIOBASE_AGGREGATED_OUTPUT, "EXIOBASE AGGREGATED", SMALL_EXTEND), (EXIOBASE_A_FILE, EXIOBASE_S_FILE, SELECTED_EXIOBASE, EXIOBASE_OUTPUT, "EXIOBASE", BIG_EXTEND)]

# COMBINED PAPRAMETERS FOR PHASE 2 EXPERIMENTS (WITH EXTRA DATA, SAME FUNCTIONAL UNIT)
SELECTED_EXIOBASE_2 = [("CN-Railway transportation services", 0), ("DE-Biodiesels", 0), ("CH-Beverages", 0), ("SE-Basic iron and steel and of ferro-alloys and first products thereof", 0), ("DE-Paraffin Waxes", 0), ("RU-Food waste for treatment: incineration", 0), ("NO-Other services (93)", 0), ("AT-Office machinery and computers (30)", 0)]
SELECTED_AGGREGATED_2 = [("RoW-Services", 0), ("EU28-Biodiesels", 0), ("EU28-Agriculture-Forestry-Fishing", 0), ("EU28-Basic iron and steel and of ferro-alloys and first products thereof", 0), ("EU28-Energy", 0), ("RoW-Waste management", 0), ("EU28-Services", 0), ("EU28-Industry", 0)]

COMBINED_PARAMETERS = [(EXIOBASE_AGGREGATED_A_FILE, EXIOBASE_AGGREGATED_S_FILE, [("exrta_column", 0)], EXIOBASE_AGGREGATED_OUTPUT, "EXIOBASE AGGREGATED", SMALL_EXTEND), (EXIOBASE_A_FILE, EXIOBASE_S_FILE, [("extra_column", 0)], EXIOBASE_OUTPUT, "EXIOBASE", BIG_EXTEND)]

# ---------- CONSTANTS FOR PLOT DRAWING ----------
FOLDER_PATH = os.path.join(os.getcwd(), "Exiobase_folders_for_Ning")
DB_SIZE = ["small", "big"]
PLT_COMP_TYPE = ["cases", "sectors"]
PLT_SAVE_PATH = os.path.join(FOLDER_PATH, "plt_output")


# ---------- CONSTANTS FOR PLOT DRAWING V2.0 (Always big then small, to prevent human mistake.)----------
LABEL_FONT = 16
TICK_FONT = 14

DATABASE_NAME = ["EXIOBASE", "EXIOBASE_aggregate"]
UNCERTAINTY_NAME = ["static_MC", "log-normal_1.106", "log-normal_1.225", "log-normal_1.363", "uniform_0.1", "uniform_0.2", "uniform_0.3"]
SECTOR_NAME_BIG = ["DE-Paraffin_Waxes", "RU-Food_waste_for_treatment__incineration", "NO-Other_services_(93)", "AT-Office_machinery_and_computers_(30)"]
SECTOR_NAME_SMALL = ["EU28-Energy", "RoW-Waste_management", "EU28-Services", "EU28-Industry"]
COMPARE_TYPE = ["sectors", "cases"]
DATABASE_SECTORS = [(COMPARE_TYPE[0], DATABASE_NAME[0], SECTOR_NAME_BIG), (COMPARE_TYPE[0], DATABASE_NAME[1], SECTOR_NAME_SMALL)]
DATABASE_UNCERTAINTIES = [(COMPARE_TYPE[1], DATABASE_NAME[0], UNCERTAINTY_NAME), (COMPARE_TYPE[1], DATABASE_NAME[1], UNCERTAINTY_NAME)]
