import os

# ---------- PARAMETERS FOR SIMULATION ---------- 

# BIG DATASET FILE PATH
BIG_DIR_INPUT = os.path.join(os.getcwd(), "IOT_2022_pxp")
BIG_DIR_OUTPUT = os.path.join(os.getcwd(), "big_output")
BIG_A_FILE = os.path.join(os.getcwd(), "IOT_2022_pxp", "A.txt")
BIG_S_FILE = os.path.join(os.getcwd(), "IOT_2022_pxp", "satellite", "S.txt")

# SMALL DATASET FILE PATH
SMALL_DIR_INPUT = os.path.join(os.getcwd(), "exiobase_2022_small")
SMALL_DIR_OUTPUT = os.path.join(os.getcwd(), "small_output")
SMALL_A_FILE = os.path.join(os.getcwd(), "exiobase_2022_small", "A.txt")
SMALL_S_FILE = os.path.join(os.getcwd(), "exiobase_2022_small", "satellite", "S.txt")

# CHOSEN ACTIVITIES
SMALL_CHOSEN_ACT = [("RoW-Services", 68), ("EU28-Biodiesels", 11), ("EU28-Agriculture-Forestry-Fishing", 0), ("EU28-Basic iron and steel and of ferro-alloys and first products thereof", 13)]
BIG_CHOSEN_ACT = [("CN-Railway transportation services", 6156), ("DE-Biodiesels", 1093), ("CH-Beverages", 7651), ("SE-Basic iron and steel and of ferro-alloys and first products thereof", 4903)]  # [6156, 1093, 7651, 4903]

# UNCERTAINTY DEFINE
# DIST_TYPE = ["baseline", "uniform", "log-normal"] # Define the types of distribution
DIST_TYPE = ["log-normal"]
U_UNIFORM = [0.1, 0.2, 0.3] # Define the uncertainty for uniform distribution
# U_LOG = [1.1, 1.2, 1.3] # Define the uncertainty for log distribution
U_LOG = [1.106, 1.225, 1.363]
AMOUNT = 4 # This is the amount of activities for 1 CASE

# ---------- PARAMETERS FOR PLOT DRAWING ---------- 

FOLDER_PATH = os.path.join(os.getcwd(), "Exiobase_folders_for_Ning")
DB_SIZE = ["small", "big"]
PLT_COMP_TYPE = ["cases", "sectors"]
PLT_SAVE_PATH = os.path.join(FOLDER_PATH, "plt_output")
