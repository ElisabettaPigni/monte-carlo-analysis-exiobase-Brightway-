# Import IO database from original table

# doing the same as here: https://github.com/brightway-lca/from-the-ground-up/blob/main/2%20-%20Building%20and%20using%20matrices%20in%20bw2calc.ipynb

# Brightway documentation for uncertainty https://stats-arrays.readthedocs.io/en/latest/#mapping-parameter-array-columns-to-uncertainty-distributions


#%% importing packages
import bw_processing as bwp #this is needed for brightway 2.5 package
import bw2calc as bc #this is needed for brightway 2.5 package
#import bw2data as bd #this is needed for brightway 2.5 package
import pandas as pd
from scipy import sparse #This is necessary to create the sparse matrix, which is a lighter matrix in which zero values have been removed.
import numpy as np
import os
import seaborn as sb #This is needed for the graphical rapresentation of the monte carlo results
import matplotlib.pyplot as plt #This is needed for the graphical rapresentation of the monte carlo results
from matplotlib.ticker import FuncFormatter #This is needed for the graphical rapresentation of the monte carlo results

#%% path to folders
os.getcwd()
os.chdir(r"C:\Users\Elisabetta\OneDrive - Alma Mater Studiorum Università di Bologna\Documents\UNIVERSITA\Progetto di dottorato\Periodo all'estero 2024\UA_Brightway\exiobase_2022_small")

#%% Importing EXIOBASE matrices (A and S) and modifying them for Brightway compatibility

A_raw = pd.read_table('A.txt')

countries = A_raw['region'].drop_duplicates().iloc[2:].tolist()
sectors = list(A_raw.iloc[2:,1].drop_duplicates())
activities = [ x + '-' + y for x in countries for y in sectors]
units = (['MEUR' for i in range(0,len(activities))])
extensions = ["349b29d1-3e58-4c66-98b9-9d1a076efd2e", "0795345f-c7ae-410c-ad25-1845784c75f5", "20185046-64bb-4c09-a8e7-e8a9e144ca98", "0795345f-c7ae-410c-ad25-1845784c75f5",
"0795345f-c7ae-410c-ad25-1845784c75f5", "0795345f-c7ae-410c-ad25-1845784c75f5", "0795345f-c7ae-410c-ad25-1845784c75f5", "0795345f-c7ae-410c-ad25-1845784c75f5",
"0795345f-c7ae-410c-ad25-1845784c75f5", "0795345f-c7ae-410c-ad25-1845784c75f5", "0795345f-c7ae-410c-ad25-1845784c75f5", "349b29d1-3e58-4c66-98b9-9d1a076efd2e",
"349b29d1-3e58-4c66-98b9-9d1a076efd2e", "35d1dff5-b535-4628-9826-4a8fce08a1f2", "e5ba3517-a93f-422e-9ce2-f1007fcf6c06", "0795345f-c7ae-410c-ad25-1845784c75f5",
"349b29d1-3e58-4c66-98b9-9d1a076efd2e", "20185046-64bb-4c09-a8e7-e8a9e144ca98", "0795345f-c7ae-410c-ad25-1845784c75f5", "349b29d1-3e58-4c66-98b9-9d1a076efd2e", "349b29d1-3e58-4c66-98b9-9d1a076efd2e"
]

# to obtain ecoinvent code directly from ecoinvent, use this code: https://github.com/brightway-lca/brightway2-io/blob/main/bw2io/data/lci/EXIOBASE-ecoinvent-biosphere.csv?plain=1

# quick check
len(sectors)
len(activities) == len(units)


A_IO = A_raw.iloc[2:,2:].astype('float').values
I = np.identity(len(A_IO))
A_ = I - A_IO
A = -A_
np.fill_diagonal(A, -A.diagonal()) # then change back again, but only the diagonal
#print(A[0:,:5])

len(A[:,1])
len(A[1,:])
len([i for i in A[:,1] if i != 0]) # non-zero values
len([i for i in A[:,1] if i == 0]) # zero values

# technology matrix A as sparse object and then coordinates
Asparse = sparse.coo_array(A)
a_data = Asparse.data # amounts or values
a_indices = np.array([tuple(coord) for coord in np.transpose(Asparse.nonzero())], dtype=bwp.INDICES_DTYPE) # indices of each exchange
a_flip = np.array([False if i[0] == i[1] else True for i in a_indices ]) # Numerical sign of the inputs needs to be flipped negative

# import environemntla extensions

S_raw = pd.read_table("satellite/S.txt", header=[0,1], index_col=[0])

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

len(B[:,1])
len(B[1,:]) == len(A[1,:])
len([i for i in B[:,1] if i != 0]) # non-zero values
len([i for i in B[:,1] if i == 0]) # zero values

# Intervention matrix B as sparse object and then coordinates

Bsparse = sparse.coo_array(B)
b_data = Bsparse.data # amounts or values
b_indices_remap = [[i[0] + len(activities),i[1]] for i in np.transpose(Bsparse.nonzero())] # need to make sure biosphere indices are different from technosphere
b_indices = np.array([tuple(coord) for coord in b_indices_remap], dtype=bwp.INDICES_DTYPE)

# matrix of charachterisation factors (method E.F. 3.1 - Climate change)

CFs = [1., 29.8, 273., 29.8, 29.8, 29.8, 29.8, 29.8, 29.8, 29.8, 29.8, 1., 1., 25200., 14600., 29.8, 1., 273., 29.8, 1., 1.]
len(extensions) == len(CFs)

C = np.matrix(np.zeros((len(CFs), len(CFs))))
C_diag = np.matrix(CFs)
np.fill_diagonal(C, C_diag)
#C

len(C[:,1])
len(C[1,:])

len([i for i in C[:,1] if i != 0]) # non-zero values
len([i for i in C[:,1] if i == 0]) # zero values


# Sparse C  matrix of characterisation factors
Csparse = sparse.coo_array(C)
c_data =  Csparse.data 
c_indices_remap = [[i[1] + len(activities),i[1]+ len(activities)] for i in np.transpose(Csparse.nonzero())] # same indices as in B
c_indices = np.array([tuple(coord) for coord in c_indices_remap], dtype=bwp.INDICES_DTYPE) 


# Creating a datapackage
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

#%% LIFE CYCLE ASSESSMENT WITH BRIGHTWAY 2.5 (Sector ??? - method E.F. 3.1 - Climate change)

 
myact = activities[1] # to find the desired code for the functional unit
print(myact)

lca = bc.LCA(
    demand={activities.index(myact): 1}, # using index because the argument of FU is an INTEGER
    data_objs=[dp_static],
)
lca.lci()
lca.lcia()
lca.score


# check with normal matrix operation
f = np.zeros(len(A))
f[1] = 1 # functional unit

print('ioscore',np.sum(C.dot(B.dot((np.linalg.inv(A_)).dot(f))))) # matrix operation
print('lca score',lca.score)  # similar values, when rounded. All is good


#%% STOCHASTIC LCA

# Description:
    # CASE 1: uniform distribution with 10% uncertainty
    # CASE 2: uniform distribution with 20% uncertainty
    # CASE 3: uniform distribution with 30% uncertainty
    # CASE 4: log-normal distribution with 1.01 uncertainty
    # CASE 5: log-normal distribution with 1.1 uncertainty 
    # CASE 6: log-normal distribution with 2 uncertainty 


#%% CASE 1: uniform distribution with 10% uncertainty - set the uncertainty and plot the results
u = [0.1, 0.2, 0.3]

for k in range(1,4):
    for j in range(u):
    #set the uncertainty for matrix A
        min_val_a = a_data - (a_data * u[j])
        max_val_a = a_data + (a_data * u[j])
    
        results_a = []
    
    # Iterate over indices and check corresponding values in a_flip
        for i in range(len(a_data)):
            if not a_flip[i]:  # If a_flip at index i is False
                parameters_a = (0, a_data[i], np.NaN, np.NaN, np.NaN, np.NaN, False)
            else:  # If a_flip at index i is True
                parameters_a = (4, np.NaN, np.NaN, np.NaN, min_val_a[i], max_val_a[i], False)
            results_a.append(parameters_a)
        
        a_uncertainty = np.array(results_a, 
            dtype=[('uncertainty_type', 'i4'), ('loc', 'f4'), ('scale', 'f4'), 
                   ('shape', 'f4'), ('minimum', 'f4'), ('maximum', 'f4'), ('negative', 'b')]
        )
    
    #The strings 'i4', 'f4', and 'b' are type codes used in NumPy to specify the data types of array elements. 
    #'i4': This stands for a 4-byte (32-bit) integer.
    #'f4': This stands for a 4-byte (32-bit) float.
    #'b': This stands for a boolean value.
    
    
    #set th uncertainty for matrix B
        min_val_b = b_data - (b_data * u[j])
        max_val_b = b_data + (b_data * u[j])
    
        results_b= []
    
    # Iterate over indices and check corresponding values in a_flip
        for i in range(len(b_data)):
            parameters_b = (4, np.NaN, np.NaN, np.NaN, min_val_b[i], max_val_b[i], False)
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
        
        myact = activities[1] # to find the desired code for the functional unit
        print(myact)
        
        lca = bc.LCA(
            demand={activities.index(myact): 1},
            data_objs=[dp_stochastic],
            use_distributions=True,
        )
        lca.lci()
        lca.lcia()


# Perform Monte Carlo simulation for CASE 1

# Define the directory for saving results
        directory = r"C:\Users\Elisabetta\OneDrive - Alma Mater Studiorum Università di Bologna\Documents\UNIVERSITA\Progetto di dottorato\Periodo all'estero 2024\UA_Brightway\MonteCarloResults_uniform_distrubution"
        os.makedirs(directory, exist_ok=True)  # Create the directory if it does not exist
        
        # Define simulation parameters
        batch_size = 50
        num_batches = 10
        
        # List to store the results
        cumulative_results = []

# Run Monte Carlo simulations and save results
        for p in range(num_batches):
            # Run simulations for the current batch
            batch_results = [lca.score for _ in zip(range(batch_size), lca)]
            
            # Accumulate results
            cumulative_results.extend(batch_results)
            
            # Convert to DataFrame
            df_cumulative = pd.DataFrame(cumulative_results, columns=["kg CO2eq"])
            
            # Define filename for saving results
            filename = os.path.join(directory, f"CASE_{k}_MC_simulations.csv")
            
            # Save to CSV
            df_cumulative.to_csv(filename, index=False)
            
            # Print progress
            print("Saved MC_simulations")

print("All simulations completed.")



# Display the results


# Custom formatter function
def scientific_format(x, pos):
    return f'{x:.2e}'

for k in range(1,4):
# Create a formatter object
    formatter = FuncFormatter(scientific_format)
    
    os.chdir(r"C:\Users\Elisabetta\OneDrive - Alma Mater Studiorum Università di Bologna\Documents\UNIVERSITA\Progetto di dottorato\Periodo all'estero 2024\UA_Brightway\MonteCarloResults_uniform_distrubution")
    
    # Load the Monte Carlo results
    mc_results = pd.read_csv(f"CASE_{k}_MC_simulations.csv", header=0)
    
    # Rename the column for clarity
    mc_results.columns = ["kg CO2eq"]
    
    # Save summary statistics to a CSV file
    summary_stats = mc_results.describe()
    summary_stats.to_csv("CASE_{k}_summary_statistics.csv")
    
    # Create and save the boxplot
    plt.figure(figsize=(10, 6))
    sb.histplot(data=mc_results, x="kg CO2eq", kde=True)
    plt.title('Climate change (kg CO2eq)')
    plt.gca().yaxis.set_major_formatter(formatter)
    plt.savefig("CASE_{k}_barplot.png")
    plt.close()


    # Create and save the violin plot
    plt.figure(figsize=(10, 6))
    sb.violinplot(data=mc_results, y="kg CO2eq")
    plt.title('Climate change (kg CO2eq)')
    plt.gca().yaxis.set_major_formatter(formatter)
    plt.savefig("CASE_{k}_violinplot.png")
    plt.close()
    
    # Create and save the scatter plot
    plt.figure(figsize=(10, 6))
    sb.scatterplot(data=mc_results, y="kg CO2eq", x=mc_results.index)
    plt.title('Climate change (kg CO2eq)')
    plt.xlabel('MC Simulations')
    plt.gca().yaxis.set_major_formatter(formatter)
    plt.savefig("CASE_{k}_scatterplot.png")
    plt.close()

print("Summary statistics and plots saved successfully.")


#%%
# Set the uncertainty for CASE 5: log-normal distribution with 1.01 uncertainty
u = 1.01

#add uncertainty for A
mu_a = np.log(a_data)
#sigma_a = np.full(len(a_data), np.log(u))
sigma = np.log(u)


results_a = []

# Iterate over indices and check corresponding values in a_flip
for i in range(len(a_data)):
    if not a_flip[i]:  # If a_flip at index i is False
        parameters_a = (0, a_data[i], np.NaN, np.NaN, np.NaN, np.NaN, False)
    else:  # If a_flip at index i is True
        parameters_a = (2, mu_a[i], sigma, np.NaN, np.NaN, np.NaN, False)
    results_a.append(parameters_a)

a_uncertainty = np.array(results_a, 
    dtype=[('uncertainty_type', 'i4'), ('loc', 'f4'), ('scale', 'f4'), 
           ('shape', 'f4'), ('minimum', 'f4'), ('maximum', 'f4'), ('negative', 'b')]
)

#The strings 'i4', 'f4', and 'b' are type codes used in NumPy to specify the data types of array elements. 
#'i4': This stands for a 4-byte (32-bit) integer.
#'f4': This stands for a 4-byte (32-bit) float.
#'b': This stands for a boolean value.


#add uncertainty for B
mu_b = np.log(b_data)
#sigma_b = np.full(len(b_data), np.log(u))


results_b = []

# Iterate over indices and check corresponding values in a_flip
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

myact = activities[1] # to find the desired code for the functional unit
print(myact)

lca = bc.LCA(
    demand={activities.index(myact): 1},
    data_objs=[dp_stochastic],
    use_distributions=True,
)
lca.lci()
lca.lcia()


# Perform Monte Carlo simulation for CASE 5

# Define the directory for saving results
directory = r"C:\Users\Elisabetta\OneDrive - Alma Mater Studiorum Università di Bologna\Documents\UNIVERSITA\Progetto di dottorato\Periodo all'estero 2024\UA_Brightway\MonteCarloResults_lognormal_distrubution"
os.makedirs(directory, exist_ok=True)  # Create the directory if it does not exist

# Define simulation parameters
batch_size = 50
num_batches = 20

# List to store the results
cumulative_results = []

# Run Monte Carlo simulations and save results
for i in range(num_batches):
    # Run simulations for the current batch
    batch_results = [lca.score for _ in zip(range(batch_size), lca)]
    
    # Accumulate results
    cumulative_results.extend(batch_results)
    
    # Convert to DataFrame
    df_cumulative = pd.DataFrame(cumulative_results, columns=["kg CO2eq"])
    
    # Define filename for saving results
    filename = os.path.join(directory, f"CASE_{j}5_cumulative_results_up_to_batch_{50+(i*50)}.csv")
    
    # Save to CSV
    df_cumulative.to_csv(filename, index=False)
    
    # Print progress
    print(f"Saved cumulative results up to batch {50+(i*50)}")

print("All simulations completed.")



# Display the results for CASE 5


# Custom formatter function
def scientific_format(x, pos):
    return f'{x:.2e}'

# Create a formatter object
formatter = FuncFormatter(scientific_format)

os.chdir(r"C:\Users\Elisabetta\OneDrive - Alma Mater Studiorum Università di Bologna\Documents\UNIVERSITA\Progetto di dottorato\Periodo all'estero 2024\UA_Brightway\MonteCarloResults_lognormal_distrubution")

# Load the Monte Carlo results
mc_results = pd.read_csv("CASE_5_cumulative_results_up_to_batch_1000.csv", header=0)

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



#%%
# Matrix Diagrams

import seaborn as sns
import numpy as np

matrix_IO = A_IO[:8,:8]
matrix_p = A[:8,:8]
I_1 = np.eye(8)
leontief_inverse = np.linalg.inv(I_1 - matrix_IO)
inverse_matrix = np.linalg.inv(matrix_p)
matrix_B = B[:8,:8] 

fig, axs = plt.subplots(2, 3, figsize=(15, 10))

sns.heatmap(matrix_IO, ax=axs[0, 0], annot=True, cbar=False)
axs[0, 0].set_title('Input-Output Matrix (A_IO)')
sns.heatmap(leontief_inverse, ax=axs[0, 2], annot=True, cbar=False)
axs[0, 2].set_title('Leontief Inverse ((I-A_IO)^-1)')
sns.heatmap(matrix_B, ax=axs[0, 1], annot=True, cbar=False)
axs[0, 1].set_title('Direct stressor/impact coefficients (S)')

sns.heatmap(matrix_p, ax=axs[1, 0], annot=True, cbar=False)
axs[1, 0].set_title('Technological Matrix (A_p)')
sns.heatmap(inverse_matrix, ax=axs[1, 2], annot=True, cbar=False)
axs[1, 2].set_title('Inverse Matrix (A_p^-1)')
sns.heatmap(matrix_B, ax=axs[1, 1], annot=True, cbar=False)
axs[1, 1].set_title('Emission Factors (B)')
sns






