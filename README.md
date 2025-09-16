# Quick start
## Run remotely
1. Clone this repository.
    ```
    git clone <repo_name>
    cd <repo_folder_name>
    ```
2. Upload your data to remote machine.
4. Install dependencies.
    ```
    pip install -r requirements.txt
    ```
5. Set up your parameters in `constants.py`.
6. Run script.
    ```
    python monte_carlo_V2.0.py
    ```
7. If you want to run cases in parallel, use this command.
   ```
    python parallel.py
   ```

## Run locally
1. Check current directory
   ```
    pwd
   ```
2. Move to the desire directory
   ```
    cd "<relative path>"
   ```
3. Make a new folder to store scripts
    ```
    mkdir <folder name>
    ```
4. Move to the new folder
   ```
    cd <folder name>
   ```
5. Clone the repository
    ```
    git clone <HTTPS URL>
    ```
6. Run the simulation (non-parallel)
    ```
    python monte_carlo_V2.0.py
    ```
7. Run plot script
    ```
    python monte_carlo_V2.0.py
    ```

## Description
- `20240805-Import-exio3big-and-monte-carlo-bw25.py`: This is the origin script from ElisabettaPigni.
- `monte_carlo_utility.py`: All functionality is implemented here.
- `monte_carlo_V2.0.py`: This is the entry file for run simulations or draw plots.
- `constants.py`: Include all the constants.
- `parallel.py`: This is the entry file for run simulations in parallel.
- `requirements.txt`: Include all the dependencies.
- `activities_big_dataset.csv` All activities from big exiobase.

### Reference:  
- [Brightway Installation](https://docs.brightway.dev/en/latest/content/installation/index.html)  
- [from-the-ground-up/2 - Building and using matrices in bw2calc.ipynb](https://github.com/brightway-lca/from-the-ground-up/blob/main/2%20-%20Building%20and%20using%20matrices%20in%20bw2calc.ipynb)  
- [Correspondance table exiobase-ecoinvent (biosphere)](https://github.com/brightway-lca/brightway2-io/blob/main/bw2io/data/lci/EXIOBASE-ecoinvent-biosphere.csv?plain=1)
- [Brightway documentation for uncertainty](https://stats-arrays.readthedocs.io/en/latest/#mapping-parameter-array-columns-to-uncertainty-distributions)  
- [EXIOBASE-ecoinvent-biosphere.csv](https://github.com/brightway-lca/brightway2-io/blob/main/bw2io/data/lci/EXIOBASE-ecoinvent-biosphere.csv?plain=1)
- [EXIOBASE Description](https://zenodo.org/records/5589597)

[![DOI](https://zenodo.org/badge/838864323.svg)](https://doi.org/10.5281/zenodo.17130112)
