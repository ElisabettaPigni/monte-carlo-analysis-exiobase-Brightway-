# Quick start
## Run remotely
1. Clone this repository
    ```
    git clone <repo_name>
    cd <repo_folder_name>
    ```
2. Upload your data to remote machine
3. Install dependencies
    ```
    pip install -r requirements.txt
    ```
4. Run script
    ```
    python monte_carlo_V2.0.py
    ```

## Run script
```
# Run simulations
python monte_carlo_V2.0.py

# Run graphical analysis
python statistic_analysis.py
```

## Description
- 20240805-Import-exio3big-and-monte-carlo-bw25.py: This is the origin script from ElisabettaPigni.
- monte_carlo_V2.0.py: This script is for simulation and result saving.
- statistic_analysis.py: This script is used to perform graphical analysis of the data.

### Reference:  
- [Brightway Installation](https://docs.brightway.dev/en/latest/content/installation/index.html)  
- [from-the-ground-up/2 - Building and using matrices in bw2calc.ipynb](https://github.com/brightway-lca/from-the-ground-up/blob/main/2%20-%20Building%20and%20using%20matrices%20in%20bw2calc.ipynb)  
- [Correspondance table exiobase-ecoinvent (biosphere)](https://github.com/brightway-lca/brightway2-io/blob/main/bw2io/data/lci/EXIOBASE-ecoinvent-biosphere.csv?plain=1)
- [Brightway documentation for uncertainty](https://stats-arrays.readthedocs.io/en/latest/#mapping-parameter-array-columns-to-uncertainty-distributions)  
- [EXIOBASE-ecoinvent-biosphere.csv](https://github.com/brightway-lca/brightway2-io/blob/main/bw2io/data/lci/EXIOBASE-ecoinvent-biosphere.csv?plain=1)  
