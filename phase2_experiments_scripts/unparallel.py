import os
import sys
sys.path.append(os.path.abspath(".."))
from constants import *
from monte_carlo_utility import *
from monte_carlo_main import *


if __name__ == "__main__":
    simu = SimulationScript()

    for param in COMBINED_PARAMETERS:
        print(f"{param[4]} simulation is running...")
        tech_matrix, bio_matrix, cf_matrix, activities = prepare_datapackage_matrices(param[0], param[1], param[5])
        parallel_params = []
        k = 0
        for t in DIST_TYPE:
            if t == "static":
                datapackage = create_static_datapackage(tech_matrix, bio_matrix, cf_matrix, activities)
                print("Datapackage are formatted.")
                for myact, index in param[2]:
                    k += 1
                    manual_tech_Data = -tech_matrix
                    np.fill_diagonal(manual_tech_Data, -manual_tech_Data.diagonal())
                    lca_score_manual = simu.manual_lca(manual_tech_Data, bio_matrix, cf_matrix, index + 1)
                    print(f"Manually calculated lca score: {lca_score_manual, myact}")
                    static_lca(index+1, datapackage, param[3], k, t, myact)
            elif t == "pedigree":
                datapackage = create_stochastic_datapackage(tech_matrix, bio_matrix, cf_matrix, activities, param[4], param[5])
                print("Datapackage are formatted.")
                for myact, index in param[2]:
                    k += 1
                    stochastic_lca(index+1, datapackage, param[3], k, t, myact)