import numpy as np

from mhn.optimizers import OmegaOptimizer, StateSpaceOptimizer


# adjust these paths so that they are correct on your machine
PATH_TO_COAD_CSV = "data/COAD_n12.csv"
PATH_TO_LUAD_CSV = "data/LUAD_n12.csv"


def compute_mhns_for_data(path: str, output_name: str):
    np.random.seed(0)
    print(f"Learn oMHN for {output_name}")
    omega_opt = OmegaOptimizer()
    omega_opt.set_penalty(omega_opt.Penalty.SYM_SPARSE)
    omega_opt.load_data_from_csv(path)
    optimal_lambda = omega_opt.find_lambda()
    print(f"optimal lambda: {optimal_lambda}")
    omega_opt.train(lam=optimal_lambda)
    omega_opt.result.save(f"results/{output_name}_oMHN")

    print(f"Learn cMHN for {output_name}")
    classical_opt = StateSpaceOptimizer()
    classical_opt.load_data_from_csv(path)
    classical_opt.set_penalty(classical_opt.Penalty.SYM_SPARSE)
    optimal_lambda = classical_opt.find_lambda()
    print(f"optimal lambda: {optimal_lambda}")
    classical_opt.train(lam=optimal_lambda)
    classical_opt.result.save(f"results/{output_name}_cMHN")


def main():
    # compute COAD MHNs
    compute_mhns_for_data(PATH_TO_COAD_CSV, "COAD")

    # compute LUAD MHNs
    compute_mhns_for_data(PATH_TO_LUAD_CSV, "LUAD")


if __name__ == '__main__':
    main()
