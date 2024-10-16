import numpy as np

from mhn.optimizers import oMHNOptimizer, cMHNOptimizer


# adjust these paths so that they are correct on your machine
PATH_TO_COAD_CSV = "data/COAD_n12.csv"
PATH_TO_LUAD_CSV = "data/LUAD_n12.csv"


def compute_mhns_for_data(path: str, output_name: str):
    np.random.seed(0)
    print(f"Learn oMHN for {output_name}")
    omega_opt = oMHNOptimizer()
    omega_opt.set_penalty(omega_opt.Penalty.SYM_SPARSE)
    omega_opt.load_data_from_csv(path)
    optimal_lambda = omega_opt.lambda_from_cv(lambda_min=0.0001, lambda_max=0.1, steps=9, nfolds=5, show_progressbar=True)
    print(f"optimal lambda: {optimal_lambda}")
    omega_opt.train(lam=optimal_lambda)
    omega_opt.result.save(f"results/{output_name}_oMHN.csv")

    print(f"Learn cMHN for {output_name}")
    classical_opt = cMHNOptimizer()
    classical_opt.load_data_from_csv(path)
    classical_opt.set_penalty(classical_opt.Penalty.SYM_SPARSE)
    optimal_lambda = classical_opt.lambda_from_cv(lambda_min=0.0001, lambda_max=0.1, steps=9, nfolds=5, show_progressbar=True)
    print(f"optimal lambda: {optimal_lambda}")
    classical_opt.train(lam=optimal_lambda)
    classical_opt.result.save(f"results/{output_name}_cMHN.csv")


def main():
    # compute COAD MHNs
    compute_mhns_for_data(PATH_TO_COAD_CSV, "COAD")

    # compute LUAD MHNs
    compute_mhns_for_data(PATH_TO_LUAD_CSV, "LUAD")


if __name__ == '__main__':
    main()
