old_configs = {
    # default config with no changes
    "default": {},
    # first big config used
    "oct_2024": {
        "p11": [0.9, 0.5, 0.1],
        "p_aru11": [0.9, 0.5, 0.1],
        "p_aru01": [0.05, 0],
        "mu": [(-2, -1.75), (-2, 0)],
        "sigma": [(0.25, 1)],
        "threshold": [-1, 0, 1],
        "nsites": [80, 200],
        "include_covar_model": [True],
        "covar_continuous": [True],
        "beta0": [-1],
        "beta1": [-1, 1],
        "nsurveys_aru": [24],
        "nsurveys_scores": [24, 8],
        "nsurveys_pc": [3],
        "aru_scores_independent_model": [True, False],
    },
}
configs = {
    # default config with no changes
    "default": {},
    # testing with only a few combs
    "testing_few_combs": {
        "model_params": {
            "include_aru_model": [True, False],
        }
    },
    # new with different model on same data
    "july_2025": {
        "data_params": {
            "p11": [0.9, 0.5, 0.1],
            "p_aru11": [0.9, 0.5, 0.1],
            "p_aru01": [0.05, 0],
            "mu": [(-2, -1.75), (-2, 0)],
            "sigma": [(0.25, 1)],
            "threshold": [-1, 0, 1],
            "nsites": [80, 200],
            "beta0": [-1],
            "beta1": [0, 1],
            "nsurveys_aru": [24],
            "nsurveys_scores": [24],
            "nsurveys_pc": [3],
        },
        "model_params": {
            "include_aru_model": [True, False],
            "include_pc_model": [True, False],
            "include_scores_model": [True, False],
            "include_covar_model": [True],
        },
    },
}
