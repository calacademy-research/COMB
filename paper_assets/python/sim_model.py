import numpy as np
from dataclasses import dataclass
from sim_data import SimParams, SimData
import pyjags.model
from typing import Union, List


@dataclass
class SimModelParams:
    """
    Dataclass to store simulation model parameters
    """

    psi_prior: str = "dbeta(2,2)"
    beta0_prior: str = "dnorm(0, 0.5)"
    beta1_prior: str = "dnorm(0, 0.5)"
    p11_prior: str = "dbeta(2, 2)"
    mu1_prior: str = "dnorm(-2, 0.2)"
    mu2_prior: str = "dnorm(-2, 0.2)"
    sigma1_prior: str = "dunif(0.1, 5)"
    sigma2_prior: str = "dunif(0.1, 5)"
    p_aru11_prior: str = "dbeta(2, 2)"
    p_aru01_prior: str = "dbeta(1, 3)I(0, 1 - p_aru11)"
    na: int = 1000  # number of iterations to adapt
    ni: int = 8000  # number of iterations
    nt: int = 1  # thinning rate
    nb: int = 1000  # number of burn-in iterations
    nc: int = 6  # number of chains
    parallel: bool = True  # whether to run each chain in parallel

    include_aru_model: bool = True
    include_pc_model: bool = True
    include_scores_model: bool = True
    include_covar_model: bool = True
    aru_scores_independent_model: bool = True

    sim_data: Union[SimData, None] = None


class SimModel(pyjags.model.Model):
    """
    SimModel class to make a model for simulation based on parameters
    """

    model_text: str = ""

    def __init__(self, params: SimModelParams):
        self.params = params
        self.model_text = self.gen_jags_model_text()
        with open("model_text.txt", "w") as f:
            f.write(self.model_text)

        self.sim_data = params.sim_data.__dict__

        # Find the variables in the model that are in the simulation data
        # Depending on the model, not all of the simulation data will be used
        # For example, if the ARU model is not included, the ARU data will not be used
        sim_vars = self.find_variables(list(self.sim_data.keys()))

        self.sim_data = {var: self.sim_data[var] for var in sim_vars}

        print(self.sim_data)

        self.inits = self.create_inits()

        self.monitored_vars = self.get_monitored_vars()

        super().__init__(
            code=self.model_text,
            data=self.sim_data,
            chains=params.nc,
            adapt=params.na,
            init=self.inits,
        )

    def sample(self):
        """
        Sample from the model
        """

        return super().sample(
            iterations=self.params.ni,
            thin=self.params.nt,
            vars=self.monitored_vars,
        )

    def gen_jags_model_text(self):
        """
        Generate JAGS model text based on the simulation model parameters
        """
        params = self.params

        string_init = "model {\n"
        priors_string = "# Priors\n"
        if not (params.include_aru_model or params.include_pc_model or params.include_scores_model):
            raise ValueError("At least one of the models must be included")

        if params.include_covar_model:
            reg_priors = f"beta0 ~ {params.beta0_prior}\n beta1 ~ {params.beta1_prior}\n"
            priors_string += reg_priors

            occ_lik = """
            # Likelihood for occupancy
            for (i in 1:nsites) {
                logit(psi[i]) <- beta0 + beta1 * covar[i]
                z[i] ~ dbern(psi[i])
            """
        else:
            priors_string += f"psi ~ {params.psi_prior}\n"
            occ_lik = """
            # Likelihood for occupancy
            for (i in 1:nsites) {
                z[i] ~ dbern(psi)
            """

        if params.include_pc_model:
            pc_lik = """
            # Point count
            p[i] <- p11 * z[i]
            for (j in 1:nsurveys_pc) {
                y_pc[i, j] ~ dbern(p[i])
            }
            """
            p11_prior_string = f"p11 ~ {params.p11_prior}\n"
            priors_string += p11_prior_string
        else:
            pc_lik = ""

        if params.include_aru_model and params.include_scores_model:
            if params.aru_scores_independent_model:
                aru_lik = """
                # ARU
                p_aru[i] <- z[i]*p_aru11 + p_aru01
                for(j in 1:nsurveys_aru) {
                    y_aru[i,j] ~ dbern(p_aru[i])
                }
                for(j in 1:nsurveys_scores) {
                    scores[i,j] ~ dnorm(mu[z[siteid[i]] + 1], tau[z[siteid[i]] + 1])
                }
                }"""
            else:
                aru_lik = """
                # ARU
                p_aru[i] <- z[i]*p_aru11 + p_aru01
                for(j in 1:nsurveys_aru) {
                    y_aru[i,j] ~ dbern(p_aru[i])
                    scores[i,j] ~ dnorm(mu[z[siteid[i]] + 1], tau[z[siteid[i]] + 1]) T(ifelse(y_aru[i,j] == 1, threshold, -10), ifelse(y_aru[i,j] == 1, 10, threshold))
                }
                }"""

            p_aru11_prior_string = f"p_aru11 ~ {params.p_aru11_prior}\n"
            p_aru01_prior_string = f"p_aru01 ~ {params.p_aru01_prior}\n"
            scores_prior = f"""
                # Priors of the observation model for the scores
                mu[1] ~ {params.mu1_prior}
                mu[2] ~ {params.mu2_prior}
                sigma[1] ~ {params.sigma1_prior}
                tau[1] <- 1 / (sigma[1] * sigma[1])
                sigma[2] ~ {params.sigma2_prior}
                tau[2] <- 1 / (sigma[2] * sigma[2])
            """
            priors_string += p_aru11_prior_string + p_aru01_prior_string + scores_prior
        elif params.include_aru_model:
            aru_lik = """
                # ARU - binomial
                p_aru[i] <- z[i]*p_aru11 + p_aru01
                for(j in 1:nsurveys_aru) {
                    y_aru[i,j] ~ dbern(p_aru[i])
                    }
                }
                """
            p_aru11_prior_string = f"p_aru11 ~ {params.p_aru11_prior}\n"
            p_aru01_prior_string = f"p_aru01 ~ {params.p_aru01_prior}\n"
            priors_string += p_aru11_prior_string + p_aru01_prior_string
        elif params.include_scores_model:
            aru_lik = """
                for(j in 1:nsurveys_scores) {
                    scores[i,j] ~ dnorm(mu[z[siteid[i]] + 1], tau[z[siteid[i]] + 1])
                }
                }
                """
            scores_prior = f"""
                # Priors of the observation model for the scores
                mu[1] ~ {params.mu1_prior}
                mu[2] ~ {params.mu2_prior}
                sigma[1] ~ {params.sigma1_prior}
                tau[1] <- 1 / (sigma[1] * sigma[1])
                sigma[2] ~ {params.sigma2_prior}
                tau[2] <- 1 / (sigma[2] * sigma[2])
            """
            priors_string += scores_prior
        else:
            aru_lik = "\n}\n"

        close_string = """
        mean_psi <- mean(psi)
        NOcc <- sum(z[])
        PropOcc <- NOcc / nsites
        }
        """

        return f"""
        {string_init}
        {priors_string}
        {occ_lik}
        {pc_lik}
        {aru_lik}
        {close_string}
        """

    def find_variables(self, search_vars: List) -> List[str]:
        """
        Find the variables in the model that are in `search_vars`

        Args:
            search_vars (List): The variables to search for in the model text

        Returns:
            List[str]: The variables in the model from `search_vars`
        """

        return [var for var in search_vars if var in self.model_text]

    def create_inits(self) -> dict:
        """
        Create initial values for the model
        """
        mu = np.random.normal(0, 1, 2)
        mu.sort()

        inits_full = {
            "mu": mu,
            "sigma": np.random.uniform(0.5, 2.5, 2),
            # "psi": np.repeat(0.5, self.sim_data["nsites"]),
            "z": np.repeat(1, self.sim_data["nsites"]),
            "beta0": 0,
            "beta1": 0,
            "p_aru11": 0.2,
            "p_aru01": 0.1,
            "p11": 0.2,
        }

        init_keys_in_model = self.find_variables(list(inits_full.keys()))
        inits_in_model = {var: inits_full[var] for var in init_keys_in_model}

        return inits_in_model

    def get_monitored_vars(self):
        """
        Get the monitored variables based on the model type
        """
        monitored_list = ["PropOcc"]
        if self.params.include_covar_model:
            monitored_list += ["beta0", "beta1", "mean_psi"]
        else:
            monitored_list += ["psi"]

        if self.params.include_pc_model:
            monitored_list += ["p11"]

        if self.params.include_aru_model:
            monitored_list += ["p_aru11", "p_aru01"]

        if self.params.include_scores_model:
            monitored_list += ["mu", "sigma"]

        return monitored_list
