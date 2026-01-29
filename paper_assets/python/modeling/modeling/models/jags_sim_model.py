import numpy as np
from .sim_data import DataParams, SimParams
import pyjags.model
from typing import List, Any
import os

# temp directories for model and data text, useful for debugging
TMP_MODEL = "tmp/model_text.txt"
TMP_DATA = "tmp/data_text.txt"
os.makedirs("tmp", exist_ok=True)


class SimData(pyjags.model.Model):
    """
    SimData class generates data for a simulation
    """

    def __init__(self, params: DataParams):
        self.params = params

        # generate the covariates
        # they are technically not parameters, but we need them to generate the data
        self.covars = self._gen_covars()

        # first generate the simulation data
        # generate jags data model text
        # this "model" will only generate data to be used in the simulation

        self.data_text = self.gen_jags_data_text()
        self.raw_sim_data = self.gen_simulated_data()

    def get_data(self):
        """
        Returns the generated data
        """
        return self.raw_sim_data, self.covars

    def gen_jags_data_text(self) -> str:
        """
        Uses JAGS to generate a "model" string
        that only generates data

        We use this to generate data from the model
        for the purposes of simulation
        """
        string_init = "data {\n"

        if self.params.include_covar_data:
            occ_lik = """
            for (i in 1:nsites) {
                logit(psi[i]) <- beta0 + beta1*covar[i]
                z[i] ~ dbern(psi[i])
            """
        else:
            occ_lik = """
            for (i in 1:nsites) {
                z[i] ~ dbern(psi)
            """

        if self.params.include_pc_data:
            pc_lik = """
            p[i] <- p11 * z[i]
            for (j in 1:nsurveys_pc) {
                y_pc[i, j] ~ dbern(p[i])
            }
            """
        else:
            pc_lik = ""

        if self.params.include_aru_data and self.params.include_scores_data:
            # independent model for aru likelihood
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

        elif self.params.include_aru_data:
            aru_lik = """
                # ARU - binomial
                p_aru[i] <- z[i]*p_aru11 + p_aru01
                for(j in 1:nsurveys_aru) {
                    y_aru[i,j] ~ dbern(p_aru[i])
                    }
                }
                """
        elif self.params.include_scores_data:
            aru_lik = """
                for(j in 1:nsurveys_scores) {
                    scores[i,j] ~ dnorm(mu[z[siteid[i]] + 1], tau[z[siteid[i]] + 1])
                }
                }
                """
        else:
            aru_lik = "\n}\n"

        close_string = """
        }
        model{
            fake <- 0
        }
        """

        data_model_text = f"""
        {string_init}
        {occ_lik}
        {pc_lik}
        {aru_lik}
        {close_string}
        """
        with open(TMP_DATA, "w") as f:
            f.write(data_model_text)

        return data_model_text

    def gen_simulated_data(self):
        """
        Generate simulated data using the JAGS data model
        """
        # find the params that are actually used in the data model
        data_vars = SimModel.find_variables(
            list(self.params.__dict__.keys()), self.data_text
        )
        data = {var: self.params.__dict__[var] for var in data_vars}
        if self.params.include_covar_data:
            data["covar"] = self.covars
            del data["psi"]

        data_model = pyjags.Model(code=self.data_text, data=data, chains=1)
        monitored = {"y_aru", "y_pc", "scores"}
        samples = data_model.sample(iterations=1, vars=monitored)

        sim_data = {var: samples[var][:, :, 0, 0] for var in samples.keys()}
        return sim_data

    def _gen_covars(self) -> np.ndarray:
        """
        Generate covariates
        """
        if self.params.covar_continuous:
            covars = np.random.normal(0, 1, self.params.nsites)
        else:
            covars = np.random.binomial(1, self.params.covar_prob, self.params.nsites)

        return covars


class SimModel(pyjags.model.Model):
    """
    SimModel class to make a model for simulation based on parameters

    This class takes the simulation parameters and generates data using JAGS
    and then runs the model using the generated data.
    """

    def __init__(self, params: SimParams, generated_data: dict, generated_covars: Any):
        self.model_params = params.model_params
        self.data_params = params.data_params
        self.generated_data = generated_data

        # generate the jags model text for the actual model
        self.model_text = self.gen_jags_model_text()

        # Find the variables in the model that are in the simulation data
        # Depending on the model, not all of the simulation data will be used
        all_vars = {
            "nsites",
            "nsurveys_aru",
            "nsurveys_pc",
            "nsurveys_scores",
            "covar",
            "siteid",
            "threshold",
        }

        sim_params = self.model_params.__dict__.copy()

        # add the data params to inits, we don't use them to generate data here though
        # because we reuse the data across models
        for k, v in self.data_params.__dict__.items():
            sim_params[k] = v

        # we add covar after adding the data params because we want the actual generated
        # covar data, not the generic one in the params
        sim_params["covar"] = generated_covars

        model_vars = set(self.find_variables(list(sim_params.keys()), self.model_text))
        model_vars = model_vars.intersection(all_vars)

        data_vars = self.find_variables(list(generated_data.keys()), self.model_text)

        # sim_data holds the data that will be fed to JAGS
        sim_data = {}

        for var in model_vars:
            sim_data[var] = sim_params[var]
        for var in data_vars:
            sim_data[var] = generated_data[var]

        self.inits = self.create_inits(len(generated_covars))

        self.monitored_vars = self.get_monitored_vars()

        super().__init__(
            code=self.model_text,
            data=sim_data,
            chains=self.model_params.nc,
            adapt=params.model_params.na,
            init=self.inits,
            threads=params.model_params.nc,
        )

    def sample(self, iterations=None, vars=None, thin=None, monitor_type="trace"):
        """
        Sample from the model
        """
        # Use instance parameters if not provided
        if iterations is None:
            iterations = self.model_params.ni
        if thin is None:
            thin = self.model_params.nt
        if vars is None:
            vars = self.monitored_vars

        return super().sample(
            iterations=iterations,
            vars=vars,
            thin=thin,
            monitor_type=monitor_type,
        )

    def gen_jags_model_text(self):
        """
        Generate JAGS model text based on the simulation model parameters
        """
        params = self.model_params

        string_init = "model {\n"
        priors_string = "# Priors\n"
        if not (
            params.include_aru_model
            or params.include_pc_model
            or params.include_scores_model
        ):
            raise ValueError("At least one of the models must be included")

        if params.include_covar_model:
            reg_priors = (
                f"beta0 ~ {params.beta0_prior}\n beta1 ~ {params.beta1_prior}\n"
            )
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

        model_text = f"""
        {string_init}
        {priors_string}
        {occ_lik}
        {pc_lik}
        {aru_lik}
        {close_string}
        """

        with open(TMP_MODEL, "w") as f:
            f.write(model_text)

        return model_text

    @staticmethod
    def find_variables(search_vars: List, model_text: str) -> List[str]:
        """
        Find the variables in the model that are in `search_vars`

        Args:
            search_vars (List): The variables to search for in the model text
            model_text (str): The model text to search

        Returns:
            List[str]: The variables in the model from `search_vars`
        """

        return [var for var in search_vars if var in model_text]

    def create_inits(self, nsites: int) -> dict:
        """
        Create initial values for the model

        Args:
            nsites (int): The number of sites in data
        """
        mu = np.random.normal(0, 1, 2)
        mu.sort()

        inits_full = {
            "mu": mu,
            "sigma": np.random.uniform(0.5, 2.5, 2),
            "psi": 0.5,
            "z": np.repeat(1, nsites),
            "beta0": 0,
            "beta1": 0,
            "p_aru11": 0.2,
            "p_aru01": 0.1,
            "p11": 0.2,
        }

        if self.model_params.include_covar_model:
            del inits_full["psi"]
        if not self.model_params.include_covar_model:
            del inits_full["beta0"]
            del inits_full["beta1"]

        init_keys_in_model = self.find_variables(
            list(inits_full.keys()), self.model_text
        )
        inits_in_model = {var: inits_full[var] for var in init_keys_in_model}

        return inits_in_model

    def get_monitored_vars(self):
        """
        Get the monitored variables based on the model type
        """
        monitored_list = ["PropOcc"]
        if self.model_params.include_covar_model:
            monitored_list += ["beta0", "beta1", "mean_psi"]
        else:
            monitored_list += ["psi"]

        if self.model_params.include_pc_model:
            monitored_list += ["p11"]

        if self.model_params.include_aru_model:
            monitored_list += ["p_aru11", "p_aru01"]

        if self.model_params.include_scores_model:
            monitored_list += ["mu", "sigma"]

        return monitored_list
