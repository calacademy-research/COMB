import numpy as np
from modeling.models.model_iterface import SimulationParams
from .model_iterface import CombinedModelInterface
from arviz import InferenceData
from caples_data import COMBData
import arviz as az
import pyjags


class JagsModelIndependent(CombinedModelInterface):
    @classmethod
    def run_model(cls, data: COMBData) -> InferenceData:
        return JagsModel.run_model(data, aru_scores_independent_model=True)

    @classmethod
    def simulate_data(cls, sim_params: SimulationParams) -> COMBData:
        return JagsModel.simulate_data(sim_params)


class JagsModelDependent(CombinedModelInterface):
    @classmethod
    def run_model(cls, data: COMBData) -> InferenceData:
        return JagsModel.run_model(data, aru_scores_independent_model=False)

    @classmethod
    def simulate_data(cls, sim_params: SimulationParams) -> COMBData:
        return JagsModel.simulate_data(sim_params)


class JagsModel:
    @classmethod
    def run_model(
        cls, data: COMBData, aru_scores_independent_model: bool
    ) -> InferenceData:
        model_text = gen_jags_model_text(aru_scores_independent_model)
        model_data = {
            "y_pc": data.y_pc,
            "y_aru": data.y_aru,
            "nsurveys_aru": data.n_surveys_aru,
            "nsurveys_pc": data.n_surveys_pc,
            "nsurveys_scores": data.n_surveys_scores,
            "scores": data.scores,
            "covar": data.covariates["covar"],
        }

        mu = np.random.normal(0, 1, 2)
        mu.sort()
        inits = {
            "mu": mu,
            "sigma": np.random.uniform(0.5, 2.5, 2),
            "psi": 0.5,
            "z": np.repeat(1, data.n_sites),
            "beta0": 0,
            "beta1": 0,
            "p_aru11": 0.2,
            "p_aru01": 0.1,
            "p11": 0.2,
        }
        model = pyjags.Model(
            code=model_text,
            data=model_data,
            init=inits,
            chains=4,
            adapt=2000,
            threads=4,
        )

        monitored = [
            "mu",
            "sigma",
            "p_aru11",
            "p_aru01",
            "p11",
            "beta0",
            "beta1",
            "mean_psi",
            "PropOcc",
        ]

        samples = model.sample(iterations=8000, vars=monitored)
        return az.from_pyjags(samples)

    @classmethod
    def simulate_data(cls, sim_params: SimulationParams) -> COMBData:
        covars = np.random.normal(0, 1, sim_params.nsites)

        jags_text = gen_jags_data_text(sim_params)

        data = {
            "beta0": sim_params.beta0,
            "beta1": sim_params.beta1,
            "p11": sim_params.p11,
            "p_aru11": sim_params.p_aru11,
            "p_aru01": sim_params.p_aru01,
            "threshold": sim_params.threshold,
            "covar": covars,
        }
        monitored = {"y_aru", "y_pc", "scores"}
        data_model = pyjags.Model(jags_text, data, chains=1)

        samples = data_model.sample(iterations=1, vars=monitored)

        return COMBData(
            n_sites=sim_params.nsites,
            n_years=1,
            n_species=1,
            y_aru=samples["y_aru"][:, :, 0, 0],
            y_index=samples["y_pc"][:, :, 0, 0],
            y_pc=samples["y_pc"][:, :, 0, 0],
            n_surveys_aru=sim_params.nsurveys_aru,
            n_surveys_pc=sim_params.nsurveys_pc,
            date_aru=np.zeros(0),
            date_pc=np.zeros(0),
            scores=samples["scores"][:, :, 0, 0],
            n_surveys_scores=sim_params.nsurveys_scores,
            covariates={"covar": covars},
            time_aru=np.zeros(0),
            time_pc=np.zeros(0),
        )


def gen_jags_data_text(sim_params: SimulationParams):
    string_init = "data {\n"
    occ_lik = """
    for (i in 1:nsites) {
        logit(psi[i]) <- beta0 + beta1*covar[i]
        z[i] ~ dbern(psi[i])
    """

    pc_lik = """
        p[i] <- p11 * z[i]
        for (j in 1:nsurveys_pc) {
            y_pc[i, j] ~ dbern(p[i])
        }
        """
    if sim_params.aru_data_independent_model:
        aru_lik = """
            # ARU detection and scores for independent model
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
            # ARU: detection + scores for dependent model
            p_aru[i] <- z[i]*p_aru11 + p_aru01
            for(j in 1:nsurveys_aru) {
            y_aru[i,j] ~ dbern(p_aru[i])
            scores[i,j] ~ dnorm(mu[z[siteid[i]] + 1], tau[z[siteid[i]] + 1]) T(ifelse(y_aru[i,j] == 1, threshold, -10), ifelse(y_aru[i,j] == 1, 10, threshold))
            }
        }"""

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

    return data_model_text


def gen_jags_model_text(aru_scores_independent_model: bool):
    string_init = "model {\n"
    priors_string = "# Priors\n"

    reg_priors = "beta0 ~ dnorm(0, 0.5)\n beta1 ~ dnorm(0, 0.5)\n"
    priors_string += reg_priors

    occ_lik = """
    # Likelihood for occupancy
    for (i in 1:nsites) {
        logit(psi[i]) <- beta0 + beta1 * covar[i]
        z[i] ~ dbern(psi[i])
    """

    pc_lik = """
    # Point count
    p[i] <- p11 * z[i]
    for (j in 1:nsurveys_pc) {
        y_pc[i, j] ~ dbern(p[i])
    }
    """
    p11_prior_string = "p11 ~ dbeta(2, 2)\n"
    priors_string += p11_prior_string

    if aru_scores_independent_model:
        aru_lik = """
        # ARU independent model
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
        # ARU dependent model
        p_aru[i] <- z[i]*p_aru11 + p_aru01
        for(j in 1:nsurveys_aru) {
            y_aru[i,j] ~ dbern(p_aru[i])
            scores[i,j] ~ dnorm(mu[z[siteid[i]] + 1], tau[z[siteid[i]] + 1]) T(ifelse(y_aru[i,j] == 1, threshold, -10), ifelse(y_aru[i,j] == 1, 10, threshold))
        }
        }"""

    p_aru11_prior_string = "p_aru11 ~ dbeta(2, 2)\n"
    p_aru01_prior_string = "p_aru01 ~ dbeta(1, 3)I(0, 1 - p_aru11)\n"
    scores_prior = """
        # Priors of the observation model for the scores
        mu[1] ~ dnorm(-2, 0.2)
        mu[2] ~ dnorm(-2, 0.2)
        sigma[1] ~ dunif(0.1, 5)
        tau[1] <- 1 / (sigma[1] * sigma[1])
        sigma[2] ~ dunif(0.1, 5)
        tau[2] <- 1 / (sigma[2] * sigma[2])
    """
    priors_string += p_aru11_prior_string + p_aru01_prior_string + scores_prior

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

    return model_text
