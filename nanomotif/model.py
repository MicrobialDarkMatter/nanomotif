"""
Models
"""
import numpy as np
from scipy.stats import beta
from scipy.special import psi

DEFAULT_PRIOR_BETA = 5
DEFAULT_PRIOR_ALPHA = 5

class BetaBernoulliModel():
    def __init__(self, alpha = DEFAULT_PRIOR_ALPHA, beta = DEFAULT_PRIOR_BETA):
        self._alpha = alpha
        self._beta = beta
        self._alpha_prior = alpha
        self._beta_prior = beta

    def __getstate__(self):
        # Return a serializable representation of the state
        return {
            "_alpha": self._alpha,
            "_beta": self._beta,
            "_alpha_prior": self._alpha_prior,
            "_beta_prior": self._beta_prior,
        }

    def __setstate__(self, state):
        # Restore from serialized state
        self._alpha = state["_alpha"]
        self._beta = state["_beta"]
        self._alpha_prior = state["_alpha_prior"]
        self._beta_prior = state["_beta_prior"]

    def get_raw_counts(self):
        return self._alpha - self._alpha_prior, self._beta - self._beta_prior

    def update(self, n_positives, n_negatives):
        self._alpha += n_positives
        self._beta += n_negatives

    def reset(self):
        self._alpha = self._alpha_prior
        self._beta = self._beta_prior

    def sample(self):
        return np.random.beta(self._alpha, self._beta)

    def mean(self):
        return self._alpha / (self._alpha + self._beta)

    def variance(self):
        a = self._alpha
        b = self._beta
        return (a * b) / ((a + b)**2 * (a + b + 1))
    
    def standard_deviation(self):
        return np.sqrt(self.variance())
    
    def log_likelihood(self):
        return (self._alpha * np.log(self.mean()) +
                self._beta * np.log(1 - self.mean()))
    
    def log_likelihood_ratio(self, other):
        return self.log_likelihood() - other.log_likelihood()

    def log_likelihood_per_obs(self):
        n_obs = self._alpha + self._beta
        if n_obs == 0:
            return 0.0
        return self.log_likelihood() / n_obs

    def log_likelihood_ratio_per_obs(self, other):
        n_obs = (self._alpha + self._beta) - (other._alpha + other._beta)
        if n_obs == 0:
            return 0.0
        return self.log_likelihood_ratio(other) / n_obs
    
    def posterior_predictive(self, n_positives, n_negatives):
        n_new = n_positives + n_negatives
        if n_new == 0:
            return 0.0
        E_log_p = psi(self._alpha) - psi(self._alpha + self._beta)
        E_log_1mp = psi(self._beta) - psi(self._alpha + self._beta)
        ll_pred = n_positives * E_log_p + n_negatives * E_log_1mp
        return ll_pred

    def posterior_predictive_per_obs(self, n_positives, n_negatives):
        n_new = n_positives + n_negatives
        if n_new == 0:
            return 0.0
        posterior_preditive = self.posterior_predictive(n_positives, n_negatives)
        return posterior_preditive / n_new
    
    def posterior_predictive_check(self, n_positives, n_negatives, num_simulations=1000):
        n_new = n_positives + n_negatives

        simulated_lls = []
        for _ in range(num_simulations):
            p_sim = self.sample()
            # simulate number of positives in n_new observations
            y_sim_positive = np.random.binomial(n_new, p_sim)
            y_sim_negative = n_new - y_sim_positive

            # log-likelihood of simulated data
            ll_sim = y_sim_positive * np.log(p_sim + 1e-12) + y_sim_negative * np.log(1 - p_sim + 1e-12)
            simulated_lls.append(ll_sim)

        simulated_lls = np.array(simulated_lls)

        # log-likelihood of actual new data under model
        p_mean = self.mean()
        observed_ll = n_positives * np.log(p_mean + 1e-12) + n_negatives * np.log(1 - p_mean + 1e-12)

        # posterior predictive p-value
        p_ppc = np.mean(simulated_lls <= observed_ll)

        return p_ppc

    def cdf(self, x):
        return beta.cdf(x, self._alpha, self._beta)

    def __str__(self):
        return f'BetaBernoulliModel(alpha={self._alpha}, beta={self._beta})'

    def __repr__(self):
        return str(self)


