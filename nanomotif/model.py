"""
Models
"""
import numpy as np
from scipy.stats import beta

DEFAULT_PRIOR_BETA = 10
DEFAULT_PRIOR_ALPHA = 10

class BetaBernoulliModel():
    def __init__(self, alpha = DEFAULT_PRIOR_ALPHA, beta = DEFAULT_PRIOR_BETA):
        self.alpha = alpha
        self.beta = beta
        self._alpha = alpha
        self._beta = beta

    def get_raw_counts(self):
        return self._alpha - DEFAULT_PRIOR_ALPHA, self._beta - DEFAULT_PRIOR_BETA

    def update(self, n_positives, n_negatives):
        self._alpha += n_positives
        self._beta += n_negatives

    def reset(self):
        self._alpha = self.alpha
        self._beta = self.beta

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
        p = self.mean()
        log_likelihood = n_positives * np.log(p + 1e-12) + n_negatives * np.log(1 - p + 1e-12)
        return log_likelihood

    def posterior_predictive_per_obs(self, n_positives, n_negatives):
        n_new = n_positives + n_negatives
        if n_new == 0:
            return 0.0
        log_likelihood = self.posterior_predictive(n_positives, n_negatives)
        return log_likelihood / n_new
    
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


