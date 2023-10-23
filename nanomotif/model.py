"""
Models
"""
import numpy as np
from scipy.stats import beta
import matplotlib.pyplot as plt

class BetaBernoulliModel():
    def __init__(self, alpha = 0, beta = 1):
        self.alpha = alpha
        self.beta = beta
        self._alpha = alpha
        self._beta = beta

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
    
    def cdf(self, x):
        return beta.cdf(x, self._alpha, self._beta)

    def __str__(self):
        return f'BetaBernoulliModel(alpha={self._alpha}, beta={self._beta})'

    def __repr__(self):
        return str(self)
    
    def plot(self, ax = None):
        if ax is None:
            fig, ax = plt.subplots()
        x = np.linspace(0, 1, 1000)
        ax.plot(x, beta.pdf(x, self._alpha, self._beta), 'r-')
        ax.set_title('Beta Distribution', fontsize='15')
        ax.set_xlabel('Values of Random Variable X (0, 1)', fontsize='15')
        ax.set_ylabel('Probability', fontsize='15')
        return ax

