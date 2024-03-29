import polars as pl
import numpy as np
import random

def set_seed(seed=42):
    """Set the seed for reproducibility."""
    random.seed(seed)
    np.random.seed(seed)
    pl.set_random_seed(seed)
