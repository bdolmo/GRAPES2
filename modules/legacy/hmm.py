import os
import sys
import re
from hmmlearn import hmm, base
import numpy as np
import string
import pandas as pd
from collections import defaultdict
from scipy.spatial import distance
from scipy.stats import nbinom, poisson, norm
from scipy.special import logsumexp
import math 

def calculate_positional_mean_variance(sample_list, analysis_dict):
    """Build per-chromosome observation dictionaries for HMM emissions."""
    normalized_bed = analysis_dict["normalized_depth"]
    if analysis_dict.get("offtarget", False):
        normalized_bed = analysis_dict["normalized_all"]

    rows = pd.read_csv(normalized_bed, sep="\t").to_dict(orient="records")

    sample_baselines = {}
    for sample in sample_list:
        baseline_samples = []
        for idx, control in enumerate(sample.references):
            if idx > 10:
                break
            baseline_samples.append(f"{control[0]}_normalized_final")
        sample_baselines[sample.name] = baseline_samples

    observations_dict = defaultdict(list)
    for row in rows:
        chromosome = row["chr"]
        exon = row["exon"]
        coord = f"{row['chr']}:{row['start']}-{row['end']}_{exon}"
        is_offtarget = "pwindow" in exon

        region_dict = {}
        for sample_name, baseline_cols in sample_baselines.items():
            sample_depth_tag = f"{sample_name}_normalized_final"
            sample_depth = row[sample_depth_tag]
            background_depth = [row[col] for col in baseline_cols]
            background_mean = round(np.median(background_depth), 6)
            background_std = round(np.std(background_depth), 6)

            region_dict[sample_name] = {
                "bg_mean": background_mean,
                "bg_std": background_std,
                "normalized_depth": sample_depth,
                "coordinate": coord,
                "is_offtarget": is_offtarget,
            }

        observations_dict[chromosome].append(region_dict)
    return observations_dict

def convert_params(mu, alpha):
    """ 
    Convert mean/dispersion parameterization of a negative binomial to the ones scipy supports

    Parameters
    ----------
    mu : float 
       Mean of NB distribution.
    alpha : float
       Overdispersion parameter used for variance calculation.

    See https://en.wikipedia.org/wiki/Negative_binomial_distribution#Alternative_formulations
    """
    var = mu + alpha * mu ** 2
    p = (var - mu) / var
    r = mu ** 2 / (var - mu)
    return r, p

class CustomHMM:
    """ """

    def __init__(
        self,
        obs_dict,
        sample,
        chr,
        n_components=5,
        transitions=None,
        emissions=None,
        start_prob=None,
    ):
        self._obs_dict = obs_dict
        self._sample = sample
        self._chr = chr
        self._n_components = n_components
        self._transitions = np.array(
            [
                [0.5, 0, 0.5, 0, 0],
                [0, 0.5, 0.5, 0, 0],
                [0.005, 0.045, 0.90, 0.045, 0.005],
                [0, 0, 0.5, 0.5, 0],
                [0, 0, 0.5, 0, 0.5],
            ]
        )

        self._emissions = emissions
        self._start_prob = np.array([0.005, 0.02, 0.95, 0.02, 0.005])
        self._emissions = self.compute_log_likelihood()

    @property
    def observations(self):
        """ """
        obs = []
        for region in self._obs_dict[self._chr]:          
            x = region[self._sample]["normalized_depth"]
            obs.append(x)
        return obs

    @property
    def list_of_rois(self):
        """ """
        regions = []
        for region in self._obs_dict[self._chr]:
           regions.append(region)
        return regions

    def compute_log_likelihood(self):
        """ """
        logp_dict = defaultdict(dict)
        emissions = []
        emissions_probs = []
        idx = 0
        epsilon = 1e-8

        mean_bg_std_list = []
        mean_bg_depth_list = []

        mean_bg_std_offtarget_list = []

        for chr in self._obs_dict:
            for region in self._obs_dict[chr]:
                sample_depth = region[self._sample]["normalized_depth"]
                bg_depth = region[self._sample]["bg_mean"]
                bg_std = region[self._sample]["bg_std"]
                if math.isnan(bg_std):
                    bg_std = 0.2

                if region[self._sample]["is_offtarget"] == False:
                    mean_bg_std_list.append(bg_std)
                else:
                    mean_bg_std_offtarget_list.append(bg_std)
        mean_bg_std = np.mean(mean_bg_std_list)
        mean_bg_std_offtarget = np.mean(mean_bg_std_offtarget_list)

        for region in self._obs_dict[self._chr]:
            tmp = re.split(':|-|_', region[self._sample]['coordinate'])
            pos = tmp[1]
            end = tmp[2]

            sample_depth = region[self._sample]["normalized_depth"]
            bg_depth = region[self._sample]["bg_mean"]
            bg_std = region[self._sample]["bg_std"]
            coord = region[self._sample]["coordinate"]
            if bg_std < epsilon:
                bg_std = epsilon

            logp_dict[idx] = defaultdict(dict)
            state_list = []
            depth_list = []
            probs_list = []
            for state in range(self._n_components):

                if state == 0:
                    state = 0.01

                ratio = state / 2
                x = sample_depth
                mean_state = bg_depth * ratio

                if region[self._sample]["is_offtarget"] == False:
                    combined_error_estimate = np.sqrt(bg_std**2 + mean_bg_std**2)
                else:
                    combined_error_estimate = np.sqrt(bg_std**2 + mean_bg_std_offtarget**2)

                if x == 0:
                    x = 0.01
                logp_dict[idx][state] = round(
                    np.log(norm.pdf(x, loc=mean_state, scale=combined_error_estimate)+epsilon), 6
                )
                probs_list.append(norm.pdf(x, loc=bg_depth, scale=bg_std+mean_bg_std))

                state_list.append(logp_dict[idx][state])
                depth_list.append(str(state) + ":" + str(x))
            emissions.append(state_list)
            idx += 1

        return emissions

    def forward(self):
        """
        Forward algorithm.
        """
        T = len(self.observations)
        M = self._n_components

        alpha = np.zeros((T, M))

        # Initialization
        alpha[0, :] = np.log(self._start_prob) + np.array(self._emissions[0])
        epsilon = 1e-8
        for t in range(1, T):
            for j in range(M):
                alpha[t, j] = logsumexp(alpha[t - 1] + np.log(self._transitions[:, j]+epsilon)) + self._emissions[t][j]

        return alpha

    def backward(self):
        """
        Backward algorithm.
        """
        T = len(self.observations)
        M = self._n_components

        beta = np.zeros((T, M))
        epsilon = 1e-8

        # Initialization
        beta[T - 1, :] = 0  # in log space, log(1) = 0

        for t in range(T - 2, -1, -1):
            for j in range(M):
                beta[t, j] = logsumexp(beta[t + 1] + np.log(self._transitions[j, :]+epsilon) + np.array(self._emissions[t + 1]))

        return beta

    def posterior_decoding(self):
        """
        Posterior decoding (forward-backward algorithm).
        """
        alpha = self.forward()
        beta = self.backward()

        # Compute posterior probabilities in log space
        log_posterior_probs = alpha + beta - logsumexp(alpha + beta, axis=1, keepdims=True)
        log_posterior_probs = alpha + beta- logsumexp(alpha + beta, axis=1, keepdims=True)

        # Convert back to normal space if necessary
        posterior_probs = np.exp(log_posterior_probs)
        return posterior_probs


    def calculate_map(self):
        """
        Calculate the Maximum A Posteriori (MAP) estimate for each region.

        Returns
        -------
        map_states : list
            A list of states representing the most probable state for each region.
        map_probabilities : list
            A list of the maximum posterior probabilities for each region.
        """
        posterior_probs = self.posterior_decoding()
        map_states = np.argmax(posterior_probs, axis=1)  # Most probable state for each region
        map_probabilities = np.max(posterior_probs, axis=1)  # Maximum posterior probability for each region

        return map_states, map_probabilities


    def decode(self):
        """
        O: Observations
        S: States (copy number: 0, 1, 2, 3, 4)
        X: Probability path
        A: State transition matrix
        B: Emission matrix
        """

        start_p = np.array([0.005, 0.02, 0.95, 0.02, 0.005])
        O = np.array(self.observations)
        S = np.array([0, 1, 2, 3, 4])
        A = np.array(
            [
                [0.5, 0, 0.5, 0, 0],
                [0, 0.5, 0.5, 0, 0],
                [0.005, 0.045, 0.90, 0.045, 0.005],
                # [0.005, 0.02, 0.95, 0.02, 0.005],
                [0, 0, 0.5, 0.5, 0],
                [0, 0, 0.5, 0, 0.5],
            ]
        )
        B = np.array(self._emissions)
        T = O.shape[0]
        M = S.shape[0]

        # Initialization
        epsilon = 1e-8  # a small constant
        omega = np.zeros((T, M))
        omega[0, :] = np.log(start_p + epsilon) + B[0, :]

        prev = np.zeros((T - 1, M))
        phred_scores = np.zeros((T, M))
        prob_scores = np.zeros((T, M))
        for t in range(1, T):
            list_phreds = []
            list_probabilites = []
            for j in range(M):
                probability = omega[t - 1] + np.log(A[:, j] + epsilon) + B[t, j]
                max_log_prob = np.max(probability)

                # Subtract max log probability for numerical stability, exponentiate and normalize
                probs = np.exp(probability - max_log_prob)
                probs /= np.sum(probs)
                # Calculate error probabilities
                error_probs = 1 - probs

                # Calculate Phred scores
                Q = -10 * np.log10(error_probs)
                # Round to nearest integer
                Q_rounded = np.round(Q)

                # Cap at 60
                Q_capped = np.clip(Q_rounded, 0, 60)
                phred_scores[t, j] = np.max(Q_capped)
                prob_scores[t, j] = np.max(probs)

                # This is our most probable state given previous state at time t (1)
                prev[t - 1, j] = np.argmax(probability)

                # This is the probability of the most probable state (2)
                omega[t, j] = np.max(probability)

                # Workaround when all state probs are inifite, this caused to immediately
                # assign a 0 state
                if np.isinf(omega[t, j]):
                    omega[t, j] = -10000000

        # Path Array
        X = np.zeros(T)

        # Find the most probable last hidden state
        last_state = np.argmax(omega[T - 1, :])

        X[0] = last_state

        backtrack_index = 1
        for i in range(T - 2, -1, -1):
            X[backtrack_index] = int(prev[i, int(last_state)])
            last_state = prev[i, int(last_state)]
            backtrack_index += 1

        # Flip the path array since we were backtracking
        X = np.flip(X, axis=0)
        # Convert numeric values to actual hidden states
        result = []
        for x in X:
            result.append(str(int(x)))

        return result, prob_scores


if __name__ == "__main__":

    model = CustomHMM()
    obs = np.random.normal(1, 0.2, 1000)
    # mean = np.mean(obs)
    model.compute_log_likelihood(obs)
