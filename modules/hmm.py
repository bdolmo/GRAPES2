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

def calculate_positional_mean_variance(sample_list, analysis_dict):
    """ """
    sample_baselines = {}
    df_dict = pd.read_csv(analysis_dict["normalized_depth"], sep="\t").to_dict(
        orient="index"
    )
    for sample in sample_list:
        baseline_samples = []
        num = 0
        sample_depth_tag = ("{}_normalized_final").format(sample.name)
        for control in sample.references:
            if num > 10:
                break
            control_depth_tag = ("{}_normalized_final").format(control[0])
            baseline_samples.append(control_depth_tag)
            num += 1
        sample_baselines[sample.name] = baseline_samples

    observations_dict = defaultdict(dict)
    # here calculate distribution params
    for region in df_dict:
        chromosome = df_dict[region]["chr"]
        coord = ("{}:{}-{}_{}").format(
            df_dict[region]["chr"],
            df_dict[region]["start"],
            df_dict[region]["end"],
            df_dict[region]["exon"],
        )

        if not chromosome in observations_dict:
            observations_dict[chromosome] = []
        region_dict = defaultdict(dict)
        for sample in sample_baselines:

            region_dict[sample] = defaultdict(dict)
            # observations_dict[coord][sample] = defaultdict(dict)
            sample_depth_tag = ("{}_normalized_final").format(sample)
            sample_depth = df_dict[region][sample_depth_tag]
            background_depth = []
            for control in sample_baselines[sample]:
                background_depth.append(df_dict[region][control])

            background_mean = round(np.median(background_depth), 6)
            background_std = round(np.std(background_depth), 6)

            region_dict[sample]["bg_mean"] = background_mean
            region_dict[sample]["bg_std"] = background_std
            region_dict[sample]["normalized_depth"] = sample_depth
            region_dict[sample]["coordinate"] = coord
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
                [0.01, 0.01, 0.96, 0.01, 0.01],
                [0, 0, 0.5, 0.5, 0],
                [0, 0, 0.5, 0, 0.5],
            ]
        )

        self._emissions = emissions
        self._start_prob = np.array([0.20, 0.20, 0.20, 0.20, 0.20])
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
        for chr in self._obs_dict:
            for region in self._obs_dict[chr]:
                # print(region[self._sample]['coordinate'])
                sample_depth = region[self._sample]["normalized_depth"]
                bg_depth = region[self._sample]["bg_mean"]
                bg_std = region[self._sample]["bg_std"]
                mean_bg_std_list.append(bg_std)
        
        mean_bg_std = np.mean(mean_bg_std_list)

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
                combined_error_estimate = np.sqrt(bg_std**2 + mean_bg_std**2)

                if x == 0:
                    x = 1
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
        alpha[0, :] = self._start_prob * np.array(self._emissions[0])

        for t in range(1, T):
            for j in range(M):
                alpha[t, j] = alpha[t - 1].dot(self._transitions[:, j]) * self._emissions[t][j]

        return alpha

    def backward(self):
        """
        Backward algorithm.
        """
        T = len(self.observations)
        M = self._n_components

        beta = np.zeros((T, M))

        # Initialization
        beta[T - 1] = np.ones((M))

        for t in range(T - 2, -1, -1):
            for j in range(M):
                beta[t, j] = (beta[t + 1] * np.array(self._emissions[t + 1]) * self._transitions[j, :]).sum()

        return beta

    def posterior_decoding(self):
        """
        Posterior decoding (forward-backward algorithm).
        """
        alpha = self.forward()
        beta = self.backward()

        posterior_probs = np.multiply(alpha, beta) / np.sum(np.multiply(alpha, beta), axis=1)[:, np.newaxis]

        # print(self._emissions[0])
        # print(posterior_probs[0])
        sys.exit()

        return posterior_probs


    def decode(self):
        """
        O: Observations
        S: States (copy number: 0, 1, 2, 3, 4)
        X: Probability path
        A: State transition matrix
        B: Emission matrix
        """

        start_p = np.array([0.25, 0.25, 0.25, 0.25, 0.25])
        O = np.array(self.observations)
        S = np.array([0, 1, 2, 3, 4])
        A = np.array(
            [
                [0.5, 0, 0.5, 0, 0],
                [0, 0.5, 0.5, 0, 0],
                [0.005, 0.02, 0.95, 0.02, 0.005],
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

                # This is our most probable state given previous state at time t (1)
                prev[t - 1, j] = np.argmax(probability)

                # This is the probability of the most probable state (2)
                omega[t, j] = np.max(probability)

                # Workaround when all state probs are inifite, this caused to immediately
                # assign a 0 state
                if np.isinf(omega[t, j]):
                    omega[t, j] = -10000000
            # print(self.list_of_rois[t]['coordinate'], phred_scores[t] )

        # print(omega)

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

        return result, phred_scores


if __name__ == "__main__":

    model = CustomHMM()
    obs = np.random.normal(1, 0.2, 1000)
    # mean = np.mean(obs)
    model.compute_log_likelihood(obs)
