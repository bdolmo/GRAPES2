import os
import sys
import re
import numpy as np
import pandas as pd
from collections import defaultdict
from scipy.special import logsumexp, gammaln
import math
import logging

########################################################################
# Original utility functions
########################################################################

def calculate_positional_mean_variance(sample_list, analysis_dict):
    """
    Calculate positional means and variances using normalized depth data.
    This function reads a BED‐like file (tab–separated) that contains
    normalized depth information per target/region.
    """
    sample_baselines = {}
    normalized_bed = analysis_dict["normalized_depth"]
    if analysis_dict.get("offtarget", False):
        normalized_bed = analysis_dict["normalized_all"]
    df_dict = pd.read_csv(normalized_bed, sep="\t").to_dict(orient="index")
    for sample in sample_list:
        baseline_samples = []
        num = 0
        sample_depth_tag = f"{sample.name}_normalized_final"
        for control in sample.references:
            if num > 10:
                break
            control_depth_tag = f"{control[0]}_normalized_final"
            baseline_samples.append(control_depth_tag)
            num += 1
        sample_baselines[sample.name] = baseline_samples

    observations_dict = defaultdict(dict)
    for region in df_dict:
        chromosome = df_dict[region]["chr"]
        coord = ("{}:{}-{}_{}").format(
            df_dict[region]["chr"],
            df_dict[region]["start"],
            df_dict[region]["end"],
            df_dict[region]["exon"],
        )
        if chromosome not in observations_dict:
            observations_dict[chromosome] = []
        region_dict = defaultdict(dict)
        for sample in sample_baselines:
            region_dict[sample] = defaultdict(dict)
            sample_depth_tag = f"{sample}_normalized_final"
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
            region_dict[sample]["is_offtarget"] = False
            if "pwindow" in df_dict[region]["exon"]:
                region_dict[sample]["is_offtarget"] = True
        observations_dict[chromosome].append(region_dict)
    return observations_dict

def convert_params(mu, alpha):
    """
    Given a negative binomial with mean mu and variance mu + alpha * mu^2,
    use the following parameterization:
      r = 1/alpha
      p = 1 / (1 + alpha * mu)
    """
    r = 1.0 / alpha
    p = 1.0 / (1.0 + alpha * mu)
    return r, p

########################################################################
# Helper function: update_state_dispersion
########################################################################

def update_state_dispersion(obs_dict, sample, chr, state_assignments, scale_factor=100, epsilon=1e-8):
    """
    For each state, compute a dispersion estimate by pooling the background
    parameters (bg_mean and bg_std) across exons assigned to that state.
    
    For each exon, we use a method‐of–moments estimator:
       alpha_i = max(((bg_std*scale_factor)^2 - (bg_mean*scale_factor)) / ((bg_mean*scale_factor)^2), epsilon)
    
    Then, for each state, return the median of the alpha_i values.
    """
    dispersions = {}
    n_states = max(state_assignments) + 1
    for s in range(n_states):
        estimates = []
        for region, state in zip(obs_dict[chr], state_assignments):
            if state == s:
                bg_mean = region[sample]["bg_mean"]
                bg_std = region[sample]["bg_std"]
                bg_mean_scaled = bg_mean * scale_factor
                bg_std_scaled = bg_std * scale_factor
                if bg_mean_scaled > epsilon:
                    alpha_i = (bg_std_scaled**2 - bg_mean_scaled) / (bg_mean_scaled**2)
                    alpha_i = max(alpha_i, epsilon)
                    estimates.append(alpha_i)
        if estimates:
            dispersions[s] = np.median(estimates)
        else:
            dispersions[s] = 0.01
    disp_array = np.zeros(n_states)
    for s in range(n_states):
        disp_array[s] = dispersions.get(s, 0.01)
    return disp_array



from scipy.optimize import minimize_scalar

def estimate_global_dispersion_MLE(obs_dict, sample, chr, scale_factor=100, epsilon=1e-8):
    """
    Estimate a single global dispersion parameter (alpha) for a given sample and chromosome
    by maximizing the overall likelihood of the NB model across all exons.
    
    For each exon, let:
      x  = normalized_depth * scale_factor
      μ  = bg_mean * scale_factor
    and assume the NB likelihood is:
      log P(x|μ,α) = gammaln(x + r) - gammaln(r) - gammaln(x+1) +
                     r * log(p) + x * log(1-p)
      where r = 1/α and p = 1/(1 + α * μ).
    
    This function defines the negative log–likelihood over all exons and minimizes it.
    
    Returns:
      A global dispersion value (α) that is the MLE across exons.
    """
    def total_neg_log_likelihood(alpha):
        total_nll = 0
        for region in obs_dict[chr]:
            # Scale the observed and expected counts
            x = region[sample]["normalized_depth"] * scale_factor
            mu = region[sample]["bg_mean"] * scale_factor
            if mu < epsilon:
                continue  # skip regions with essentially no signal
            if alpha < epsilon:
                return np.inf
            r = 1.0 / alpha
            p = 1.0 / (1.0 + alpha * mu)
            # NB negative log-likelihood for this exon
            nll = -(gammaln(x + r) - gammaln(r) - gammaln(x + 1) +
                    r * np.log(p + epsilon) + x * np.log(1 - p + epsilon))
            total_nll += nll
        return total_nll

    res = minimize_scalar(total_neg_log_likelihood, bounds=(epsilon, 1.0), method='bounded')
    if res.success:
        # You can also enforce a minimum value if desired (e.g. 0.05)
        return max(res.x, 0.01)
    else:
        return 0.01


########################################################################

def update_global_dispersion(obs_dict, sample, chr, scale_factor=100, epsilon=1e-8):
    """
    Compute a global dispersion parameter by pooling background estimates from all exons.
    
    For each exon, we estimate dispersion as:
         alpha_i = max( ( (bg_std*scale_factor)^2 - (bg_mean*scale_factor) ) / ((bg_mean*scale_factor)^2), epsilon )
    and then return the median of these estimates.
    """
    estimates = []
    min_alpha = 0.05
    
    for region in obs_dict[chr]:
        bg_mean = region[sample]["bg_mean"]
        bg_std = region[sample]["bg_std"]
        bg_mean_scaled = bg_mean * scale_factor
        bg_std_scaled = bg_std * scale_factor

        print(bg_mean_scaled)

        if bg_mean_scaled > epsilon:
            alpha_i = abs((bg_std_scaled**2 - bg_mean_scaled) / (bg_mean_scaled**2))
            print("alpha_i", alpha_i)
            alpha_i = max(alpha_i, epsilon)

            estimates.append(alpha_i)
    print("before global_dispersion", np.median(estimates))
    if estimates:
        median_estimates = np.median(estimates)
        return max(median_estimates, min_alpha)
        # if median_estimates < 0.01:
        #     return 0.01
        # else:
        #     return median_estimates
    else:
        return min_alpha



########################################################################
# CustomHMM class (Frequentist Hierarchical Approach)
########################################################################

class CustomHMM:
    """Hidden Markov Model with Negative Binomial emission probabilities.
    
    This version scales float normalized depths into count-like data and
    can learn a (global or state-specific) dispersion parameter via an EM-like procedure.
    """
    def __init__(
        self,
        obs_dict,
        sample,
        chr,
        n_components=5,
        transitions=None,
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
                [0.025, 0.025, 0.99, 0.025, 0.025],
                [0, 0, 0.5, 0.5, 0],
                [0, 0, 0.5, 0, 0.5],
            ]
        ) if transitions is None else transitions
        self._start_prob = np.array([0.025, 0.025, 0.99, 0.025, 0.025]) if start_prob is None else start_prob
        self._scale_factor = 100
        self._emissions = self.compute_log_likelihood()
        self._state_dispersion = None

    @property
    def observations(self):
        return [region[self._sample]["normalized_depth"] for region in self._obs_dict[self._chr]]

    @property
    def list_of_rois(self):
        return [region for region in self._obs_dict[self._chr]]

    def compute_log_likelihood(self, state_dispersion=None):
        logp_dict = defaultdict(dict)
        emissions = []
        epsilon = 1e-8
        sf = self._scale_factor

        # Global background std for smoothing (scaled)
        mean_bg_std_list = []
        mean_bg_std_offtarget_list = []
        for chr_key in self._obs_dict:
            for region in self._obs_dict[chr_key]:
                bg_std = region[self._sample]["bg_std"]
                if math.isnan(bg_std):
                    bg_std = 0.2
                bg_std_scaled = bg_std * sf
                if not region[self._sample]["is_offtarget"]:
                    mean_bg_std_list.append(bg_std_scaled)
                else:
                    mean_bg_std_offtarget_list.append(bg_std_scaled)
        mean_bg_std = np.mean(mean_bg_std_list)

        if mean_bg_std_offtarget_list:
            mean_bg_std_offtarget = np.mean(mean_bg_std_offtarget_list)

        for region in self._obs_dict[self._chr]:
            sample_depth = region[self._sample]["normalized_depth"] * sf
            bg_depth = region[self._sample]["bg_mean"] * sf
            bg_std = region[self._sample]["bg_std"] * sf
            if bg_std < epsilon:
                bg_std = epsilon

            state_list = []
            for state in range(self._n_components):
                effective_ratio = 0.01 if state == 0 else state / 2.0
                mu_state = bg_depth * effective_ratio
                if mu_state < epsilon:
                    mu_state = epsilon

                if not region[self._sample]["is_offtarget"]:
                    effective_std = np.sqrt(bg_std**2 + mean_bg_std**2)
                else:
                    effective_std = np.sqrt(bg_std**2 + mean_bg_std_offtarget**2)

                if state_dispersion is not None:
                    alpha_used = state_dispersion[state]
                else:
                    if bg_depth > epsilon:
                        alpha_used = max((bg_std**2 - bg_depth) / (bg_depth**2), 0.001)
                    else:
                        alpha_used = 0.01
                # print("dispersion_alpha:", alpha_used)
                r, p_param = convert_params(mu_state, alpha_used)
                x_obs = sample_depth
                if x_obs < 0:
                    x_obs = epsilon
                logp = (gammaln(x_obs + r) - gammaln(r) - gammaln(x_obs + 1) +
                        r * np.log(p_param + epsilon) + x_obs * np.log(1 - p_param + epsilon))
                logp = round(logp, 6)
                state_list.append(logp)
            emissions.append(state_list)
        return emissions

    def forward(self):
        T = len(self.observations)
        M = self._n_components
        alpha = np.zeros((T, M))
        epsilon = 1e-8
        alpha[0, :] = np.log(self._start_prob + epsilon) + np.array(self._emissions[0])
        for t in range(1, T):
            for j in range(M):
                alpha[t, j] = logsumexp(alpha[t - 1] + np.log(self._transitions[:, j] + epsilon)) + self._emissions[t][j]
        return alpha

    def backward(self):
        T = len(self.observations)
        M = self._n_components
        beta = np.zeros((T, M))
        epsilon = 1e-8
        beta[T - 1, :] = 0
        for t in range(T - 2, -1, -1):
            for j in range(M):
                beta[t, j] = logsumexp(beta[t + 1] + np.log(self._transitions[j, :] + epsilon) + np.array(self._emissions[t + 1]))
        return beta

    def posterior_decoding(self):
        alpha = self.forward()
        beta = self.backward()
        log_post = alpha + beta - logsumexp(alpha + beta, axis=1, keepdims=True)
        return np.exp(log_post)

    def calculate_map(self):
        posterior_probs = self.posterior_decoding()
        map_states = np.argmax(posterior_probs, axis=1)
        map_probs = np.max(posterior_probs, axis=1)
        return map_states, map_probs

    def decode(self):
        """
        Viterbi decoding algorithm to infer the most probable state path.
        """
        start_p = self._start_prob
        O = np.array(self.observations)
        S = np.arange(self._n_components)
        A = np.array(self._transitions)
        B = np.array(self._emissions)
        T = O.shape[0]
        M = S.shape[0]
        epsilon = 1e-16

        omega = np.zeros((T, M))
        omega[0, :] = np.log(start_p + epsilon) + B[0, :]

        prev = np.zeros((T - 1, M))
        phred_scores = np.zeros((T, M))
        prob_scores = np.zeros((T, M))
        for t in range(1, T):
            for j in range(M):
                probability = omega[t - 1] + np.log(A[:, j] + epsilon) + B[t, j]
                max_log_prob = np.max(probability)
                probs = np.exp(probability - max_log_prob)
                probs /= np.sum(probs)
                error_probs = 1 - probs
                Q = -10 * np.log10(error_probs + epsilon)
                Q_rounded = np.round(Q)
                Q_capped = np.clip(Q_rounded, 0, 60)
                phred_scores[t, j] = np.max(Q_capped)
                prob_scores[t, j] = np.max(probs)
                prev[t - 1, j] = np.argmax(probability)
                omega[t, j] = np.max(probability)
                if np.isinf(omega[t, j]):
                    omega[t, j] = -10000000

        X = np.zeros(T, dtype=int)
        X[T - 1] = np.argmax(omega[T - 1, :])
        for t in range(T - 2, -1, -1):
            X[t] = int(prev[t, X[t + 1]])
        result = [str(x) for x in X]
        return result, prob_scores

    def fit_dispersion(self, max_iter=10, tol=1e-3):
        # new_disp = update_global_dispersion(self._obs_dict, self._sample, self._chr, scale_factor=self._scale_factor)
        new_disp = estimate_global_dispersion_MLE(self._obs_dict, self._sample, self._chr, scale_factor=100, epsilon=1e-8)
        global_disp_array = np.ones(self._n_components) * new_disp
        self._emissions = self.compute_log_likelihood(state_dispersion=global_disp_array)
        self._state_dispersion = global_disp_array

        return new_disp


########################################################################
# Main
########################################################################

if __name__ == "__main__":
    dummy_obs = {"chr1": []}
    for i in range(100):
        region = {
            "sample1": {
                "normalized_depth": np.random.uniform(0.8, 1.2),
                "bg_mean": 1.0,
                "bg_std": 0.2,
                "coordinate": f"chr1:{1000 + i*100}-{1100 + i*100}_exon{i}",
                "is_offtarget": False
            }
        }
        dummy_obs["chr1"].append(region)
        
    sample_name = "sample1"
    print("=== Running CustomHMM (global dispersion approach) ===")
    custom_model = CustomHMM(obs_dict=dummy_obs, sample=sample_name, chr="chr1", n_components=5)
    custom_model.fit_dispersion(max_iter=10, tol=1e-3)
    alpha = custom_model.forward()
    print("Forward probabilities (log-scale):")
    print(alpha)
    posterior = custom_model.posterior_decoding()
    print("Posterior probabilities:")
    print(posterior)
    state_path, prob_scores = custom_model.decode()
    print("Viterbi state path:")
    print(state_path)
    if custom_model._state_dispersion is not None:
        print("Global dispersion used for all states:")
        print(custom_model._state_dispersion)
