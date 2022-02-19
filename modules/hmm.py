import os
import sys
from hmmlearn import hmm, base
import numpy as np
from scipy.stats import nbinom, poisson, norm
import string
import pandas as pd
from collections import defaultdict

def calculate_positional_mean_variance(sample_list, analysis_dict):
    '''
    '''
    sample_baselines = {}
    df_dict = pd.read_csv(analysis_dict['normalized_depth'], sep="\t").to_dict(orient="index")
    for sample in sample_list:
        baseline_samples = []
        num = 0
        sample_depth_tag = ("{}_normalized_final").format(sample.name)
        for control in sample.references:
            if num > 15:
                break
            control_depth_tag = ("{}_normalized_final").format(control[0])
            baseline_samples.append(control_depth_tag)
            num+=1
        sample_baselines[sample.name] = baseline_samples

    observations_dict = defaultdict(dict)
    # here calculate distribution params
    for region in df_dict:
        chromosome = df_dict[region]['chr']
        coord = ("{}:{}-{}_{}").format(df_dict[region]['chr'], df_dict[region]['start'],
            df_dict[region]['end'], df_dict[region]['exon'])
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

            background_mean = round(np.median(background_depth),3)
            background_std  = round(np.std(background_depth),3)

            region_dict[sample]['bg_mean'] = background_mean
            region_dict[sample]['bg_std'] = background_std
            region_dict[sample]['normalized_depth'] = sample_depth
            region_dict[sample]['coordinate'] = coord
        observations_dict[chromosome].append(region_dict)
    # sys.exit()
    return observations_dict

class CustomHMM():
    '''
    '''
    def __init__(self, obs_dict, sample, chr, n_components=5, transitions=None, emissions=None, start_prob=None):
        self._obs_dict    = obs_dict
        self._sample      = sample
        self._chr         = chr
        self._n_components= n_components
        self._transitions = transitions
        self._emissions   = emissions
        self._start_prob  = start_prob
        self._emissions   = self.compute_log_likelihood()

    @property
    def observations(self):
        '''
        '''
        obs = []
        for region in self._obs_dict[self._chr]:
            x = region[self._sample]['normalized_depth']
            obs.append(x)
        return obs

    def compute_log_likelihood(self):
        '''
        '''
        logp_dict = defaultdict(dict)
        emissions = []
        idx = 0
        for region in self._obs_dict[self._chr]:
            sample_depth = region[self._sample]['normalized_depth']
            bg_depth = region[self._sample]['bg_mean']
            bg_std   = region[self._sample]['bg_std']
            coord    = region[self._sample]['coordinate']
            logp_dict[idx]= defaultdict(dict)
            state_list = []
            depth_list = []
            probs_list = []
            for state in range(self._n_components):
                ratio = state/2
                x = sample_depth

                # x = sample_depth*ratio
                mean_state =bg_depth*ratio
                std_state = bg_std
                if x == 0:
                    x = 1
                logp_dict[idx][state] = round(np.log(norm.pdf(x, loc=mean_state, scale=bg_std)),3)
                probs_list.append(norm.pdf(x, loc=bg_depth, scale=bg_std))
                state_list.append(logp_dict[idx][state])
                depth_list.append(str(state) + ":" +str(x))
            emissions.append(state_list)
            idx+=1
        # sys.exit()
        return emissions

    def decode(self):
        '''
            O: Observations
            S: States (copy number: 0, 1, 2, 3, 4)
            X: Probability path
            A: State transition matrix
            B: Emission matrix
        '''

        start_p = np.array([.20, .20, .20, .20, .20])
        O = np.array(self.observations)
        S = np.array([0, 1, 2, 3, 4])
        A = np.array([[0.5, 0, 0.5, 0, 0],
                     [0, 0.5, 0.5, 0, 0],
                     [0.01, 0.01, 0.96, 0.01, 0.01],
                     [0, 0, 0.5, 0.5, 0],
                     [0, 0, 0.5, 0, 0.5]])
        B = np.array(self._emissions)
        T = O.shape[0]
        M = S.shape[0]

        # Initialization
        omega = np.zeros((T, M))
        omega[0,:] = np.log(start_p) + B[0,:]

        prev = np.zeros((T - 1, M))

        for t in range(1, T):
            for j in range(M):
                probability = omega[t - 1] + np.log(A[:, j]) + B[t, j]

                # This is our most probable state given previous state at time t (1)
                prev[t - 1, j] = np.argmax(probability)

                # This is the probability of the most probable state (2)
                omega[t, j] = np.max(probability)

                # Workaround when all state probs are inifite, this caused to immediately
                # assign 0 state
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

        return result

if __name__ == '__main__':

    model = CustomHMM()
    obs = np.random.normal(1, 0.2, 1000)
    # mean = np.mean(obs)
    model.compute_log_likelihood(obs)
