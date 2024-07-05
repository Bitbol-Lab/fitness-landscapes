#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from itertools import product
from decimal import *

def transition_probability_star(f_wt, f_m, D, carrying_cap, alpha, g_m, g_wt, L, walk):
    s = f_m / f_wt - 1
    s2 = f_wt / f_m - 1
    N_w = carrying_cap * (1 - g_wt / f_wt)
    N_m = carrying_cap * (1 - g_m / f_m)


    if walk == 'w-f':
        p_m = (1 - np.exp(- 2 * s)) / (1 - np.exp(- 2 * s * N_w))
        p_w = (1 - np.exp(- 2 * s2)) / (1 - np.exp(- 2 * s2 * N_m))
    elif walk == 'moran':
        p_m = (1 - (s + 1) ** (-1)) / (1 - ( s + 1) ** (-N_w))
        p_w = (1 - (s2 + 1) ** (-1)) / (1 - (s2 + 1) ** (-N_m))

    gamma = N_w * p_w / (N_m * p_m)
    numerator = p_m * (1 - gamma ** 2) * (gamma + alpha * D + gamma * alpha ** 2 * (D - 1))
    denominator = D * (alpha + gamma) * (1 + alpha * gamma - gamma * (alpha + gamma) * (gamma * (1 + alpha * gamma) / (alpha + gamma)) ** (D - 1))
    return numerator / denominator / L


def star_walk(landscape, genotype_space, D, alpha, C, g_m, g_wt, L, rep, walk):
    """
    Parameters
    ----------
    landscape: dictionary
        Keys are tuples denoting genotype sequences, e.g. (0, 0, 0), and the corresponding values are fitnesses

    genotype_space: list of tuples
        Each element is a tuple denoting a genotype sequence, e.g. (0, 0, 0)

    D: integer
        D is the number of demes

    alpha: float
        alpha is the migration asymmetry

    C: integer
        C is the carrying capacity

    g_m: float
        g_m denotes mutant death rate

    g_wt: float
        g_wt denotes wild-type death rate

    rep: integer
        rep is the number of walks performed per starting node

    walk: chain of character
        walk = 'moran' or 'w-f', it determines if the code simulates the Moran or Wright-Fisher star walk

    Returns
    -------
    mean_l: float
        average length of the walks (number of successful fixations) before reaching a fitness peak

    mean_h: float
        mean fitness of the first fitness peak reached

    SD_h: float
         Standard deviation of the fitness of the first fitness peak reached
    """

    l_list = np.zeros(rep*2**L) # lengths of the walks
    h_list = np.zeros(rep*2**L) # fitness values of the first encountered peak


    for initial_node in range(2**L):
        for r in range(rep): # We perform the same number of walks for each starting node
            wt = initial_node # index of the current wild-type sequence
            f_wt = landscape[genotype_space[wt]] # fitness of the wild type

            trajectory = [wt]

            while True:
                neighbor_fitnesses = []  # list of the fitnesses of the neighbors
                neighbors = []  # list of the indices of the neighbors
                transition_proba = [] # transition probabilities towards the neighbors
                counter = 0  # counter of fitter neighbors to detect fitness peaks

                for i in range(L):  # We check all possible mutations
                    mutant = list(genotype_space[wt])
                    mutant[i] = (-1) ** genotype_space[wt][i] + genotype_space[wt][i]
                    mutant = tuple(mutant)
                    neighbors.append(genotype_space.index(mutant))
                    f_m = landscape[mutant]
                    neighbor_fitnesses.append(f_m)

                    transition_proba.append(transition_probability_star(f_wt, f_m, D, C, alpha, g_m, g_wt, L, walk))

                    if f_m > f_wt: # We check if there is still at least one beneficial mutation or not
                        counter += 1

                # If we have reached a peak, we stop the walk:
                if counter == 0:
                    break

                # Here, we normalize the transition probabilities towards the neighbors.
                # To compute the number of mutation events (fixing or not), denoted by t in the paper,
                # instead of the number of fixations l, the line below should be commented out and replaced by:
                # neighbors.append(wt)
                # transition_proba.append(1-np.sum(transition_proba))
                # neighbor_fitnesses.append(f_wt)
                transition_proba = np.array(transition_proba) / np.sum(transition_proba)

                cumulative_proba = [np.sum(transition_proba[:i]) for i in range(1, L + 1)]
                draw = np.searchsorted(cumulative_proba, np.random.uniform(0, 1))

                trajectory.append(neighbors[draw]) # The mutation is done
                wt = neighbors[draw]
                f_wt = neighbor_fitnesses[draw]

            l_list[rep * initial_node + r] = len(trajectory) - 1
            h_list[rep * initial_node + r] = f_wt


    mean_l = np.mean(l_list)
    mean_h = np.mean(h_list)
    SD_h = np.std(h_list)

    return mean_l, mean_h, SD_h


if __name__ == "__main__":
    # Working example:

    # Arguments
    L = 3
    landscape = {(0, 0, 0): 1.554, (0, 0, 1): 1.674, (0, 1, 0): 2.051, (0, 1, 1): 1.181,
                 (1, 0, 0): 2.332, (1, 0, 1): 2.452, (1, 1, 0): 2.004, (1, 1, 1): 1.134}
    genotype_space = list(product(range(2), repeat=L))
    D = 10
    alpha = 10**(-1)
    g_wt = 0.01
    g_m = 0.01
    C = 20

    rep = 10000
    moran_or_wf = 'moran' # 'moran' or 'wf' version of the star walk

    # Simulation
    moran_star_simulation = star_walk(landscape, genotype_space, D, alpha, C, g_m, g_wt, L, rep, 'moran')
    wf_star_simulation = star_walk(landscape, genotype_space, D, alpha, C, g_m, g_wt, L, rep, 'w-f')

    # Outcome
    print('Moran star walk simulation:')
    print('Mean length = ' + str(moran_star_simulation[0]))
    print('Mean height = ' + str(moran_star_simulation[1]))
    print('SD height = ' + str(moran_star_simulation[2]))

    print('')

    print('Wright-Fisher star walk simulation:')
    print('Mean length = ' + str(wf_star_simulation[0]))
    print('Mean height = ' + str(wf_star_simulation[1]))
    print('SD height = ' + str(wf_star_simulation[2]))
