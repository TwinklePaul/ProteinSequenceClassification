import numpy as np
from . import aminos


def encodeDistances(sequence):
    encoded_distances = {}

    for amino in aminos.amino_acids:
        encoded_distances[f'{amino}'] = []

    for index, amino in enumerate(sequence):
        # dict[amino] =
        #               0th index of list = prev occurance of the amino
        #               from 1st index -> store index - dict[amino][0]
        if len(encoded_distances[f'{amino}']) == 0:
            encoded_distances[f'{amino}'].append(index+1)
        else:
            prev_distance = encoded_distances[f'{amino}'][0] - 1
            encoded_distances[f'{amino}'].append(index - prev_distance)

    for amino in aminos.amino_acids:
        if len(encoded_distances[f'{amino}']) > 1:
            encoded_distances[f'{amino}'].pop(0)

    return encoded_distances


def applyMeanNormalization(encoded_list, mean):
    min_distance = min(encoded_list)
    max_distance = max(encoded_list)

    if (len(encoded_list) == 1):
        range_distance = 1
    else:
        range_distance = max_distance - min_distance

    for index, distance in enumerate(encoded_list):
        encoded_list[index] = (distance - min_distance)/range_distance

    return encoded_list


def performDistanceBasedEncoding(sequence):

    n_mean = {}
    n_SD = {}
    n_CV = {}

    distances = encodeDistances(sequence)

    for amino in aminos.amino_acids:
        # find mean
        mean = np.mean(distances[f'{amino}'])

        distances[f'{amino}'] = applyMeanNormalization(
            distances[f'{amino}'], mean)

        # find the normalized values
        n_mean[f'{amino}'] = np.mean(distances[f'{amino}'])
        n_SD[f'{amino}'] = np.std(distances[f'{amino}'])

        if n_mean[f'{amino}'] != 0:
            n_CV[f'{amino}'] = n_SD[f'{amino}'] / n_mean[f'{amino}']
        else:
            n_CV[f'{amino}'] = 0

    return (n_mean, n_SD, n_CV)
