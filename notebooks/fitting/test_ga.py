import numpy as np

import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader

from ga import GeneticAlgorithm, Individual

import pytest


class PlaceHolderNet(nn.Module):
    def __init__(self):
        super(self).__init__()

    def forward(self, x):
        return x


class TestGeneticAlgorithm:
    def test_initialize(self):
        pass

    def test_train(self):
        pass

    def test_evaluate(self):
        pass

    def test_mutate(self):
        pass

    def test_discard(self):
        pass


class TestIndividual:
    def __init__(self):
        self.Individual = Individual()

    def test_initialize(self):
        pass

    def test_get_fitness(self):
        assert self.Individual.get_fitness() == 0

    def test_set_fitness(self):
        self.Individual.set_fitness(1)
        assert self.Individual.fitness == 1

    def test_get_feature_set(self):
        pass

    def test_set_feature_set(self):
        pass

    def test_train(self):
        # make fake dataset
        pass
