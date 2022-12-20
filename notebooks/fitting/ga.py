from typing import Sequence, Callable, Any, Optional

from math import floor
import numpy as np

import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader


class GeneticAlgorithm:
    """Object which handles genetic algorithm for feature selection.
    Creates many instances of a neural network and selects varying
    combinations of input features.

    Parameters
    ----------
    model : nn.Module
        the neural network model to use
    model_args : Sequence[Any]
        arguments to create the NN
    model_kwargs : dict
        keyword arguments
    fitness_fn : Callable
        Function with which to evaluate the fitness of the networks
    feature_list : Sequence[int]
    n_gen : optional, int
        Number of generations to run GA
    n_models : optional, int
        Number of individual models
    frac_keep : optional, float
        Fraction of best individual models to keep at each iteration
    mut_rate : optional, float
    """

    def __init__(
        self,
        model_constructor: Callable,
        model_args: Sequence[Any],
        model_kwargs: dict,
        training_fn: Callable,
        train_params: dict,
        data: DataLoader,
        optimizer,
        fitness_fn: Callable,
        feature_list: list[int],
        n_gen: Optional[int] = 1000,
        n_models: Optional[int] = 1000,
        frac_keep: Optional[float] = 0.1,
        mut_rate: Optional[float] = 0.05,
    ):
        self.model_constructor = model_constructor
        self.model_args = model_args
        self.n_input_feat = model_args[0]
        self.model_kwargs = model_kwargs

        self.training_fn = training_fn
        self.train_params = train_params
        self.data = data
        self.optimizer = optimizer
        self.individuals = []
        self.fitnesses = None

        self.fitness_fn = fitness_fn
        self.feature_list = feature_list
        self.n_gen = n_gen
        self.n_models = n_models
        self.frac_keep = frac_keep
        self.mut_rate = mut_rate

        self.history = []

        self._initialize()

    def _initialize(self):
        for i in range(self.n_models):
            # sample indices of features randomly from feature_list
            feature_set = np.random.choice(
                self.feature_list, self.n_input_feat, replace=False
            )
            self.individuals.append(
                Individual(
                    self.model_constructor,
                    self.model_args,
                    self.model_kwargs,
                    feature_set,
                    self.training_fn,
                    self.optimizer,
                    self.fitness_fn,
                )
            )

    def run(self):
        for i in range(self.n_gen):
            self._run_iteration()
            step_history = {}
            step_history["fitnesses"] = self.fitnesses
            step_history["model_features"]
            self.history.append(step_history)

    def _run_iteration(self):
        # train each model separately
        self.train()
        # evaluate individuals
        self.evaluate()
        # choose best `frac_keep` individuals
        n_keep = floor(self.n_gen * self.frac_keep)
        n_discard = self.n_models - n_keep
        top_individuals_ids = np.argsort(fitnesses)[:, n_keep]
        top_individuals = self.individuals[top_individuals_ids]
        # discard

        # generate new children by mutation
        self.mutate(top_individuals)

    def train(self, dataset):
        for ind in self.individuals:
            ind.train(dataset, **self.train_params)

    def evaluate(self, test_data):
        """Evaluate each of the individual models using the fitness metric"""
        if self.fitnesses is None:
            self.fitnesses = []
        for ind in self.individuals:
            f = ind.get_fitness()
            self.fitnesses.append(f)
        return


class Individual:
    """Class containing information for individuals in genetic algorithm

    Parameters
    ----------
    model: nn.Module,
    model_args: Sequence[Any],
    model_kwargs: dict,
    feature_set: list[int],
    training_fn: Callable[[nn.Module, nn.DataLoader, Callable, torch.optim.optimizer.Optimizer], None]
    tirness_fn: Callable
    """

    def __init__(
        self,
        model_constructor: Callable,
        model_args: Sequence[Any],
        model_kwargs: dict,
        feature_set: list[int],
        training_fn: Callable,
        optimizer,
        fitness_fn: Callable,
    ):
        self.fitness = None
        self.feature_set = feature_set
        self.model_args = model_args
        self.n_input_feat = self.model_args[0]  # input dimension to model
        self.model_kwargs = model_kwargs
        self.model = model_constructor(*self.model_args, **self.model_kwargs)

        self.training_fn = training_fn
        self.optimizer = optimizer(self.model.parameters())
        self.fitness_fn = fitness_fn
        self.is_trained = False

    def __repr__(self):
        print(
            f"""Model: {self.model}
Model Args: {self.model_args}
Model Kwargs: {self.model_kwargs}
Feature set: {self.feature_set}
Training function: {self.training_fn}"""
        )

    def get_fitness(self):
        if self.fitness is None:
            # return 0
            self.fitness = self.fitness_fn()
        else:
            return self.fitness

    def set_fitness(self, fitness):
        self.fitness = fitness

    # def update_fitness(self, train_params):
    # self.training_fn(train_params)

    def get_feature_set(self):
        return self.feature_set

    def set_feature_set(self, features):
        self.feature_set = features

    def train(self, dataset, **train_params):
        self.training_fn(
            self.model,
            dataset,
            self.fitness_fn,
            self.optimizer,
            self.feature_set,
            **train_params,
        )


class SingleLayerNet(nn.Module):
    """Neural network with a single hidden layer
    and sigmoid activation.
    """

    def __init__(self, input_dim, hidden_dim=30, output_dim=1):
        super(SingleLayerNet, self).__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, output_dim),
            nn.Sigmoid(),
        )

    def forward(self, x):
        return self.net(x)


class MultiLayerNet(nn.Module):
    """Neural network with a multiple hidden layers
    and sigmoid activation.

    Parameters
    ----------
    """

    def __init__(self, input_dim, hidden_dim=30, output_dim=1, n_hidden=5):
        super(MultiLayerNet, self).__init__()
        layers = []
        layers.append(nn.Linear(input_dim, hidden_dim))
        layers.append(nn.ReLU())
        for _ in range(n_hidden):
            layers.append(nn.Linear(hidden_dim, hidden_dim))
            layers.append(nn.ReLU())
        layers.append(nn.Linear(hidden_dim, output_dim))
        layers.append(nn.Sigmoid())
        self.net = nn.Sequential(*layers)

    def forward(self, x):
        return self.net(x)


class CommittorDataset(Dataset):
    def __init__(self, data, q):
        self.data = data
        self.q = q

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        return (self.data[idx], self.q[idx])


def train_loop(dataloader, model, loss_fn, optimizer, verbose=False):
    size = len(dataloader.dataset)
    for batch, (X, y) in enumerate(dataloader):
        # Compute prediction and loss
        pred = model(X)
        loss = loss_fn(pred, y)

        # Backpropagation
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        if verbose:
            if batch % 100 == 0:
                loss, current = loss.item(), batch * len(X)
                print(f"loss: {loss:>7f}  [{current:>5d}/{size:>5d}]")
    return loss.item()


def test_loop(dataloader, model, loss_fn, verbose=False):
    test_loss = 0
    n_batches = len(dataloader)

    with torch.no_grad():
        for (X, y) in iter(dataloader):
            pred = model(X)
            test_loss += loss_fn(pred, y).item()

    test_loss /= n_batches
    if verbose:
        print(
            f"Test Error: \n RMSE: {np.sqrt(test_loss):>7f}%, Avg loss: {test_loss:>8f} \n"
        )
    return test_loss
