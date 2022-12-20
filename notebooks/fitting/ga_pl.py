import copy
from typing import Sequence, Callable, Any, Optional

from math import floor
import numpy as np
from sklearn.model_selection import train_test_split

import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
import pytorch_lightning as pl


class GeneticAlgorithmPL:
    """Object which handles genetic algorithm for feature selection.
    Creates many instances of a neural network and selects varying
    combinations of input features.

    Parameters
    ----------
    model_constructor : Callable[Any] -> pl.LightningModule
        a method which constructs a new PyTorch Lightning module
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
        model_constructor: Callable[Any, pl.LightningModule],
        model_args: Sequence[Any],
        model_kwargs: dict,
        data,
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

        self.X, self.y = data
        # self._make_train_test_data()

        self.individuals = []
        self.fitnesses = None

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
                )
            )

    def _make_train_test_data(self):
        # train, test split
        train_X, train_y, test_X, test_y = train_test_split(self.X, self.y, test_size=0.2)
        self.train_dataset = CommittorDataset(train_X, train_y)
        self.test_dataset = CommittorDataset(test_X, test_y)

    def run(self, verbose=True):
        for i in range(self.n_gen):
            self._run_iteration()
            step_history = {}
            step_history["fitnesses"] = self.fitnesses
            step_history["model_features"] = [ind.feature_list for ind in self.individuals]
            self.history.append(step_history)

    def _run_iteration(self):
        # train each model separately
        self.train()
        # evaluate individuals
        self.evaluate()
        # choose best `frac_keep` individuals
        top_individuals_ids = self.select()
        top_individuals = self.individuals[top_individuals_ids]
        # discard
        for i in range(self.n_models):
            if i not in top_individuals_ids:
                del self.individuals[i]
        # generate new children by mutation
        self.mutate(top_individuals, n_discard)

    def train(self, **training_kwargs):
        # TODO: make parallel dispatch using joblib
        dataloader = DataLoader(self.train_dataset, batch_size=self.batch_size)
        for ind in self.individuals:
            ind.train(dataloader, **training_kwargs)

    def evaluate(self, test_data):
        """Evaluate each of the individual models using the fitness metric"""
        if self.fitnesses is None:
            self.fitnesses = []
        dataloader = DataLoader(self.test_dataset, batch_size=self.batch_size)
        for ind in self.individuals:
            # perform on test
            f = ind.get_fitness(dataloader)
            self.fitnesses.append(f)
 
    def selection(self):
        n_keep = floor(self.n_gen * self.frac_keep)
        n_discard = self.n_models - n_keep
        top_individuals_ids = np.argsort(self.fitnesses)[:, n_keep]
        return top_individuals_ids

    def mutate(self, top_individuals_ids, n_discard):
        """Creates n_discard new individuals by mutating
        the top individuals."""
        for i in range(n_discard):
            # choose good top individual to mutate/reproduce
            old_individual = self.individuals[np.random.choice(top_individuals_ids)]
            test = np.random.random()
            if test < self.mut_rate:
                # mutate a single descriptor
                new_feature_set = copy.copy(old_individual.feature_set)
                feature_to_change = np.random.choice(old_individual.feature_set)
                new_feature = np.random.choice(self.feature_set, 1)
                new_feature_set[feature_to_change] = new_feature
                self.individuals.append(
                    Individual(
                        self.model_constructor,
                        self.model_args,
                        self.model_kwargs,
                        new_feature_set,
                    )
                )
            else:
                # else, duplicate old individual
                new_individual = copy.deep_copy(old_individual)
                self.individuals.append(new_individual)

    def cross_over(self, first_individual, second_individual):
        random_point = np.random.randint(0, self.n_input_feat)
        for i in range(self.n_input_feat):
            if i > random_point:
                first_individual.features_set[i], second_individual.feature_set[i]
                    = first_individual.feature_set[i], second_individual.feature_set[i]
        # return first_candidate, second_candidate

    def mutate(self, to_mutate):
        for i in range(len(to_mutate):
            if np.random.random() < self.mut_rate:
                new_feature_set = copy.copy(old_individual.feature_set)
                feature_to_change = np.random.choice(old_individual.feature_set)
                new_feature = np.random.choice(self.feature_set, 1)
                new_feature_set[feature_to_change] = new_feature
                self.individuals.append(
                    Individual(
                        self.model_constructor,
                        self.model_args,
                        self.model_kwargs,
                        new_feature_set,
                    )
                )
        return offspring


class Individual:
    """Class containing information for individuals in genetic algorithm

    Parameters
    ----------
    model_constructor : (Any) -> pl.LightningModule
        a method which constructs a new PyTorch Lightning module
    model_args: Sequence[Any],
    model_kwargs: dict,
    feature_set: list[int],
    """

    def __init__(
        self,
        model_constructor: Callable,
        model_args: Sequence[Any],
        model_kwargs: dict,
        feature_set: list[int],
    ):
        self.fitness = None
        self.feature_set = feature_set
        self.model_args = model_args
        self.n_input_feat = self.model_args[0]  # input dimension to model
        self.model_kwargs = model_kwargs
        self.model = model_constructor(*self.model_args, **self.model_kwargs)

        self.is_trained = False
        self.trainer = None

    def __repr__(self):
        print(
            f"""Model: {self.model}
Model Args: {self.model_args}
Model Kwargs: {self.model_kwargs}
Feature set: {self.feature_set}"""
        )

    def get_fitness(self, test_dataloader):
        if self.fitness is None:
            # return 0
            self.fitness = self.trainer.test(self.model, test_dataloader)
        return self.fitness

    def set_fitness(self, fitness):
        self.fitness = fitness

    # def update_fitness(self, train_params):
    # self.training_fn(train_params)

    def get_feature_set(self):
        return self.feature_set

    def set_feature_set(self, features):
        self.feature_set = features

    def train(self, train_dataset, val_dataset, **training_kwargs):
        if self.trainer is None:
            self.trainer = pl.Trainer(**training_kwargs)
        # train, val, test split
        self.trainer.fit(model=self.model, train_dataloaders=train_dataset, val_dataloaders=val_dataset)
        self.is_trained = True


class SingleLayerNet(pl.LightningModule):
    """Neural network with a single hidden layer
    and sigmoid activation.
    """

    def __init__(self, input_dim, hidden_dim=30, output_dim=1, verbose=False):
        super(SingleLayerNet, self).__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, output_dim),
            nn.Sigmoid(),
        )
        self.verbose = verbose

    def forward(self, x):
        return self.net(x)

    def training_step(self, batch, batch_idx):
        X, y = batch
        # Compute prediction and loss
        pred = self.net(X)
        loss = nn.MSELoss(pred, y)

        if self.verbose:
            if batch_idx % 100 == 0:
                loss, current = loss.item(), batch * len(X)
                print(f"loss: {loss:>7f}  [{current:>5d}/]")
        return loss

    def configure_optimizers(self):
        optimizer = torch.optim.Adam(self.parameters(), lr=4e-3)
        return optimizer


class MultiLayerNet(pl.LightningModule):
    """Neural network with a multiple hidden layers
    and sigmoid activation.

    Parameters
    ----------
    """

    def __init__(
        self, input_dim, hidden_dim=30, output_dim=1, n_hidden=5, verbose=True
    ):
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
        self.verbose = verbose

    def forward(self, x):
        return self.net(x)

    def training_step(self, batch, batch_idx):
        X, y = batch
        # Compute prediction and loss
        pred = self.net(X)
        loss = nn.MSELoss(pred, y)

        if self.verbose:
            if batch_idx % 100 == 0:
                self.log("train_loss", loss, on_step=True, on_epoch=True, prog_bar=True, logger=True)
        return loss

    def validation_step(self, batch, batch_idx):
        x, y = batch
        pred = self.net(x)
        loss = F.mse_loss(pred, y)
        if self.verbose:
            if batch_idx % 100 == 0:
                self.log("val_loss", loss)

    def test_step(self, batch, batch_idx):
        x, y = batch
        pred = self.net(x)
        loss = F.mse_loss(pred, y)

    def configure_optimizers(self):
        optimizer = torch.optim.Adam(self.parameters(), lr=4e-3)
        return optimizer


class CommittorDataset(Dataset):
    def __init__(self, data, q):
        self.data = data
        self.q = q

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        return (self.data[idx], self.q[idx])
