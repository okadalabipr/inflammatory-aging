from multiprocessing import Condition
import numpy as np

from biomass.dynamics.solver import *

from .name2idx import C, V
from .set_model import DifferentialEquation

class Observable(DifferentialEquation):
    """
    Correlating model simulations and experimental measurements.

    Attributes
    ----------
    obs_names : list of strings
        Names of model observables.

    t : range
        Simulation time span.

    conditions : list of strings
        Experimental conditions.

    simulations : numpy.ndarray
        The numpy array to store simulation results.

    normalization : nested dict
        * 'timepoint' : Optional[int]
            The time point at which simulated values are normalized.
            If :obj:`None`, the maximum value will be used for normalization.

        * 'condition' : list of strings
            The experimental conditions to use for normalization.
            If empty, all conditions defined in ``sim.conditions`` will be used.

    experiments : list of dict
        Time series data.

    error_bars : list of dict
        Error bars to show in figures.

    """

    def __init__(self):
        super(Observable, self).__init__(perturbation={})
        self.obs_names: list = [
            'TNFR',
            'active_TTR',
            'active_IKK',
            'inactive_IKK',
            'total_IkBa',
            'IkBa_mRNA',
            'A20_mRNA',
            'nuclear_NFkB',
            'total_nuclear_NFkB',
        ]
        self.t: range = range(0, 180*60+1)
        self.conditions: list = ['control']
        self.simulations: np.ndarray = np.empty(
            (len(self.obs_names), len(self.t), len(self.conditions))
        )
        self.normalization: dict = {}
        for observable in self.obs_names:
            self.normalization[observable] = {"timepoint": None, "condition": []}
        self.experiments: list = [None] * len(self.obs_names)
        self.error_bars: list = [None] * len(self.obs_names)

    def simulate(self, x: list, y0: list, _perturbation: dict = {}):
        if _perturbation:
            self.perturbation = _perturbation

        # unperturbed steady state
        for i, condition in enumerate(self.conditions):
            if condition == 'control':
                pass

            sol = solve_ode(self.diffeq, y0, self.t, tuple(x))
            if sol is None:
                return False
            else:
                self.simulations[self.obs_names.index('TNFR'), :, i] = (
                    sol.y[V.TNFR, :]
                )
                self.simulations[self.obs_names.index('active_TTR'), :, i] = (
                    sol.y[V.aTTR, :]
                )
                self.simulations[self.obs_names.index('active_IKK'), :, i] = (
                    sol.y[V.pIKK, :]
                )
                self.simulations[self.obs_names.index('inactive_IKK'), :, i] = (
                    sol.y[V.iIKK, :]
                )
                self.simulations[self.obs_names.index('total_IkBa'), :, i] = (
                    sol.y[V.IkBa, :]+sol.y[V.IkBa_NFkB, :]+
                    sol.y[V.pIkBa, :]+sol.y[V.pIkBa_NFkB, :]+
                    sol.y[V.IkBan, :]+sol.y[V.IkBa_NFkBn, :]
                )
                self.simulations[self.obs_names.index('IkBa_mRNA'), :, i] = (
                    sol.y[V.ikbamRNA, :]
                )
                self.simulations[self.obs_names.index('A20_mRNA'), :, i] = (
                    sol.y[V.a20mRNA, :]
                )
                self.simulations[self.obs_names.index('nuclear_NFkB'), :, i] = (
                    sol.y[V.NFkBn, :]
                )
                self.simulations[self.obs_names.index('total_nuclear_NFkB'), :, i] = (
                    sol.y[V.NFkBn, :]+sol.y[V.IkBa_NFkBn, :]
                )


    def set_data(self):
        # Experimental data for model learning
        self.experiments[self.obs_names.index('active_IKK')] = {
            'control' : [
                0.0,
                1.0,
                0.498284054,
                0.42485667,
                0.297049043,
                0.210507463,
                0.17685879,
                0.27337504,
                0.369912374,
                0.335024647,
                0.307125737,
                0.312581656,
                0.349592363
            ]
        }
        self.experiments[self.obs_names.index('total_IkBa')] = {
            'control' : [
                0.64483497,
                0.573655529,
                0.470261714,
                0.627503511,
                0.88683425,
                1.0,
                0.982390062,
                0.864739294,
                0.862507473,
                0.892113279,
                0.928203721,
                0.807215831,
                0.733951887
            ]
        }
        self.experiments[self.obs_names.index('total_nuclear_NFkB')] = {
            'control': [
                0.0,
                0.456246824,
                1.0,
                0.811670778,
                0.318807521,
                0.180022166,
                0.212843169,
                0.447420167,
                0.582405513,
                0.596044696,
                0.371012416,
                0.395512392,
                0.58077599
            ]
        }

        # Experimental data for model prediction
        self.experiments[self.obs_names.index('IkBa_mRNA')] = {
            'control' : [
                0.0,
                0.003913513,
                0.106755668,
                0.549977707,
                0.976703035,
                1.0,
                0.668329674,
                0.324153072,
                0.268227863,
                0.556378073,
                0.616180895,
                0.64597117,
                0.351993405,
            ]
        }
        self.experiments[self.obs_names.index('A20_mRNA')] = {
            'control' : [
                0.00476252,
                0.0,
                0.017948479,
                0.284928491,
                0.771260819,
                1.0,
                0.781345851,
                0.430022553,
                0.240113191,
                0.344950397,
                0.443265783,
                0.53033821,
                0.334726872,
            ]
        }

    @staticmethod
    def get_timepoint(obs_name: str):

        return [15*t*60 for t in range(13)]
