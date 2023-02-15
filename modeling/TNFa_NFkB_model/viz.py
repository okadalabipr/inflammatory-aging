from biomass.plotting import *
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

from .observable import Observable

class Visualization(Observable):
    """
    Plotting parameters for customizing figure properties.

    Attributes
    ----------
    cm : matplotlib.colors.ListedColormap (default: ``plt.cm.get_cmap('tab10')``)
        Choosing colormaps for ``cmap``.
    single_observable_options : list of SingleObservable
        Visualization options for time-course simulation (single-observable).
    multiple_observables_options : MultipleObservables
        Visualization options for time-course simulation (multi-observables).
    sensitivity_options : SensitivityOptions
        Visualization options for sensitivity analysis results.
    """

    def __init__(self):
        super().__init__()

        self.cm = plt.cm.get_cmap("tab10")
        self.single_observable_options = [
            SingleObservable(self.cm, obs_name) for obs_name in self.obs_names
        ]
        self.multiple_observables_options = MultipleObservables(self.cm)
        self.sensitivity_options = SensitivityOptions(self.cm)

    def get_single_observable_options(self):

        for i, _ in enumerate(self.obs_names):
            self.single_observable_options[i].figsize = (4, 3)
            self.single_observable_options[i].divided_by = 60
            self.single_observable_options[i].xlabel = 'TNFα treatment (min)'
            self.single_observable_options[i].ylabel = 'Intensity (a.u.)'
            self.single_observable_options[i].yticks = [0.5*i for i in range(3)]
            self.single_observable_options[i].dont_show = []

            time = 180
            if time == 180:
                self.single_observable_options[i].xlim = (-0.5, 180.5)
                self.single_observable_options[i].xticks = [60*i for i in range(3+1)]

            if time == 480:
                self.single_observable_options[i].xlim = (-0.5, 480.5)
                self.single_observable_options[i].xticks = [120*i for i in range(4+1)]

        return self.single_observable_options

    def get_multiple_observables_options(self) -> MultipleObservables:

        return self.multiple_observables_options

    def get_sensitivity_options(self) -> SensitivityOptions:

        return self.sensitivity_options

    @staticmethod
    def set_timecourse_rcParams():
        plt.rcParams["font.size"] = 24
        plt.rcParams['font.family'] = 'Arial'
        plt.rcParams['mathtext.fontset'] = 'custom'
        plt.rcParams['mathtext.it'] = 'Arial:italic'
        plt.rcParams["axes.linewidth"] = 1.5
        plt.rcParams["xtick.major.width"] = 1.5
        plt.rcParams["ytick.major.width"] = 1.5
        plt.rcParams["lines.linewidth"] = 3.0
        plt.rcParams["lines.markersize"] = 12
        plt.rcParams["savefig.dpi"] = 300
        plt.rcParams["savefig.bbox"] = "tight"

    @staticmethod
    def set_sensitivity_rcParams():
        """figure/sensitivity"""
        plt.rcParams["font.size"] = 18
        plt.rcParams['font.family'] = 'Arial'
        plt.rcParams['mathtext.fontset'] = 'custom'
        plt.rcParams['mathtext.it'] = 'Arial:italic'
        plt.rcParams["axes.linewidth"] = 1.2
        plt.rcParams["xtick.major.width"] = 1.2
        plt.rcParams["ytick.major.width"] = 1.2
        plt.rcParams["savefig.dpi"] = 300
        plt.rcParams["savefig.bbox"] = "tight"

    @staticmethod
    def convert_species_name(name):

        return name
