import numpy as np
from biomass.estimation import convert_scale, initialize_search_param

from .name2idx import C, V
from .set_model import initial_values, param_values


class SearchParam(object):
    """Specify model parameters and/or initial values to optimize"""

    def __init__(self):
        # parameters
        self.idx_params = [
            C.kf1,
            C.kf2,
            C.kf3,
            C.kf4,
            C.kf5,
            C.kf6,
            C.kf7,
            C.kf8,
            C.kf9,
            C.kr9,
            # C.kf10,
            # C.kr10,
            C.kf11,
            C.kf12,
            C.kf13,
            C.kf14,
            C.kf15,
            C.kf16,
            C.kf17,
            C.kf18,
            C.kf19,
            C.V20,
            C.K20,
            # C.n20,
            # C.kf21,
            # C.kf22,
            # C.kf23,
            # C.kf24,
            C.V25,
            C.K25,
            # C.n25,
            # C.kf26,
            # C.kf27,
            # C.kf28,
        ]

        # initial values
        self.idx_initials = [
            V.TTR,
            V.IKK,
            V.IkBa_NFkB,
        ]

    def get_region(self):
        x = param_values()
        y0 = initial_values()

        search_param = initialize_search_param(
            parameters=C.NAMES,
            species=V.NAMES,
            param_values=x,
            initial_values=y0,
            estimated_params=self.idx_params,
            estimated_initials=self.idx_initials,
        )

        search_rgn = np.zeros((2, len(x) + len(y0)))
        # Default: 0.1 ~ 10
        for i, j in enumerate(self.idx_params):
            search_rgn[0, j] = search_param[i] * 0.1  # lower bound
            search_rgn[1, j] = search_param[i] * 10   # upper bound
        # Default: 0.5 ~ 2
        for i, j in enumerate(self.idx_initials):
            search_rgn[0, j + len(x)] = search_param[i + len(self.idx_params)] * 1e-2  # lower bound
            search_rgn[1, j + len(x)] = search_param[i + len(self.idx_params)] * 1e2  # upper bound

        search_rgn[:,C.kf1] = [1e-10, 1e-7]
        search_rgn[:,C.kf2] = [1e-10, 1e-7]
        search_rgn[:,C.kf5] = [1e-4, 1e1]
        search_rgn[:,C.kf7] = [1e-4, 1e1]
        search_rgn[:,C.kf11] = [1e-7, 1e-5]

        search_rgn = convert_scale(
            region=search_rgn,
            parameters=C.NAMES,
            species=V.NAMES,
            estimated_params=self.idx_params,
            estimated_initials=self.idx_initials,
        )

        return search_rgn

    def update(self, indiv):
        x = param_values()
        y0 = initial_values()

        for i, j in enumerate(self.idx_params):
            x[j] = indiv[i]
        for i, j in enumerate(self.idx_initials):
            y0[j] = indiv[i + len(self.idx_params)]

        # parameter constraints
        x[C.kf10] = x[C.kf9]
        x[C.kr10] = x[C.kr9]
        x[C.kf24] = x[C.kf23]

        return x, y0

    def gene2val(self, indiv_gene):
        search_rgn = self.get_region()
        indiv = 10 ** (indiv_gene * (search_rgn[1, :] - search_rgn[0, :]) + search_rgn[0, :])

        return indiv

    def val2gene(self, indiv):
        search_rgn = self.get_region()
        indiv_gene = (np.log10(indiv) - search_rgn[0, :]) / (search_rgn[1, :] - search_rgn[0, :])

        return indiv_gene
