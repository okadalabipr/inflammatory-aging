from .name2idx import C, V

class DifferentialEquation(object):
    def __init__(self, perturbation):
        super(DifferentialEquation, self).__init__()
        self.perturbation = perturbation

    def diffeq(self, t, y, *x):
        """Kinetic equations"""
        # v : flux vector
        v = {}
        v[1] = x[C.kf1] * y[V.TNF]
        v[2] = x[C.kf2] * y[V.TNFR]
        v[3] = x[C.kf3] * y[V.TNFR] * y[V.TTR]
        v[4] = x[C.kf4] * y[V.pIKK] * y[V.TTR]
        v[5] = x[C.kf5] * y[V.A20] * y[V.aTTR]
        v[6] = x[C.kf6] * y[V.aTTR] * y[V.IKK]
        v[7] = x[C.kf7] * y[V.A20] * y[V.pIKK]
        v[8] = x[C.kf8] * y[V.iIKK]
        v[9] = x[C.kf9] * y[V.IkBa] * y[V.NFkB] - x[C.kr9] * y[V.IkBa_NFkB]
        v[10] = x[C.kf10] * y[V.IkBan] * y[V.NFkBn] - x[C.kr10] * y[V.IkBa_NFkBn]
        v[11] = x[C.kf11] * y[V.pIKK] * y[V.IkBa_NFkB]
        v[12] = x[C.kf12] * y[V.pIkBa_NFkB]
        v[13] = x[C.kf13] * y[V.IkBa_NFkBn]
        v[14] = x[C.kf14] * y[V.NFkB]
        v[15] = x[C.kf15] * y[V.NFkBn]
        v[16] = x[C.kf16] * y[V.pIKK] * y[V.IkBa]
        v[17] = x[C.kf17] * y[V.pIkBa]
        v[18] = x[C.kf18] * y[V.IkBa]
        v[19] = x[C.kf19] * y[V.IkBan]
        v[20] = x[C.V20] * y[V.NFkBn] ** x[C.n20] / (x[C.K20] ** x[C.n20] + y[V.NFkBn] ** x[C.n20])
        v[21] = x[C.kf21] * y[V.ikbamRNA]
        v[22] = x[C.kf22] * y[V.ikbamRNA]
        v[23] = x[C.kf23] * y[V.IkBa]
        v[24] = x[C.kf24] * y[V.IkBan]
        v[25] = x[C.V25] * y[V.NFkBn] ** x[C.n25] / (x[C.K25] ** x[C.n25] + y[V.NFkBn] ** x[C.n25])
        v[26] = x[C.kf26] * y[V.a20mRNA]
        v[27] = x[C.kf27] * y[V.a20mRNA]
        v[28] = x[C.kf28] * y[V.A20]

        if self.perturbation:
            for i, dv in self.perturbation.items():
                v[i] = v[i] * dv

        dydt = [0] * V.NUM
        dydt[V.TNFR] = + v[1] - v[2]
        dydt[V.TTR] = - v[3] - v[4] + v[5]
        dydt[V.aTTR] = + v[3] + v[4] - v[5]
        dydt[V.IKK] = - v[6] + v[8]
        dydt[V.pIKK] = + v[6] - v[7]
        dydt[V.iIKK] = + v[7] - v[8]
        dydt[V.NFkB] = - v[9] + v[12] - v[14] + v[15]
        dydt[V.NFkBn] = - v[10] + v[14] - v[15]
        dydt[V.IkBa] = - v[9] - v[16] - v[18] + v[19] + v[21] - v[23]
        dydt[V.IkBan] = - v[10] + v[18] - v[19] - v[24]
        dydt[V.IkBa_NFkB] = + v[9] - v[11] + v[13]
        dydt[V.IkBa_NFkBn] = + v[10] - v[13]
        dydt[V.pIkBa] = + v[16] - v[17]
        dydt[V.pIkBa_NFkB] = + v[11] - v[12]
        dydt[V.ikbamRNA] = + v[20] - v[22]
        dydt[V.a20mRNA] = + v[25] - v[27]
        dydt[V.A20] = + v[26] - v[28]

        return dydt


def param_values():
    """Parameter values"""
    x = [1] * C.NUM
    x[C.kf6] = 0.004
    x[C.kf7] = 0.003
    x[C.kf8] = 0.0006
    x[C.kf9] = 0.5
    x[C.kr9] = 0.0005
    x[C.kf10] = x[C.kf9]
    x[C.kr10] = x[C.kr9]
    x[C.kf11] = 0.185
    x[C.kf12] = 0.1
    x[C.kf13] = 0.01
    x[C.kf14] = 0.0026
    x[C.kf15] = 0.00052
    x[C.kf16] = 0.074
    x[C.kf17] = 0.1
    x[C.kf18] = 0.00067
    x[C.kf19] = 0.000335
    x[C.n20] = 2
    x[C.V20] = 1.4E-7
    x[C.K20] = 0.065
    x[C.kf21] = 0.5
    x[C.kf22] = 3E-4
    x[C.kf23] = 5E-4
    x[C.kf24] = x[C.kf23]
    x[C.n25] = 2
    x[C.V25] = 1.4E-7
    x[C.K25] = 0.065
    x[C.kf26] = 0.5
    x[C.kf27] = 4.8E-4
    x[C.kf28] = 4.5E-3

    return x


def initial_values():
    """Values of the initial condition"""
    y0 = [0] * V.NUM
    y0[V.TNF] = 1
    y0[V.TTR] = 1
    y0[V.IKK] = 1
    y0[V.IkBa_NFkB] = 1

    return y0
