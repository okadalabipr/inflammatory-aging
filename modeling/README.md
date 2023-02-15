# Mathematical modeling of NF-κB dynamics on inflammatory aging

We used [`biomass==0.5.5`](https://github.com/biomass-dev/biomass) for mathematical modeling.
The package can be installed via pip:

```
$ pip install biomass==0.5.5
```

This requires Python 3.7 or later.

## Description

| Name                                                           | Content                                                                                                  |
| -------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------- |
| [`name2idx/`](./TNFa_NFkB_model/name2idx/)                     | Names of model parameters and species                                                                    |
| [`out/`](./TNFa_NFkB_model/out/)                               | Parameter values that are estimated from experimental data                                               |
| [`reaction_network.py`](./TNFa_NFkB_model/reaction_network.py) | Reaction indices grouped according to biological processes                                               |
| [`set_model.py`](./TNFa_NFkB_model/set_model.py)               | Differential equation, parameters and initial condition                                                  |
| [`observalbe.py`](./TNFa_NFkB_model/observable.py)             | Observables, simulations and experimental data                                                           |
| [`viz.py`](./TNFa_NFkB_model/viz.py)                           | Plotting parameters for customizing figure properties                                                    |
| [`set_search_param.py`](./TNFa_NFkB_model/set_search_param.py) | Lower and upper bounds of model parameters to be estimated                                               |
| [`fitness.py`](./TNFa_NFkB_model/fitness.py)                   | An objective function to be minimized, i.e., the distance between model simulation and experimental data |

## Usage

To get simulation results, please run the following code.

```python
>>> from biomass import Model, run_simulation
>>> model = Model("TNFa_NFkB_model").create()
>>> run_simulation(model, viz_type="average", show_all=True)
```

Then the time-course simulation results will be saved in `TNFa_NFkB_model/figure/simulation/average/`.
