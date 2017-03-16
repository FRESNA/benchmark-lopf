# coding: utf-8

import logging
logging.basicConfig(level=logging.WARNING)

from vresutils import timer
import gurobipy as grb
from itertools import chain
import numpy as np
import re

import pypsa
pypsa.io.logger.setLevel(logging.ERROR)

def run_timing(log_fn, lp_fn, method, network):
    nums = {'read_time': np.nan,
            'presol_time': np.nan,
            'sol_time': np.nan}

    env = grb.Env(log_fn)
    with timer(verbose=False) as t_read:
        model = grb.read(lp_fn, env)
    nums['read_time'] = t_read.usec / 1e6

    model.setParam("LogToConsole", 0)
    model.setParam("Method", method)
    model.setParam("TimeLimit", 2*60*60)  # noqa: E226
    if method == 3:
        model.setParam("BarHomogeneous", 1)
        model.setParam("Threads", 4)

    nums['presol_time'] = 0
    with timer(verbose=False) as t_solve:
        model.optimize()
    if model.Status == grb.GRB.OPTIMAL:
        nums['sol_time'] = t_solve.usec / 1e6

    # Extract info from log file
    with open(log_fn) as log_fp:
        logs = log_fp.read()
    m = re.search(r"Presolve time: ([\d\.]+)s", logs)
    if m is not None:
        nums['presol_time'] = float(m.group(1))
        nums['sol_time'] -= nums['presol_time']
    m = re.search(r"Presolved: (\d+) rows, (\d+) columns, (\d+) nonzeros", logs)
    if m is not None:
        nums['pConstrs'] = int(m.group(1))
        nums['pVars'] = int(m.group(2))
        nums['pNZs'] = int(m.group(3))
    else:
        nums['pConstrs'] = np.nan
        nums['pVars'] = np.nan
        nums['pNZs'] = np.nan

    m = re.search(r"Solved with ([a-z ]+)", logs)
    if m is not None:
        nums['solved_with'] = m.group(1)

    nums['N'] = len(network.buses)
    nums['L'] = len(network.lines) + len(network.transformers)
    nums['Constrs'] = model.NumConstrs
    nums['Vars'] = model.NumVars
    nums['NZs'] = model.NumNZs

    return nums

nums = run_timing(snakemake.params.gurobi_log,
                  snakemake.input[0], int(snakemake.wildcards.method),
                  pypsa.Network(csv_folder_name=snakemake.input[1]))
columns = dict(chain(snakemake.wildcards.items(), nums.items()))

with open(snakemake.output[0], 'w') as f:
    f.write(','.join(str(columns[k]) for k in snakemake.params.header.split(',')) + "\n")
