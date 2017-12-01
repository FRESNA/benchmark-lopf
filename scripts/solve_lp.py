# coding: utf-8

import logging
logging.basicConfig(level=logging.WARNING)

from vresutils import timer
import gurobipy as grb
from itertools import chain
import numpy as np
import re, os
import subprocess as sp
from textwrap import dedent

import pypsa
pypsa.io.logger.setLevel(logging.ERROR)

def run_timing_gurobi(log_fn, lp_fn, method, network):
    nums = {'read_time': np.nan,
            'presol_time': np.nan,
            'sol_time': np.nan}

    env = grb.Env(log_fn)
    with timer(verbose=False) as t_read:
        model = grb.read(lp_fn, env)
    nums['read_time'] = t_read.usec / 1e6

    model.setParam("LogToConsole", 0)
    model.setParam("Method", method)
    model.setParam("TimeLimit", 3*60*60)  # noqa: E226
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

    nums['Constrs'] = model.NumConstrs
    nums['Vars'] = model.NumVars
    nums['NZs'] = model.NumNZs

    nums['objective'] = getattr(model, 'ObjVal', np.nan)

    return nums


def run_timing_cplex(log_fn, lp_fn, method, network):
    nums = {'read_time': np.nan,
            'presol_time': np.nan,
            'sol_time': np.nan}

    assert method == 3

    p = sp.run(["cplex"],
               input=dedent("""
        set parallel -1
        set threads 4
        set lpmethod 6
        set timelimit {timelimit}
        set logfile {log_fn}
        read {lp_fn}
        display problem stats
        optimize
               """).format(timelimit=3*60*60, log_fn=log_fn, lp_fn=lp_fn),
               universal_newlines=True,
               stdout=sp.PIPE)

    stdout = p.stdout

    m = re.search("^Read time = ([\d\.]+) sec\.", stdout, flags=re.MULTILINE)
    if m is not None:
        nums['read_time'] = float(m.group(1))

    m = re.search(r"^Variables\s+:\s+(\d+).*\n"
                r".*\n.*\n"
                r"Linear constraints\s+:\s+(\d+).*\n"
                r"\s+Nonzeros\s+:\s+(\d+)\n", stdout, flags=re.MULTILINE)
    if m is not None:
        nums["Vars"] = int(m.group(1))
        nums["Constrs"] = int(m.group(2))
        nums["NZs"] = int(m.group(3))

    m = re.search(r"^Presolve time = ([\d\.]+) sec\.", stdout, flags=re.MULTILINE)
    if m is not None:
        nums["presol_time"] = float(m.group(1))

    m = re.search(r"^Reduced LP has (\d+) rows, (\d+) columns, and (\d+) nonzeros\.", stdout, flags=re.MULTILINE)
    if m is not None:
        nums['pConstrs'] = int(m.group(1))
        nums['pVars'] = int(m.group(2))
        nums['pNZs'] = int(m.group(3))
    else:
        nums['pConstrs'] = np.nan
        nums['pVars'] = np.nan
        nums['pNZs'] = np.nan

    m = re.search("^(\w+) solved model\.", stdout, flags=re.MULTILINE)
    if m is not None:
        nums['solved_with'] = m.group(1).lower()

    m = re.search("^Solution time =\s+([\d\.]+) sec\.", stdout, flags=re.MULTILINE)
    if m is not None:
        nums['sol_time'] = float(m.group(1))

    m = re.search("Objective =\s+([\d\.e+]+)$", stdout, flags=re.MULTILINE)
    if m is not None:
        nums['objective'] = float(m.group(1))
    else:
        nums['objective'] = np.nan

    return nums

network = pypsa.Network(csv_folder_name=snakemake.input[1])
log_fn = snakemake.params.solver_log
log_dir = os.path.dirname(log_fn)
if not os.path.isdir(log_dir):
    os.makedirs(os.path.dirname(log_fn))

if snakemake.wildcards.solver == 'gurobi':
    nums = run_timing_gurobi(log_fn,
                             snakemake.input[0], int(snakemake.wildcards.method),
                             network)
elif snakemake.wildcards.solver == 'cplex':
    nums = run_timing_cplex(log_fn,
                            snakemake.input[0], int(snakemake.wildcards.method),
                            network)
else:
    raise NotImplementedError

nums['N'] = len(network.buses)
nums['L'] = len(network.lines) + len(network.transformers)
columns = dict(chain(snakemake.wildcards.items(), nums.items()))
with open(snakemake.output[0], 'w') as f:
    f.write(','.join(str(columns.get(k)) for k in snakemake.params.header.split(',')) + "\n")
