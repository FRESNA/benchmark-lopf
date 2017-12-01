# coding: utf-8

import logging
logging.basicConfig(level=logging.WARNING)

import pypsa  # noqa: E402
pypsa.io.logger.setLevel(logging.ERROR)

import numpy as np
import pandas as pd

import networkx as nx

class AbortException(Exception):
    pass


def save_lp_and_abort(network, snapshots):
    raise AbortException()


def write_lp_file(network, formulation, fn, type, io_options={}):
    try:
        if type == "lopf":
            network.lopf(snapshots=network.snapshots, formulation=formulation,
                         extra_functionality=save_lp_and_abort, ptdf_tolerance=1e-5)
        elif type.startswith("sclopf"):
            # Check branch_outages that do not disconnect the graph
            network.lopf(snapshots=network.snapshots, formulation="kirchhoff",
                         solver_name="gurobi", solver_options={"Method": 1})

            branch_outages = []
            G = network.graph()
            for u, v, k in G.edges(keys=True):
                H = pypsa.descriptors.OrderedGraph(G)
                H.remove_edge(u, v, k)
                if nx.is_connected(H):
                    # branch_outages.append(k)
                    branch_outages.append((k.__class__.__name__, k.name))
            branch_outages = pd.MultiIndex.from_tuples(branch_outages)

            line_loading = pd.concat(dict(Line=(abs(network.lines_t.p0).mean()/network.lines.s_nom)))
            if type[len("sclopf")] == "-":
                n_outages = int(type[len("sclopf-"):])
                branch_outages = line_loading.loc[branch_outages].nlargest(n_outages).index
            else:
                branch_outages = line_loading.loc[branch_outages].loc[lambda s: s>0.5].index

            assert len(branch_outages) > 0

            # Load shedding
            network.add("Carrier", "load")
            network.import_components_from_dataframe(
                pd.DataFrame(
                    dict(marginal_cost=1e4,
                         # intersect between macroeconomic and surveybased willingness to pay
                         # http://journal.frontiersin.org/article/10.3389/fenrg.2015.00055/full
                         p_nom=1e6,
                         carrier='load',
                         bus=network.buses.index),
                    index=network.buses.index + ' load'),
                "Generator"
            )

            network.sclopf(branch_outages=branch_outages,
                           snapshots=network.snapshots,
                           formulation=formulation,
                           extra_functionality=save_lp_and_abort,
                           ptdf_tolerance=1e-5)

        elif type.startswith("lopf-ne"):
            network.generators['p_nom_extendable'] = True
            network.storage_units['p_nom_extendable'] = True

            network.lopf(snapshots=network.snapshots, formulation=formulation,
                         extra_functionality=save_lp_and_abort, ptdf_tolerance=1e-5)
    except AbortException as e:
        pass
    network.model.write(fn, io_options=io_options)


if __name__ == '__main__':
    network = pypsa.Network(snakemake.input[0])
    network.lines['type'] = np.nan
    network.transformers['type'] = np.nan

    write_lp_file(network, fn=snakemake.output[0],
                formulation=snakemake.wildcards.formulation,
                type=snakemake.wildcards.type)
