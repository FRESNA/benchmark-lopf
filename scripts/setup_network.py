# coding: utf-8

import os
import subprocess
from tempfile import mkstemp

import pandas as pd
import scipy.io
import numpy as np
import re
import pypsa

import logging
pypsa.io.logger.setLevel(logging.ERROR)

octave = "octave"
def read_matpowercase(case, matpower_dir="data/matpower_cases"):  # noqa: E302
    fd, tmpfn = mkstemp(); os.close(fd)  # noqa: E702

    proc = subprocess.Popen([octave, "--no-gui"],
                            stdin=subprocess.PIPE,
                            cwd=matpower_dir)
    proc.communicate("c = {case};\nsave -6 {tmpfn} c\n"
                     .format(case=case, tmpfn=tmpfn)
                     .encode('utf-8'))
    if proc.poll() != 0:
        raise RuntimeError("octave did not return successfully".format(octave))

    s = scipy.io.loadmat(tmpfn, squeeze_me=True, struct_as_record=False)['c']
    os.remove(tmpfn)
    return {f: getattr(s, f)
            for f in {'version', 'baseMVA', 'bus',
                      'gen', 'branch', 'gencost', 'bus_name'}
            if hasattr(s, f)}


def load_case(case, overwrite_zero_s_nom=1e3):
    ppc = read_matpowercase(case)
    network = pypsa.Network()
    network.name = "MATPOWER {}".format(case)
    pypsa.io.import_from_pypower_ppc(network, ppc,
                                     overwrite_zero_s_nom=overwrite_zero_s_nom)
    network.buses.drop('type', axis=1, errors='ignore', inplace=True)
    network.lines.drop('type', axis=1, errors='ignore', inplace=True)
    network.transformers.drop('type', axis=1, errors='ignore', inplace=True)

    return network


def extend_load_variation(network, nhours=24*2, scale=0.2, mode=''):  # noqa: E226
    if len(network.snapshots) == 1:
        network.set_snapshots(pd.RangeIndex(nhours))
        network.loads_t.p_set = pd.DataFrame(
            network.loads.p_set.values *
            (1. - abs(np.random.normal(scale=scale, size=(nhours, len(network.loads))))),
            index=network.snapshots, columns=network.loads.index
        )

    network.now = network.snapshots[0]


def modify_generators(network, nhours=24*2):  # noqa: E226
    from vresutils import costdata as vcostdata
    ocgt_costs = vcostdata.get_cost('diw').loc['OCGT']
    network.generators.marginal_cost = (
        ocgt_costs['wbi'] * (1 + 0.1 * np.random.random(len(network.generators)))
    )


def add_storage_unit(network, mode, max_hours=6, p_nom=10):
    m = re.search(r's(\d+)', mode)
    assert m is not None, "storage modes must be defined with a number"
    n_storage = int(m.group(1))
    p_nom = (network.loads_t.p_set.mean().groupby(network.loads.bus).sum() / 3.).nlargest(n_storage)
    storage_df = pd.DataFrame(
        dict(max_hours=max_hours,
             marginal_cost=1e-2 + 2e-3 * (np.random.random(len(p_nom.index)) - 0.5),
             p_nom=p_nom,
             p_max_pu=+1, p_min_pu=-1.,
             efficiency_store=0.9,
             efficiency_dispatch=0.9,
             cyclic_state_of_charge=True,
             bus=p_nom.index)
    )

    network.import_components_from_dataframe(storage_df, "StorageUnit")


def add_renewable_generators(network, mode):
    scigrid_network = pypsa.Network(csv_folder_name="data/scigrid-with-load-gen-trafos-96")

    wind_gens_b = scigrid_network.generators.carrier == 'Wind Onshore'
    wind_gens = np.random.choice(scigrid_network.generators.index[wind_gens_b],
                                 len(network.buses))
    wind_index = 'WO' + network.buses.index

    p_max_pu = pd.DataFrame(scigrid_network.generators_t.p_max_pu[wind_gens]
                            .values[:len(network.snapshots)],
                            index=network.snapshots, columns=wind_index)
    p_max_pu.where(lambda df: df > 1e-2, other=0., inplace=True)

    p_nom = network.loads_t.p_set.mean().mean() / p_max_pu.mean().mean()

    network.import_components_from_dataframe(
        pd.DataFrame(
            dict(carrier='windon',
                 p_nom=p_nom,
                 marginal_cost=(1e-2 + 2e-3 * (np.random.random(len(network.buses.index)) - 0.5)),
                 bus=network.buses.index),
            index=wind_index),
        "Generator"
    )

    network.import_series_from_dataframe(
        p_max_pu,
        "Generator",
        "p_max_pu"
    )


def modify_network_according_to_metadata(network, meta, extend_load=True):
    if 's' in meta['mode']:
        add_storage_unit(network, mode=meta['mode'])

    if extend_load:
        extend_load_variation(network, nhours=meta['nhours'], mode=meta['mode'])

    modify_generators(network, nhours=meta['nhours'])

    if 'r' in meta['mode']:
        add_renewable_generators(network, meta['mode'])


def load_network(meta):
    if meta.get('case', 'scigrid') != 'scigrid':
        network = load_case(meta['case'], overwrite_zero_s_nom=meta.get('overwrite_zero_s_nom', 1e3))
        modify_network_according_to_metadata(network, meta, extend_load=True)
    else:
        network = pypsa.Network(csv_folder_name="data/scigrid-with-load-gen-trafos-96")

        network.generators.query("(carrier != 'Wind Onshore') and "
                                 "(carrier != 'Wind Offshore') and "
                                 "(carrier != 'Solar')", inplace=True)
        network.generators_t.p_max_pu.drop(network.generators_t.p_max_pu.columns,
                                           axis=1, inplace=True)
        network.storage_units.drop(network.storage_units.index, inplace=True)

        # There are some infeasibilities at the edge of the network where loads
        # are supplied by foreign lines - just add some extra capacity in Germany
        for line_name in ["350", "583", "316", "602"]:
            network.lines.loc[line_name, "s_nom"] += 500

        network.set_snapshots(network.snapshots[:meta['nhours']])
        modify_network_according_to_metadata(network, meta, extend_load=False)

    return network


if __name__ == "__main__":
    logging.basicConfig(level=logging.WARNING)
    meta = dict(snakemake.wildcards.items())
    meta['nhours'] = int(meta['nhours'])
    meta['overwrite_zero_s_nom'] = snakemake.params.overwrite_zero_s_nom

    network = load_network(meta)
    network.export_to_csv_folder(snakemake.output[0])
