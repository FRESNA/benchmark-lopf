# coding: utf-8

import logging
logging.basicConfig(level=logging.WARNING)

import pypsa  # noqa: E402
pypsa.io.logger.setLevel(logging.ERROR)


class AbortException(Exception):
    pass


def save_lp_and_abort(network, snapshots):
    raise AbortException()


def write_lp_file(network, formulation, fn, io_options={}):
    try:
        network.lopf(snapshots=network.snapshots, formulation=formulation,
                     extra_functionality=save_lp_and_abort, ptdf_tolerance=1e-5)
    except AbortException as e:
        pass
    network.model.write(fn, io_options=io_options)


network = pypsa.Network(snakemake.input[0])
write_lp_file(network, fn=snakemake.output[0],
              formulation=snakemake.wildcards.formulation)
