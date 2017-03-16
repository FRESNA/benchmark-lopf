from tempfile import mkstemp
import os

configfile: "config.yaml"

localrules: all, setup_network, combine_timing

NOS = {which: range(config['nos_max'][which]) for which in ['ptdf', 'no-ptdf']}

def mem_requirements(wildcards):
    factor = 1.5 if hasattr(wildcards, 'method') else 1.
    if wildcards.formulation == 'ptdf':
        factor *= 3
    elif wildcards.formulation == 'ptdf-flows':
        factor *= 4

    case = wildcards.case

    if (case in {'case2383wp', 'case2869pegase'} and
        wildcards.formulation == 'ptdf-flows' and
        hasattr(wildcards, 'method')):
        factor *= 1.7
    if (case in {'case2383wp', 'case2869pegase'} and
        wildcards.formulation in {'ptdf-flows', 'ptdf'}):
        factor *= 1.3

    if case in {'case118', 'case300', 'scigrid'}:
        return int(factor * 3000)
    elif case in {'case1354pegase', 'case1951rte'}:
        return int(factor * 6000)
    elif case in {'case2383wp'}:
        return int(factor * 7000)
    elif case in {'case2869pegase'}:
        return int(factor * 8000)

rule all:
    input:
        'timings.csv'

rule setup_network:
    output: 'networks/{case}_{mode}_{nhours}_{no}'
    params:
        overwrite_zero_s_nom=lambda wildcards: (2000 if (wildcards.case == 'case1354pegase' and wildcards.mode == 'p')
                                           else 1000)
    script: 'scripts/setup_network.py'

rule write_lp:
    input: 'networks/{case}_{mode}_{nhours}_{no}'
    resources: mem=mem_requirements
    output: temp('lps/{case}_{mode}_{nhours}_{no}_{formulation}.lp')
    script: 'scripts/write_lp.py'

rule solve_lp:
    input:
        'lps/{case}_{mode}_{nhours}_{no}_{formulation}.lp',
        'networks/{case}_{mode}_{nhours}_{no}'
    params:
        gurobi_log='logs/gurobi/{case}_{mode}_{nhours}_{no}_{formulation}_{method}.log',
        header=config['header'],
    threads: 4
    resources: mem=mem_requirements
    output: 'timings/{case}_{mode}_{nhours}_{no}_{formulation}_{method}'
    script: 'scripts/solve_lp.py'

rule combine_timing:
    input:
        expand('timings/{case}_{mode}_{nhours}_{no}_{formulation}_{method}',
               case=config['cases'],
               mode=config['modes'],
               formulation=config['formulations']['no-ptdf'],
               nhours=config['nhours'],
               no=NOS['no-ptdf'],
               method=config['method']),
        expand('timings/{case}_{mode}_{nhours}_{no}_{formulation}_{method}',
               case=config['cases'],
               mode=config['modes'],
               formulation=config['formulations']['ptdf'],
               nhours=config['nhours'],
               no=NOS['ptdf'],
               method=config['method']),
    output: 'timings.csv'
    params: header=config['header']
    run:
        fd, inputfile = mkstemp(text=True)
        with os.fdopen(fd, 'w') as f:
            f.writelines(fn + "\n" for fn in input)
        shell('''
            echo {params.header} > {output}
            xargs -a {inputfile} cat >> {output}
        ''')
        os.unlink(inputfile)
