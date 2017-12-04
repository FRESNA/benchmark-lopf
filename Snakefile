from tempfile import mkstemp
import os

configfile: "config.yaml"

localrules: all, setup_network, combine_timing

wildcard_constraints:
    type="[a-zA-Z0-9]+",
    solver="[a-zA-Z0-9]+",

def mem_requirements(wildcards):
    if wildcards.type == 'sclopf':
        return 12000 # 58000 # 3000
    elif wildcards.type == 'lopf-ne':
        return 1000
    elif wildcards.type == 'lopf-ne-year':
        return 12000

    factor = 1.

    if wildcards.type == 'lopf-ne':
        factor *= 2
    if hasattr(wildcards, 'method'):
        factor *= 1.5
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
        expand('timings-{type}-{solver}.csv', type=config['types'], solver=config['solvers'])

rule setup_network:
    output: 'networks/{case}_{mode}_{nhours}_{no}'
    params:
        overwrite_zero_s_nom=lambda wildcards: (2000 if (wildcards.case == 'case1354pegase' and wildcards.mode == 'p')
                                           else 1000)
    script: 'scripts/setup_network.py'

rule write_lp:
    input: 'networks/{case}_{mode}_{nhours}_{no}'
    resources: mem=mem_requirements
    output: temp('lps/{type}/{case}_{mode}_{nhours}_{no}_{formulation}.lp')
    script: 'scripts/write_lp.py'

rule solve_lp:
    input:
        'lps/{type}/{case}_{mode}_{nhours}_{no}_{formulation}.lp',
        'networks/{case}_{mode}_{nhours}_{no}'
    output: 'timings/{type}_{solver}/{case}_{mode}_{nhours}_{no}_{formulation}_{method}'
    params:
        solver_log='logs/{type}_{solver}/{case}_{mode}_{nhours}_{no}_{formulation}_{method}.log',
        header=config['header'],
    threads: 4
    resources: mem=mem_requirements
    script: 'scripts/solve_lp.py'

def combine_timing_input(wildcards):
    c = config[wildcards.type]
    return sum((expand('timings/{type}_{solver}/{case}_{mode}_{nhours}_{no}_{formulation}_{method}',
                       type=wildcards.type,
                       solver=wildcards.solver,
                       case=c['cases'],
                       mode=c['modes'],
                       formulation=c['formulations'][yesnoptdf],
                       nhours=c['nhours'],
                       no=range(c['nos_max'][yesnoptdf]),
                       method=c['method'])
                for yesnoptdf in ('ptdf', 'no-ptdf')),
               [])

rule combine_timing:
    input: combine_timing_input
    output: 'results/timings-{type}-{solver}.csv'
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
