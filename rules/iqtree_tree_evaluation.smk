rule reevaluate_iqtree_pars_tree:
    """
    Rule that re-evaluates the given parsimony search tree.
    """
    input:
        best_tree_of_run    = f"{iqtree_tree_inference_prefix_pars}.treefile"
    output:
        log         = f"{iqtree_tree_eval_prefix_pars}.log",
        best_tree   = f"{iqtree_tree_eval_prefix_pars}.treefile",
        eval_log    = f"{iqtree_tree_eval_prefix_pars}.iqtree",
    params:
        prefix  = iqtree_tree_eval_prefix_pars,
        msa     = lambda wildcards: msas[wildcards.msa],
        model   = lambda wildcards: iqtree_models[wildcards.msa],
        threads = config["software"]["iqtree"]["threads"]
    log:
        f"{iqtree_tree_eval_prefix_pars}.snakelog"
    shell:
        "{iqtree_command} "
        "-s {params.msa} "
        "-m {params.model} "
        "-te {input.best_tree_of_run} "
        "-seed 0 "
        "-pre {params.prefix} "
        "-nt 2 "
         "> {log} 2>&1"


rule reevaluate_iqtree_rand_tree:
    """
    Rule that re-evaluates the given random search tree.
    """
    input:
        best_tree_of_run    = f"{iqtree_tree_inference_prefix_rand}.treefile"
    output:
        log         = f"{iqtree_tree_eval_prefix_rand}.log",
        best_tree   = f"{iqtree_tree_eval_prefix_rand}.treefile",
        eval_log    = f"{iqtree_tree_eval_prefix_rand}.iqtree",
    params:
        prefix  = iqtree_tree_eval_prefix_rand,
        msa     = lambda wildcards: msas[wildcards.msa],
        model   = lambda wildcards: iqtree_models[wildcards.msa],
        threads = config["software"]["iqtree"]["threads"]
    log:
        f"{iqtree_tree_eval_prefix_rand}.snakelog"

    shell:
        "{iqtree_command} "
        "-s {params.msa} "
        "-te {input.best_tree_of_run} "
        "-m {params.model} "
        "-pre {params.prefix} "
        "-nt {params.threads} "
        "-seed 0 "
        "> {log} 2>&1"
