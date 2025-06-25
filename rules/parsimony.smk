rule parsimony_tree:
    output:
        parsimony_tree  = parsimony_tree_file_name,
        log             = parsimony_log_file_name,
    params:
        msa     = lambda wildcards: msas[wildcards.msa],
        prefix  = output_files_parsimony_trees + "seed_{seed}",
        model   = lambda wildcards: iqtree_models[wildcards.msa]
    log:
        "results/{msa}/output_files/parsimony/seed_{seed}.snakelog",
    run:
        # Use iqtree
        cmd = [
            "{iqtree_command} ",
            "-s {params.msa} ",
            "-nt 1 ", 
            "-m {params.model} ",
            "-seed {wildcards.seed} ",
            "-pre {params.prefix} ",
            "> {log} 2>&1"
        ]
        shell("".join(cmd))
