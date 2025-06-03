rule iqtree_rf_distance_search_trees:
    """
    Rule that computes the RF-Distances between all search trees using IQ-TREE.
    """
    input:
        all_search_trees = rules.collect_search_trees.output.all_search_trees
    output:
        rfDist      = f"{iqtree_tree_inference_dir}inference.iqtree.rfdist",
        rfDist_log  = f"{iqtree_tree_inference_dir}inference.iqtree.rfDistances.log",
    params:
        prefix = f"{iqtree_tree_inference_dir}inference"
    log:
        f"{iqtree_tree_inference_dir}inference.iqtree.rfDistances.snakelog",
    shell:
        "{iqtree_command} "
        "-rf_all " 
        "-t {input.all_search_trees} "
        "-pre {params.prefix} " 
        ">> {output.rfDist_log} "

rule iqtree_rfdistance_eval_trees:
    """
    Rule that computes the RF-Distances between all eval trees using iqtree-NG.
    """
    input:
        all_eval_trees = rules.collect_eval_trees.output.all_eval_trees
    output:
        rfDist      = f"{iqtree_tree_eval_dir}eval.iqtree.rfdist",
        rfDist_log  = f"{iqtree_tree_eval_dir}eval.iqtree.rfDistances.log",
    params:
        prefix = f"{iqtree_tree_eval_dir}eval"
    log:
        f"{iqtree_tree_eval_dir}eval.iqtree.rfDistances.snakelog",
    shell:
        "{iqtree_command} "
        "-rf_all "
        "-t {input.all_eval_trees} "
        "-pre {params.prefix} "
        ">> {output.rfDist_log} "


rule iqtree_rfdistance_plausible_trees:
    """
    Rule that computes the RF-Distances between all plausible trees using iqtree-NG.
    """
    input:
        all_plausible_trees = rules.collect_plausible_trees.output.all_plausible_trees
    output:
        rfDist      = f"{iqtree_tree_eval_dir}plausible.iqtree.rfdist",
        rfDist_log  = f"{iqtree_tree_eval_dir}plausible.iqtree.rfDistances.log",
    params:
        prefix = f"{iqtree_tree_eval_dir}plausible"
    log:
        f"{iqtree_tree_eval_dir}plausible.iqtree.rfDistances.snakelog",
    run:
        num_plausible = len(open(input.all_plausible_trees).readlines())
        # we need this distinction because iqtree-NG requires more than one tree in the input file
        # in order to compute the RF Distance
        # but there might be no plausible trees for a given dataset
        if num_plausible <= 1:
            # write 0.0 as RF-Distance in a dummy log
            with open(output.rfDist_log, "w") as f:
                f.write("""
                Number of unique topologies in this tree set: 1
                Average absolute RF distance in this tree set: 0.0
                Average relative RF distance in this tree set: 0.0
                """)

            with open(output.rfDist, "w") as f:
                f.write("0 1 0.0 0.0")
        else:
            shell("{iqtree_command} -rf_all -t {input.all_plausible_trees} -pre {params.prefix} >> {output.rfDist_log}")


rule iqtree_rfdistance_parsimony_trees:
    """
    Rule that computes the RF-Distances between all parsimony trees inferred with Parsimonator using iqtree-NG.
    """
    input:
        all_parsimony_trees = f"{output_files_parsimony_trees}AllParsimonyTrees.trees",
    output:
        rfDist      = f"{output_files_parsimony_trees}parsimony.iqtree.rfdist",
        rfDist_log  = f"{output_files_parsimony_trees}parsimony.iqtree.rfDistances.log",
    params:
        prefix = f"{output_files_parsimony_trees}parsimony"
    log:
        f"{output_files_parsimony_trees}parsimony.iqtree.rfDistances.snakelog",
    shell:
        "{iqtree_command} "
        "-rf_all "
        "-t {input.all_parsimony_trees} "
        "-pre {params.prefix} "
        ">> {output.rfDist_log} "

