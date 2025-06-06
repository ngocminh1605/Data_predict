rule iqtree_filter_unique_tree_topologies:
    """
    The statistical tests can be biased by the number of trees in the candidate set. 
    Therefore, before perfoming the tests on the eval trees we filter duplicate topologies.
    This rule produces two output:
    - A .trees file containing the unique topologies
    - A .pkl file (pickle file) that contains a list of sets with the clusters, so we can later match 
        IQ-Tree test results to newick tree strings 
    Tránh thiên lệch thống kê do có quá nhiều cây giống nhau (về topology) trong tập candidate trees.
    Chúng ta sẽ lọc ra các cây có topology duy nhất.
    input tat ca cay tu eval_tree, file log chua rf cac cay 
    output 
    .trees: chứa các cây duy nhất (1 đại diện cho mỗi nhóm topology).
    .pkl: lưu danh sách các cluster (mỗi cluster là tập các cây giống nhau).
        Dùng để ánh xạ lại kết quả IQ-TREE sau này.
    """
    input:
        all_eval_trees              = rules.collect_eval_trees.output.all_eval_trees,
        eval_trees_rfdistances_log  = rules.iqtree_rfdistance_eval_trees.output.rfDist_log,
    output:
        filtered_trees  = f"{output_files_iqtree_dir}filteredEvalTrees.trees",
        clusters        = f"{output_files_iqtree_dir}filteredEvalTrees.clusters.pkl",
    script:
        "scripts/filter_tree_topologies.py"


rule iqtree_significance_tests_on_eval_trees:
    """
    Perfoms all significance tests as implemented in IQ-Tree on the set of filtered trees.
    As reference tree for estimating the model parameters, we pass the best tree 
    (i.e. with the highest log-likelihood) of the dataset
    """
    input:
        filtered_trees  = rules.iqtree_filter_unique_tree_topologies.output.filtered_trees,
        best_tree       = rules.save_best_eval_tree.output.best_eval_tree
    output:
        summary     = f"{output_files_iqtree_dir}significance.iqtree",
        iqtree_log  = f"{output_files_iqtree_dir}significance.iqtree.log",
    params:
        msa         = lambda wildcards: msas[wildcards.msa],
        data_type   = lambda wildcards: data_types[wildcards.msa],
        prefix      = f"{output_files_iqtree_dir}significance",
        model       = lambda wildcards: iqtree_models[wildcards.msa],
        model_str   = "-p" if partitioned else "-m",
        threads     = config["software"]["iqtree"]["threads"]
    log:
        f"{output_files_iqtree_dir}significance.iqtree.snakelog",
    run:
        morph = "-st MORPH " if params.data_type == "MORPH" else ""
        shell("{iqtree_command} "
        "-s {params.msa} "
        "{morph} "
        "{params.model_str} {params.model} "
        "-pre {params.prefix} "
        "-z {input.filtered_trees} "
        "-te {input.best_tree} "
        "-n 0 "
        "-zb 10000 "
        "-zw "
        "-au "
        "-nt {params.threads} "
        "-seed 0 "
        "> {output.iqtree_log} ")
