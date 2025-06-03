rule save_data:
    input:
        # một danh sách của các tên file hoặc đường dẫn bằng cách thay thế các tham số trong chuỗi format
        # Tree seach tree files and logs
        # Lưu trữ và thu thập tất cả các dữ liệu cần thiết cho việc phân tích vào một cơ sở dữ liệu SQLite.
        pars_search_trees   = expand(iqtree_tree_inference_prefix_pars + ".treefile",seed=pars_seeds,allow_missing=True),
        pars_starting_trees = expand(iqtree_tree_inference_prefix_pars + ".log",seed=pars_seeds,allow_missing=True),
        pars_search_logs    = expand(iqtree_tree_inference_prefix_pars + ".iqtree",seed=pars_seeds,allow_missing=True),
        
        rand_search_trees   = expand(iqtree_tree_inference_prefix_rand + ".treefile", seed=rand_seeds, allow_missing=True),
        rand_search_logs    = expand(iqtree_tree_inference_prefix_rand + ".iqtree",seed=rand_seeds,allow_missing=True),
        
        search_logs_collected = f"{iqtree_tree_inference_dir}AllSearchLogs.log",

        # Tree search tree RFDistance logs
        search_rfdistance = f"{iqtree_tree_inference_dir}inference.iqtree.rfDistances.log",

        # Eval tree files and logs
        pars_eval_trees = expand(iqtree_tree_eval_prefix_pars + ".treefile", seed=pars_seeds, allow_missing=True),
        pars_eval_logs  = expand(iqtree_tree_eval_prefix_pars + ".iqtree",seed=pars_seeds,allow_missing=True),
        rand_eval_trees = expand(iqtree_tree_eval_prefix_rand + ".treefile", seed=rand_seeds, allow_missing=True),
        rand_eval_logs  = expand(iqtree_tree_eval_prefix_rand + ".iqtree",seed=rand_seeds,allow_missing=True),
        eval_logs_collected = f"{iqtree_tree_eval_dir}AllEvalLogs.log",

        # Eval tree RFDistance logs
        eval_rfdistance = f"{iqtree_tree_eval_dir}eval.iqtree.rfDistances.log",

        # Plausible tree RFDistance logs
        plausible_rfdistance = f"{iqtree_tree_eval_dir}plausible.iqtree.rfDistances.log",
        plausible_trees_collected = f"{iqtree_tree_eval_dir}AllPlausibleTrees.trees",

        # IQ-Tree significance test results and clusters
        iqtree_results  = f"{output_files_iqtree_dir}significance.iqtree",
        clusters        = f"{output_files_iqtree_dir}filteredEvalTrees.clusters.pkl",

        # MSA Features
        msa_features = f"{output_files_dir}msa_features.json",

        # Parsimony Trees and logs
        parsimony_trees = f"{output_files_parsimony_trees}AllParsimonyTrees.trees",
        parsimony_logs = f"{output_files_parsimony_trees}AllParsimonyLogs.log",
        parsimony_rfdistance = f"{output_files_parsimony_trees}parsimony.iqtree.rfDistances.log",
    output:
        database = "{msa}_data.sqlite3"
    params:
        iqtree_command = iqtree_command,
        msa             = lambda wildcards: msas[wildcards.msa],
    script:
        "scripts/save_data.py"


rule move_db:
    # due to an issue with our lab webservers, I cannot directly create the database on the mounted fs
    # therefore I creat it in the current workdir and then move it to the mounted fs
    input:
        "{msa}_data.sqlite3"
    output:
        database = f"{db_path}data.sqlite3"
    shell:
        "mv {input} {output}"


rule database_to_training_dataframe:
    input:
        database = rules.move_db.output.database,
    output:
        dataframe = f"{db_path}training_data.parquet"
    params:
        num_pars_trees = num_pars_trees,
        num_rand_trees = num_rand_trees,
        num_parsimony_trees = num_parsimony_trees
    script:
        "scripts/database_to_dataframe.py"