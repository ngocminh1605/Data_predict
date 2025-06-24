rule iqtree_pars_tree:
    """
    Rule that infers a single tree based on a parsimony starting tree using iqtree-NG.
    # file output là cây tối ưu, cây khởi tạo, mô hình tiến hóa tốt nhất, log
    """
    output:
        iqtree_best_tree     = f"{iqtree_tree_inference_prefix_pars}.treefile",
        iqtree_starting_tree = f"{iqtree_tree_inference_prefix_pars}.log",
        iqtree_best_model    = f"{iqtree_tree_inference_prefix_pars}.ckp.gz",
        iqtree_log           = f"{iqtree_tree_inference_prefix_pars}.iqtree",
    # truyền các giá trị trung gian vào lệnh shell, lambda là hàm ẩn danh truy cập các giá trị đại diện
    params:
        prefix  = iqtree_tree_inference_prefix_pars,
        msa     = lambda wildcards: msas[wildcards.msa],
        model   = lambda wildcards: iqtree_models[wildcards.msa],
        threads = config["software"]["iqtree"]["threads"]
    #  file log riêng của Snakemake, ghi lại toàn bộ thông tin quá trình thực thi lệnh shell
    log:
        f"{iqtree_tree_inference_prefix_pars}.snakelog",
    shell:
        "{iqtree_command} "
        "-s {params.msa} "
        "-m {params.model} "
        "-ninit 1 "
        "-pre {params.prefix} "
        "-nt {params.threads} "
        "-seed {wildcards.seed} "
         "> {log} 2>&1"

rule iqtree_rand_tree:
    """
    Rule that infers a single tree based on a random starting tree using iqtree-NG.
    """
    output:
        iqtree_best_tree     = f"{iqtree_tree_inference_prefix_rand}.treefile",
        iqtree_starting_tree = f"{iqtree_tree_inference_prefix_rand}.log",
        iqtree_best_model    = f"{iqtree_tree_inference_prefix_rand}.ckp.gz",
        iqtree_log           = f"{iqtree_tree_inference_prefix_rand}.iqtree",
    params:
        prefix  = iqtree_tree_inference_prefix_rand,
        msa     = lambda wildcards: msas[wildcards.msa],
        model   = lambda wildcards: iqtree_models[wildcards.msa],
        threads = config["software"]["iqtree"]["threads"]
    log:
        f"{iqtree_tree_inference_prefix_rand}.snakelog",
    shell:
        "{iqtree_command} "
        "-s {params.msa} "
        "-m {params.model} "
        "-pre {params.prefix} "
        "-seed {wildcards.seed} "
        "-nt {params.threads} "
        "-t RANDOM "
         "> {log} 2>&1"
