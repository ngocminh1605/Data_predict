import os
import sys
from pathlib import Path

sys.path.append("rules/scripts")

from predict.msa import MSA

configfile: "config.yaml"

# Lấy thông tin từ config
raxmlng_command = config["software"]["raxml-ng"]["command"]
iqtree_command = config["software"]["iqtree"]["command"]

num_pars_trees = config["_debug"]["_num_pars_trees"]
num_rand_trees  = config["_debug"]["_num_rand_trees"]
num_parsimony_trees = config["_debug"]["_num_parsimony_trees"]

# Số giá trị khởi tạo
pars_seeds = range(num_pars_trees)
rand_seeds = range(num_pars_trees, num_pars_trees + num_rand_trees)
parsimony_seeds = range(1, num_parsimony_trees + 1)

# Định nghĩa đường dẫn đến file MSA (chỉ có một MSA)
msa_paths = config["msa_paths"]
part_paths_raxmlng = []
part_paths_iqtree = []
partitioned = False

# kiểm tra tên file trùng lặp
if isinstance(msa_paths[0], list):
    # in this case the MSAs are partitioned
    msa_paths, part_paths_raxmlng, part_paths_iqtree = zip(*msa_paths)
    partitioned = True

#  tạo dict name_msa : path
# This assumes, that each msa
msa_names = [os.path.split(pth)[1] for pth in msa_paths]
msas = dict(zip(msa_names, msa_paths))


# xác định kiểu dữ liệu DNA
data_types = {}
for msa, name in zip(msa_paths, msa_names):
    msa = MSA(msa)
    data_types[name] = msa.data_type

# Lựa chọn mô hình tiến hóa phù hợp
if partitioned:
    raxmlng_models = dict(list(zip(msa_names, part_paths_raxmlng)))
    iqtree_models = dict(list(zip(msa_names, part_paths_iqtree)))
else:
    # infer the data type for each MSA
    raxmlng_models = []
    iqtree_models = []
    for name, msa in msas.items():
        msa = MSA(msa)
        raxmlng_model = msa.get_raxmlng_model()
        raxmlng_models.append((name, raxmlng_model))

        if msa.data_type == "MORPH":
            iqtree_models.append((name, "MK"))
        else:
            iqtree_models.append((name, f"{raxmlng_model}4+FO"))

    raxmlng_models = dict(raxmlng_models)
    iqtree_models = dict(iqtree_models)


# Thiết lập đường dẫn output
outdir = config["outdir"]
db_path = outdir + "{msa}/"  # Thư mục riêng cho từng MSA
output_files_dir = outdir + "{msa}/output_files/"  # Nơi chứa file output cho mỗi MSA

# File paths cho IQ-Tree
output_files_iqtree_dir = output_files_dir + "iqtree/"
iqtree_tree_inference_dir = output_files_iqtree_dir + "inference/"
iqtree_tree_inference_prefix_pars = iqtree_tree_inference_dir + "pars_{seed}"
iqtree_tree_inference_prefix_rand = iqtree_tree_inference_dir + "rand_{seed}"

# Tree evaluation
iqtree_tree_eval_dir = output_files_iqtree_dir + "evaluation/"
iqtree_tree_eval_prefix_pars = iqtree_tree_eval_dir + "pars_{seed}"
iqtree_tree_eval_prefix_rand = iqtree_tree_eval_dir + "rand_{seed}"

# Parsimony trees
output_files_parsimony_trees = output_files_dir + "parsimony/"
parsimony_tree_file_name = output_files_parsimony_trees + "seed_{seed}.treefile"
parsimony_log_file_name = output_files_parsimony_trees + "seed_{seed}.log"


# rule chính tạo training data. parquet cho từng MSA 
rule all:
    input:
        expand(f"{db_path}training_data.parquet", msa = msa_names)

include: "rules/iqtree_tree_inference.smk"
include: "rules/iqtree_tree_evaluation.smk"
include: "rules/collect_data.smk"
include: "rules/iqtree_rfdistance.smk"
include: "rules/iqtree_significance_tests.smk"
include: "rules/msa_features.smk"
include: "rules/parsimony.smk"
include: "rules/save_data.smk"



