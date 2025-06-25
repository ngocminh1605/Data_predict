import json

from predict.msa import MSA
from predict.iqtree import IQTree

msa_file = snakemake.params.msa
model = snakemake.params.model

msa = MSA(msa_file)

# the Biopython DistanceCalculator does not support morphological data
# so for morphological data we cannot compute the treelikeness at the moment
# compute_treelikeness = msa.data_type != "MORPH"

# iqtree = IQTree(snakemake.params.iqtree_command)

# patterns, gaps, invariant = iqtree.get_patterns_gaps_invariant(msa_file, model)

msa_features = {
    "taxa": msa.n_taxa,
    "sites": msa.n_sites,
    "patterns": msa.n_patterns,  # cached_property → không cần dấu ()
    "gaps": msa.proportion_gaps,
    "invariant": msa.proportion_invariant,
    "entropy": msa.entropy(),  # hàm thường → cần gọi với ()
    # "column_entropies": msa.column_entropy(),  # nếu có hàm này
    "bollback": msa.bollback_multinomial(),
    # "treelikeness": msa.treelikeness_score() if compute_treelikeness else None,
}

# lưu kết quả vào file json
with open(snakemake.output.msa_features, "w") as f:
    json.dump(msa_features, f)
