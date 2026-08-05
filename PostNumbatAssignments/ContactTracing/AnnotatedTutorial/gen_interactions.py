import pandas as pd

BASE = (
    "https://raw.githubusercontent.com/"
    "ventolab/cellphonedb-data/v4.1.0/data"
)

interactions = pd.read_csv(f"{BASE}/interaction_input.csv")
genes = pd.read_csv(f"{BASE}/gene_input.csv")

# Some UniProt entries occur more than once because they map to multiple
# Ensembl transcripts, but their HGNC symbol is the same.
uniprot_to_hgnc = (
    genes[["uniprot", "hgnc_symbol"]]
    .dropna()
    .drop_duplicates(subset="uniprot")
    .set_index("uniprot")["hgnc_symbol"]
)

# Retain direct protein-protein interactions where both partners are
# individual proteins rather than CellPhoneDB complexes or non-protein entities.
keep = (
    interactions["is_ppi"].eq(True)
    & interactions["partner_a"].isin(uniprot_to_hgnc.index)
    & interactions["partner_b"].isin(uniprot_to_hgnc.index)
)

output = interactions.loc[keep, ["partner_a", "partner_b"]].copy()

output["ligand"] = output["partner_a"].map(uniprot_to_hgnc)
output["receptor"] = output["partner_b"].map(uniprot_to_hgnc)

output = output[["ligand", "receptor"]]

assert len(output) == 980
assert output.duplicated().sum() == 0

output.to_csv(
    "interactions_cellphonedb.txt",
    sep="\t",
    index=False,
)
