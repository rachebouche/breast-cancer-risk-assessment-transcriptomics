import pandas as pd
from cmapPy.pandasGEXpress.parse import parse

def extract_l1000_subset(raw_data, cell_line, metadata_genes, metadata_compounds, metadata_signatures, landmarks=True):
    
    # Select subset on metadata
    subset_sig_ids = metadata_signatures[(metadata_signatures["cell_iname"] == cell_line)]["sig_id"]

    # Parse raw data with cmapPy and select only subset of data
    subset_only_gctoo = parse(raw_data, make_multiindex= True ) # Parse all data, need lot of RAM

    if landmarks==True:
        landmark_metadata = metadata_genes[metadata_genes["feature_space"].isin(["landmark","best inferred"]) ]
        subset_only_gctoo = subset_only_gctoo.data_df.loc[subset_only_gctoo.data_df.index.isin(landmark_metadata["gene_id"]),subset_only_gctoo.data_df.columns.isin(subset_sig_ids)] # Select only subset of data based on metadata and landmark genes if True
    else :
        subset_only_gctoo = subset_only_gctoo.data_df.loc[:,subset_only_gctoo.data_df.columns.isin(subset_sig_ids)]
    
    
    subset_sig = metadata_signatures[metadata_signatures["sig_id"].isin(subset_only_gctoo.columns.tolist())]

    doses = subset_sig["pert_idose"]
    times = subset_sig["pert_itime"]
    cells = subset_sig["cell_iname"]
    compounds_brd = subset_sig["pert_id"]
    inchikeys = [metadata_compounds[metadata_compounds["pert_id"] == brd]["inchi_key"].tolist()[0] for brd in compounds_brd] # This step take about 2 min depending on the size of subset. Can be optimized by using a dictionary ?

    if len(doses)== len(times) == len(cells) == len(compounds_brd) == len(inchikeys):
        print("All metadata have the good length")
    else:
        print(f"Metadata length are not the same({len(doses)}, {len(times)}, {len(cells)}, {len(compounds_brd)}, {len(inchikeys)})")
        
    allcond = [f"{compound}_{cell}_{dose}_{time}" for compound, cell, dose, time in zip(inchikeys, cells, doses, times)]
    meta_subset =  pd.DataFrame({"sig_id":subset_only_gctoo.columns.tolist(), "brd_id": compounds_brd, "inchikey":inchikeys, "dose":doses, "time":times, "cell":cells, "allcond":allcond}).set_index('sig_id')
    sig_subset = subset_only_gctoo.T.rename_axis(None, axis=1).rename_axis("sig_id", axis=0)
    L1000_subset = pd.concat([meta_subset, sig_subset], axis = 1,)

    return L1000_subset

def main():

    L1000_level5 = "../../data/external/L1000/level5_beta_trt_cp_n720216x12328.gctx"
    cmpd_info = pd.read_csv("../../data/lincs-metadata/compoundinfo_beta.txt", sep="\t")
    gene_info = pd.read_csv("../../data/lincs-metadata/geneinfo_beta.txt", sep="\t", dtype=str)
    sig_info = pd.read_csv( "../../data/lincs-metadata/siginfo_beta.txt", sep="\t")

    L1000_MCF7_l_bi = extract_l1000_subset(raw_data=L1000_level5, cell_line="MCF7", metadata_genes=gene_info, metadata_compounds=cmpd_info, metadata_signatures=sig_info, landmarks=True)

    # Keep, for each compound-dose-time condition, the replicate with the highest
    # transcriptional activity score (TAS) 
    keep_sig_ids = []
    for cond_group in L1000_MCF7_l_bi.groupby("allcond"):
        best_tas_id = sig_info[sig_info["sig_id"].isin(cond_group[1].index.tolist())]["tas"].idxmax()
        sig_selected = sig_info[sig_info["sig_id"].isin(cond_group[1].index.tolist())].loc[best_tas_id]["sig_id"]
        keep_sig_ids.append(sig_selected)

    L1000_MCF7_l_bi_selected = L1000_MCF7_l_bi.loc[keep_sig_ids]

    L1000_MCF7_l_bi_selected.dropna(subset=["inchikey"]).to_csv("../../data/external/L1000/L1000_MCF7_l_bi_bestTAS.csv")

if __name__ == "__main__":
    main()
