"""Functions shared by the analysis notebooks.

Imported as:

    import sys; sys.path.insert(0, "..")
    from common import gene_ids_to_symbols, g_scores, enrich
"""

import numpy as np
import pandas as pd

Z_CUTOFF = 2  # default z-score threshold for calling a gene deregulated


def gene_ids_to_symbols(gene_info, gene_ids):
    """Map L1000 numeric gene ids to HGNC symbols."""
    return gene_info[gene_info["gene_id"].isin(gene_ids)]["gene_symbol"].tolist()


def g_scores(data, z_cutoff=Z_CUTOFF, direction="up", return_scores=False):
    """Genes consistently deregulated across the signatures of a compound class.

    For each gene, the signatures are split into those passing the z-score
    threshold (s) and those not (ns); the score weights the mean of each group
    by its share of the signatures:

        G = |mean(s) * (1 + |s|/N)| - |mean(ns) * (1 + |ns|/N)|

    (the absolute values apply to the "down" direction only, so that strong
    negative z-scores score high). Genes above an overall cut-off derived from
    the same means are returned.

    `data` holds signatures as rows and genes as columns (z-scores).
    Returns the selected gene ids, or `(ids, scores, cutoff)` if
    `return_scores` is set.
    """
    if direction not in ("up", "down"):
        raise ValueError(f"direction must be 'up' or 'down', got {direction!r}")

    scores, all_s, all_ns = {}, [], []
    for gene in data:
        values = data[gene]
        if direction == "up":
            s, ns = values[values > z_cutoff], values[values <= z_cutoff]
        else:
            s, ns = values[values < -z_cutoff], values[values >= -z_cutoff]
        all_s += s.tolist()
        all_ns += ns.tolist()

        s_bar, len_s = (s.mean(), len(s)) if len(s) else (0, 0)
        ns_bar, len_ns = (ns.mean(), len(ns)) if len(ns) else (0, 0)

        n_sig = len(values)
        weighted_s = s_bar * (1 + len_s / n_sig)
        weighted_ns = ns_bar * (1 + len_ns / n_sig)
        scores[int(gene)] = (
            weighted_s - weighted_ns if direction == "up"
            else abs(weighted_s) - abs(weighted_ns)
        )

    score_cutoff = abs(np.mean(all_s) * 1.5) - abs(np.mean(all_ns) * 1.5)
    ranked = pd.Series(scores).sort_values()
    selected = ranked[ranked > score_cutoff]
    if return_scores:
        return selected.index.tolist(), selected, score_cutoff
    return selected.index.tolist()


def enrich(gene_symbols, gene_set, cutoff=0.05):
    """Enrichr over-representation test on one gene-set library.

    `gseapy.enrichr` changed signature between major versions: it required a
    `description` argument up to 0.x, dropped it in 1.0, and may drop others
    later. The call is therefore built from the signature actually available, so
    the notebooks run on either generation of the library. `description` was
    only a label and never reached the results, so omitting it changes nothing.
    """
    import inspect

    import gseapy

    wanted = {
        "gene_list": gene_symbols,
        "gene_sets": gene_set,
        "cutoff": cutoff,
        "outdir": None,
        "description": "pathway",   # pre-1.0 only
    }
    accepted = inspect.signature(gseapy.enrichr).parameters
    return gseapy.enrichr(**{k: v for k, v in wanted.items() if k in accepted})


def dose_window(doses, low=0.1, high=1.3):
    """Doses within the micromolar window used throughout the study.

    `doses` are LINCS strings such as "0.37 uM".
    """
    return [d for d in doses if low < float(d.split()[0]) < high]


def simulate_null_cosine(n_genes, n_simulations=100_000, seed=2025):
    """Mean and sd of the cosine similarity between pairs of random profiles.

    The spread of cosine similarity between unrelated profiles shrinks as the
    number of genes grows, so a similarity threshold has to be calibrated per
    signature length rather than fixed. Vectorised over simulations.
    """
    rng = np.random.default_rng(seed)
    a = rng.normal(size=(n_simulations, n_genes))
    b = rng.normal(size=(n_simulations, n_genes))
    sims = np.einsum("ij,ij->i", a, b) / (
        np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1)
    )
    return sims.mean(), sims.std()


def flatten(profiles, cutoff=Z_CUTOFF):
    """Set to zero every value whose magnitude is below `cutoff`.

    Keeps the similarity between two profiles driven by the genes at least one
    of them actually moves, instead of by the inert majority.
    """
    flattened = np.array(profiles, dtype=float)
    flattened[np.abs(flattened) < cutoff] = 0
    return flattened
