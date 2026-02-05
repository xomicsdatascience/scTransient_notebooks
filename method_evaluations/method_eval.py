from fontTools.t1Lib import font_dictionary_keys

import wavelet_pseudotime.synthetic
from time import time
from warnings import warn
import pandas as pd
import numpy as np
import GPfates.GPfates
from contextlib import redirect_stdout
import io
import anndata as ad
from typing import Tuple
import os
import pickle as pkl
from typing import List, Tuple, Optional, Collection, Literal
from matplotlib import pyplot as plt

import rpy2.robjects as robjects
from rpy2.robjects import r as R
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import default_converter
from rpy2.robjects import numpy2ri


def compute_precision(score_list: list,
                      n_signals: int):
    """Computes how many of the top-scoring entries should be top-scoring. Assumes that all signals
    are at the start of the list (i.e., score_list[0:n_signals] are all actual signals.)"""
    # Get pos of highest scores
    scores = np.array(score_list)
    score_buffer = np.array(score_list)

    score_buffer[np.isnan(score_buffer)] = -np.inf
    best_idx = np.argsort(score_buffer)[-n_signals:]
    # If NaN is in the first n_signals, remove it from the count
    num_nan = np.sum(np.isnan(scores[best_idx[best_idx < n_signals]]))

    numerator = np.sum(best_idx < n_signals) - num_nan
    return numerator / n_signals


def evaluate_synth_gpfates(synth_params: dict,
                           data_dist: Literal["gaussian", "nb"] = "nb") -> Tuple[float, float, dict]:
    """
    Evaluate GPfates using synthetic data producd by the wavelet_pseudotime.synthetic module.
    Parameters
    ----------
    synth_params : dict
        Dict containing parameters for generating synthetic data.

    Returns
    -------
    score : float
    time : float
    exp_info : dict
    """
    if "seed" not in synth_params:
        warn("No random seed has been specified; results won't be reproducible.")

    if data_dist == "nb":
        synth_data, pseudotime, _ = wavelet_pseudotime.synthetic.generate_synthetic_data_nb_spikes(**synth_params)
    elif data_dist == "gaussian":
        synth_data, pseudotime, _ = wavelet_pseudotime.synthetic.generate_synthetic_data(**synth_params)
    else:
        raise ValueError(f"Invalid data distribution: {data_dist}. Supported distributions are 'nb' and 'gaussian'.")

    # Convert data to format needed by GPfates
    data_df = pd.DataFrame(synth_data, index=[f"cell_{i}" for i in range(synth_data.shape[0])],
                           columns=[f"gene_{i}" for i in range(synth_data.shape[1])])
    metadata_np = np.random.randn(synth_data.shape[0], 1)
    metadata = pd.DataFrame(metadata_np, index=[f"cell_{i}" for i in range(metadata_np.shape[0])],
                            columns=["pseudotime"])
    metadata["pseudotime"] = pseudotime

    t_start = time()
    # Run GPfates
    m = GPfates.GPfates.GPfates(sample_info=metadata,
                                expression_matrix=data_df.T)
    m.dimensionality_reduction()
    m.store_dr()

    # Lots of output to stdout; suppress it.
    with io.StringIO() as buf, redirect_stdout(buf):
        m.model_fates(maxiter=1000)
    stats = m.calculate_bifurcation_statistics()

    # Evaluate results
    t_end = time()
    dstats = stats["D"].values
    n_signal_genes = synth_params["n_signal_genes"]
    # "How many of the genes with the top D-stat are signal genes?"
    score = compute_precision(dstats, n_signals=n_signal_genes)

    results = synth_params.copy()
    del results["pseudotime_density"]
    results["method"] = "GPfates"
    results["score"] = score
    results["time"] = t_end - t_start
    return score, t_end - t_start, results


def evaluate_synth_sctransient(synth_params: dict,
                               data_dist: Literal["gaussian", "nb"] = "nb") -> Tuple[float, float, dict]:
    if "seed" not in synth_params:
        warn("No random seed has been specified; results won't be reproducible.")

    if data_dist == "nb":
        synth_data, pseudotime, _ = wavelet_pseudotime.synthetic.generate_synthetic_data_nb_spikes(**synth_params)
    elif data_dist == "gaussian":
        synth_data, pseudotime, _ = wavelet_pseudotime.synthetic.generate_synthetic_data(**synth_params)

    t_start = time()
    # De-mean data
    synth_data -= np.mean(synth_data, axis=0)
    adata = ad.AnnData(synth_data)
    adata.obs["psupertime"] = pseudotime
    adata.obs["phase"] = 1  # for assumptions

    waves, scores, pt_sigs_dict, adata2 = wavelet_pseudotime.process.pipeline3(adata,
                                                                               scoring_threshold=0.2)
    t_end = time()
    score_list = []
    n_genes = synth_params["n_genes"]
    n_signal_genes = synth_params["n_signal_genes"]
    for i in range(n_genes):
        score_list.append(scores[str(i)])
    s = compute_precision(score_list, n_signals=n_signal_genes)

    results = synth_params.copy()
    del results["pseudotime_density"]
    results["method"] = "scTransient"
    results["score"] = s
    results["time"] = t_end - t_start

    return s, t_end - t_start, results


def evaluate_synth_tradeseq(synth_params: dict,
                            tradeSeq_R_params: dict,
                            data_dist: Literal["gaussian", "nb"] = "nb") -> Tuple[float, float, dict]:
    if "seed" not in synth_params:
        warn("No random seed has been specified; results won't be reproducible.")

    if data_dist == "nb":
        synth_data, pseudotime, _ = wavelet_pseudotime.synthetic.generate_synthetic_data_nb_spikes(**synth_params)
    elif data_dist == "gaussian":
        synth_data, pseudotime, _ = wavelet_pseudotime.synthetic.generate_synthetic_data(**synth_params)
        synth_data -= np.min(synth_data, axis=0)

    pseudotime = np.expand_dims(pseudotime, 1)
    counts = synth_data.T
    cell_weights = np.ones(pseudotime.shape)

    # Set up R environment
    R("library(tradeSeq)")
    R("library(RColorBrewer)")
    R("library(SingleCellExperiment)")
    R("library(slingshot)")
    R('options(device=NULL)')
    R('library(R.devices)')
    # Export to R
    with localconverter(default_converter + numpy2ri.converter):
        robjects.globalenv["counts"] = counts
        robjects.globalenv["pseudotime"] = pseudotime
        robjects.globalenv["cell_weights"] = cell_weights
        robjects.globalenv["knots"] = tradeSeq_R_params["knots"]
        R(f"knots <- {tradeSeq_R_params['knots']}")
        t_start = time()

        R(f"icMat <- suppressGraphics(evaluateK(counts=counts, pseudotime=pseudotime, cellWeights=cell_weights, k=knots))")
        R(f"nknots <- knots[order(colMeans(icMat))[1]]")
        R(f"sce <- fitGAM(counts=counts, pseudotime=pseudotime, cellWeights=cell_weights, nknots=nknots)")
        R(f"pthresh <- 0.2")  # not really needed
        R(f"assoc <- associationTest(sce, inverse='eigen')")
        R(f"assocp <- assoc$pvalue")
        t_end = time()
        pvals = np.array(robjects.globalenv["assocp"])
        # pvals[np.isnan(pvals)] = 1.0
        assoc = robjects.globalenv["assoc"]
        score = compute_precision(-pvals, n_signals=synth_params["n_signal_genes"])
        # R(f"R_score <- (sum(assoc_argmin[1:100] <= 100) - sum( is.na(assoc$pvalue[assoc_argmin[1:100]]))) / {synth_data.shape[1]}")
    results = synth_params.copy()
    del results["pseudotime_density"]
    results["method"] = "tradeSeq"
    results["score"] = score
    results["time"] = t_end - t_start

    return score, t_end - t_start, results


def gather_results(results_dir) -> list:
    """Load pickled results from a directory and gather the contents into a list."""
    results = []
    for file in os.listdir(results_dir):
        if file.endswith(".pkl"):
            with open(os.path.join(results_dir, file), "rb") as f:
                results.append(pkl.load(f))
    return results


def extract_grouped_values(data: list,
                           independent_key: str,
                           dependent_key: str,
                           filter: Optional[dict] = None):
    """
    Extracts and groups values from a list of dictionaries based on given conditions and keys.
    Parameters
    ----------
    data : List[dict]
        List of dicts to parse.
    independent_key : str
        Key to use for the independent variable.
    dependent_key : str
        Key to use for the dependent variable.
    filter : dict, optional
        Restrict the data to only those matching the given conditions.

    Returns
    -------

    """

    # Dictionary to store independent value as key and list of dependent values as value
    grouped_data = {}
    if filter is None:
        filter = {}
    # Iterate over each record in the data list
    for record in data:
        skip_record = False
        for filt_key, filt_value in filter.items():
            record_value = record.get(filt_key)
            if isinstance(record_value, list) and len(record_value) >= 1:
                record_value = record_value[0]
            if record_value != filt_value:
                skip_record = True
                break  # skip if record doesn't match the requested filter value
        if skip_record:
            continue
        indep_val = record.get(independent_key)
        if isinstance(indep_val, Collection) and len(indep_val) >= 1:
            indep_val = indep_val[0]
        dep_val = record.get(dependent_key)
        if indep_val is not None:
            # Append the dependent value to the list for the corresponding independent value
            if indep_val in grouped_data:
                grouped_data[indep_val].append(dep_val)
            else:
                grouped_data[indep_val] = [dep_val]

    # Extract the unique independent values and the corresponding dependent values lists.
    # Sorting the keys is optional; if order matters, you might choose to preserve the order
    # in which they were first encountered.
    independent_values = list(grouped_data.keys())
    dependent_values = [grouped_data[key] for key in independent_values]

    return independent_values, dependent_values


def determine_unique_values(data: List[dict],
                            key: str,
                            filter: Optional[dict] = None) -> List:
    values = set()
    for dat in data:
        if check_retain_data(dat, filter):
            value = dat.get(key)
            if isinstance(value, Collection) and len(value) == 1:
                values.add(value[0])
            else:
                values.add(value)
    return sorted(list(values))


def check_retain_data(data: dict,
                      filter: dict = None) -> bool:
    """Checks whether a record should be retained. For entries in `filter` that are list, checks that the matching
    entry is in the list."""
    if filter is None:
        return True
    for k, v in filter.items():
        if isinstance(v, list) or isinstance(v, tuple) or isinstance(v, set):
            if data.get(k) not in v:
                return False
        else:
            if data.get(k) != v:
                return False
    return True


def plot_results(results_dir: str,
                 x_axis_key: str,
                 y_axis_key: str,
                 horizontal_plot_split_key: str,
                 vertical_plot_split_key: str,
                 multiplot_key: str,
                 multiplot_value_order: Optional[list] = None,
                 horizontal_split_values: Optional[list] = None,
                 vertical_split_values: Optional[list] = None,
                 axis_horizontal_label: Optional[str] = None,
                 axis_vertical_label: Optional[str] = None,
                 ylim: Optional[tuple] = None,
                 title_numerical_precision: int = 2,
                 filter: Optional[dict] = None,
                 save_path: Optional[str] = None,
                 figheight: int = 10,
                 figwidth: int = 20,
                 suplabels: dict = None):
    """
    Produces plots from the results. The figure is divided horizontally by the horizontal_plot_split_key (so that each
    column is a separate value for that key), and similarly vertically. The `_split_values` allow for a subset to be
    selected by value; if None, all values are used.
    In each plot, the x_axis_key and y_axis_key determine what is plotted.

    Parameters
    ----------
    results_dir
    x_axis_key
    y_axis_key
    horizontal_plot_split_key
    vertical_plot_split_key
    multiplot_key: str
    horizontal_split_values
    vertical_split_values
    ylim : tuple, optional
    title_numerical_precision : int, optional
    filter

    Returns
    -------

    """
    results = gather_results(results_dir)
    if horizontal_split_values is None:
        horizontal_split_values = determine_unique_values(results, horizontal_plot_split_key, filter)
    if vertical_split_values is None:
        vertical_split_values = determine_unique_values(results, vertical_plot_split_key, filter)
    multiplot_values = determine_unique_values(results, multiplot_key, filter)
    if multiplot_value_order is not None:
        multiplot_values = reorder_to_reference(multiplot_values, multiplot_value_order)

    local_filter = filter.copy() if filter is not None else {}
    fig, axs = plt.subplots(len(vertical_split_values), len(horizontal_split_values), constrained_layout=True)
    for h_idx, hv in enumerate(horizontal_split_values):
        local_filter[horizontal_plot_split_key] = hv
        for v_idx, vv in enumerate(vertical_split_values):
            # ax = axs[h_idx, v_idx]
            ax = axs[v_idx, h_idx]
            local_filter[vertical_plot_split_key] = vv
            for m_idx, mv in enumerate(multiplot_values):
                local_filter[multiplot_key] = mv

                x_values, y_values = extract_grouped_values(results, x_axis_key, y_axis_key, local_filter)
                y_means = [np.mean(y) for y in y_values]
                y_std = [np.std(y) for y in y_values]
                y_stderr = [np.std(y) / np.sqrt(len(y)) for y in y_values]
                x_arg_idx = np.argsort(x_values)
                x_values = np.array(x_values)[x_arg_idx]
                y_means = np.array(y_means)[x_arg_idx]
                y_std = np.array(y_std)[x_arg_idx]
                y_stderr = np.array(y_stderr)[x_arg_idx]
                # ax.errorbar(x_values, y_means, yerr=y_std, label=mv, marker="o")
                ax.errorbar(x_values, y_means, yerr=y_stderr, label=mv, marker="o")
                if ylim is not None:
                    ax.set_ylim(ylim)

                if vv is int or np.floor(vv) == vv:
                    vv_str = f"{vv}"
                else:
                    vv_str = f"{vv:.{title_numerical_precision}f}"
                if axis_vertical_label is not None:
                    vv_str = f"{axis_vertical_label}={vv_str}"
                else:
                    vv_str = f"{vertical_plot_split_key}={vv_str}"

                if hv is int or np.floor(hv) == hv:
                    hv_str = f"{hv}"
                else:
                    hv_str = f"{hv:.{title_numerical_precision}f}"
                if axis_horizontal_label is not None:
                    hv_str = f"{axis_horizontal_label}={hv_str}"
                else:
                    hv_str = f"{horizontal_plot_split_key}={hv_str}"
                title_str = f"{hv_str}, {vv_str}"
                ax.set_title(title_str)
                # Remove x-axis ticks and ticklabels unless in the last row
                if v_idx != len(vertical_split_values)-1:
                    ax.set_xticks([])
                    ax.set_xticklabels([])
                # Same for y-axis
                if h_idx != 0:
                    ax.set_yticks([])
                    ax.set_yticklabels([])
    ax.legend(fontsize=12)  # apply only to last axis
    apply_suplabels(fig, suplabels)
    # fig.tight_layout()
    fig.set_figheight(figheight)
    fig.set_figwidth(figwidth)
    if save_path is not None:
        fig.savefig(save_path)
        with open(save_path + ".pkl", "wb") as f:
            pkl.dump(fig, f)
    else:
        return fig


def reorder_to_reference(to_reorder: list,
                         reference: list) -> list:
    """Reorders a list to match the order of a reference list; entries not in the reference are placed at the end."""
    to_reorder_set = set(to_reorder)
    reference_set = set(reference)
    to_reorder_diff = to_reorder_set.difference(reference_set)
    new_order = []

    for ref in reference:
        if ref in to_reorder_set:
            new_order.append(ref)
    for reorder in to_reorder:
        if reorder not in reference_set:
            new_order.append(reorder)
    return new_order


def apply_suplabels(fig,
                    suplabels: dict):
    if suplabels is None:
        return
    keys = suplabels.keys()
    if "fontsize" in keys:
        fontsize = suplabels["fontsize"]
    else:
        fontsize = 18
    if "supxlabel" in keys:
        fig.supxlabel(suplabels["supxlabel"], fontsize=fontsize)
    if "supylabel" in keys:
        fig.supylabel(suplabels["supylabel"], fontsize=fontsize)
    if "suptitle" in keys:
        fig.suptitle(suplabels["suptitle"], fontsize=fontsize)
    return