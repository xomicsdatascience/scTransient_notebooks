import argparse
from method_eval import evaluate_synth_tradeseq, evaluate_synth_gpfates, evaluate_synth_sctransient
import numpy as np
import wavelet_pseudotime.synthetic
import pickle
import os

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Execute method evaluation with specific index')
    parser.add_argument('--index', type=int, required=True, help='Index for parameter list.')
    parser.add_argument('--seed_offset', type=int, default=0, help='Random seed offset.')
    parser.add_argument('--count', action="store_true", help='Dry run; counts the number of possible jobs.')
    parser.add_argument("--output_dir", type=str, default="results", help='Output directory.')
    # parser.add_argument('--method', type=str, default="sctransient",
    #                     choices=["sctransient", "tradeseq", "gpfates"], help='Method to evaluate.')
    parser.add_argument("--parameter_file", type=str, required=False, default=None,
                        help="File containing list of parameters.")
    args = parser.parse_args()

    if args.parameter_file is None:
        # Default parameters
        n_cells = [20, 50, 100, 200, 500]
        amplitude_values = [[i] for i in np.linspace(0.01, 80, 10)]
        num_repeats = 20
        methods = ["sctransient", "tradeseq", "gpfates"]

        # Data parameters
        n_genes = 2000
        spike_widths = [[0.1], [0.05], [0.01], [0.005], [0.001]]
        pdens = wavelet_pseudotime.synthetic.pseudotime_density("uniform")
        n_signal_genes = 100
        spike_times = [0.5]
        data_nb_r = 10
        data_nb_p = 0.5
        autocorr = 0.5
        correct_genes = [str(i) for i in range(n_signal_genes)]
        all_genes = [str(i) for i in range(n_genes)]
        data_dist = "nb"
    else:
        import importlib.util

        spec = importlib.util.spec_from_file_location("parameters", args.parameter_file)
        parameters = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(parameters)

        n_cells = parameters.n_cells
        amplitude_values = parameters.amplitude_values
        num_repeats = parameters.num_repeats
        methods = parameters.methods
        n_genes = parameters.n_genes
        spike_widths = parameters.spike_widths
        spike_times = parameters.spike_times
        pdens = parameters.pdens
        n_signal_genes = parameters.n_signal_genes
        data_dist = parameters.data_dist
        if data_dist == "nb":
            data_nb_r = parameters.data_nb_r
            data_nb_p = parameters.data_nb_p
        elif data_dist == "gaussian":
            noise_gene_mean = parameters.noise_gene_mean
            noise_gene_sd = parameters.noise_gene_sd

        autocorr = parameters.autocorr
        correct_genes = parameters.correct_genes
        all_genes = parameters.all_genes


    idx_tuple = np.unravel_index(args.index, (len(n_cells), len(amplitude_values), len(spike_widths), len(methods), num_repeats))
    if args.count:
        print(f"Total jobs: {len(n_cells) * len(amplitude_values) * len(spike_widths) * len(methods) * num_repeats}")
        print(f"Respective indices: {idx_tuple}")
        exit(0)
    n_cells = n_cells[idx_tuple[0]]
    amplitude_values = amplitude_values[idx_tuple[1]]
    spike_widths = spike_widths[idx_tuple[2]]
    method = methods[idx_tuple[3]]
    seed = args.index  # so the different methods see the same data

    if data_dist == "nb":
        params = {"n_cells": n_cells,
                  "n_genes": n_genes,
                  "pseudotime_density": pdens,
                  "n_signal_genes": n_signal_genes,
                  "spike_amplitudes": amplitude_values,
                  "spike_widths": spike_widths,
                  "spike_times": spike_times,
                  "data_nb_r": data_nb_r,
                  "data_nb_p": data_nb_p,
                  "autocorr": autocorr,
                  "seed": seed + args.seed_offset}
    elif data_dist == "gaussian":
        params = {"n_cells": n_cells,
                  "n_genes": n_genes,
                  "pseudotime_density": pdens,
                  "n_signal_genes": n_signal_genes,
                  "spike_amplitudes": amplitude_values,
                  "spike_widths": spike_widths,
                  "spike_times": spike_times,
                  "noise_gene_mean": noise_gene_mean,
                  "noise_gene_sd": noise_gene_sd,
                  "seed": seed
        }

    if method.lower() == "sctransient":
        _, _, details = evaluate_synth_sctransient(params, data_dist=data_dist)
    elif method.lower() == "tradeseq":
        _, _, details = evaluate_synth_tradeseq(params, {"knots": "3:14"}, data_dist=data_dist)
    elif method.lower() == "gpfates":
        _, _, details = evaluate_synth_gpfates(params, data_dist=data_dist)
    else:
        raise ValueError(f"Unrecognized method: {method}")

    os.makedirs(args.output_dir, exist_ok=True)
    output_file = f"{args.output_dir}/{args.index}_{method}_{n_cells}_{amplitude_values[0]}_{spike_widths[0]}.pkl"
    with open(output_file, "wb") as f:
        pickle.dump(details, f)

    # with open(f"{args.output_dir}/{args.index}_{method}_{n_cells}_{amplitude_values}_{spike_widths}.txt", "wb") as f:
    #     f.write(details)
