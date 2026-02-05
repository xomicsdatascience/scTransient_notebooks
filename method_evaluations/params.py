import numpy as np
import wavelet_pseudotime.synthetic

data_dist = "nb"
n_cells = [20, 50, 100, 200, 500]
num_repeats = 20
methods = ["sctransient", "tradeseq", "gpfates"]
#methods = ["tradeseq"]

# Data parameters
n_genes = 2000
spike_widths = [[0.1], [0.05], [0.01], [0.005], [0.001]]
pdens = wavelet_pseudotime.synthetic.pseudotime_density("uniform")
n_signal_genes = 100
spike_times = [0.5]

autocorr = 0.5
correct_genes = [str(i) for i in range(n_signal_genes)]
all_genes = [str(i) for i in range(n_genes)]

if data_dist == "nb":
    data_nb_r = 10
    data_nb_p = 0.5
    amplitude_values = [[i] for i in np.linspace(0.01, 80, 10)]
elif data_dist == "gaussian":
    noise_gene_mean = 0
    noise_gene_sd = 0.1
    amplitude_values = [[i] for i in np.linspace(0.01, 6, 10)]
