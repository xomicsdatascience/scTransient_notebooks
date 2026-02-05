# Method Comparison

This directory contains the materials needed to obtain performance comparison of the three methods described in the scTransient paper: scTransient, tradeSeq, and GPFates. It uses an early development version of scTransient ("wavelet_pseudotime"), but the methods can be inspected by opening the .whl file to see that they correspond to the scTransient package.

There are thousands of evaluations being performed, and since tradeSeq and GPFates can take several hours per experiment, the evaluation was done in parallel using an HPC managed by SLURM. The SLURMArrayJobScript.sh file can be submitted to run the experiments. Note that the "--array=0-xxx" parameter should match the number of parameters params.py. The evaluation script (method_execution.py) and parameter file (params.py) are kept outside the container so that they can be modified to investigate different parameter settings.

Results for the execution time of the methods can be found in `results_exectime.zip` and plotted using the `method_exectie_plot.ipynb` notebook.
