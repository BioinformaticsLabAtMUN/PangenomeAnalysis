#!/usr/bin/env python3
import matplotlib
matplotlib.use('Agg')  # Use headless backend to avoid Qt errors

import numpy as np
import pandas as pd
import scipy.sparse
import scipy.optimize
from pathlib import Path
import matplotlib.pyplot as plt  # For plotting

def estimate_pan_core_size(gene_matrix, num_iter=100, log_batch=10):
    """
    Computes pan/core genome size curves for many randomizations of genome order.
    """
    num_genes, num_strains = gene_matrix.shape
    print(f"Matrix shape: {gene_matrix.shape}")
    
    gene_data = gene_matrix.tocsr()
    pan_genomes = np.zeros((num_iter, num_strains))
    core_genomes = np.zeros((num_iter, num_strains))
    
    print('Generating pan/core curves from shuffled strains')
    for i in range(num_iter):
        if log_batch > 0 and ((i+1) % log_batch) == 0:
            print(f'\tIteration {i+1} of {num_iter}')
        shuffle_indices = np.arange(num_strains)
        np.random.shuffle(shuffle_indices)
        gene_incidence = np.zeros(num_genes, dtype=int)
        for j, shuffle_col in enumerate(shuffle_indices):
            current_col = gene_data[:, shuffle_col].toarray().flatten()
            gene_incidence += current_col
            pan_genomes[i,j] = np.sum(gene_incidence > 0)
            core_genomes[i,j] = np.sum(gene_incidence == j+1)
            if j % 10 == 0:
                print(f"\t\tProcessed genome {j+1}/{num_strains}")
    
    iter_index = [f'Iter{x}' for x in range(1, num_iter+1)]
    pan_cols = [f'Pan{x}' for x in range(1, num_strains+1)]
    core_cols = [f'Core{x}' for x in range(1, num_strains+1)]
    
    print("Creating output DataFrame...")
    df_pan_core = pd.DataFrame(index=iter_index, 
                               columns=pan_cols + core_cols,
                               data=np.hstack([pan_genomes, core_genomes]))
    
    print(f"Pan-genome sizes range: {np.min(pan_genomes)} - {np.max(pan_genomes)}")
    print(f"Core-genome sizes range: {np.min(core_genomes)} - {np.max(core_genomes)}")
    
    return df_pan_core

def __fit_heaps_single__(df_freqs):
    """Fits a single iteration to Heaps Law: PG size = kappa * (genes)^alpha"""
    heaps = lambda x, alpha, kappa: kappa * np.power(x, alpha)
    n_strains = df_freqs.shape[0]
    p0 = [0.5, float(min(df_freqs.values))]
    
    try:
        popt, _ = scipy.optimize.curve_fit(heaps, 
            np.arange(1, n_strains+1), df_freqs.values, p0=p0)
        return popt
    except Exception as e:
        print(f"Warning: Curve fit failed: {str(e)}")
        return [np.nan, np.nan]

def fit_heaps_law(pan_core_df):
    """Fits Heaps Law to each iteration."""
    print("Starting Heaps law fitting...")
    df = pan_core_df.iloc[:,:int(pan_core_df.shape[1]/2)].T
    n_samples, n_iters = df.shape
    print(f"Fitting {n_iters} iterations...")
    
    heaps_fits = {}
    for i, iter_label in enumerate(df.columns):
        try:
            alpha, kappa = __fit_heaps_single__(df.iloc[:,i])
            heaps_fits[iter_label] = {'alpha': alpha, 'kappa': kappa}
            if i % 10 == 0:
                print(f"\tFitted iteration {i+1}/{n_iters}")
        except Exception as e:
            print(f"Warning: Fitting failed for iteration {iter_label}: {str(e)}")
            continue
    
    fits_df = pd.DataFrame.from_dict(heaps_fits, orient='index')
    print("Calculating summary statistics...")
    alpha_median = fits_df['alpha'].median()
    alpha_ci = np.percentile(fits_df['alpha'].dropna(), [2.5, 97.5])
    kappa_median = fits_df['kappa'].median()
    
    results = {
        'alpha': alpha_median,
        'kappa': kappa_median,
        'ci_alpha_lower': alpha_ci[0],
        'ci_alpha_upper': alpha_ci[1],
        'pangenome_status': 'open' if alpha_median < 1 else 'closed',
        'iteration_fits': fits_df
    }
    
    return results

def analyze_pangenome(matrix_file, labels_file, output_dir, num_iter=100):
    print(f"Loading matrix from {matrix_file}")
    matrix = scipy.sparse.load_npz(matrix_file)
    
    print(f"Loading labels from {labels_file}")
    with open(labels_file, 'r') as f:
        labels = [line.strip() for line in f]
    n_rows = matrix.shape[0]
    row_labels = labels[:n_rows]
    col_labels = labels[n_rows:]
    
    print(f"Matrix dimensions: {matrix.shape}")
    print(f"Number of row labels: {len(row_labels)}")
    print(f"Number of column labels: {len(col_labels)}")
    
    print("Computing pan/core genome sizes...")
    df_pan_core = estimate_pan_core_size(matrix, num_iter)
    
    print("Fitting Heaps law...")
    results = fit_heaps_law(df_pan_core)
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    print("Saving results...")
    results_df = pd.DataFrame({
        'alpha': [results['alpha']],
        'kappa': [results['kappa']],
        'ci_alpha_lower': [results['ci_alpha_lower']],
        'ci_alpha_upper': [results['ci_alpha_upper']],
        'pangenome_status': [results['pangenome_status']]
    })
    results_file = output_path / 'heaps_law_results.csv'
    results_df.to_csv(results_file, index=False)
    fits_file = output_path / 'heaps_law_iterations.csv'
    results['iteration_fits'].to_csv(fits_file)
    data_file = output_path / 'heaps_law_data.npz'
    np.savez(data_file,
             pan_core_curves=df_pan_core.values,
             genome_numbers=np.arange(1, matrix.shape[1] + 1),
             row_labels=row_labels,
             col_labels=col_labels)
    
    # Generate a sample plot
    sample_curve = df_pan_core.iloc[0, :int(matrix.shape[1])].values
    genomes = np.arange(1, matrix.shape[1] + 1)
    p0 = [0.5, float(min(sample_curve))]
    try:
        popt, _ = scipy.optimize.curve_fit(lambda x, alpha, kappa: kappa * x**alpha, genomes, sample_curve, p0=p0)
        fitted_curve = popt[1] * genomes**popt[0]
        plt.figure(figsize=(8,6))
        plt.plot(genomes, sample_curve, 'o', label='Observed pan-genome curve')
        plt.plot(genomes, fitted_curve, '-', label=f'Fitted: alpha={popt[0]:.3f}, kappa={popt[1]:.3f}')
        plt.xlabel("Number of genomes")
        plt.ylabel("Unique gene count")
        plt.legend()
        plt.title("Heaps Law Fit")
        plot_path = output_path / 'heaps_law_plot.png'
        plt.savefig(plot_path)
        plt.close()
        print(f"Saved plot to {plot_path}")
    except Exception as e:
        print(f"Warning: Could not generate plot: {str(e)}")
    
    # Debug: check if expected files exist
    for fname in [results_file, fits_file, data_file, plot_path]:
        if not fname.exists():
            print(f"Warning: {fname} was not created!")
    
    return results

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Perform Heaps Law analysis on pangenome data')
    parser.add_argument('--matrix', required=True, help='Path to presence/absence matrix (.npz)')
    parser.add_argument('--labels', required=True, help='Path to genome labels file')
    parser.add_argument('--output', required=True, help='Output directory')
    parser.add_argument('--iterations', type=int, default=100, help='Number of random iterations')
    
    args = parser.parse_args()
    results = analyze_pangenome(args.matrix, args.labels, args.output, args.iterations)
    print("\nHeaps Law Analysis Results:")
    print(f"Alpha: {results['alpha']:.3f} (CI: {results['ci_alpha_lower']:.3f} - {results['ci_alpha_upper']:.3f})")
    print(f"Kappa: {results['kappa']:.3f}")
    print(f"Pangenome status: {results['pangenome_status']}")

