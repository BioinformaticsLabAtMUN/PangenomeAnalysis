#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd
import scipy.sparse
import scipy.optimize
import scipy.stats
from statsmodels.stats.stattools import durbin_watson
from scipy.special import betaln
import sys
import traceback
import logging
import matplotlib.pyplot as plt
import sparse_utils   

########################################
# Beta-binomial helper functions  
########################################

def betabin_logpmf(x, n, a, b):
    """
    Beta-binomial log-PMF as in Jason's code.
    """
    k = np.floor(x)
    combiln = -np.log(n + 1) - betaln(n - k + 1, k + 1)
    return combiln + betaln(k + a, n - k + b) - betaln(a, b)

def ecdf_from_counts(vals, counts, limit):
    """
    Computes the empirical CDF from unique values and counts for x = np.arange(limit).
    """
    pmf = np.zeros(limit)
    for i in range(len(vals)):
        pmf[vals[i]] += counts[i]
    return np.cumsum(pmf) / pmf.sum()

def draw_bbn(n, a, b, size, sim_limit=1000):
    """
    Draws random samples from a beta-binomial BBN(n,a,b) over the range(sim_limit).
    """
    Xs = np.arange(sim_limit)
    probs = np.exp(betabin_logpmf(Xs, n, a, b))
    probs /= probs.sum()
    return np.random.choice(Xs, size=size, p=probs)

def ks_montecarlo_bbn(observed_series, n, a, b, iterations=1000, sim_limit=1000):
    """
    Monte-Carlo KS test for a beta-binomial model.
    observed_series is a pandas Series with index=miss counts and values=number of genes.
    """
    Xrange = np.arange(sim_limit)
    model_pmf = np.exp(betabin_logpmf(Xrange, n, a, b))
    model_cdf = np.cumsum(model_pmf)
    
    obs_vals = observed_series.index.to_numpy()
    obs_counts = observed_series.values
    observed_ecdf = ecdf_from_counts(obs_vals, obs_counts, sim_limit)
    ks_stat = np.max(np.abs(observed_ecdf - model_cdf))
    
    total_count = observed_series.sum()
    ks_sim = np.zeros(iterations)
    for i in range(iterations):
        draws = draw_bbn(n, a, b, size=total_count, sim_limit=sim_limit)
        unique, counts = np.unique(draws, return_counts=True)
        sim_ecdf = ecdf_from_counts(unique, counts, sim_limit)
        ks_sim[i] = np.max(np.abs(sim_ecdf - model_cdf))
    ks_pvalue = np.sum(ks_sim > ks_stat) / float(iterations)
    return ks_pvalue, ks_stat, ks_sim

########################################
# Bernoulli grid functions (sparse version)
########################################
'''
def __bernoulli_grid_loglikelihood__(X, P, Q):
    """
    Computes the log-likelihood for the gene presence/absence matrix X,
    given gene true frequencies P and genome capture rates Q.
    
    X is a sparse CSR matrix.
    The likelihood is computed as:
      LL = sum_{i,j} [ x_ij * log(P_i*Q_j) + (1-x_ij)*log(1-P_i*Q_j) ]
    """
    epsilon = 1e-10
    P = np.maximum(P, epsilon)
    Q = np.maximum(Q, epsilon)
    PQ = np.outer(P, Q)
    PQ = np.minimum(PQ, 1 - epsilon)
    
    # Full likelihood over all (i,j): start with log(1-P_i*Q_j)
    ll_all = np.sum(np.log(1 - PQ))
    # For nonzero entries in X, add the correction [log(PQ)-log(1-PQ)]
    X = X.tocsr()
    row, col = X.nonzero()
    ll_nonzero = np.sum(np.log(PQ[row, col]) - np.log(1 - PQ[row, col]))
    return ll_all + ll_nonzero
'''    
    
def __bernoulli_grid_loglikelihood__(X, P, Q):
    """Memory-efficient sparse version"""
    epsilon = 1e-10
    P = np.maximum(P, epsilon)
    Q = np.maximum(Q, epsilon)
    
    # Calculate complete log-likelihood for all zeros efficiently
    # This is mathematically equivalent to np.sum(np.log(1 - np.outer(P,Q)))
    ll_all = 0
    for i in range(len(P)):
        ll_all += np.sum(np.log(1 - P[i] * Q))
    
    # Handle non-zero entries efficiently (only process actual 1s in the matrix)
    X = X.tocsr()
    row, col = X.nonzero()
    pq_vals = P[row] * Q[col]  # Only compute values we need
    ll_nonzero = np.sum(np.log(pq_vals) - np.log(1 - pq_vals))
    
    return ll_all + ll_nonzero
    
def __bernoulli_grid_loglikelihood_gradient__(X, P, Q):
    """
    Memory-optimized gradient function that avoids creating full matrices.
    Computes the gradient of the log-likelihood with respect to P and Q.
    X is a sparse CSR matrix.
    Returns a concatenated vector [dL/dP, dL/dQ].
    """
    m, n = X.shape
    epsilon = 1e-10
    P = np.maximum(P, epsilon)
    Q = np.maximum(Q, epsilon)
    
    # Initialize gradient vectors
    dLdp = np.zeros(m)
    dLdq = np.zeros(n)
    
    # Get row and column sums (efficiently calculated from sparse X)
    X_sum_rows = np.array(X.sum(axis=1)).flatten()  # shape (m,)
    X_sum_cols = np.array(X.sum(axis=0)).flatten()  # shape (n,)
    
    # First term of gradient: sum_j x_ij / p_i for each gene i
    dLdp = X_sum_rows / P
    
    # First term of gradient: sum_i x_ij / q_j for each genome j
    dLdq = X_sum_cols / Q
    
    # Convert X to CSR format for efficient row access
    X_csr = X.tocsr()
    
    # Second term of gradient for P: sum_j (1-x_ij) * q_j / (1-p_i*q_j) for each gene i
    for i in range(m):
        # Find which columns have 1s in this row
        row_start = X_csr.indptr[i]
        row_end = X_csr.indptr[i+1]
        nonzero_cols = X_csr.indices[row_start:row_end]
        
        # Create a mask where we have zeros (complement of nonzero_cols)
        zero_mask = np.ones(n, dtype=bool)
        zero_mask[nonzero_cols] = False
        
        # Only compute for columns where x_ij = 0
        q_values = Q[zero_mask]
        denominator = 1.0 - P[i] * q_values
        dLdp[i] -= np.sum(q_values / denominator)
    
    # Convert X to CSC format for efficient column access
    X_csc = X.tocsc()
    
    # Second term of gradient for Q: sum_i (1-x_ij) * p_i / (1-p_i*q_j) for each genome j
    for j in range(n):
        # Find which rows have 1s in this column
        col_start = X_csc.indptr[j]
        col_end = X_csc.indptr[j+1]
        nonzero_rows = X_csc.indices[col_start:col_end]
        
        # Create a mask where we have zeros (complement of nonzero_rows)
        zero_mask = np.ones(m, dtype=bool)
        zero_mask[nonzero_rows] = False
        
        # Only compute for rows where x_ij = 0
        p_values = P[zero_mask]
        denominator = 1.0 - p_values * Q[j]
        dLdq[j] -= np.sum(p_values / denominator)
    
    return np.concatenate((dLdp, dLdq))

########################################
# Core Genome Analysis: Bernoulli Grid Model
########################################

def compute_bernoulli_grid_core_genome(df_genes_dense, 
                                       prob_bounds=(1e-8, 1-1e-8), 
                                       init_capture_prob=0.99, 
                                       init_gene_freqs=None):
    """
    Models gene presence/absence as Bernoulli random variables and estimates the
    true gene frequencies (p) and genome capture rates (q) by maximizing the 
    log-likelihood.
    
    Only genes with observed frequency > 10% are optimized.
    Core genes are defined as those with estimated frequency > 99.99%.
    """
    print("Starting Bernoulli grid analysis...")
    print(f"Matrix shape: {df_genes_dense.shape}")
    
    # Get the underlying sparse matrix.
    gene_data = df_genes_dense.data
    if gene_data.shape[1] > gene_data.shape[0]:
        gene_data = gene_data.T
    X = gene_data.tocsr()
    n_genes, n_genomes = X.shape
    print(f"Matrix dimensions: {n_genes} genes x {n_genomes} genomes")
    
    # Compute observed gene frequencies using the LightSparseDataFrame's sum method.
    raw_counts = df_genes_dense.sum(axis='index')
    raw_gene_counts = np.array(raw_counts) / float(n_genomes)
    
    # Select genes with observed frequency > 10%
    valid_idx = raw_gene_counts > 0.1  
    print(f"Number of genes used in optimization: {np.sum(valid_idx)} out of {n_genes}")
    
    if init_gene_freqs is None:
        P_guess_valid = raw_gene_counts[valid_idx]
    else:
        P_guess_valid = np.array(init_gene_freqs)[valid_idx]
    
    Q_guess = init_capture_prob * np.ones(n_genomes)
    PQ_guess = np.concatenate((P_guess_valid, Q_guess))
    PQ_guess = np.clip(PQ_guess, prob_bounds[0], prob_bounds[1])
    
    # Compute initial log-likelihood using full raw frequencies and Q_guess.
    init_ll = __bernoulli_grid_loglikelihood__(X, raw_gene_counts, Q_guess)
    print(f"Initial loglikelihood (using raw frequencies): {init_ll}")
    
    print("Running optimization on genes with observed frequency > 10%...")
    try:
        res = scipy.optimize.minimize(
            lambda PQ: -__bernoulli_grid_loglikelihood__(X[valid_idx, :],
                                                         PQ[:np.sum(valid_idx)],
                                                         PQ[np.sum(valid_idx):]),
            PQ_guess,
            jac=lambda PQ: -__bernoulli_grid_loglikelihood_gradient__(X[valid_idx, :],
                                                                      PQ[:np.sum(valid_idx)],
                                                                      PQ[np.sum(valid_idx):]),
            method='L-BFGS-B',
            bounds=[prob_bounds] * len(PQ_guess),
            options={'disp': True, 'maxiter': 1000}
        )
        if not res.success:
            print("Optimization failed for valid genes! Using observed frequencies for them.")
            optimized_P_valid = P_guess_valid
        else:
            optimized_P_valid = res.x[:np.sum(valid_idx)]
    except Exception as e:
        print(f"Optimization error: {str(e)}")
        print("Using observed frequencies for valid genes.")
        optimized_P_valid = P_guess_valid
    
    print(f"Optimization complete. Success: {res.success if 'res' in locals() else False}")
    print(f"Optimized frequencies (for valid genes) shape: {optimized_P_valid.shape}")
    
    gene_freqs_values = np.copy(raw_gene_counts)
    gene_freqs_values[valid_idx] = optimized_P_valid
    gene_freqs = pd.Series(data=gene_freqs_values, index=df_genes_dense.index, name='estimated_frequency')
    
    # Identify core genes (estimated frequency > 99.99%)
    core_threshold = 0.99 
    core_genes = gene_freqs[gene_freqs > core_threshold].index.tolist()
    print(f"Identified {len(core_genes)} core genes with threshold {core_threshold}")
    
    return gene_freqs_values, core_genes, gene_freqs

########################################
# Core Genome Analysis: Beta-binomial Model
########################################

def compute_beta_binomial_core_genome(df_genes, frac_recovered=0.999, num_points=100, ks_iter=1000):
    """
    Estimates the core gene frequency threshold using a beta-binomial error model.
    Also computes several fit quality statistics and prints information about the fitting.
    
    This version restricts the analysis to potential core genes—that is, genes with an 
    observed frequency > 90% (as described in the paper).
    
    Parameters:
      df_genes : LightSparseDataFrame or pd.DataFrame
          Sparse binary gene x genome DataFrame.
      frac_recovered : float
          The fraction of the core genome to be recovered.
      num_points : int
          Number of points to include from the high-frequency (core) end.
      ks_iter : int
          Number of iterations for the Monte Carlo KS test.
          
    Returns:
      A dictionary with the fitted parameters (alpha, beta), cutoff, loglikelihood, and
      quality metrics (MAE, normalized MAE, max-norm MAE, relative MAE, Shapiro-Wilk p-value,
      Durbin-Watson statistic, and Monte Carlo KS-test p-value).
    """
    # Use the underlying sparse matrix if available.
    if hasattr(df_genes, 'data'):
        gene_mat = df_genes.data
    else:
        gene_mat = sparse_utils.sparse_arrays_to_spmatrix(df_genes)
    n_genes, n_genomes = gene_mat.shape

    # Compute raw gene presence counts per gene.
    gene_presence_counts = np.array(gene_mat.sum(axis=1)).flatten()
    # Compute observed frequency (fraction of genomes in which each gene is present).
    observed_freq = gene_presence_counts / float(n_genomes)

    # Restrict to potential core genes: those with observed frequency > 90%
    potential_core = observed_freq > 0.9
    if potential_core.sum() == 0:
        print("No potential core genes found (observed frequency > 90%).")
        return {'alpha': np.nan, 'beta': np.nan, 'cutoff': np.nan, 'loglikelihood': np.nan,
                'mae': np.nan, 'ks_pvalue': np.nan, 'sw_pvalue': np.nan, 'dwstat': np.nan}
    print(f"Potential core genes: {potential_core.sum()} out of {n_genes}")

    # Use only the counts for potential core genes.
    filtered_counts = gene_presence_counts[potential_core]

    # Create a frequency distribution (number of genes with a given frequency)
    # The index will be the raw frequency (e.g., 46, 47, …, up to n_genomes).
    df_counts = pd.Series(filtered_counts).value_counts().sort_index()

    # Reindex to include all possible frequency values (1 to n_genomes).
    df_counts = df_counts.reindex(np.arange(1, n_genomes + 1), fill_value=0)
    # Select the bottom (i.e. highest frequency) num_points rows.
    df_core = df_counts.iloc[-min(num_points, len(df_counts)):]
    # Convert frequencies to "miss counts": a gene present in f genomes is missing in n_genomes - f genomes.
    df_core.index = n_genomes - df_core.index
    # Sort by miss count in descending order.
    df_core = df_core.sort_index(ascending=False)

    if len(df_core.index.unique()) == 1:
        print("Only one unique miss count present; cannot fit beta-binomial model.")
        return {'alpha': np.nan, 'beta': np.nan, 'cutoff': np.nan, 'loglikelihood': np.nan,
                'mae': np.nan, 'ks_pvalue': np.nan, 'sw_pvalue': np.nan, 'dwstat': np.nan}

    # Let X = miss counts, Y = counts.
    X = df_core.index.values
    Y = df_core.values

    def loglikelihood(ab):
        a, b = ab
        if a <= 0 or b <= 0:
            return -1e12
        return np.dot(Y, betabin_logpmf(X, n_genomes, a, b))

    print("Starting beta-binomial optimization with initial guess (1, 100) on potential core genes...")
    res = scipy.optimize.minimize(
        lambda ab: -loglikelihood(ab),
        x0=(1, 100),
        method='L-BFGS-B',
        bounds=[(1e-8, None), (1e-8, None)],
        options={'maxiter': ks_iter}
    )

    if res is None or not res.success:
        print("Beta-binomial optimization failed!")
        return {'alpha': np.nan, 'beta': np.nan, 'cutoff': np.nan, 'loglikelihood': np.nan,
                'mae': np.nan, 'ks_pvalue': np.nan, 'sw_pvalue': np.nan, 'dwstat': np.nan}
    else:
        a, b = res.x
        print(f"Beta-binomial optimization succeeded: alpha = {a:.4f}, beta = {b:.4f}")

        cutoff = 0
        cdf = np.exp(betabin_logpmf(cutoff, n_genomes, a, b))
        while cdf < frac_recovered and cutoff <= n_genomes:
            cutoff += 1
            cdf += np.exp(betabin_logpmf(cutoff, n_genomes, a, b))
        cutoff = max(0, cutoff - 1)
        loglike = -res.fun

        total_genes = Y.sum()
        model_prob = np.exp(betabin_logpmf(X, n_genomes, a, b))
        Yhat = total_genes * model_prob
        residuals = Y - Yhat
        mae = np.mean(np.abs(residuals))
        mae_normalized = mae / len(residuals)
        mae_max = mae / np.max(Y) if np.max(Y) > 0 else np.nan
        mae_relative = np.mean(np.abs(residuals) / (Y + 1e-10))

        try:
            stat, sw_pvalue = scipy.stats.shapiro(residuals)
        except Exception as e:
            sw_pvalue = np.nan
            print("Shapiro-Wilk test failed:", e)
        # Use the globally imported durbin_watson function.
        dwstat = durbin_watson(residuals)

        print(f"Beta-binomial fitting results:")
        print(f"  Cutoff (max misses allowed): {cutoff}")
        print(f"  Log-likelihood: {loglike:.4f}")
        print(f"  MAE: {mae:.4f} (normalized: {mae_normalized:.4f}, max norm: {mae_max:.4f}, relative: {mae_relative:.4f})")
        print(f"  Shapiro-Wilk p-value: {sw_pvalue:.4f}")
        print(f"  Durbin-Watson statistic: {dwstat:.4f}")

        observed_series = df_core  # already indexed by miss count, with counts as values
        Xrange = np.arange(n_genomes + 1)
        model_pmf_full = np.exp(betabin_logpmf(Xrange, n_genomes, a, b))
        model_cdf_full = np.cumsum(model_pmf_full)
        err = 1 - model_cdf_full
        sim_limit_candidates = np.where(err < 1e-8)[0]
        sim_limit = sim_limit_candidates[0] if len(sim_limit_candidates) > 0 else n_genomes + 1
        # Ensure sim_limit covers the maximum observed miss count.
        if sim_limit < (np.max(observed_series.index) + 1):
            sim_limit = int(np.max(observed_series.index) + 1)
        if sim_limit > 0:
            ks_pvalue, ks_stat, ks_sim = ks_montecarlo_bbn(observed_series, n_genomes, a, b,
                                                             iterations=ks_iter, sim_limit=sim_limit)
        else:
            ks_pvalue = np.nan

        print(f"  Monte Carlo KS-test p-value: {ks_pvalue:.4f}")

        return {'alpha': a, 'beta': b, 'cutoff': cutoff, 'loglikelihood': loglike,
                'mae': mae, 'mae_normalized': mae_normalized, 'mae_max': mae_max, 'mae_relative': mae_relative,
                'ks_pvalue': ks_pvalue, 'sw_pvalue': sw_pvalue, 'dwstat': dwstat}






########################################
# Visualization
########################################

def visualize_beta_binomial_fit(df_genes, a, b, output_file):
    """
    Plots the observed distribution of miss counts for potential core genes (observed frequency > 90%)
    along with the fitted beta-binomial probability mass function (PMF).
    
    Parameters:
      df_genes : LightSparseDataFrame or pd.DataFrame
          The gene presence/absence matrix.
      a, b : float
          Fitted beta-binomial parameters.
      output_file : str
          The file name to save the plot.
    """
    # Use the underlying sparse matrix if available.
    if hasattr(df_genes, 'data'):
        gene_mat = df_genes.data
    else:
        gene_mat = sparse_utils.sparse_arrays_to_spmatrix(df_genes)
        
    n_genes, n_genomes = gene_mat.shape
    # Compute gene presence counts.
    gene_presence_counts = np.array(gene_mat.sum(axis=1)).flatten()
    observed_freq = gene_presence_counts / float(n_genomes)
    
    # Filter to potential core genes (observed frequency > 90%).
    potential_core = observed_freq > 0.9
    if potential_core.sum() == 0:
        print("No potential core genes (frequency > 90%) for plotting.")
        return
    filtered_counts = gene_presence_counts[potential_core]
    
    # Build a frequency distribution for potential core genes.
    # Keys: gene count (number of genomes in which gene is present)
    # Values: number of genes with that count.
    df_counts = pd.Series(filtered_counts).value_counts().sort_index()
    # Ensure that we have counts for every possible value from 1 to n_genomes.
    df_counts = df_counts.reindex(np.arange(1, n_genomes+1), fill_value=0)
    
    # Because potential core genes should have high counts, we take only the bottom portion:
    # For example, select the bottom num_points (largest counts).
    num_points = min(100, len(df_counts))
    df_core = df_counts.iloc[-num_points:]
    
    # Convert counts to miss counts: a gene present in f genomes is missing in n_genomes - f genomes.
    # (For example, if n_genomes=50 and f=48 then miss count = 2.)
    # Set the index to be the miss counts.
    df_core.index = n_genomes - df_core.index
    # For plotting, sort in ascending order so that miss=0 is at left.
    df_core = df_core.sort_index()
    
    # Observed miss counts and frequencies.
    x_obs = df_core.index.values.astype(int)
    y_obs = df_core.values.astype(float)
    # Normalize observed frequencies so that the total area equals 1.
    y_obs_norm = y_obs / y_obs.sum()
    
    # Define an x-range for the fitted PMF. We use the range 0 to max observed miss count.
    x_range = np.arange(0, np.max(x_obs)+1)
    # Compute the fitted PMF at each miss count.
    fitted_pmf = np.exp(betabin_logpmf(x_range, n_genomes, a, b))
    fitted_pmf /= fitted_pmf.sum()  # normalize to sum to 1

    # Plot the observed distribution and fitted PMF.
    import matplotlib.pyplot as plt
    plt.figure(figsize=(8,6))
    plt.bar(x_obs, y_obs_norm, width=0.8, alpha=0.6, color='skyblue', edgecolor='black',
            label='Observed miss distribution')
    plt.plot(x_range, fitted_pmf, 'ro-', linewidth=2, markersize=6,
             label=f'Fitted beta-binomial\n(a={a:.2f}, b={b:.2f})')
    plt.xlabel("Number of misses (genomes missing the gene)")
    plt.ylabel("Normalized frequency")
    plt.title("Beta-binomial fit to potential core gene miss counts (observed frequency > 90%)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()



########################################
# Main
########################################

def main():
    try:
        parser = argparse.ArgumentParser(description='Compute core genome using statistical models')
        parser.add_argument('--matrix', required=True, help='Path to gene presence/absence matrix (NPZ)')
        parser.add_argument('--labels', required=True, help='Path to labels file')
        parser.add_argument('--output-prefix', required=True, help='Prefix for output files')
        args = parser.parse_args()
        
        logging.basicConfig(filename=f"{args.output_prefix}.log", level=logging.DEBUG,
                            format='%(asctime)s - %(levelname)s - %(message)s')
        print("\n=== Starting Core Genome Analysis ===")
        
        print('\nLoading gene presence/absence matrix...')
        df_genes = sparse_utils.read_lsdf(args.matrix, args.labels)
        print(f"Loaded matrix with shape: {df_genes.shape}")
        print(f"Number of index labels: {len(df_genes.index)}")
        print(f"Number of column labels: {len(df_genes.columns)}")
        
        print('\nEstimating gene frequencies (Bernoulli grid model)...')
        opt_values, core_genes, gene_freqs = compute_bernoulli_grid_core_genome(df_genes)
        print(f"Number of core genes identified: {len(core_genes)}")
        
        # Save core gene list and frequency estimates.
        with open(f"{args.output_prefix}_core_genes.txt", 'w') as f:
            for gene in core_genes:
                f.write(f"{gene}\n")
        gene_freqs.to_frame().to_csv(f"{args.output_prefix}_frequency_estimates.csv")
        
        print('\nFitting beta-binomial error model...')
        bb_results = compute_beta_binomial_core_genome(df_genes)
        pd.DataFrame([bb_results]).to_csv(f"{args.output_prefix}_beta_binomial_results.csv", index=False)
        
        if np.isnan(bb_results.get('alpha', np.nan)) or np.isnan(bb_results.get('beta', np.nan)):
            print("Beta-binomial fitting failed; skipping visualization.")
        else:
            visualize_beta_binomial_fit(df_genes, bb_results['alpha'], bb_results['beta'],
                                        f"{args.output_prefix}_beta_binomial_fit.png")
        
        print("Saving core genome matrix as NPZ...")
        try:
            core_df = df_genes.labelslice(indices=core_genes)
            core_df.to_npz(f"{args.output_prefix}_core_matrix.npz")
        except Exception as e:
            print(f"Error saving core genome matrix: {str(e)}", file=sys.stderr)
        
        print("\nCore genome analysis complete!")
    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        traceback.print_exc()
        sys.exit(1)

if __name__ == '__main__':
    main()
