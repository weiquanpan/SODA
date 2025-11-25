# -*- coding: utf-8 -*-
"""
SODA: Sinkhorn Optimal Transport with Divergence-Based Adaptation
for Well-Calibrated Single-Cell Integration

This implementation provides SODA batch correction along with benchmarking
against standard methods (Harmony, ComBat, BBKNN, scVI).
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import cdist
from scipy import sparse
from sklearn.metrics import adjusted_rand_score
from sklearn.cluster import KMeans
import warnings
warnings.filterwarnings('ignore')

try:
    import ot
    HAS_OT = True
except ImportError:
    HAS_OT = False
    print("Warning: POT library not installed. SODA will be skipped.")
    print("Install with: pip install POT")

# Global settings
plt.rcParams['figure.dpi'] = 100
plt.rcParams['savefig.dpi'] = 300
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")
sc.settings.verbosity = 0


class SODABatchCorrector:
    """
    SODA: Sinkhorn Optimal Transport with Divergence-Based Adaptation
    
    SODA adopts Harmony's iterative clustering paradigm but replaces linear
    centroid-based correction with entropy-regularized optimal transport.
    It uses Sinkhorn divergence SD(P‖Q) = W²₂(P,Q) − ½W²₂(P,P) − ½W²₂(Q,Q)
    for proper batch effect quantification and adaptive cluster-wise correction.
    """
    
    def __init__(self, n_clusters=50, max_iter=10, epsilon=0.05, 
                 divergence_threshold=0.05):
        """
        Initialize SODA batch corrector.
        
        Parameters:
        -----------
        n_clusters : int
            Number of clusters for iterative clustering (Harmony-style)
        max_iter : int
            Maximum iterations for convergence
        epsilon : float
            Entropic regularization parameter for Sinkhorn algorithm
        divergence_threshold : float
            Sinkhorn divergence threshold for adaptive correction
        """
        self.n_clusters = n_clusters
        self.max_iter = max_iter
        self.epsilon = epsilon
        self.divergence_threshold = divergence_threshold
        self.correction_history = []
    
    def compute_sinkhorn_divergence(self, X_i, X_j):
        """
        Compute Sinkhorn divergence: SD(P‖Q) = W²₂(P,Q) − ½W²₂(P,P) − ½W²₂(Q,Q)
        
        This provides proper batch effect quantification by accounting for
        within-batch variation.
        """
        try:
            # Cross-batch Wasserstein distance W²₂(P,Q)
            M_ij = ot.dist(X_i, X_j, metric='sqeuclidean')
            a = np.ones(len(X_i)) / len(X_i)
            b = np.ones(len(X_j)) / len(X_j)
            W2_ij = ot.sinkhorn2(a, b, M_ij, self.epsilon)
            
            # Within-batch distances W²₂(P,P) and W²₂(Q,Q)
            M_ii = ot.dist(X_i, X_i, metric='sqeuclidean')
            M_jj = ot.dist(X_j, X_j, metric='sqeuclidean')
            W2_ii = ot.sinkhorn2(a, a, M_ii, self.epsilon)
            W2_jj = ot.sinkhorn2(b, b, M_jj, self.epsilon)
            
            # Sinkhorn divergence
            sinkhorn_div = W2_ij - 0.5 * W2_ii - 0.5 * W2_jj
            return sinkhorn_div
        except:
            return 0.0
    
    def detect_global_batch_effect(self, adata, batch_key='batch', use_rep='X_pca'):
        """
        Quantify global batch effect using Sinkhorn divergence.
        
        Computes pairwise Sinkhorn divergence between all batches to
        assess overall batch effect magnitude.
        """
        X = adata.obsm[use_rep] if use_rep in adata.obsm else (
            adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X)
        
        divergences = []
        batches = adata.obs[batch_key].unique()
        
        if len(batches) < 2:
            return 0.0
        
        for i, batch_i in enumerate(batches):
            for j, batch_j in enumerate(batches):
                if i < j:
                    mask_i = adata.obs[batch_key] == batch_i
                    mask_j = adata.obs[batch_key] == batch_j
                    X_i, X_j = X[mask_i], X[mask_j]
                    
                    if len(X_i) < 5 or len(X_j) < 5:
                        continue
                    
                    sinkhorn_div = self.compute_sinkhorn_divergence(X_i, X_j)
                    if sinkhorn_div > 0:
                        divergences.append(sinkhorn_div)
        
        return np.mean(divergences) if divergences else 0.0
    
    def compute_cluster_divergence(self, X_cluster, batches):
        """
        Compute Sinkhorn divergence within a cluster across batches.
        
        Used for divergence-based adaptive correction: only clusters
        exceeding the threshold undergo transport-based alignment.
        """
        unique_batches = batches.unique()
        if len(unique_batches) < 2:
            return 0.0
        
        divergences = []
        for i, batch_i in enumerate(unique_batches):
            for j, batch_j in enumerate(unique_batches):
                if i < j:
                    X_i = X_cluster[batches == batch_i]
                    X_j = X_cluster[batches == batch_j]
                    
                    if len(X_i) < 5 or len(X_j) < 5:
                        continue
                    
                    sinkhorn_div = self.compute_sinkhorn_divergence(X_i, X_j)
                    if sinkhorn_div > 0:
                        divergences.append(sinkhorn_div)
        
        return np.mean(divergences) if divergences else 0.0
    
    def correct(self, adata, batch_key='batch', use_rep='X_pca'):
        """
        Apply SODA batch correction.
        
        Iterative process:
        1. Cluster cells (Harmony paradigm)
        2. For each cluster, compute Sinkhorn divergence
        3. If divergence > threshold, apply OT-based correction
        4. Repeat until convergence
        
        Returns corrected embedding.
        """
        if use_rep not in adata.obsm:
            raise ValueError(f"{use_rep} not found in adata.obsm")
        
        X = adata.obsm[use_rep].copy()
        
        for iteration in range(self.max_iter):
            # Harmony-style iterative clustering
            kmeans = KMeans(n_clusters=self.n_clusters, random_state=42, n_init=10)
            clusters = kmeans.fit_predict(X)
            
            # Cluster-wise divergence-based correction
            for cluster_id in range(self.n_clusters):
                mask_cluster = clusters == cluster_id
                if mask_cluster.sum() < 10:
                    continue
                
                # Compute cluster-specific Sinkhorn divergence
                cluster_divergence = self.compute_cluster_divergence(
                    X[mask_cluster], adata.obs.loc[mask_cluster, batch_key])
                
                # Adaptive correction: only correct if divergence exceeds threshold
                if cluster_divergence < self.divergence_threshold:
                    continue
                
                # Apply transport-based alignment
                X_cluster = X[mask_cluster]
                batches_cluster = adata.obs.loc[mask_cluster, batch_key]
                
                # Select reference batch (largest batch as anchor)
                ref_batch = batches_cluster.value_counts().idxmax()
                mask_ref = (batches_cluster == ref_batch).values
                X_ref = X_cluster[mask_ref]
                
                # Align other batches to reference via optimal transport
                for batch in batches_cluster.unique():
                    if batch == ref_batch:
                        continue
                    
                    mask_batch = (batches_cluster == batch).values
                    X_batch = X_cluster[mask_batch]
                    
                    if len(X_batch) == 0 or len(X_ref) == 0:
                        continue
                    
                    try:
                        # Compute entropy-regularized transport plan
                        M = ot.dist(X_batch, X_ref, metric='sqeuclidean')
                        a = np.ones(len(X_batch)) / len(X_batch)
                        b = np.ones(len(X_ref)) / len(X_ref)
                        T = ot.sinkhorn(a, b, M, reg=self.epsilon)
                        
                        # Apply probabilistic transport-based correction
                        X_corrected_batch = (T / T.sum(axis=1, keepdims=True)) @ X_ref
                        
                        # Soft correction to preserve biological structure
                        batch_indices = np.where(mask_cluster)[0][mask_batch]
                        X[batch_indices] = 0.7 * X_corrected_batch + 0.3 * X_batch
                    except:
                        continue
            
            # Check convergence using global Sinkhorn divergence
            if iteration >= 3:
                adata_temp = adata.copy()
                adata_temp.obsm[use_rep] = X
                global_divergence = self.detect_global_batch_effect(
                    adata_temp, batch_key=batch_key, use_rep=use_rep)
                
                if global_divergence < self.divergence_threshold:
                    break
        
        return X


class SingleCellBatchBenchmark:
    """
    Comprehensive benchmarking framework for single-cell batch correction.
    
    Evaluates SODA against baseline methods:
    - Harmony (iterative clustering baseline)
    - ComBat (linear method)
    - BBKNN (graph-based)
    - scVI (deep learning)
    """
    
    def __init__(self, h5ad_path="data/pbmc4k-0.h5ad", output_dir='results'):
        self.h5ad_path = h5ad_path
        self.output_dir = output_dir
        self.adata = None
        self.adata_dict = {}
        self.metrics = {}
        self._check_dependencies()
    
    def _check_dependencies(self):
        """Check availability of batch correction methods."""
        try:
            import bbknn
            self.has_bbknn = True
        except ImportError:
            self.has_bbknn = False
            print("Warning: bbknn not installed. Install with: pip install bbknn")
        
        try:
            import scvi
            self.has_scvi = True
        except ImportError:
            self.has_scvi = False
            print("Warning: scvi-tools not installed. Install with: pip install scvi-tools")
        
        self.has_harmony = hasattr(sc.external.pp, 'harmony_integrate')
        if not self.has_harmony:
            print("Warning: harmony not available. Update scanpy to get harmony support.")
        
        self.has_soda = HAS_OT
    
    def load_preprocessed_data(self):
        """Load preprocessed single-cell data from H5AD file."""
        
        print(f"Loading data from {self.h5ad_path}...")
        
        if not os.path.exists(self.h5ad_path):
            raise FileNotFoundError(
                f"Data file not found: {self.h5ad_path}\n"
                f"Please ensure the preprocessed data is available."
            )
        
        adata = sc.read_h5ad(self.h5ad_path)
        
        print(f"Data loaded successfully!")
        print(f"  - Cells: {adata.n_obs:,}")
        print(f"  - Genes (HVG): {adata.n_vars:,}")
        print(f"  - Batches: {adata.obs['batch'].unique()}")
        print(f"  - Batch 0: {(adata.obs['batch'] == 0).sum():,} cells")
        print(f"  - Batch 1: {(adata.obs['batch'] == 1).sum():,} cells")
        print(f"  - Original clusters: {adata.obs['leiden'].nunique()}")
        
        # Compute UMAP for original data if not present
        if 'X_umap' not in adata.obsm:
            print("  - Computing UMAP for original data...")
            if 'X_pca' not in adata.obsm:
                sc.tl.pca(adata, n_comps=50, random_state=42)
            sc.pp.neighbors(adata, n_pcs=30, random_state=42)
            sc.tl.umap(adata, random_state=42)
        
        self.adata = adata.copy()
        self.adata_dict['Original'] = adata.copy()
        
        return adata
    
    def apply_combat(self):
        """Apply ComBat batch correction (linear baseline)."""
        print("\nApplying ComBat...")
        try:
            adata_combat = self.adata.copy()
            sc.pp.combat(adata_combat, key='batch')
            sc.pp.neighbors(adata_combat, n_pcs=30, random_state=42)
            sc.tl.leiden(adata_combat, resolution=0.5, random_state=42)
            sc.tl.umap(adata_combat, random_state=42)
            self.adata_dict['ComBat'] = adata_combat
            print("  ✓ ComBat completed")
        except Exception as e:
            print(f"  ✗ ComBat failed: {e}")
    
    def apply_bbknn(self):
        """Apply BBKNN batch correction (graph-based)."""
        if not self.has_bbknn:
            print("\nSkipping BBKNN (not installed)")
            return
        
        print("\nApplying BBKNN...")
        try:
            import bbknn
            adata_bbknn = self.adata.copy()
            sc.tl.pca(adata_bbknn, n_comps=50, random_state=42)
            bbknn.bbknn(adata_bbknn, batch_key='batch', neighbors_within_batch=3)
            sc.tl.leiden(adata_bbknn, resolution=0.5, random_state=42)
            sc.tl.umap(adata_bbknn, random_state=42)
            self.adata_dict['BBKNN'] = adata_bbknn
            print("  ✓ BBKNN completed")
        except Exception as e:
            print(f"  ✗ BBKNN failed: {e}")
    
    def apply_scvi(self, max_epochs=100):
        """Apply scVI batch correction (deep learning)."""
        if not self.has_scvi:
            print("\nSkipping scVI (not installed)")
            return
        
        print("\nApplying scVI...")
        try:
            import scvi
            adata_scvi = self.adata.copy()
            
            adata_scvi.obs['batch'] = adata_scvi.obs['batch'].astype(str)
            
            # Prepare count data
            if adata_scvi.raw is not None:
                adata_full = adata_scvi.raw.to_adata()
                hvg_genes = adata_scvi.var_names
                if 'counts' in adata_full.layers:
                    adata_scvi.layers['counts'] = adata_full[:, hvg_genes].layers['counts']
                else:
                    adata_scvi.layers['counts'] = adata_full[:, hvg_genes].X.copy()
            elif 'counts' in adata_scvi.layers:
                pass
            else:
                print("  Warning: No raw counts available, using approximation")
                adata_scvi.layers['counts'] = adata_scvi.X.copy()
            
            if not sparse.issparse(adata_scvi.layers['counts']):
                adata_scvi.layers['counts'] = sparse.csr_matrix(adata_scvi.layers['counts'])
            
            scvi.settings.seed = 42
            scvi.model.SCVI.setup_anndata(adata_scvi, layer='counts', batch_key='batch')
            model = scvi.model.SCVI(adata_scvi, n_latent=30, n_layers=2)
            
            try:
                model.train(max_epochs=max_epochs, train_size=0.9, early_stopping=True, 
                          check_val_every_n_epoch=10, accelerator='cpu')
            except:
                model.train(max_epochs=max_epochs, train_size=0.9, early_stopping=True, 
                          check_val_every_n_epoch=10, use_gpu=False)
            
            adata_scvi.obsm['X_scVI'] = model.get_latent_representation()
            sc.pp.neighbors(adata_scvi, use_rep='X_scVI', random_state=42)
            sc.tl.leiden(adata_scvi, resolution=0.5, random_state=42)
            sc.tl.umap(adata_scvi, random_state=42)
            self.adata_dict['scVI'] = adata_scvi
            print("  ✓ scVI completed")
        except Exception as e:
            print(f"  ✗ scVI failed: {e}")
    
    def apply_harmony(self):
        """Apply Harmony batch correction (iterative clustering baseline)."""
        if not self.has_harmony:
            print("\nSkipping Harmony (not available)")
            return
        
        print("\nApplying Harmony...")
        try:
            adata_harmony = self.adata.copy()
            if 'X_pca' not in adata_harmony.obsm:
                sc.tl.pca(adata_harmony, n_comps=50, random_state=42)
            
            if not isinstance(adata_harmony.obs['batch'].dtype, pd.CategoricalDtype):
                adata_harmony.obs['batch'] = adata_harmony.obs['batch'].astype(str).astype('category')
            
            sc.external.pp.harmony_integrate(adata_harmony, key='batch', basis='X_pca',
                                            adjusted_basis='X_harmony', max_iter_harmony=10)
            sc.pp.neighbors(adata_harmony, use_rep='X_harmony', random_state=42)
            sc.tl.leiden(adata_harmony, resolution=0.5, random_state=42)
            sc.tl.umap(adata_harmony, random_state=42)
            self.adata_dict['Harmony'] = adata_harmony
            print("  ✓ Harmony completed")
        except Exception as e:
            print(f"  ✗ Harmony failed: {e}")
    
    def apply_soda(self, n_clusters=50, max_iter=10, epsilon=0.05, divergence_threshold=0.05):
        """Apply SODA batch correction (proposed method)."""
        if not self.has_soda:
            print("\nSkipping SODA (POT not installed)")
            return
        
        print("\nApplying SODA...")
        try:
            adata_soda = self.adata.copy()
            if 'X_pca' not in adata_soda.obsm:
                sc.tl.pca(adata_soda, n_comps=50, random_state=42)
            
            corrector = SODABatchCorrector(
                n_clusters=n_clusters, 
                max_iter=max_iter, 
                epsilon=epsilon, 
                divergence_threshold=divergence_threshold
            )
            X_corrected = corrector.correct(adata_soda, batch_key='batch', use_rep='X_pca')
            adata_soda.obsm['X_soda'] = X_corrected
            
            sc.pp.neighbors(adata_soda, use_rep='X_soda', random_state=42)
            sc.tl.leiden(adata_soda, resolution=0.5, random_state=42)
            sc.tl.umap(adata_soda, random_state=42)
            self.adata_dict['SODA'] = adata_soda
            print("  ✓ SODA completed")
        except Exception as e:
            print(f"  ✗ SODA failed: {e}")
    
    def apply_all_methods(self, scvi_epochs=100, soda_n_clusters=50, soda_max_iter=10):
        """Apply all available batch correction methods."""
        print("\n" + "="*80)
        print("APPLYING BATCH CORRECTION METHODS")
        print("="*80)
        
        self.apply_combat()
        self.apply_bbknn()
        self.apply_scvi(max_epochs=scvi_epochs)
        self.apply_harmony()
        self.apply_soda(n_clusters=soda_n_clusters, max_iter=soda_max_iter)
        
        print(f"\n✓ Applied {len(self.adata_dict)-1} batch correction methods")
    
    def calculate_neighborhood_preservation(self, k=30):
        """
        Calculate nearest neighbor rank change (neighborhood preservation metric).
        
        Lower values indicate better preservation of local structure.
        """
        print("\nCalculating neighborhood preservation...")
        
        X_orig = self.adata.X.toarray() if sparse.issparse(self.adata.X) else self.adata.X
        dist_orig = cdist(X_orig, X_orig, metric='euclidean')
        
        def get_top_k_neighbors(dist_matrix, k):
            neighbors = np.zeros((dist_matrix.shape[0], k), dtype=int)
            for i in range(dist_matrix.shape[0]):
                neighbors[i, :] = np.argsort(dist_matrix[i, :])[1:k+1]
            return neighbors
        
        top_k_orig = get_top_k_neighbors(dist_orig, k)
        results = {}
        
        for method_name, adata in self.adata_dict.items():
            if method_name == 'Original':
                continue
            
            X_corrected = adata.X.toarray() if sparse.issparse(adata.X) else adata.X
            dist_corrected = cdist(X_corrected, X_corrected, metric='euclidean')
            
            rank_changes = []
            for i in range(dist_orig.shape[0]):
                cell_rank_changes = []
                for j, neighbor_idx in enumerate(top_k_orig[i]):
                    new_rank = np.where(np.argsort(dist_corrected[i, :]) == neighbor_idx)[0][0]
                    rank_change = abs(new_rank - (j + 1))
                    cell_rank_changes.append(rank_change)
                rank_changes.append(np.median(cell_rank_changes))
            
            results[method_name] = rank_changes
        
        self.metrics['neighborhood_preservation'] = pd.DataFrame(results)
        print("  ✓ Neighborhood preservation calculated")
        return self.metrics['neighborhood_preservation']
    
    def calculate_cluster_concordance(self):
        """
        Calculate cluster concordance score.
        
        Lower values indicate better preservation of original biological structure.
        """
        print("\nCalculating cluster concordance...")
        
        orig_clusters = self.adata.obs['leiden'].astype(int).values
        results = {}
        
        for method_name, adata in self.adata_dict.items():
            if method_name == 'Original':
                continue
            
            corrected_clusters = adata.obs['leiden'].astype(int).values
            n_orig_clusters = len(np.unique(orig_clusters))
            n_corrected_clusters = len(np.unique(corrected_clusters))
            
            confusion = np.zeros((n_orig_clusters, n_corrected_clusters))
            for i in range(len(orig_clusters)):
                confusion[orig_clusters[i], corrected_clusters[i]] += 1
            
            cluster_changes = []
            for i in range(n_orig_clusters):
                if confusion[i].sum() > 0:
                    cluster_changes.append((confusion[i].sum() - confusion[i].max()) / confusion[i].sum())
            
            cluster_changes_reverse = []
            for j in range(n_corrected_clusters):
                if confusion[:, j].sum() > 0:
                    cluster_changes_reverse.append((confusion[:, j].sum() - confusion[:, j].max()) / confusion[:, j].sum())
            
            results[method_name] = (np.mean(cluster_changes) + np.mean(cluster_changes_reverse)) / 2
        
        self.metrics['cluster_concordance'] = pd.Series(results)
        print("  ✓ Cluster concordance calculated")
        return self.metrics['cluster_concordance']
    
    def calculate_ari(self):
        """
        Calculate Adjusted Rand Index (ARI).
        
        Higher values indicate better preservation of biological clusters.
        """
        print("\nCalculating ARI...")
        
        orig_clusters = self.adata.obs['leiden'].values
        results = {}
        
        for method_name, adata in self.adata_dict.items():
            if method_name == 'Original':
                continue
            results[method_name] = adjusted_rand_score(orig_clusters, adata.obs['leiden'].values)
        
        self.metrics['ari'] = pd.Series(results)
        print("  ✓ ARI calculated")
        return self.metrics['ari']
    
    def calculate_all_metrics(self):
        """Calculate all evaluation metrics."""
        print("\n" + "="*80)
        print("CALCULATING EVALUATION METRICS")
        print("="*80)
        
        self.calculate_neighborhood_preservation(k=30)
        self.calculate_cluster_concordance()
        self.calculate_ari()
        
        print("\n✓ All metrics calculated")
        return self.metrics
    
    def plot_umap_comparison(self, save_path=None):
        """Plot UMAP visualizations colored by batch."""
        if save_path is None:
            save_path = os.path.join(self.output_dir, 'umap_comparison.png')
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        
        valid_methods = {name: adata for name, adata in self.adata_dict.items() 
                        if 'X_umap' in adata.obsm}
        
        if not valid_methods:
            print("  ✗ No methods with UMAP available for plotting")
            return
        
        n_methods = len(valid_methods)
        ncols, nrows = 3, (n_methods + 2) // 3
        
        fig, axes = plt.subplots(nrows, ncols, figsize=(5*ncols, 5*nrows))
        axes = axes.flatten() if n_methods > 1 else [axes]
        
        for idx, (method_name, adata) in enumerate(valid_methods.items()):
            ax = axes[idx]
            batch_categories = adata.obs['batch'].cat.categories if hasattr(adata.obs['batch'], 'cat') else adata.obs['batch'].unique()
            
            for batch in batch_categories:
                mask = adata.obs['batch'] == batch
                batch_int = int(float(batch)) if isinstance(batch, str) else int(batch)
                ax.scatter(adata.obsm['X_umap'][mask, 0], adata.obsm['X_umap'][mask, 1],
                          c=f'C{batch_int}', label=f'Batch {batch_int}', s=8, alpha=0.6, edgecolors='none')
            
            ax.set_title(f'{method_name}', fontsize=14, fontweight='bold', pad=10)
            ax.set_xlabel('UMAP 1', fontsize=11)
            ax.set_ylabel('UMAP 2', fontsize=11)
            ax.legend(markerscale=2, loc='upper right', fontsize=9)
            ax.grid(True, alpha=0.2)
            ax.set_aspect('equal', 'box')
        
        for idx in range(n_methods, len(axes)):
            axes[idx].axis('off')
        
        plt.suptitle('UMAP Visualization - Batch Integration', fontsize=16, fontweight='bold', y=1.00)
        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  ✓ Saved: {save_path}")
    
    def plot_umap_clusters(self, save_path=None):
        """Plot UMAP visualizations colored by biological cluster."""
        if save_path is None:
            save_path = os.path.join(self.output_dir, 'umap_clusters.png')
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        valid_methods = {name: adata for name, adata in self.adata_dict.items() 
                        if 'X_umap' in adata.obsm}
        
        if not valid_methods:
            print("  ✗ No methods with UMAP available for plotting")
            return
        
        n_methods = len(valid_methods)
        ncols, nrows = 3, (n_methods + 2) // 3
        
        fig, axes = plt.subplots(nrows, ncols, figsize=(5*ncols, 5*nrows))
        axes = axes.flatten() if n_methods > 1 else [axes]
        
        for idx, (method_name, adata) in enumerate(valid_methods.items()):
            ax = axes[idx]
            ax.scatter(adata.obsm['X_umap'][:, 0], adata.obsm['X_umap'][:, 1],
                      c=adata.obs['leiden'].astype(int), cmap='tab20', s=8, alpha=0.6, edgecolors='none')
            
            n_clusters = len(adata.obs['leiden'].unique())
            ax.set_title(f'{method_name}\n({n_clusters} clusters)', fontsize=12, fontweight='bold', pad=10)
            ax.set_xlabel('UMAP 1', fontsize=11)
            ax.set_ylabel('UMAP 2', fontsize=11)
            ax.grid(True, alpha=0.2)
            ax.set_aspect('equal', 'box')
        
        for idx in range(n_methods, len(axes)):
            axes[idx].axis('off')
        
        plt.suptitle('UMAP Visualization - Biological Clusters', fontsize=16, fontweight='bold', y=1.00)
        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  ✓ Saved: {save_path}")
    
    def plot_metrics_summary(self, save_path=None):
        """Plot comprehensive metrics summary."""
        if save_path is None:
            save_path = os.path.join(self.output_dir, 'metrics_summary.png')
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        
        if not self.metrics:
            return
        
        fig, axes = plt.subplots(1, 3, figsize=(20, 5))
        
        # Neighborhood Preservation
        if 'neighborhood_preservation' in self.metrics:
            data = self.metrics['neighborhood_preservation']
            ax = axes[0]
            bp = ax.boxplot([data[col].values for col in data.columns], labels=data.columns,
                           patch_artist=True, showfliers=True, widths=0.6)
            
            colors = sns.color_palette("Set2", len(data.columns))
            for patch, color in zip(bp['boxes'], colors):
                patch.set_facecolor(color)
                patch.set_alpha(0.7)
            
            ax.set_title('Neighborhood Preservation\n(Lower = Better)', fontsize=13, fontweight='bold', pad=10)
            ax.set_ylabel('Median Rank Change', fontsize=11)
            ax.grid(True, alpha=0.3, axis='y')
            plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')
            
            medians = data.median().values
            for i, median in enumerate(medians):
                ax.text(i+1, median, f'{median:.1f}', ha='center', va='bottom', fontweight='bold', fontsize=9)
        
        # Cluster Concordance
        if 'cluster_concordance' in self.metrics:
            data = self.metrics['cluster_concordance']
            ax = axes[1]
            colors = sns.color_palette("Set2", len(data))
            bars = ax.bar(range(len(data)), data.values, color=colors, alpha=0.8)
            
            ax.set_xticks(range(len(data)))
            ax.set_xticklabels(data.index, rotation=45, ha='right')
            ax.set_title('Cluster Concordance\n(Lower = Better)', fontsize=13, fontweight='bold', pad=10)
            ax.set_ylabel('Concordance Score', fontsize=11)
            ax.grid(True, alpha=0.3, axis='y')
            
            for bar, value in zip(bars, data.values):
                ax.text(bar.get_x() + bar.get_width()/2, value, f'{value:.4f}', 
                       ha='center', va='bottom', fontsize=9)
        
        # ARI
        if 'ari' in self.metrics:
            data = self.metrics['ari']
            ax = axes[2]
            colors = sns.color_palette("Set2", len(data))
            bars = ax.bar(range(len(data)), data.values, color=colors, alpha=0.8)
            
            ax.set_xticks(range(len(data)))
            ax.set_xticklabels(data.index, rotation=45, ha='right')
            ax.set_title('Adjusted Rand Index\n(Higher = Better)', fontsize=13, fontweight='bold', pad=10)
            ax.set_ylabel('ARI Score', fontsize=11)
            ax.set_ylim([0, 1.1])
            ax.axhline(y=0.8, color='red', linestyle='--', alpha=0.5, label='Excellent (>0.8)')
            ax.grid(True, alpha=0.3, axis='y')
            ax.legend(fontsize=9)
            
            for bar, value in zip(bars, data.values):
                ax.text(bar.get_x() + bar.get_width()/2, value, f'{value:.4f}', 
                       ha='center', va='bottom', fontsize=9)
        
        plt.suptitle('SODA Benchmark: Performance Metrics', fontsize=16, fontweight='bold', y=1.00)
        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  ✓ Saved: {save_path}")
    
    def plot_neighborhood_distribution(self, save_path=None):
        """Plot distribution of neighborhood preservation scores."""
        if save_path is None:
            save_path = os.path.join(self.output_dir, 'neighborhood_distribution.png')
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        
        if 'neighborhood_preservation' not in self.metrics:
            return
        
        data = self.metrics['neighborhood_preservation']
        fig, ax = plt.subplots(figsize=(14, 6))
        
        parts = ax.violinplot([data[col].values for col in data.columns],
                             positions=range(len(data.columns)), showmeans=True, showmedians=True, widths=0.7)
        
        colors = sns.color_palette("Set2", len(data.columns))
        for i, pc in enumerate(parts['bodies']):
            pc.set_facecolor(colors[i])
            pc.set_alpha(0.7)
        
        ax.set_xticks(range(len(data.columns)))
        ax.set_xticklabels(data.columns, rotation=45, ha='right')
        ax.set_title('Neighborhood Preservation Distribution', fontsize=14, fontweight='bold', pad=15)
        ax.set_ylabel('Rank Change (log scale)', fontsize=12)
        ax.set_yscale('log')
        ax.grid(True, alpha=0.3, axis='y', which='both')
        
        medians = data.median().values
        ax.plot(range(len(medians)), medians, 'r-', linewidth=2, label='Median', marker='o', markersize=8)
        ax.legend(fontsize=11)
        
        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  ✓ Saved: {save_path}")
    
    def generate_all_plots(self):
        """Generate all visualization plots."""
        print("\n" + "="*80)
        print("GENERATING VISUALIZATIONS")
        print("="*80)
        
        umap_comp_path = os.path.join(self.output_dir, 'umap_comparison.png')
        umap_clust_path = os.path.join(self.output_dir, 'umap_clusters.png')
        metrics_path = os.path.join(self.output_dir, 'metrics_summary.png')
        neighborhood_path = os.path.join(self.output_dir, 'neighborhood_distribution.png')
        
        self.plot_umap_comparison(save_path=umap_comp_path)
        self.plot_umap_clusters(save_path=umap_clust_path)
        self.plot_metrics_summary(save_path=metrics_path)
        self.plot_neighborhood_distribution(save_path=neighborhood_path)
        
        print("\n✓ All visualizations generated")
    
    def generate_report(self, save_path=None):
        """Generate comprehensive benchmark report."""
        if save_path is None:
            save_path = os.path.join(self.output_dir, 'benchmark_report.txt')
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        
        with open(save_path, 'w', encoding='utf-8') as f:
            f.write("=" * 80 + "\n")
            f.write(" " * 15 + "SODA: SINGLE-CELL BATCH CORRECTION BENCHMARK\n")
            f.write("=" * 80 + "\n\n")
            
            # Dataset info
            f.write("1. DATASET INFORMATION\n" + "-" * 80 + "\n")
            f.write(f"Total cells: {self.adata.n_obs:,}\n")
            f.write(f"Total genes (HVG): {self.adata.n_vars:,}\n")
            
            batch_col = self.adata.obs['batch']
            n_batch0 = ((batch_col == 0) | (batch_col == '0')).sum()
            n_batch1 = ((batch_col == 1) | (batch_col == '1')).sum()
            f.write(f"Batch 0 cells: {n_batch0:,}\n")
            f.write(f"Batch 1 cells: {n_batch1:,}\n")
            f.write(f"Original clusters: {len(self.adata.obs['leiden'].unique())}\n\n")
            
            # Methods
            f.write("2. EVALUATED METHODS\n" + "-" * 80 + "\n")
            for method in self.adata_dict.keys():
                if method != 'Original':
                    n_clusters = len(self.adata_dict[method].obs['leiden'].unique())
                    f.write(f"  ✓ {method:20s}: {n_clusters} clusters\n")
            f.write(f"\nTotal methods: {len(self.adata_dict)-1}\n\n")
            
            # Metrics
            if self.metrics:
                f.write("3. PERFORMANCE METRICS\n" + "-" * 80 + "\n\n")
                
                # Neighborhood Preservation
                if 'neighborhood_preservation' in self.metrics:
                    f.write("3.1 Neighborhood Preservation (Lower = Better)\n")
                    medians = self.metrics['neighborhood_preservation'].median().sort_values()
                    for rank, (method, value) in enumerate(medians.items(), 1):
                        grade = "Excellent" if value < 10 else "Good" if value < 30 else "Fair" if value < 50 else "Poor"
                        f.write(f"    {rank}. {method:20s}: {value:8.2f}  [{grade}]\n")
                    f.write(f"\n    ★ Best: {medians.idxmin()}\n\n")
                
                # Cluster Concordance
                if 'cluster_concordance' in self.metrics:
                    f.write("3.2 Cluster Concordance (Lower = Better)\n")
                    cc = self.metrics['cluster_concordance'].sort_values()
                    for rank, (method, value) in enumerate(cc.items(), 1):
                        grade = "Excellent" if value < 0.1 else "Good" if value < 0.3 else "Fair" if value < 0.5 else "Poor"
                        f.write(f"    {rank}. {method:20s}: {value:8.4f}  [{grade}]\n")
                    f.write(f"\n    ★ Best: {cc.idxmin()}\n\n")
                
                # ARI
                if 'ari' in self.metrics:
                    f.write("3.3 Adjusted Rand Index (Higher = Better)\n")
                    ari = self.metrics['ari'].sort_values(ascending=False)
                    for rank, (method, value) in enumerate(ari.items(), 1):
                        grade = "Excellent" if value > 0.8 else "Good" if value > 0.5 else "Fair" if value > 0.3 else "Poor"
                        f.write(f"    {rank}. {method:20s}: {value:8.4f}  [{grade}]\n")
                    f.write(f"\n    ★ Best: {ari.idxmax()}\n\n")
                
                # Overall ranking
                f.write("4. OVERALL RANKING\n" + "-" * 80 + "\n")
                scores = {}
                for method in self.adata_dict.keys():
                    if method == 'Original':
                        continue
                    
                    score, n_metrics = 0, 0
                    
                    if 'neighborhood_preservation' in self.metrics:
                        nn = self.metrics['neighborhood_preservation'].median()
                        nn_range = nn.max() - nn.min()
                        if nn_range > 0:
                            score += (1 - (nn[method] - nn.min()) / nn_range)
                            n_metrics += 1
                    
                    if 'cluster_concordance' in self.metrics:
                        cc = self.metrics['cluster_concordance']
                        cc_range = cc.max() - cc.min()
                        if cc_range > 0:
                            score += (1 - (cc[method] - cc.min()) / cc_range)
                            n_metrics += 1
                    
                    if 'ari' in self.metrics:
                        ari = self.metrics['ari']
                        ari_range = ari.max() - ari.min()
                        if ari_range > 0:
                            score += (ari[method] - ari.min()) / ari_range
                            n_metrics += 1
                    
                    if n_metrics > 0:
                        scores[method] = score / n_metrics
                
                sorted_methods = sorted(scores.items(), key=lambda x: x[1], reverse=True)
                for rank, (method, score) in enumerate(sorted_methods, 1):
                    stars = "★" * min(5, int(score * 5) + 1)
                    f.write(f"  {rank}. {method:20s}: {score:.4f}  {stars}\n")
                
                if sorted_methods:
                    f.write(f"\n  ★★★ TOP PERFORMER: {sorted_methods[0][0]} ★★★\n\n")
            
            f.write("=" * 80 + "\n")
            f.write("SODA: Sinkhorn Optimal Transport with Divergence-Based Adaptation\n")
            f.write("=" * 80 + "\n")
        
        # Print to console
        with open(save_path, 'r', encoding='utf-8') as f:
            print("\n" + f.read())
        
        print(f"✓ Report saved to: {save_path}")
    
    def run_full_benchmark(self, scvi_epochs=100, output_dir='results', 
                          soda_n_clusters=50, soda_max_iter=10):
        """Run complete batch correction benchmark."""
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        
        print("\n" + "="*80)
        print("SODA: SINGLE-CELL BATCH CORRECTION BENCHMARK")
        print("="*80)
        
        self.load_preprocessed_data()
        self.apply_all_methods(scvi_epochs=scvi_epochs, 
                              soda_n_clusters=soda_n_clusters, 
                              soda_max_iter=soda_max_iter)
        self.calculate_all_metrics()
        self.generate_all_plots()
        
        report_path = os.path.join(output_dir, 'benchmark_report.txt')
        self.generate_report(save_path=report_path)
        
        print("\n" + "="*80)
        print(f"BENCHMARK COMPLETE! Results saved to {output_dir}/")
        print("="*80)


def main():
    """Main execution function with command-line argument support."""
    parser = argparse.ArgumentParser(
        description='SODA: Sinkhorn Optimal Transport with Divergence-Based Adaptation for Single-Cell Integration',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python soda_benchmark.py --input data/simulated_data.h5ad
  python soda_benchmark.py --input data/pbmc4k-0.h5ad --epochs 200 --output results_pbmc
  python soda_benchmark.py --input data/my_data.h5ad --batch-key batch_id
        """
    )
    
    parser.add_argument(
        '--input', '-i',
        type=str,
        required=True,
        help='Path to input H5AD file (required)'
    )
    
    parser.add_argument(
        '--output', '-o',
        type=str,
        default='results',
        help='Output directory for results (default: results)'
    )
    
    parser.add_argument(
        '--epochs',
        type=int,
        default=100,
        help='Maximum epochs for scVI training (default: 100)'
    )
    
    parser.add_argument(
        '--batch-key',
        type=str,
        default='batch',
        help='Key in adata.obs for batch information (default: batch)'
    )
    
    parser.add_argument(
        '--n-clusters',
        type=int,
        default=50,
        help='Number of clusters for SODA (default: 50)'
    )
    
    parser.add_argument(
        '--max-iter',
        type=int,
        default=10,
        help='Maximum iterations for SODA (default: 10)'
    )
    
    args = parser.parse_args()
    
    h5ad_path = args.input
    
    if not os.path.exists(h5ad_path):
        print(f"Error: Data file not found: {h5ad_path}")
        print("\nPlease ensure you have:")
        print("1. Downloaded the required single-cell dataset")
        print("2. Run preprocessing to generate the H5AD file")
        print(f"\nUsage: python soda_benchmark.py --input <path_to_h5ad_file>")
        sys.exit(1)
    
    print(f"Input file: {h5ad_path}")
    print(f"Output directory: {args.output}")
    print(f"scVI epochs: {args.epochs}")
    print(f"Batch key: {args.batch_key}")
    print(f"SODA n_clusters: {args.n_clusters}")
    print(f"SODA max_iter: {args.max_iter}")
    
    # Run SODA benchmark
    benchmark = SingleCellBatchBenchmark(h5ad_path=h5ad_path, output_dir=args.output)
    benchmark.run_full_benchmark(
        scvi_epochs=args.epochs,
        output_dir=args.output,
        soda_n_clusters=args.n_clusters,
        soda_max_iter=args.max_iter
    )


if __name__ == '__main__':
    main()
    