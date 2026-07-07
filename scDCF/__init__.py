__version__ = '0.1.16'

# Import commonly used functions for convenience
from .utils import read_gene_symbols
from .control_genes import generate_control_genes
from .analysis import monte_carlo_comparison
from .post_analysis import (
    load_monte_carlo_results,
    combine_p_values_across_iterations,
    visualize_combined_p_values,
    export_final_celltype_summary,
    organize_results,
    apply_dataset_level_cell_fdr,
    compute_celltype_enrichment
)
from .trait_association import get_trait_association_scores

# Cell metadata helpers (new in v0.1.13)
from .cell_metadata import add_cell_metadata, create_cell_id_mapping, validate_cell_ids, get_cell_info_by_id

# Parallel processing (new in v0.1.12)
try:
    from .parallel import parallel_monte_carlo_comparison, auto_monte_carlo
    __parallel_available__ = True
except ImportError:
    __parallel_available__ = False
    # Fallback if multiprocessing issues
    def parallel_monte_carlo_comparison(*args, **kwargs):
        import warnings
        warnings.warn("Parallel processing not available, falling back to serial")
        return monte_carlo_comparison(*args, **kwargs)
    auto_monte_carlo = monte_carlo_comparison

# Public output-organization aliases
organize_output = organize_results
organize_final_output = organize_results
