__version__ = '0.1.14'

# Import commonly used functions for convenience
from .utils import read_gene_symbols
from .control_genes import generate_control_genes
from .analysis import monte_carlo_comparison
from .post_analysis import (
    load_monte_carlo_results,
    combine_p_values_across_iterations,
    visualize_combined_p_values,
    export_final_celltype_summary
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

# Import organize_final_output from post_trait_test
import sys
import os
import importlib.util

# Add the parent directory to sys.path to import post_trait_test
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

# Import organize_final_output from post_trait_test
try:
    from post_trait_test import organize_final_output
    # Create an alias for backward compatibility
    organize_output = organize_final_output
except ImportError:
    # Define fallback if import fails
    def organize_output(source_dir, dest_dir):
        """Fallback function if post_trait_test can't be imported."""
        print(f"Organizing output from {source_dir} to {dest_dir}")
        print("Warning: Using fallback function - post_trait_test module not found")
        return dest_dir
    
    organize_final_output = organize_output

