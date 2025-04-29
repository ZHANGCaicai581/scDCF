import os
import logging
import json
import numpy as np
import matplotlib.pyplot as plt

def read_gene_symbols(file_path, zstat_column=None):
    """
    Read gene symbols from a file.
    
    Args:
        file_path: Path to the file containing gene symbols
        zstat_column: Column name containing z-statistic values (optional)
        
    Returns:
        DataFrame containing gene symbols and z-statistics
    """
    import pandas as pd
    
    file_ext = os.path.splitext(file_path)[1].lower()
    
    # Handle CSV files
    if file_ext == '.csv':
        gene_df = pd.read_csv(file_path)
        
        # Check if 'gene_name' column exists
        if 'gene_name' not in gene_df.columns:
            # Try to use the first column as gene names
            if gene_df.shape[1] >= 1:
                gene_df.rename(columns={gene_df.columns[0]: 'gene_name'}, inplace=True)
                logging.warning(f"Renamed first column to 'gene_name'")
            else:
                raise ValueError("CSV file must contain at least one column for gene names")
        
        # Check if 'zstat' column exists or was specified
        if zstat_column is not None and zstat_column in gene_df.columns:
            gene_df.rename(columns={zstat_column: 'zstat'}, inplace=True)
        elif 'zstat' not in gene_df.columns:
            # If there's at least a second column, use it as z-statistic
            if gene_df.shape[1] >= 2:
                gene_df.rename(columns={gene_df.columns[1]: 'zstat'}, inplace=True)
                logging.warning(f"Renamed second column to 'zstat'")
            else:
                # If no z-statistic column, add a default value of 1.0
                gene_df['zstat'] = 1.0
                logging.warning("No z-statistic column found. Using default value of 1.0")
    
    # Handle text files
    elif file_ext in ['.txt', '.tsv']:
        with open(file_path, 'r') as f:
            lines = f.readlines()
        
        # Parse the text file
        genes = []
        for line in lines:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
                
            parts = line.split('\t' if file_ext == '.tsv' else None)
            if len(parts) >= 2:
                gene_name = parts[0]
                try:
                    zstat = float(parts[1])
                except ValueError:
                    zstat = 1.0  # Default value if not a valid number
                genes.append({'gene_name': gene_name, 'zstat': zstat})
            else:
                # If only gene name is provided, use default zstat of 1.0
                genes.append({'gene_name': parts[0], 'zstat': 1.0})
        
        gene_df = pd.DataFrame(genes)
    else:
        raise ValueError(f"Unsupported file extension: {file_ext}. Use .csv, .txt, or .tsv")
    
    return gene_df 

def filter_valid_genes(gene_list, available_genes):
    """
    Filter a list of genes to include only those available in the given dataset.
    
    Args:
        gene_list: List of gene symbols to filter
        available_genes: List of gene symbols available in the dataset
        
    Returns:
        List of valid gene symbols
    """
    valid_genes = [gene for gene in gene_list if gene in available_genes]
    if len(valid_genes) < len(gene_list):
        logging.warning(f"Filtered out {len(gene_list) - len(valid_genes)} genes not found in the dataset")
    return valid_genes


def load_control_genes(file_path):
    """
    Load control genes from a JSON file.
    
    Args:
        file_path: Path to the JSON file containing control genes
        
    Returns:
        Tuple of (disease_control_genes, healthy_control_genes)
    """
    if not os.path.exists(file_path):
        logging.warning(f"Control genes file not found: {file_path}")
        return None, None
    
    try:
        with open(file_path, 'r') as f:
            control_genes = json.load(f)
        
        disease_control_genes = control_genes.get('disease_control_genes', {})
        healthy_control_genes = control_genes.get('healthy_control_genes', {})
        
        logging.info(f"Loaded control genes for {len(disease_control_genes)} disease genes and {len(healthy_control_genes)} healthy genes")
        
        return disease_control_genes, healthy_control_genes
    except Exception as e:
        logging.error(f"Error loading control genes file: {e}")
        return None, None 

def normalize_weights(weights):
    """Normalize weights to sum to 1."""
    if not isinstance(weights, (np.ndarray, list)):
        logging.error("Weights should be provided as a NumPy array or list.")
        return np.array([])

    weights = np.abs(np.array(weights))  # Ensure weights are positive and convert to NumPy array if needed
    total = np.sum(weights)
    if total == 0:
        logging.warning("Sum of weights is zero. Normalizing to equal weights.")
        return np.ones_like(weights) / len(weights)
    return weights / total 

def save_results(df, output_file):
    """Save the results to a CSV file."""
    output_dir = os.path.dirname(output_file)
    try:
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        df.to_csv(output_file, index=False)
        logging.info(f"Results saved to {output_file}, containing {len(df)} entries.")
    except OSError as e:
        logging.error(f"Failed to create directory {output_dir}: {e}")
    except Exception as e:
        logging.error(f"Error saving results to {output_file}: {e}") 

plt.xlabel('Combined P-value (Fisher\'s Method)') 