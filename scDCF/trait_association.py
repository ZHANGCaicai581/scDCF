def get_trait_association_scores(output_dir, cell_type):
    """
    Calculate trait association scores for a specific cell type.
    
    Args:
        output_dir: Base output directory
        cell_type: Cell type to analyze
    
    Returns:
        DataFrame with trait association scores
    """
    import os
    import pandas as pd
    import numpy as np
    import logging
    from scipy.stats import rankdata
    
    logging.info(f"Calculating trait association scores for {cell_type}")
    
    # Note: Combined p-values are saved per cell-type directory by post_analysis.
    # We'll search multiple filename patterns to maximize compatibility.
    
    # Define contexts to analyze
    contexts = ['disease', 'healthy']
    
    # Define score calculation methods
    methods = ['negative_log10', 'z_score']
    
    # Get cell type directory
    cell_type_dir = os.path.join(output_dir, cell_type)
    
    # Check if directory exists
    if not os.path.exists(cell_type_dir):
        logging.warning(f"Directory for {cell_type} not found at {cell_type_dir}")
        return None
    
    # Store results for each context
    context_results = {}
    
    # Process each context (disease/healthy)
    for context in contexts:
        # Candidate file patterns (ordered by preference)
        candidate_files = [
            os.path.join(cell_type_dir, f"{context}_combined_p_values.csv"),
            os.path.join(cell_type_dir, f"{cell_type}_{context}_combined.csv"),
            os.path.join(output_dir, f"{cell_type}_{context}_combined.csv"),  # legacy top-level
        ]

        combined_file = None
        for path in candidate_files:
            if os.path.exists(path):
                combined_file = path
                break

        if combined_file is None:
            logging.warning(f"Combined p-values file not found for {context} in {cell_type}")
            logging.info(f"Searched candidates: {[os.path.abspath(p) for p in candidate_files]}")
            continue
        
        # Load data
        df = pd.read_csv(combined_file)
        
        # Rename 'gene' to 'cell_id' if needed
        if 'gene' in df.columns and 'cell_id' not in df.columns:
            df = df.rename(columns={'gene': 'cell_id'})
            logging.info("Renamed 'gene' column to 'cell_id'")
        
        if df.empty:
            logging.warning(f"No data found for {context} in {cell_type}")
            continue
        
        # Determine p-value column name - handle different possible names
        p_value_col = None
        for possible_col in ['combined_p_value_fisher', 'combined_p_value_stouffer', 'combined_p_value', 
                             'p_value', 'adjusted_p_value', 'p_value_adj', 'pvalue', 'min_p_value']:
            if possible_col in df.columns:
                p_value_col = possible_col
                logging.info(f"Found p-value column: {p_value_col}")
                break
                
        if p_value_col is None:
            logging.warning(f"Could not find p-value column in {context} file for {cell_type}")
            logging.info(f"Available columns: {df.columns.tolist()}")
            continue
        
        # Calculate scores using different methods
        results = {}
        
        # Method 1: -log10(p-value)
        if 'negative_log10' in methods:
            # Avoid log(0) by replacing 0 with smallest non-zero value
            min_non_zero = df[p_value_col][df[p_value_col] > 0].min() if any(df[p_value_col] > 0) else 1e-300
            df['p_value_used'] = df[p_value_col].replace(0, min_non_zero)
            df['neg_log10_p'] = -np.log10(df['p_value_used'])
            
            # Normalize to 0-1 range
            max_val = df['neg_log10_p'].max()
            df['neg_log10_score'] = df['neg_log10_p'] / max_val if max_val > 0 else 0
            
            results['negative_log10'] = df[['cell_id', 'neg_log10_score']].rename(
                columns={'neg_log10_score': f'{context}_neg_log10_score'}
            )
        
        # Method 2: Z-score of p-values
        if 'z_score' in methods:
            df['z_score'] = (df[p_value_col] - df[p_value_col].mean()) / df[p_value_col].std() if df[p_value_col].std() > 0 else 0
            # Invert so that lower p-values have higher Z-scores
            df['z_score'] = -df['z_score']
            
            # Normalize to 0-1 range
            min_val = df['z_score'].min()
            max_val = df['z_score'].max()
            if max_val > min_val:
                df['z_score_norm'] = (df['z_score'] - min_val) / (max_val - min_val)
            else:
                df['z_score_norm'] = 0.5  # Default value if all z-scores are identical
            
            results['z_score'] = df[['cell_id', 'z_score_norm']].rename(
                columns={'z_score_norm': f'{context}_z_score'}
            )
        
        # Save context results
        context_results[context] = results
    
    # Combine results across contexts and methods
    combined_scores = None
    
    for method in methods:
        method_dfs = []
        
        for context in contexts:
            if context in context_results and method in context_results[context]:
                method_dfs.append(context_results[context][method])
        
        if method_dfs:
            # Merge all dataframes for this method
            method_df = method_dfs[0]
            for df in method_dfs[1:]:
                method_df = pd.merge(method_df, df, on='cell_id', how='outer')
            
            # Fill NaN with 0
            method_df.fillna(0, inplace=True)
            
            # If this is the first method, initialize combined_scores
            if combined_scores is None:
                combined_scores = method_df
            else:
                # Merge with existing scores
                combined_scores = pd.merge(combined_scores, method_df, on='cell_id', how='outer')
    
    # If no results were generated, return None
    if combined_scores is None:
        logging.warning(f"No trait association scores generated for {cell_type}")
        return None
    
    # Save combined scores
    output_file = os.path.join(cell_type_dir, "trait_association_scores.csv")
    combined_scores.to_csv(output_file, index=False)
    logging.info(f"Trait association scores saved to {output_file}")
    
    return combined_scores 