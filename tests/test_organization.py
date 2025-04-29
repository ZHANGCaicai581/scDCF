#!/usr/bin/env python
import os
import shutil
import pandas as pd
import json
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def create_test_data():
    """Create test data for organization testing"""
    # Create test directory
    if os.path.exists("test_output"):
        shutil.rmtree("test_output")
    os.makedirs("test_output")
    
    # Create cell metrics files for each cell type
    for cell_type in ["T_cell", "B_cell", "NK_cell"]:
        # Create cell metrics
        cell_metrics = pd.DataFrame({
            "cell_id": [f"cell{i}" for i in range(1, 11)],
            "cell_type": [cell_type] * 10,
            "t_stat": [2.5, 1.8, 2.1, 1.9, 2.3, 2.0, 1.7, 2.2, 2.4, 1.6],
            "p_value": [0.01, 0.05, 0.03, 0.04, 0.02, 0.03, 0.06, 0.02, 0.01, 0.07],
            # Add z-score columns
            "disease_z_score": [0.8, 0.5, 0.7, 0.6, 0.7, 0.6, 0.5, 0.7, 0.8, 0.4],
            "healthy_z_score": [0.2, 0.3, 0.2, 0.3, 0.2, 0.3, 0.4, 0.2, 0.2, 0.5],
            # Add -log10 columns
            "neg_log10_p_disease": [2.0, 1.3, 1.5, 1.4, 1.7, 1.5, 1.2, 1.7, 2.0, 1.1],
            "neg_log10_p_healthy": [0.7, 0.5, 0.6, 0.5, 0.6, 0.5, 0.4, 0.6, 0.7, 0.3],
            # Add trait score columns
            "trait_age_score": [0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7],
            "trait_gender_score": [0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3],
            "trait_severity_score": [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
        })
        
        # Save cell metrics
        os.makedirs(os.path.join("test_output", cell_type), exist_ok=True)
        cell_metrics.to_csv(os.path.join("test_output", cell_type, "cell_metrics.csv"), index=False)
        
        # Create supporting data directory
        supporting_dir = os.path.join("test_output", cell_type, "supporting_data")
        os.makedirs(supporting_dir, exist_ok=True)
        
        # Create Monte Carlo iterations directory
        monte_carlo_dir = os.path.join(supporting_dir, "monte_carlo_iterations")
        os.makedirs(monte_carlo_dir, exist_ok=True)
        
        # Create control genes directory
        control_genes_dir = os.path.join(supporting_dir, "control_genes")
        os.makedirs(control_genes_dir, exist_ok=True)
        
        # Create disease monte carlo results
        monte_carlo_results = pd.DataFrame({
            "cell_id": [f"cell{i}" for i in range(1, 11)],
            "p_value": [0.01, 0.05, 0.03, 0.04, 0.02, 0.03, 0.06, 0.02, 0.01, 0.07],
            "t_stat": [2.5, 1.8, 2.1, 1.9, 2.3, 2.0, 1.7, 2.2, 2.4, 1.6]
        })
        monte_carlo_results.to_csv(os.path.join(supporting_dir, f"{cell_type}_disease_monte_carlo_results.csv"), index=False)
        
        # Create Monte Carlo iterations files (at least 2)
        for i in range(1, 3):
            monte_carlo_iter = pd.DataFrame({
                "cell_id": [f"cell{j}" for j in range(1, 11)],
                "p_value": [0.01, 0.05, 0.03, 0.04, 0.02, 0.03, 0.06, 0.02, 0.01, 0.07],
                "t_stat": [2.5, 1.8, 2.1, 1.9, 2.3, 2.0, 1.7, 2.2, 2.4, 1.6]
            })
            monte_carlo_iter.to_csv(os.path.join(monte_carlo_dir, f"{cell_type}_disease_monte_carlo_results_iteration{i}.csv"), index=False)
        
        # Create trait scores
        trait_scores = pd.DataFrame({
            "trait": ["age", "gender", "severity"],
            "score": [0.7, 0.3, 0.5]
        })
        trait_scores.to_csv(os.path.join(supporting_dir, f"{cell_type}_trait_scores.csv"), index=False)
        
        # Create healthy monte carlo results (if needed)
        healthy_monte_carlo = pd.DataFrame({
            "cell_id": [f"healthy_cell{i}" for i in range(1, 6)],
            "p_value": [0.02, 0.04, 0.03, 0.05, 0.01],
            "t_stat": [1.8, 1.5, 1.7, 1.4, 2.0]
        })
        healthy_monte_carlo.to_csv(os.path.join(supporting_dir, f"{cell_type}_healthy_monte_carlo_results.csv"), index=False)
        
        # Create trait association scores
        trait_assoc = pd.DataFrame({
            "trait": ["age", "gender", "severity"],
            "association": [0.6, 0.4, 0.5]
        })
        trait_assoc.to_csv(os.path.join(supporting_dir, "trait_association_scores.csv"), index=False)
        
        # Create control genes JSON file
        control_genes = {
            "disease": ["GENE1", "GENE2", "GENE3"],
            "healthy": ["GENE4", "GENE5", "GENE6"]
        }
        
        # Save in both locations (subfolder and root for backward compatibility)
        with open(os.path.join(control_genes_dir, f"{cell_type}_control_genes.json"), "w") as f:
            json.dump(control_genes, f, indent=2)
        
        with open(os.path.join("test_output", f"{cell_type}_control_genes.json"), "w") as f:
            json.dump(control_genes, f, indent=2)
    
    logging.info("Created test data in test_output directory")
    return True

def test_organization():
    """Test basic organization functionality"""
    # Create test data
    create_test_data()
    
    # Use the organize_final_output function to organize the data
    try:
        from post_trait_test import organize_final_output
        organize_final_output(source_dir="test_output", dest_dir="test_organized")
    except ImportError:
        logging.error("Could not import organize_final_output from post_trait_test")
        return False
    
    # Run check_cell_metrics.py to verify the output
    logging.info("Running check_cell_metrics.py...")
    os.system("python check_cell_metrics.py --output_dir test_organized")
    
    return True

if __name__ == "__main__":
    success = test_organization()
    if success:
        logging.info("✅ Basic organization test completed successfully")
    else:
        logging.error("❌ Basic organization test failed")