# 🔖 Cell ID Guide for Researchers

## 🎯 **Your Question: What Are These Cell IDs?**

You're seeing cell_ids like: `2, 3, 9, 11, 13, 14, 19...`

**Answer**: These **ARE** your original cell IDs from your AnnData object!

---

## 🔍 **What's Happening**

### **Your AnnData Structure:**

```python
import scanpy as sc

adata = sc.read_h5ad("your_data.h5ad")
print(adata.obs_names[:20])

# Output: ['0', '1', '2', '3', '9', '11', ...]
# or: [0, 1, 2, 3, 9, 11, ...]  (as integers)
```

**Your AnnData has integer row indices as cell names!**

This is common when:
- Data was subset/filtered (creating non-sequential indices)
- Original processing used integer indices
- Cells were selected from larger dataset

### **The Good News:**

✅ **Cell IDs in results MATCH your AnnData exactly!**

```python
# Result cell_id '2' corresponds to:
adata[adata.obs_names == '2']  # or adata.obs_names[2]

# You can directly map back!
```

---

## 🔧 **Solution 1: Add Meaningful Metadata (RECOMMENDED)**

Create enhanced results with all cell information:

```python
import scanpy as sc
import pandas as pd
import scDCF

# Load your data and results
adata = sc.read_h5ad("data.h5ad")
results = pd.read_csv("c_Mo_disease_monte_carlo_results.csv")

# Create cell metadata mapping
cell_info = adata.obs.copy()
cell_info['original_cell_id'] = adata.obs_names

# Merge with results
results['cell_id_str'] = results['cell_id'].astype(str)
cell_info['cell_id_str'] = cell_info.index.astype(str)

enhanced_results = pd.merge(
    results,
    cell_info,
    on='cell_id_str',
    how='left'
)

# Now you have everything!
enhanced_results.to_csv("results_with_metadata.csv", index=False)

print(enhanced_results[['cell_id', 'p_value', 'celltype', 'sample', 'batch']].head())
```

**Result:**
```csv
cell_id, p_value, celltype, sample, batch, patient_id, ...
2,       0.051,   c_Mo,     S01,    B1,   P001, ...
3,       0.823,   c_Mo,     S01,    B1,   P001, ...
9,       0.756,   c_Mo,     S02,    B1,   P001, ...
```

---

## 🔧 **Solution 2: Use Better Cell Names in AnnData**

If you want more meaningful cell IDs, set them BEFORE running scDCF:

```python
import scanpy as sc

# Load data
adata = sc.read_h5ad("data.h5ad")

# Create meaningful cell names
# Option A: Use barcode (if available)
if 'barcode' in adata.obs.columns:
    adata.obs_names = adata.obs['barcode'].values

# Option B: Combine sample + barcode
adata.obs_names = [
    f"{sample}_{barcode}" 
    for sample, barcode in zip(adata.obs['sample'], adata.obs['barcode'])
]

# Option C: Use existing ID column
if 'cell_id' in adata.obs.columns:
    adata.obs_names = adata.obs['cell_id'].values

# Now run scDCF - results will have meaningful cell_ids!
results = scDCF.monte_carlo_comparison(...)
```

---

## 🔧 **Solution 3: Helper Function (I'll Create This)**

Add a built-in function to scDCF to automatically enhance results:

```python
import scDCF

# Run analysis
results = scDCF.monte_carlo_comparison(...)

# Enhance with metadata
enhanced = scDCF.add_cell_metadata(
    results, 
    adata, 
    metadata_columns=['celltype', 'sample', 'batch', 'patient']
)

# Now enhanced results have all info
enhanced.to_csv("results_enhanced.csv")
```

---

## 📊 **What Cell IDs Actually Are**

### **In Your Case:**

```python
# Your AnnData obs_names:
['0', '1', '2', '3', ..., '345641']  # String integers

# After filtering to c_Mo disease:
['2', '3', '9', '11', ...]  # Subset of above (non-sequential)

# In results:
cell_id: 2, 3, 9, 11, ...  # Same as filtered obs_names!
```

**They ARE consistent!** ✅

The IDs match your `adata.obs_names` exactly. They just happen to be integers because your data uses integer indices.

---

## 🎯 **For Researchers: Best Practice**

### **Recommended Workflow:**

```python
import scDCF
import scanpy as sc
import pandas as pd

# 1. Load data
adata = sc.read_h5ad("data.h5ad")

# 2. OPTIONAL: Set meaningful cell names
if 'barcode' in adata.obs.columns:
    adata.obs_names = adata.obs['barcode']

# 3. Run scDCF
results = scDCF.monte_carlo_comparison(
    adata, "T_cell", "celltype", genes, 
    disease_ctrl, healthy_ctrl, "results/",
    iterations=10
)

# 4. Add metadata for easier interpretation
results_enhanced = results.merge(
    adata.obs,
    left_on='cell_id',
    right_index=True,
    how='left'
)

# 5. Save enhanced results
results_enhanced.to_csv("results_with_metadata.csv")
```

**Result: Results with full cell information!**

---

## 📝 **What I'll Add to Package**

### **New Helper Function:**

```python
def add_cell_metadata(results_df, adata, cell_id_column='cell_id', 
                     metadata_columns=None):
    """
    Add cell metadata from AnnData to scDCF results.
    
    Args:
        results_df: DataFrame from monte_carlo_comparison
        adata: Original AnnData object
        cell_id_column: Column in results with cell IDs
        metadata_columns: Specific columns to add (None = all)
        
    Returns:
        DataFrame with cell metadata added
        
    Example:
        >>> results = scDCF.monte_carlo_comparison(...)
        >>> enhanced = scDCF.add_cell_metadata(results, adata)
        >>> enhanced.to_csv("results_enhanced.csv")
    """
    # Implementation...
```

This will make it **super easy** for researchers to get full cell information!

---

## ✅ **Summary**

### **Your Question:**
> "Can they get cell_id consistent with adata? Currently I see 1, 2, 3, 4, 5"

### **Answer:**

**YES, they already are consistent!** ✅

```
cell_id in results = adata.obs_names (exactly!)
```

The numbers you see (2, 3, 9...) ARE your original obs_names - they just happen to be integers in your dataset.

### **For Researchers:**

**To get more information:**
```python
# Map results back to original data
results['cell_id']  # These match adata.obs_names

# Get full cell info:
cell_info = adata.obs.loc[results['cell_id'].astype(str)]
```

**Or use the helper function I'll create:**
```python
enhanced = scDCF.add_cell_metadata(results, adata)
# Now has: cell_id, p_value, celltype, sample, batch, etc.
```

---

**Would you like me to:**
1. **Create the helper function** `add_cell_metadata()` (10 minutes)?
2. **Add documentation** explaining cell ID mapping?
3. **Both**?

This will make it crystal clear for researchers! 🎯

