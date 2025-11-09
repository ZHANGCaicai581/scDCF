# Methods (detailed)

This document provides an academic description of the scDCF methodology, aligned with the manuscript.

## GWAS-to-gene prioritization

Let G be the set of genes prioritized by MAGMA (or TWAS), and E the set of genes expressed in the scRNA-seq dataset. Downstream analyses use only
\[ G^\* = G \cap E. \]
MAGMA gene-level statistics (e.g., Z_g) are used as association weights.

## Healthy reference pools matched by library size

For each target cell c, define the absolute library-size difference to a healthy cell r ∈ 𝓗 as
\[ \Delta_L(c, r) = |L_c - L_r|. \]
Rank all healthy cells by \(\Delta_L(c,r)\) and retain the 1000 nearest to form \(\mathcal{H}_{1000}(c)\).
In each Monte Carlo iteration b = 1,…,B, draw without replacement a subset of 100 healthy cells
\[ \mathcal{H}_b(c) \subset \mathcal{H}_{1000}(c),\quad |\mathcal{H}_b(c)|=100, \]
which provides a bootstrap reference for background expression.

## Control-gene matching within cell types

For each g ∈ G*, select a pool 𝓒(g) of 10 control genes matched on expression properties within the same annotated cell type, performed separately for disease and healthy groups.
Let (μ_g^(D), σ_g^{2(D)}) and (μ_g^(H), σ_g^{2(H)}) denote mean and variance of g in disease and healthy groups, respectively. For candidate control gene c, define distances
\[ d^{(D)}(g,c) = \sqrt{(\mu_g^{(D)} - \mu_c^{(D)})^2 + (\sigma_{g}^{2(D)} - \sigma_{c}^{2(D)})^2}, \]
\[ d^{(H)}(g,c) = \sqrt{(\mu_g^{(H)} - \mu_c^{(H)})^2 + (\sigma_{g}^{2(H)} - \sigma_{c}^{2(H)})^2}. \]
When evaluating a disease target cell, d^(D) is used; when evaluating a healthy target cell, d^(H) is used. The 10 nearest controls by the relevant distance are retained in 𝓒(g).

## Per-cell expression deviations (Monte Carlo framework)

For target cell c and prioritized gene g, let x_{g,c} denote expression of g in c, and x_{g,r} the expression in healthy reference cell r ∈ 𝓗_b(c). In iteration b, define the difference
\[ \delta^{(b)}_{g,c} = x_{g,c} - \frac{1}{|\mathcal{H}_b(c)|} \sum_{r\in\mathcal{H}_b(c)} x_{g,r}. \]
For the matched control gene c_b(g) ∈ 𝓒(g) sampled at iteration b, define analogously
\[ \delta^{(b)}_{c_b(g),c} = x_{c_b(g),c} - \frac{1}{|\mathcal{H}_b(c)|} \sum_{r\in\mathcal{H}_b(c)} x_{c_b(g),r}. \]

## Difference-of-differences and MAGMA weighting

The disease-relevant signal is isolated via a difference-of-differences statistic
\[ \Delta^{(b)}_{g,c} = \delta^{(b)}_{g,c} - \delta^{(b)}_{c_b(g),c}, \quad g\in G^\*. \]
Each gene-wise difference is weighted by the MAGMA Z-score:
\[ \widetilde{\Delta}^{(b)}_{g,c} = Z_g\, \Delta^{(b)}_{g,c}. \]
Aggregate across prioritized genes to obtain a per-cell statistic in iteration b:
\[ S^{(b)}_{c} = \frac{1}{|G^\*|} \sum_{g\in G^\*} \widetilde{\Delta}^{(b)}_{g,c}. \]

## Meta-analysis across iterations and cell-level calling

Within iteration b, compute the gene-wise variance about S^{(b)}_c:
\[ \widehat{\sigma}^{2(b)}_c = \frac{1}{|G^\*|-1} \sum_{g\in G^\*} ( \widetilde{\Delta}^{(b)}_{g,c} - S^{(b)}_c )^2, \quad \mathrm{SE}^{(b)}_c = \sqrt{\frac{\widehat{\sigma}^{2(b)}_c}{|G^\*|}}. \]
Standardize to form Z-statistics and one-sided p-values:
\[ Z^{(b)}_c = \frac{S^{(b)}_c}{\mathrm{SE}^{(b)}_c}, \qquad p^{(b)}_c = 1 - \Phi(Z^{(b)}_c). \]
Combine evidence across iterations using Fisher’s method:
\[ X_c = -2 \sum_{b=1}^{B} \log p^{(b)}_c \;\sim\; \chi^2_{2B} \quad (H_0). \]
Adjust meta-analytic p-values across all tested cells via Benjamini–Hochberg to control FDR; call cells with FDR < α (e.g., 0.05) as disease-associated.

## Cell-type-level enrichment

For each annotated cell type, form a 2×2 contingency table of (disease-associated vs. not) × (patient vs. control). Enrichment is assessed using Fisher’s exact test, declaring enrichment when the patient proportion exceeds the control proportion with FDR-controlled significance.

## Practical defaults and notes

- Reference pool size: 1000 nearest healthy cells; per-iteration reference subset: 100 healthy cells.
- Iteration count: user-controlled via CLI flag `--iterations`.
- Control-gene matching: 10 controls per prioritized gene, matched within cell type and group.
- Approximate independence across iterations via independent resampling of references and control selections.
