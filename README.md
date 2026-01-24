# MRDAG
**MR-DAG** is an R framework for inferring directed acyclic graphs (DAGs) among phenotypes using genetic variants as instrumental variables. It combines Two-Stage Least Squares (2SLS) regressions with the **NOTEARS** continuous optimization constraint to enforce acyclicity and recover causal structures.

## Prerequisites

Ensure the following R packages are installed. This workflow relies on `BiocManager` for graph-theoretic dependencies.

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("graph", "RBGL", "Rgraphviz"))
install.packages(c("Matrix", "pcalg", "AER", "expm", "MASS"))
```

## Data Generation

The included simulation engine generates synthetic datasets that mirror the assumptions of Mendelian Randomization (MR). The data generation process is encapsulated in `simulate_mr_data` and follows these steps:

1.  **Instrument Assignment (cis-SNPs):**
    Each phenotype node $X_j$ is assigned a specific set of valid genetic instruments (SNPs) denoted as $S_j$. These represent cis-acting variants.

2.  **Genotype Simulation ($Z$):**
    Genotypes are simulated as binomial variables (0, 1, 2) based on random Minor Allele Frequencies (MAF), representing the genetic anchors.

3.  **Topology Construction ($B_{true}$):**
    A ground-truth weighted adjacency matrix $B_{true}$ is generated. A random permutation ensures the graph is a DAG (no cycles), with edge weights sampled uniformly.

4.  **Structural Equation Modeling:**
    The phenotype matrix $X$ is generated using the following structural equation:
    
    $$X = (Z\Gamma^T + U\Pi + E)(I - B_{true}^T)^{-1}$$
    
    Where:
    * **$Z\Gamma^T$**: The genetic component (valid instruments).
    * **$U\Pi$**: Unobserved confounding effects affecting multiple phenotypes.
    * **$E$**: Random noise.
    * **$(I - B_{true}^T)^{-1}$**: The mixing matrix that propagates causal effects through the network.

### Example Generation

```r
source("run_simulation.R")

# Generate a dataset with 5000 samples and 10 phenotype nodes
sim_data <- simulate_mr_data(
  n = 5000,          # Sample size
  p = 10,            # Number of phenotypes (nodes)
  q_min = 3,         # Min SNPs per phenotype
  q_max = 8,         # Max SNPs per phenotype
  edge_prob = 0.3    # Probability of an edge existing
)

# Extract components
X <- sim_data$X      # Phenotype Matrix (n x p)
Z <- sim_data$Z      # Genotype Matrix (n x q)
S <- sim_data$S      # List of valid instrument indices for each node
B_true <- sim_data$B_true # Ground truth causal matrix
```

## Usage

To estimate the causal graph, use the `MR_DAG` function. This function performs a two-step estimation:
1.  **IV Regression:** Estimates pairwise effects using the provided genetic instruments ($S$).
2.  **NOTEARS Projection:** Projects the initial estimate onto the space of DAGs to resolve cycles and enforce sparsity.

```r
# Run MR-DAG estimation
# lam: Sparsity penalty (L1 regularization)
# w_threshold: Cutoff for small edge weights
B_est <- MR_DAG(X, Z, S, lam = 0.01, w_threshold = 0.1)

# View the estimated adjacency matrix
print(round(B_est, 3))
```

## Interpretation

The output `B_est` is a weighted adjacency matrix of dimensions $p \times p$.

* **Rows vs Columns:** In this implementation, the matrix is oriented as **Target $\leftarrow$ Source**.
    * Row $j$, Column $i$ ($B_{ji}$) represents the causal effect of **Node $i$ on Node $j$**.
* **Values:**
    * **0:** No direct causal link.
    * **Non-zero:** Represents the estimated strength and direction (positive/negative) of the causal effect.

### Comparing Truth vs Estimate

You can validate the model by comparing the ground truth matrix used in generation against the estimated matrix:

```r
# Check a specific edge (e.g., Effect of Node 1 on Node 2)
true_effect <- B_true[2, 1] 
est_effect  <- B_est[2, 1]

cat(sprintf("True Effect (1->2): %.3f\n", true_effect))
cat(sprintf("Est. Effect (1->2): %.3f\n", est_effect))

# Calculate Frobenius Norm (Distance between matrices)
frob_dist <- norm(B_true - B_est, type = "F")
cat(sprintf("Total Matrix Error (Frobenius): %.3f\n", frob_dist))
```
