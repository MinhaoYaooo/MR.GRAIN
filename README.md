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

## Data Generation

In this section, the data generation process is explicitly defined below. We simulate a network Mendelian randomization setting where:
* **n = 5000**: Sample size.
* **p = 10**: Number of phenotypes (e.g., protein levels).
* **q = 40**: Total genetic variants (instruments), with exactly **4 distinct SNPs** acting as valid instruments for each phenotype.

The simulation accounts for pleiotropy via unobserved confounders ($U$) and generates phenotypes ($X$) based on the structural equation $X = (Z\Gamma^T + U\Pi + E)(I - B^T)^{-1}$.

### Example Code

You can run this snippet directly to generate the simulation data components (`X`, `Z`, `S`, `B_true`) and visualize the true causal graph.

```r
library(Rgraphviz)
library(graph)

# --- 1. Parameters ---
set.seed(123)      # For reproducibility
n <- 5000          # Sample size
p <- 10            # Number of nodes
snps_per_node <- 4 # Fixed SNPs per node
q <- p * snps_per_node # Total SNPs = 40

# --- 2. Instruments (S) & Gamma ---
# S is a list where S[[j]] contains indices of valid instruments for node j
S <- vector("list", p)
Gamma <- matrix(0, p, q) # Effect of SNPs on phenotypes (p x q)

curr_idx <- 1
for (j in 1:p) {
  # Assign block of 4 SNPs to each node
  S[[j]] <- curr_idx:(curr_idx + snps_per_node - 1)
  
  # Set instrument strength (Gamma) for valid SNPs
  Gamma[j, S[[j]]] <- runif(snps_per_node, 0.2, 0.4)
  curr_idx <- curr_idx + snps_per_node
}

# --- 3. Genotypes (Z) ---
# Simulate genotypes as Binomial(2, MAF)
Z <- matrix(rbinom(n * q, 2, 0.3), nrow = n, ncol = q)

# --- 4. Causal Topology (B_true) ---
# Fixed DAG structure (Target, Source)
# B[j, i] != 0 implies edge i -> j
B_true <- matrix(0, p, p)

# Define specific causal edges (Cascade: 1->2->3->... and some branches)
B_true[2, 1] <-  0.5   # Node 1 -> Node 2
B_true[3, 2] <- -0.4   # Node 2 -> Node 3
B_true[4, 2] <-  0.3   # Node 2 -> Node 4
B_true[5, 3] <-  0.4   # Node 3 -> Node 5
B_true[6, 1] <- -0.3   # Node 1 -> Node 6
B_true[7, 5] <-  0.5   # Node 5 -> Node 7
B_true[8, 6] <-  0.4   # Node 6 -> Node 8
B_true[9, 8] <- -0.4   # Node 8 -> Node 9
B_true[10,9] <-  0.3   # Node 9 -> Node 10

# Visualize the Ground Truth DAG
g <- new("graphAM", adjMat = t(B_true) != 0, edgemode = "directed")
plot(g, attrs = list(node = list(fillcolor = "lightblue", shape = "circle")))

# --- 5. Phenotypes (X) ---
# Model: X = (Z*Gamma' + U*Pi + E) * (I - B')^-1
# Add unobserved confounders (U) and noise (E)
r <- 3 # Number of hidden confounders
U <- matrix(rnorm(n * r), n, r)
Pi <- matrix(runif(r * p, 0.3, 0.6), r, p) # Confounder effects
E <- matrix(rnorm(n * p), n, p)

I_minus_B_T <- diag(p) - t(B_true)
RHS <- (Z %*% t(Gamma)) + (U %*% Pi) + E

# Solve structural equation
X <- RHS %*% solve(I_minus_B_T)

# Data is now ready: X (Phenotypes), Z (Genotypes), S (Instrument Map)
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
