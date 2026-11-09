# Attribute likelihood weights and non-Gaussian attributes

This experiment reproduces Supplementary Sections S4 and S8. Each condition uses ten independent data replicates. Comparisons across conditions sharing graphs or attribute draws are paired by replicate.

## Common settings

The undirected, loop-free SBM has 1,000 nodes, five communities of 200 nodes, and expected average degree 20. The structural ratio is $r=p_{\mathrm{in}}/p_{\mathrm{out}}$, with $p_{\mathrm{out}}=20/(199r+800)$. The attributes have 50 dimensions. The coordinates of each community mean are independent standard normal draws.

For each input, Louvain estimates $K$ at resolution 1. Attribute K-means uses that count (at least two), `nstart=5`, and `iter.max=10`. The K-means partition is the estimated-$K$ reference and initializes aSBM-Louvain and aSBM. Both fits start with the same post-initialization RNG state. The Louvain reference is the partition used to estimate $K$.

aSBM-Louvain uses initial annealing strength 1, cooling factor 0.8, relative objective tolerance 0.001, absolute modularity tolerance 0.0001, and an unchanged community count for three successive updates. The outer-iteration budget is 1,000 and the per-fit soft time limit is 300 seconds. The aSBM baseline uses objective tolerance 1 and the compiled correction routine. Timings are diagnostic.

## Likelihood weights

The weighted variational objective is $L_A+\beta L_X-\max(\tau,10^{-5})\operatorname{KL}(q\Vert\pi)$. The structural likelihood counts each undirected pair once. The conditional update for node $i$ uses the full sum of messages from its other nodes and scales its attribute log-likelihood by $\beta$. Both updates use the same attribute weight. The $\beta=0$ arm still uses attributes in initialization.

The grid crosses $r\in\{1.5,4\}$, RMS noise $\sigma\in\{1.5,4\}$, high-to-low community noise ratio $h\in\{1,3\}$, and $\beta\in\{0,0.25,0.5,1,2,4\}$. Two communities are designated high-noise in each replicate. Their scale multipliers $a_k$ are $h$ and the other three are 1. Each community receives noise SD $\sigma a_k/\sqrt{\operatorname{mean}(a_k^2)}$. This grid gives 480 aSBM-Louvain fits and 240 reference results.

## Non-Gaussian attributes

The grid crosses $r\in\{1.5,4\}$ and $\sigma\in\{1.5,4\}$ at $h=1$. The fitted attribute model stays Gaussian. For data generation, add $\sigma E$ to the planted community mean using one of four residual distributions, each with population covariance $I$.

- Gaussian uses a standard normal vector $Z$.
- Student-$t_5$ uses $Z/\sqrt{U/3}$, where $U\sim\chi^2_5$ is shared across the coordinates of each node.
- Laplace uses independent Laplace coordinates with scale $1/\sqrt{2}$, obtained by mapping each Gaussian coordinate through its CDF and the inverse Laplace CDF.
- Contaminated Gaussian uses $SZ/\sqrt{3.4}$, with a node-level multiplier $S$ equal to 5 with probability 0.1 and 1 otherwise.

The four methods are aSBM-Louvain at $\beta=1$, aSBM, Louvain, and K-means. The 640 method records share 160 Gaussian reference records with the likelihood-weight study. The combined run yields 1,200 unique records across 200 condition-replicate pairs. The simulated graphs are shared across settings with the same $r$ and replicate. Community means and base normal and tail draws are shared across conditions within a replicate.

## Seeds and output

For replicate `b` and structural ratio index `i` in `c(1.5, 4)`, the seeds are:

| Draw | Seed |
| --- | --- |
| Graph | `92210000 + 100000*i + b` |
| Means | `92211000 + b` |
| Base residuals | `92212000 + b` |
| High-noise communities | `92213000 + b` |
| Student-t scales | `92214000 + b` |
| Contamination scales | `92215000 + b` |
| Louvain | `92700000 + 1000*i + b` |
| K-means | `92800000 + b` |

`run.R` writes per-fit metrics, inputs, and checkpoints under `outputs/weight_misspecification`. `analyze.py` summarizes the complete grid, including mean and sample SD of ARI and paired differences. Use `WM_OUTPUT` to override the output directory for both scripts. `WM_WORKERS` sets the number of fitting workers (default 4). `WM_SMOKE=1` runs one replicate of a reduced grid and should use a separate output directory.
