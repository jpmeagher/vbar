
<!-- README.md is generated from README.Rmd. Please edit that file -->

# vbar

<!-- badges: start -->
<!-- badges: end -->

The goal of `vbar` is to implement Variational Bayes Ancestral
Reconstruction for collections of discrete and continuous phenotypic
traits given a phylogeny of evolutionary relationships between species.
Discrete traits may be ordinal or nominal, while continuous traits may
be scalar- or function-valued. Ancestral Reconstruction is based on the
Phylogenetic Latent Variable Model (PLVM) for trait evolution. This
model accommodates repeated measurements for extant species.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jpmeagher/vbar")
```

## Synthetic Example

Consider a synthetic example. Its analysis requires the following
packages.

``` r
library(vbar)
library(ape)
library(ggplot2)
```

We will reconstruct
![P = 4](https://latex.codecogs.com/png.latex?P%20%3D%204 "P = 4")
manifest traits for the common ancestors of
![S = 128](https://latex.codecogs.com/png.latex?S%20%3D%20128 "S = 128")
extant taxa, given
![N\_s = 3](https://latex.codecogs.com/png.latex?N_s%20%3D%203 "N_s = 3")
recordings of these traits for each extant species
![s = 1, \\dots S](https://latex.codecogs.com/png.latex?s%20%3D%201%2C%20%5Cdots%20S "s = 1, \dots S")
such that we have recordings for
![N = 384](https://latex.codecogs.com/png.latex?N%20%3D%20384 "N = 384")
individuals in total. These traits will consist of one ordinal, one
nominal, one scalar-valued, and one function-valued trait. Our first
task is to simulate these traits from the PLVM generative model.

### Phylogenetic Tree

The shared evolutionary history between taxa is modeled as a known,
fixed phylogeny
![\\mathcal T = \\left\\{ \\mathcal V, \\mathcal B \\right\\}](https://latex.codecogs.com/png.latex?%5Cmathcal%20T%20%3D%20%5Cleft%5C%7B%20%5Cmathcal%20V%2C%20%5Cmathcal%20B%20%5Cright%5C%7D "\mathcal T = \left\{ \mathcal V, \mathcal B \right\}").
![\\mathcal T](https://latex.codecogs.com/png.latex?%5Cmathcal%20T "\mathcal T")
is a bifurcating, directed, acyclic graph tree with nodes
![\\mathcal V = \\left\\{ v\_1, \\dots, v\_{2S - 1} \\right\\}](https://latex.codecogs.com/png.latex?%5Cmathcal%20V%20%3D%20%5Cleft%5C%7B%20v_1%2C%20%5Cdots%2C%20v_%7B2S%20-%201%7D%20%5Cright%5C%7D "\mathcal V = \left\{ v_1, \dots, v_{2S - 1} \right\}")
linked by branches
![\\mathcal B = \\left\\{ b\_1, \\dots, b\_{2S - 2} \\right\\}](https://latex.codecogs.com/png.latex?%5Cmathcal%20B%20%3D%20%5Cleft%5C%7B%20b_1%2C%20%5Cdots%2C%20b_%7B2S%20-%202%7D%20%5Cright%5C%7D "\mathcal B = \left\{ b_1, \dots, b_{2S - 2} \right\}").
The graph originates at the degree-2 root node
![v\_{2S - 1}](https://latex.codecogs.com/png.latex?v_%7B2S%20-%201%7D "v_{2S - 1}")
and terminates at the degree-1 terminal nodes
![v\_1. \\dots, v\_{S}](https://latex.codecogs.com/png.latex?v_1.%20%5Cdots%2C%20v_%7BS%7D "v_1. \dots, v_{S}")
corresponding to the extant taxa. Node
![v\_s](https://latex.codecogs.com/png.latex?v_s "v_s") is linked to its
parent
![v\_{\\operatorname{pa} \\left( s \\right)}](https://latex.codecogs.com/png.latex?v_%7B%5Coperatorname%7Bpa%7D%20%5Cleft%28%20s%20%5Cright%29%7D "v_{\operatorname{pa} \left( s \right)}")
by the branch of length
![b\_s](https://latex.codecogs.com/png.latex?b_s "b_s") for all
![s = 1, \\dots, 2S - 2](https://latex.codecogs.com/png.latex?s%20%3D%201%2C%20%5Cdots%2C%202S%20-%202 "s = 1, \dots, 2S - 2").
For notational convenience, we let
![\\boldsymbol t \\in \\mathcal T](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20t%20%5Cin%20%5Cmathcal%20T "\boldsymbol t \in \mathcal T")
denote a continuous position along a branch of
![\\mathcal T](https://latex.codecogs.com/png.latex?%5Cmathcal%20T "\mathcal T"),
where
![\\boldsymbol t\_s](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20t_s "\boldsymbol t_s")
corresponds to node
![v\_s](https://latex.codecogs.com/png.latex?v_s "v_s") for all
![s = 1, \\dots, 2S - 1](https://latex.codecogs.com/png.latex?s%20%3D%201%2C%20%5Cdots%2C%202S%20-%201 "s = 1, \dots, 2S - 1").
Furthermore, we define the patristic distance operator
![d\_{\\mathcal T} \\left( \\boldsymbol t, \\boldsymbol t' \\right)](https://latex.codecogs.com/png.latex?d_%7B%5Cmathcal%20T%7D%20%5Cleft%28%20%5Cboldsymbol%20t%2C%20%5Cboldsymbol%20t%27%20%5Cright%29 "d_{\mathcal T} \left( \boldsymbol t, \boldsymbol t' \right)"),
which is the the shortest path from
![\\boldsymbol t](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20t "\boldsymbol t")
to
![\\boldsymbol t'](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20t%27 "\boldsymbol t'")
over
![\\mathcal T](https://latex.codecogs.com/png.latex?%5Cmathcal%20T "\mathcal T").

``` r
S <- 2^7
set.seed(98)
phy <- rcoal(S)
phy <- scale_phylo(phy, max_dist = 1)
#> Loading required namespace: testthat
```

<img src="man/figures/README-plot_phylogeny-1.png" width="100%" />

### Phylogenetic Latent Variables

Given
![\\mathcal T](https://latex.codecogs.com/png.latex?%5Cmathcal%20T "\mathcal T"),
we assume that the evolution of manifest traits is driven by
![L = 3](https://latex.codecogs.com/png.latex?L%20%3D%203 "L = 3")
latent traits, modelled as independent Ornstein-Uhlenbeck (OU) processes
on
![\\mathcal T](https://latex.codecogs.com/png.latex?%5Cmathcal%20T "\mathcal T").
These are the phylogenetic latent variables. For
![\\boldsymbol t \\in \\mathcal{T}](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20t%20%5Cin%20%5Cmathcal%7BT%7D "\boldsymbol t \in \mathcal{T}"),
we have

![
\\begin{aligned}
f\_l \\left( \\boldsymbol t \\right) &\\sim \\mathcal{GP} \\left( 0, k\_l \\left( \\boldsymbol t, \\boldsymbol t' \\right)\\right), \\\\
k\_l \\left( \\boldsymbol t, \\boldsymbol t' \\right) &=  h\_l^2 \\exp \\left( - \\frac{d\_{\\mathcal T} \\left( \\boldsymbol t, \\boldsymbol t' \\right) }{\\ell} \\right) + \\left( 1 - h\_l^2\\right) \\delta \\left( \\boldsymbol t \\in \\left\\{ \\boldsymbol t\_1, \\dots, \\boldsymbol t\_S\\right\\}\\right)
\\end{aligned}
](https://latex.codecogs.com/png.latex?%0A%5Cbegin%7Baligned%7D%0Af_l%20%5Cleft%28%20%5Cboldsymbol%20t%20%5Cright%29%20%26%5Csim%20%5Cmathcal%7BGP%7D%20%5Cleft%28%200%2C%20k_l%20%5Cleft%28%20%5Cboldsymbol%20t%2C%20%5Cboldsymbol%20t%27%20%5Cright%29%5Cright%29%2C%20%5C%5C%0Ak_l%20%5Cleft%28%20%5Cboldsymbol%20t%2C%20%5Cboldsymbol%20t%27%20%5Cright%29%20%26%3D%20%20h_l%5E2%20%5Cexp%20%5Cleft%28%20-%20%5Cfrac%7Bd_%7B%5Cmathcal%20T%7D%20%5Cleft%28%20%5Cboldsymbol%20t%2C%20%5Cboldsymbol%20t%27%20%5Cright%29%20%7D%7B%5Cell%7D%20%5Cright%29%20%2B%20%5Cleft%28%201%20-%20h_l%5E2%5Cright%29%20%5Cdelta%20%5Cleft%28%20%5Cboldsymbol%20t%20%5Cin%20%5Cleft%5C%7B%20%5Cboldsymbol%20t_1%2C%20%5Cdots%2C%20%5Cboldsymbol%20t_S%5Cright%5C%7D%5Cright%29%0A%5Cend%7Baligned%7D%0A "
\begin{aligned}
f_l \left( \boldsymbol t \right) &\sim \mathcal{GP} \left( 0, k_l \left( \boldsymbol t, \boldsymbol t' \right)\right), \\
k_l \left( \boldsymbol t, \boldsymbol t' \right) &=  h_l^2 \exp \left( - \frac{d_{\mathcal T} \left( \boldsymbol t, \boldsymbol t' \right) }{\ell} \right) + \left( 1 - h_l^2\right) \delta \left( \boldsymbol t \in \left\{ \boldsymbol t_1, \dots, \boldsymbol t_S\right\}\right)
\end{aligned}
")

for
![l = 1, \\dots, L](https://latex.codecogs.com/png.latex?l%20%3D%201%2C%20%5Cdots%2C%20L "l = 1, \dots, L"),
where
![h\_l^2 \\in \\left( 0, 1 \\right)](https://latex.codecogs.com/png.latex?h_l%5E2%20%5Cin%20%5Cleft%28%200%2C%201%20%5Cright%29 "h_l^2 \in \left( 0, 1 \right)")
is the heritability of the
![l^{th}](https://latex.codecogs.com/png.latex?l%5E%7Bth%7D "l^{th}")
latent trait,
![\\ell = 2](https://latex.codecogs.com/png.latex?%5Cell%20%3D%202 "\ell = 2")
is the fixed phylogenetic length-scale, and
![\\delta \\left( \\cdot \\right)](https://latex.codecogs.com/png.latex?%5Cdelta%20%5Cleft%28%20%5Ccdot%20%5Cright%29 "\delta \left( \cdot \right)")
is the indicator function.

The OU process is a Gauss-Markov process, a property we will exploit to
develop a computationally efficient inference scheme. If we denote
![f\_{s,l} = f\_l \\left( \\boldsymbol t\_s \\right)](https://latex.codecogs.com/png.latex?f_%7Bs%2Cl%7D%20%3D%20f_l%20%5Cleft%28%20%5Cboldsymbol%20t_s%20%5Cright%29 "f_{s,l} = f_l \left( \boldsymbol t_s \right)"),
then

![
\\begin{aligned}
f\_{2S - 1, l} &\\sim \\mathcal N \\left( 0, h^2 \\right)\\\\
f\_{s,l} \\mid f\_{\\operatorname{pa} \\left( s \\right),l}, h &\\sim \\mathcal N \\left( \\nu\_{s, l} f\_{\\operatorname{pa} \\left( s \\right),l}, \\eta\_{s, l}^2 \\right)
\\end{aligned}
](https://latex.codecogs.com/png.latex?%0A%5Cbegin%7Baligned%7D%0Af_%7B2S%20-%201%2C%20l%7D%20%26%5Csim%20%5Cmathcal%20N%20%5Cleft%28%200%2C%20h%5E2%20%5Cright%29%5C%5C%0Af_%7Bs%2Cl%7D%20%5Cmid%20f_%7B%5Coperatorname%7Bpa%7D%20%5Cleft%28%20s%20%5Cright%29%2Cl%7D%2C%20h%20%26%5Csim%20%5Cmathcal%20N%20%5Cleft%28%20%5Cnu_%7Bs%2C%20l%7D%20f_%7B%5Coperatorname%7Bpa%7D%20%5Cleft%28%20s%20%5Cright%29%2Cl%7D%2C%20%5Ceta_%7Bs%2C%20l%7D%5E2%20%5Cright%29%0A%5Cend%7Baligned%7D%0A "
\begin{aligned}
f_{2S - 1, l} &\sim \mathcal N \left( 0, h^2 \right)\\
f_{s,l} \mid f_{\operatorname{pa} \left( s \right),l}, h &\sim \mathcal N \left( \nu_{s, l} f_{\operatorname{pa} \left( s \right),l}, \eta_{s, l}^2 \right)
\end{aligned}
")

where, given
![k\_{s, s', l} = k\_l \\left( \\boldsymbol t\_s, \\boldsymbol t\_{s'} \\right)](https://latex.codecogs.com/png.latex?k_%7Bs%2C%20s%27%2C%20l%7D%20%3D%20k_l%20%5Cleft%28%20%5Cboldsymbol%20t_s%2C%20%5Cboldsymbol%20t_%7Bs%27%7D%20%5Cright%29 "k_{s, s', l} = k_l \left( \boldsymbol t_s, \boldsymbol t_{s'} \right)")
and
![k\_{s, l} = k\_l \\left( \\boldsymbol t\_s, \\boldsymbol t\_{s} \\right)](https://latex.codecogs.com/png.latex?k_%7Bs%2C%20l%7D%20%3D%20k_l%20%5Cleft%28%20%5Cboldsymbol%20t_s%2C%20%5Cboldsymbol%20t_%7Bs%7D%20%5Cright%29 "k_{s, l} = k_l \left( \boldsymbol t_s, \boldsymbol t_{s} \right)")

![
\\begin{aligned}
\\nu\_{s, l} &= k\_{s, \\operatorname{pa} \\left( s \\right),l} k\_{\\operatorname{pa} \\left( s \\right), l}^{-1}, \\\\
&= \\exp\\left( -\\frac{d\_{\\mathcal T} \\left( \\boldsymbol t\_s, \\boldsymbol t\_{\\operatorname{pa} \\left( s \\right)} \\right) }{\\ell}\\right), \\\\
\\eta\_{s, l}^2 &= k\_{s, l} - \\nu\_{s, l} \\, k\_{\\operatorname{pa} \\left( s \\right), s, l}, \\\\
&= \\frac{\\ell h\_l^2}{2} \\left( 1 - \\exp \\left( - 2 \\frac{d\_{\\mathcal T} \\left( \\boldsymbol t\_s, \\boldsymbol t\_{\\operatorname{pa} \\left( s \\right)} \\right) }{\\ell} \\right) \\right) + \\left( 1 - h\_l^2\\right) \\delta \\left( \\boldsymbol t \\in \\left\\{ \\boldsymbol t\_1, \\dots, \\boldsymbol t\_S\\right\\}\\right).
\\end{aligned}
](https://latex.codecogs.com/png.latex?%0A%5Cbegin%7Baligned%7D%0A%5Cnu_%7Bs%2C%20l%7D%20%26%3D%20k_%7Bs%2C%20%5Coperatorname%7Bpa%7D%20%5Cleft%28%20s%20%5Cright%29%2Cl%7D%20k_%7B%5Coperatorname%7Bpa%7D%20%5Cleft%28%20s%20%5Cright%29%2C%20l%7D%5E%7B-1%7D%2C%20%5C%5C%0A%26%3D%20%5Cexp%5Cleft%28%20-%5Cfrac%7Bd_%7B%5Cmathcal%20T%7D%20%5Cleft%28%20%5Cboldsymbol%20t_s%2C%20%5Cboldsymbol%20t_%7B%5Coperatorname%7Bpa%7D%20%5Cleft%28%20s%20%5Cright%29%7D%20%5Cright%29%20%7D%7B%5Cell%7D%5Cright%29%2C%20%5C%5C%0A%5Ceta_%7Bs%2C%20l%7D%5E2%20%26%3D%20k_%7Bs%2C%20l%7D%20-%20%5Cnu_%7Bs%2C%20l%7D%20%5C%2C%20k_%7B%5Coperatorname%7Bpa%7D%20%5Cleft%28%20s%20%5Cright%29%2C%20s%2C%20l%7D%2C%20%5C%5C%0A%26%3D%20%5Cfrac%7B%5Cell%20h_l%5E2%7D%7B2%7D%20%5Cleft%28%201%20-%20%5Cexp%20%5Cleft%28%20-%202%20%5Cfrac%7Bd_%7B%5Cmathcal%20T%7D%20%5Cleft%28%20%5Cboldsymbol%20t_s%2C%20%5Cboldsymbol%20t_%7B%5Coperatorname%7Bpa%7D%20%5Cleft%28%20s%20%5Cright%29%7D%20%5Cright%29%20%7D%7B%5Cell%7D%20%5Cright%29%20%5Cright%29%20%2B%20%5Cleft%28%201%20-%20h_l%5E2%5Cright%29%20%5Cdelta%20%5Cleft%28%20%5Cboldsymbol%20t%20%5Cin%20%5Cleft%5C%7B%20%5Cboldsymbol%20t_1%2C%20%5Cdots%2C%20%5Cboldsymbol%20t_S%5Cright%5C%7D%5Cright%29.%0A%5Cend%7Baligned%7D%0A "
\begin{aligned}
\nu_{s, l} &= k_{s, \operatorname{pa} \left( s \right),l} k_{\operatorname{pa} \left( s \right), l}^{-1}, \\
&= \exp\left( -\frac{d_{\mathcal T} \left( \boldsymbol t_s, \boldsymbol t_{\operatorname{pa} \left( s \right)} \right) }{\ell}\right), \\
\eta_{s, l}^2 &= k_{s, l} - \nu_{s, l} \, k_{\operatorname{pa} \left( s \right), s, l}, \\
&= \frac{\ell h_l^2}{2} \left( 1 - \exp \left( - 2 \frac{d_{\mathcal T} \left( \boldsymbol t_s, \boldsymbol t_{\operatorname{pa} \left( s \right)} \right) }{\ell} \right) \right) + \left( 1 - h_l^2\right) \delta \left( \boldsymbol t \in \left\{ \boldsymbol t_1, \dots, \boldsymbol t_S\right\}\right).
\end{aligned}
")

for
![s = 1, \\dots 2S - 1](https://latex.codecogs.com/png.latex?s%20%3D%201%2C%20%5Cdots%202S%20-%201 "s = 1, \dots 2S - 1")
and
![d\_{\\mathcal T} \\left( \\boldsymbol t\_s, \\boldsymbol t\_{\\operatorname{pa} \\left( s \\right)} \\right) = b\_s](https://latex.codecogs.com/png.latex?d_%7B%5Cmathcal%20T%7D%20%5Cleft%28%20%5Cboldsymbol%20t_s%2C%20%5Cboldsymbol%20t_%7B%5Coperatorname%7Bpa%7D%20%5Cleft%28%20s%20%5Cright%29%7D%20%5Cright%29%20%3D%20b_s "d_{\mathcal T} \left( \boldsymbol t_s, \boldsymbol t_{\operatorname{pa} \left( s \right)} \right) = b_s").

In this example, we let
![\\boldsymbol h = \\left( h\_1, \\dots, h\_L \\right)^\\top = \\left( \\sqrt{0.95}, \\sqrt{0.66}, \\sqrt{0.25} \\right)^\\top](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20h%20%3D%20%5Cleft%28%20h_1%2C%20%5Cdots%2C%20h_L%20%5Cright%29%5E%5Ctop%20%3D%20%5Cleft%28%20%5Csqrt%7B0.95%7D%2C%20%5Csqrt%7B0.66%7D%2C%20%5Csqrt%7B0.25%7D%20%5Cright%29%5E%5Ctop "\boldsymbol h = \left( h_1, \dots, h_L \right)^\top = \left( \sqrt{0.95}, \sqrt{0.66}, \sqrt{0.25} \right)^\top"),
such that the heritability of latent traits ranges from high to low.

``` r
h <- sqrt(c(0.95, 0.66, 0.25))
latent_taxa_traits <- sapply(
  h, function(x) {
    simulate_phylogenetic_ou(
      phy = phy, 
      heritable_amplitude = x, length_scale = 2,
      environmental_amplitude = sqrt(1 - x^2),
      internal = TRUE
    )
  }
)
```
