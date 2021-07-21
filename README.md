
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

## The Generalised Phylogenetic Latent Variable Model

Consider the
![D](https://latex.codecogs.com/png.latex?D "D")-dimensional observation
of ![P](https://latex.codecogs.com/png.latex?P "P") discrete- or
continuous-valued manifest traits
![\\boldsymbol y\_{sn} = \\left( y\_{sn 1}, \\dots, y\_{snD} \\right)^\\top](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20y_%7Bsn%7D%20%3D%20%5Cleft%28%20y_%7Bsn%201%7D%2C%20%5Cdots%2C%20y_%7BsnD%7D%20%5Cright%29%5E%5Ctop "\boldsymbol y_{sn} = \left( y_{sn 1}, \dots, y_{snD} \right)^\top"),
specific to individual
![n = 1, \\dots, N\_s](https://latex.codecogs.com/png.latex?n%20%3D%201%2C%20%5Cdots%2C%20N_s "n = 1, \dots, N_s")
within taxon
![s = 1, \\dots, S](https://latex.codecogs.com/png.latex?s%20%3D%201%2C%20%5Cdots%2C%20S "s = 1, \dots, S"),
where traits may be ordinal-, nominal-, scalar- or function-valued.

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
denotes the position of node
![v\_s](https://latex.codecogs.com/png.latex?v_s "v_s") for all
![s = 1, \\dots, 2S - 1](https://latex.codecogs.com/png.latex?s%20%3D%201%2C%20%5Cdots%2C%202S%20-%201 "s = 1, \dots, 2S - 1").
Furthermore, the patristic distance operator
![d\_{\\mathcal T} \\left( \\boldsymbol t, \\boldsymbol t' \\right)](https://latex.codecogs.com/png.latex?d_%7B%5Cmathcal%20T%7D%20%5Cleft%28%20%5Cboldsymbol%20t%2C%20%5Cboldsymbol%20t%27%20%5Cright%29 "d_{\mathcal T} \left( \boldsymbol t, \boldsymbol t' \right)")
defines the shortest path from
![\\boldsymbol t](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20t "\boldsymbol t")
to
![\\boldsymbol t'](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20t%27 "\boldsymbol t'")
over
![\\mathcal T](https://latex.codecogs.com/png.latex?%5Cmathcal%20T "\mathcal T").

For each manifest trait, we introduce auxiliary traits
![\\boldsymbol x\_{sn} = \\left( x\_{s n 1}, \\dots, x\_{sn D'} \\right)^\\top \\in \\mathbb R^{D'}](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20x_%7Bsn%7D%20%3D%20%5Cleft%28%20x_%7Bs%20n%201%7D%2C%20%5Cdots%2C%20x_%7Bsn%20D%27%7D%20%5Cright%29%5E%5Ctop%20%5Cin%20%5Cmathbb%20R%5E%7BD%27%7D "\boldsymbol x_{sn} = \left( x_{s n 1}, \dots, x_{sn D'} \right)^\top \in \mathbb R^{D'}")
such that
![D' \\geq D](https://latex.codecogs.com/png.latex?D%27%20%5Cgeq%20D "D' \geq D"),
there exists a deterministic map
![g : \\boldsymbol x\_{sn} \\to \\boldsymbol y\_{sn}](https://latex.codecogs.com/png.latex?g%20%3A%20%5Cboldsymbol%20x_%7Bsn%7D%20%5Cto%20%5Cboldsymbol%20y_%7Bsn%7D "g : \boldsymbol x_{sn} \to \boldsymbol y_{sn}"),
and auxiliary traits mapping
![y\_{snd}](https://latex.codecogs.com/png.latex?y_%7Bsnd%7D "y_{snd}")
are indexed by ![d'](https://latex.codecogs.com/png.latex?d%27 "d'").
This allows us to adopt a probit model for discrete traits.

If
![y\_{snd} = k](https://latex.codecogs.com/png.latex?y_%7Bsnd%7D%20%3D%20k "y_{snd} = k")
for
![k = 1, \\dots, K](https://latex.codecogs.com/png.latex?k%20%3D%201%2C%20%5Cdots%2C%20K "k = 1, \dots, K")
is a nominal trait belonging to one of
![K](https://latex.codecogs.com/png.latex?K "K") **unordered**
categories, we have a
![K](https://latex.codecogs.com/png.latex?K "K")-to-1 mapping such that

![
y\_{snd} = g \\left( x\_{s, n, d'} , \\dots, x\_{s, n,d' + K - 1} \\right) = c\_{k}, \\text{ if } x\_{s, n, d' + k - 1} = \\max \\left\\{ x\_{s, n, d'} , \\dots, x\_{s, n, d' + K - 1} \\right\\}.
](https://latex.codecogs.com/png.latex?%0Ay_%7Bsnd%7D%20%3D%20g%20%5Cleft%28%20x_%7Bs%2C%20n%2C%20d%27%7D%20%2C%20%5Cdots%2C%20x_%7Bs%2C%20n%2Cd%27%20%2B%20K%20-%201%7D%20%5Cright%29%20%3D%20c_%7Bk%7D%2C%20%5Ctext%7B%20if%20%7D%20x_%7Bs%2C%20n%2C%20d%27%20%2B%20k%20-%201%7D%20%3D%20%5Cmax%20%5Cleft%5C%7B%20x_%7Bs%2C%20n%2C%20d%27%7D%20%2C%20%5Cdots%2C%20x_%7Bs%2C%20n%2C%20d%27%20%2B%20K%20-%201%7D%20%5Cright%5C%7D.%0A "
y_{snd} = g \left( x_{s, n, d'} , \dots, x_{s, n,d' + K - 1} \right) = c_{k}, \text{ if } x_{s, n, d' + k - 1} = \max \left\{ x_{s, n, d'} , \dots, x_{s, n, d' + K - 1} \right\}.
")

If
![y\_{id'} = k](https://latex.codecogs.com/png.latex?y_%7Bid%27%7D%20%3D%20k "y_{id'} = k")
for
![k = 1, \\dots, K](https://latex.codecogs.com/png.latex?k%20%3D%201%2C%20%5Cdots%2C%20K "k = 1, \dots, K")
is an ordinal trait belonging to one of
![K](https://latex.codecogs.com/png.latex?K "K") **ordered** categories,
we have a 1-to-1 mapping

![
y\_{snd} = g \\left( x\_{snd} \\right) = k, \\text{if } \\gamma\_{k-1} \\leq  x\_{snd} &lt;  \\gamma\_{k},
](https://latex.codecogs.com/png.latex?%0Ay_%7Bsnd%7D%20%3D%20g%20%5Cleft%28%20x_%7Bsnd%7D%20%5Cright%29%20%3D%20k%2C%20%5Ctext%7Bif%20%7D%20%5Cgamma_%7Bk-1%7D%20%5Cleq%20%20x_%7Bsnd%7D%20%3C%20%20%5Cgamma_%7Bk%7D%2C%0A "
y_{snd} = g \left( x_{snd} \right) = k, \text{if } \gamma_{k-1} \leq  x_{snd} <  \gamma_{k},
")

given thresholds
![\\boldsymbol \\gamma = \\left( \\gamma\_{0}, \\gamma\_{1} \\dots, \\gamma\_{K} \\right)^\\top](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20%5Cgamma%20%3D%20%5Cleft%28%20%5Cgamma_%7B0%7D%2C%20%5Cgamma_%7B1%7D%20%5Cdots%2C%20%5Cgamma_%7BK%7D%20%5Cright%29%5E%5Ctop "\boldsymbol \gamma = \left( \gamma_{0}, \gamma_{1} \dots, \gamma_{K} \right)^\top")
where
![\\gamma\_{k - 1} &lt; \\gamma\_k](https://latex.codecogs.com/png.latex?%5Cgamma_%7Bk%20-%201%7D%20%3C%20%5Cgamma_k "\gamma_{k - 1} < \gamma_k"),
![\\gamma\_0 = -\\infty](https://latex.codecogs.com/png.latex?%5Cgamma_0%20%3D%20-%5Cinfty "\gamma_0 = -\infty"),
![\\gamma\_1 = 0](https://latex.codecogs.com/png.latex?%5Cgamma_1%20%3D%200 "\gamma_1 = 0"),
and
![\\gamma\_K = \\infty](https://latex.codecogs.com/png.latex?%5Cgamma_K%20%3D%20%5Cinfty "\gamma_K = \infty").
Finally, if
![y\_{snd}](https://latex.codecogs.com/png.latex?y_%7Bsnd%7D "y_{snd}")
takes continuous values, then any bijective map to the real numbers is
appropriate. We then assume that

![
\\boldsymbol x\_{sn} \\sim \\mathcal N \\left( \\boldsymbol \\mu\_{sn}, \\Lambda^{-1}\\right), 
](https://latex.codecogs.com/png.latex?%0A%5Cboldsymbol%20x_%7Bsn%7D%20%5Csim%20%5Cmathcal%20N%20%5Cleft%28%20%5Cboldsymbol%20%5Cmu_%7Bsn%7D%2C%20%5CLambda%5E%7B-1%7D%5Cright%29%2C%20%0A "
\boldsymbol x_{sn} \sim \mathcal N \left( \boldsymbol \mu_{sn}, \Lambda^{-1}\right), 
")

where
![\\boldsymbol \\mu\_{sn} = \\left( \\mu\_{sn1}, \\dots, \\mu\_{snD'}^s \\right)^\\top \\in \\mathbb R^{D'}](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20%5Cmu_%7Bsn%7D%20%3D%20%5Cleft%28%20%5Cmu_%7Bsn1%7D%2C%20%5Cdots%2C%20%5Cmu_%7BsnD%27%7D%5Es%20%5Cright%29%5E%5Ctop%20%5Cin%20%5Cmathbb%20R%5E%7BD%27%7D "\boldsymbol \mu_{sn} = \left( \mu_{sn1}, \dots, \mu_{snD'}^s \right)^\top \in \mathbb R^{D'}")
is the expected auxiliary trait and
![\\Lambda](https://latex.codecogs.com/png.latex?%5CLambda "\Lambda") is
a diagonal precision matrix parameterised by
![\\boldsymbol \\lambda = \\left( \\lambda\_1, \\dots, \\lambda\_{P} \\right) \\in \\mathbb R\_+^{P}](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20%5Clambda%20%3D%20%5Cleft%28%20%5Clambda_1%2C%20%5Cdots%2C%20%5Clambda_%7BP%7D%20%5Cright%29%20%5Cin%20%5Cmathbb%20R_%2B%5E%7BP%7D "\boldsymbol \lambda = \left( \lambda_1, \dots, \lambda_{P} \right) \in \mathbb R_+^{P}").
Note that we consider one independent precision parameter for each of
the ![P](https://latex.codecogs.com/png.latex?P "P") manifest traits,
where
![\\Lambda\_{d'd'} = \\lambda\_p](https://latex.codecogs.com/png.latex?%5CLambda_%7Bd%27d%27%7D%20%3D%20%5Clambda_p "\Lambda_{d'd'} = \lambda_p")
when
![x\_{snd'}](https://latex.codecogs.com/png.latex?x_%7Bsnd%27%7D "x_{snd'}")
maps to the
![p^{th}](https://latex.codecogs.com/png.latex?p%5E%7Bth%7D "p^{th}")
manifest trait and
![\\lambda\_p = 1](https://latex.codecogs.com/png.latex?%5Clambda_p%20%3D%201 "\lambda_p = 1")
is fixed when that trait is discrete.

In order to define the PLVM, let

![
\\boldsymbol \\mu\_{sn} = \\boldsymbol W \\boldsymbol z\_{sn}
](https://latex.codecogs.com/png.latex?%0A%5Cboldsymbol%20%5Cmu_%7Bsn%7D%20%3D%20%5Cboldsymbol%20W%20%5Cboldsymbol%20z_%7Bsn%7D%0A "
\boldsymbol \mu_{sn} = \boldsymbol W \boldsymbol z_{sn}
")

where
![\\boldsymbol W \\in \\mathbb R^{D' \\times L}](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20W%20%5Cin%20%5Cmathbb%20R%5E%7BD%27%20%5Ctimes%20L%7D "\boldsymbol W \in \mathbb R^{D' \times L}")
is the loading matrix with columns denoted
![\\boldsymbol w\_l](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20w_l "\boldsymbol w_l")
for
![l = 1, \\dots, L](https://latex.codecogs.com/png.latex?l%20%3D%201%2C%20%5Cdots%2C%20L "l = 1, \dots, L"),
rows
![\\boldsymbol w\_{d'}](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20w_%7Bd%27%7D "\boldsymbol w_{d'}")
for
![d' = 1, \\dots, D'](https://latex.codecogs.com/png.latex?d%27%20%3D%201%2C%20%5Cdots%2C%20D%27 "d' = 1, \dots, D'"),
and
![\\boldsymbol z\_{sn} = \\left( z\_{sn1}, \\dots, z\_{snL} \\right)^\\top \\in \\mathbb R^{L}](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20z_%7Bsn%7D%20%3D%20%5Cleft%28%20z_%7Bsn1%7D%2C%20%5Cdots%2C%20z_%7BsnL%7D%20%5Cright%29%5E%5Ctop%20%5Cin%20%5Cmathbb%20R%5E%7BL%7D "\boldsymbol z_{sn} = \left( z_{sn1}, \dots, z_{snL} \right)^\top \in \mathbb R^{L}")
is the individual-specific latent trait. We link
![\\boldsymbol z\_{sn}](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20z_%7Bsn%7D "\boldsymbol z_{sn}")
to the taxon-specific latent trait
![\\boldsymbol f\_s = \\left( f\_{s1}, \\dots, f\_{sL} \\right)^\\top \\in \\mathbb R^{L}](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20f_s%20%3D%20%5Cleft%28%20f_%7Bs1%7D%2C%20%5Cdots%2C%20f_%7BsL%7D%20%5Cright%29%5E%5Ctop%20%5Cin%20%5Cmathbb%20R%5E%7BL%7D "\boldsymbol f_s = \left( f_{s1}, \dots, f_{sL} \right)^\top \in \mathbb R^{L}")
by setting

![
z\_{snl} =  \\sqrt{1 - \\tau\_l^2}f\_{sl} +  \\tau\_l \\, \\epsilon\_{snl},
](https://latex.codecogs.com/png.latex?%0Az_%7Bsnl%7D%20%3D%20%20%5Csqrt%7B1%20-%20%5Ctau_l%5E2%7Df_%7Bsl%7D%20%2B%20%20%5Ctau_l%20%5C%2C%20%5Cepsilon_%7Bsnl%7D%2C%0A "
z_{snl} =  \sqrt{1 - \tau_l^2}f_{sl} +  \tau_l \, \epsilon_{snl},
")

where
![\\tau\_l \\in \\left( 0, 1 \\right)](https://latex.codecogs.com/png.latex?%5Ctau_l%20%5Cin%20%5Cleft%28%200%2C%201%20%5Cright%29 "\tau_l \in \left( 0, 1 \right)")
is the standard deviation of the inter-taxon variation driven by
![\\epsilon\_{snl} \\sim \\mathcal N \\left(0, 1 \\right)](https://latex.codecogs.com/png.latex?%5Cepsilon_%7Bsnl%7D%20%5Csim%20%5Cmathcal%20N%20%5Cleft%280%2C%201%20%5Cright%29 "\epsilon_{snl} \sim \mathcal N \left(0, 1 \right)")
and we have
![\\boldsymbol \\tau = \\left( \\tau\_1, \\dots, \\tau\_L \\right)^\\top](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20%5Ctau%20%3D%20%5Cleft%28%20%5Ctau_1%2C%20%5Cdots%2C%20%5Ctau_L%20%5Cright%29%5E%5Ctop "\boldsymbol \tau = \left( \tau_1, \dots, \tau_L \right)^\top").
Taxon-specific latent traits are modelled as independent phylogenetic
Ornstein-Uhlenbeck (OP) processes over
![\\mathcal T](https://latex.codecogs.com/png.latex?%5Cmathcal%20T "\mathcal T")
such that

![
\\begin{aligned}
f\_{sl} &= f\_l \\left( \\boldsymbol t\_s \\right), \\\\
f\_l \\left( \\boldsymbol t \\right) &\\sim \\mathcal{GP} \\left( 0, k\_l \\left( \\boldsymbol t, \\boldsymbol t'  \\right) \\right), \\\\
k\_l \\left( \\boldsymbol t, \\boldsymbol t' \\right) &=  h\_l^2 \\exp \\left( - \\frac{d\_{\\mathcal T} \\left( \\boldsymbol t, \\boldsymbol t' \\right) }{\\ell} \\right) + \\left( 1 - h\_l^2\\right) \\delta\_{\\boldsymbol t} \\, \\delta\_{\\boldsymbol t'}
\\end{aligned}
](https://latex.codecogs.com/png.latex?%0A%5Cbegin%7Baligned%7D%0Af_%7Bsl%7D%20%26%3D%20f_l%20%5Cleft%28%20%5Cboldsymbol%20t_s%20%5Cright%29%2C%20%5C%5C%0Af_l%20%5Cleft%28%20%5Cboldsymbol%20t%20%5Cright%29%20%26%5Csim%20%5Cmathcal%7BGP%7D%20%5Cleft%28%200%2C%20k_l%20%5Cleft%28%20%5Cboldsymbol%20t%2C%20%5Cboldsymbol%20t%27%20%20%5Cright%29%20%5Cright%29%2C%20%5C%5C%0Ak_l%20%5Cleft%28%20%5Cboldsymbol%20t%2C%20%5Cboldsymbol%20t%27%20%5Cright%29%20%26%3D%20%20h_l%5E2%20%5Cexp%20%5Cleft%28%20-%20%5Cfrac%7Bd_%7B%5Cmathcal%20T%7D%20%5Cleft%28%20%5Cboldsymbol%20t%2C%20%5Cboldsymbol%20t%27%20%5Cright%29%20%7D%7B%5Cell%7D%20%5Cright%29%20%2B%20%5Cleft%28%201%20-%20h_l%5E2%5Cright%29%20%5Cdelta_%7B%5Cboldsymbol%20t%7D%20%5C%2C%20%5Cdelta_%7B%5Cboldsymbol%20t%27%7D%0A%5Cend%7Baligned%7D%0A "
\begin{aligned}
f_{sl} &= f_l \left( \boldsymbol t_s \right), \\
f_l \left( \boldsymbol t \right) &\sim \mathcal{GP} \left( 0, k_l \left( \boldsymbol t, \boldsymbol t'  \right) \right), \\
k_l \left( \boldsymbol t, \boldsymbol t' \right) &=  h_l^2 \exp \left( - \frac{d_{\mathcal T} \left( \boldsymbol t, \boldsymbol t' \right) }{\ell} \right) + \left( 1 - h_l^2\right) \delta_{\boldsymbol t} \, \delta_{\boldsymbol t'}
\end{aligned}
")

where
![h\_l^2 \\in \\left( 0, 1 \\right)](https://latex.codecogs.com/png.latex?h_l%5E2%20%5Cin%20%5Cleft%28%200%2C%201%20%5Cright%29 "h_l^2 \in \left( 0, 1 \right)")
is the heritability of the
![l^{th}](https://latex.codecogs.com/png.latex?l%5E%7Bth%7D "l^{th}")
latent trait such that
![\\boldsymbol h = \\left(h\_1, \\dots, h\_L\\right)^\\top](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20h%20%3D%20%5Cleft%28h_1%2C%20%5Cdots%2C%20h_L%5Cright%29%5E%5Ctop "\boldsymbol h = \left(h_1, \dots, h_L\right)^\top"),
![\\ell = 2](https://latex.codecogs.com/png.latex?%5Cell%20%3D%202 "\ell = 2")
is the fixed phylogenetic length-scale, and

![
\\delta\_{\\mathcal t} =
\\begin{cases}
1, & \\text{if }\\boldsymbol t   \\in \\left\\{ \\boldsymbol t\_1, \\dots, \\boldsymbol t\_S\\right\\}, \\\\
0, & \\text{otherwise,}
\\end{cases}
](https://latex.codecogs.com/png.latex?%0A%5Cdelta_%7B%5Cmathcal%20t%7D%20%3D%0A%5Cbegin%7Bcases%7D%0A1%2C%20%26%20%5Ctext%7Bif%20%7D%5Cboldsymbol%20t%20%20%20%5Cin%20%5Cleft%5C%7B%20%5Cboldsymbol%20t_1%2C%20%5Cdots%2C%20%5Cboldsymbol%20t_S%5Cright%5C%7D%2C%20%5C%5C%0A0%2C%20%26%20%5Ctext%7Botherwise%2C%7D%0A%5Cend%7Bcases%7D%0A "
\delta_{\mathcal t} =
\begin{cases}
1, & \text{if }\boldsymbol t   \in \left\{ \boldsymbol t_1, \dots, \boldsymbol t_S\right\}, \\
0, & \text{otherwise,}
\end{cases}
")

identifies extant taxa. Restricting
![h\_l](https://latex.codecogs.com/png.latex?h_l "h_l") and
![\\ell](https://latex.codecogs.com/png.latex?%5Cell "\ell") in this
manner avoids scale invariance in the Gaussian likelihood, but still
provides a flexible model for latent trait evolution.

The OU process is a Gauss-Markov process. As such, it can be shown that

![
\\begin{aligned}
f\_{2S - 1, l} &\\sim \\mathcal N \\left( 0, h\_l^2 \\right)\\\\
f\_{sl} \\mid f\_{\\operatorname{pa} \\left( s \\right),l} &\\sim \\mathcal N \\left( \\nu\_{s, l} f\_{\\operatorname{pa} \\left( s \\right),l}, \\eta\_{s, l}^2 \\right)
\\end{aligned}
](https://latex.codecogs.com/png.latex?%0A%5Cbegin%7Baligned%7D%0Af_%7B2S%20-%201%2C%20l%7D%20%26%5Csim%20%5Cmathcal%20N%20%5Cleft%28%200%2C%20h_l%5E2%20%5Cright%29%5C%5C%0Af_%7Bsl%7D%20%5Cmid%20f_%7B%5Coperatorname%7Bpa%7D%20%5Cleft%28%20s%20%5Cright%29%2Cl%7D%20%26%5Csim%20%5Cmathcal%20N%20%5Cleft%28%20%5Cnu_%7Bs%2C%20l%7D%20f_%7B%5Coperatorname%7Bpa%7D%20%5Cleft%28%20s%20%5Cright%29%2Cl%7D%2C%20%5Ceta_%7Bs%2C%20l%7D%5E2%20%5Cright%29%0A%5Cend%7Baligned%7D%0A "
\begin{aligned}
f_{2S - 1, l} &\sim \mathcal N \left( 0, h_l^2 \right)\\
f_{sl} \mid f_{\operatorname{pa} \left( s \right),l} &\sim \mathcal N \left( \nu_{s, l} f_{\operatorname{pa} \left( s \right),l}, \eta_{s, l}^2 \right)
\end{aligned}
")

where, setting
![k\_{s,s', l} = k\_l \\left( \\boldsymbol t\_s, \\boldsymbol t\_{s'} \\right)](https://latex.codecogs.com/png.latex?k_%7Bs%2Cs%27%2C%20l%7D%20%3D%20k_l%20%5Cleft%28%20%5Cboldsymbol%20t_s%2C%20%5Cboldsymbol%20t_%7Bs%27%7D%20%5Cright%29 "k_{s,s', l} = k_l \left( \boldsymbol t_s, \boldsymbol t_{s'} \right)")
and
![k\_{s, l} = k\_l \\left( \\boldsymbol t\_s, \\boldsymbol t\_{s} \\right)](https://latex.codecogs.com/png.latex?k_%7Bs%2C%20l%7D%20%3D%20k_l%20%5Cleft%28%20%5Cboldsymbol%20t_s%2C%20%5Cboldsymbol%20t_%7Bs%7D%20%5Cright%29 "k_{s, l} = k_l \left( \boldsymbol t_s, \boldsymbol t_{s} \right)"),

![
\\begin{aligned}
\\nu\_{sl} &= k\_{s, \\operatorname{pa} \\left( s \\right),l} k\_{\\operatorname{pa} \\left( s \\right), l}^{-1}, \\\\
&= \\exp\\left( -\\frac{d\_{\\mathcal T} \\left( \\boldsymbol t\_s, \\boldsymbol t\_{\\operatorname{pa} \\left( s \\right)} \\right) }{\\ell}\\right), \\\\
\\eta\_{sl}^2 &= k\_{s, l} - \\nu\_{s, l} \\, k\_{\\operatorname{pa} \\left( s \\right), s, l}, \\\\
&= h\_l^2 \\left( 1 - \\exp \\left( - 2 \\frac{d\_{\\mathcal T} \\left( \\boldsymbol t\_s, \\boldsymbol t\_{\\operatorname{pa} \\left( s \\right)} \\right) }{\\ell} \\right) \\right) + \\left( 1 - h\_l^2\\right) \\delta\_{\\boldsymbol t\_s}
\\end{aligned}
](https://latex.codecogs.com/png.latex?%0A%5Cbegin%7Baligned%7D%0A%5Cnu_%7Bsl%7D%20%26%3D%20k_%7Bs%2C%20%5Coperatorname%7Bpa%7D%20%5Cleft%28%20s%20%5Cright%29%2Cl%7D%20k_%7B%5Coperatorname%7Bpa%7D%20%5Cleft%28%20s%20%5Cright%29%2C%20l%7D%5E%7B-1%7D%2C%20%5C%5C%0A%26%3D%20%5Cexp%5Cleft%28%20-%5Cfrac%7Bd_%7B%5Cmathcal%20T%7D%20%5Cleft%28%20%5Cboldsymbol%20t_s%2C%20%5Cboldsymbol%20t_%7B%5Coperatorname%7Bpa%7D%20%5Cleft%28%20s%20%5Cright%29%7D%20%5Cright%29%20%7D%7B%5Cell%7D%5Cright%29%2C%20%5C%5C%0A%5Ceta_%7Bsl%7D%5E2%20%26%3D%20k_%7Bs%2C%20l%7D%20-%20%5Cnu_%7Bs%2C%20l%7D%20%5C%2C%20k_%7B%5Coperatorname%7Bpa%7D%20%5Cleft%28%20s%20%5Cright%29%2C%20s%2C%20l%7D%2C%20%5C%5C%0A%26%3D%20h_l%5E2%20%5Cleft%28%201%20-%20%5Cexp%20%5Cleft%28%20-%202%20%5Cfrac%7Bd_%7B%5Cmathcal%20T%7D%20%5Cleft%28%20%5Cboldsymbol%20t_s%2C%20%5Cboldsymbol%20t_%7B%5Coperatorname%7Bpa%7D%20%5Cleft%28%20s%20%5Cright%29%7D%20%5Cright%29%20%7D%7B%5Cell%7D%20%5Cright%29%20%5Cright%29%20%2B%20%5Cleft%28%201%20-%20h_l%5E2%5Cright%29%20%5Cdelta_%7B%5Cboldsymbol%20t_s%7D%0A%5Cend%7Baligned%7D%0A "
\begin{aligned}
\nu_{sl} &= k_{s, \operatorname{pa} \left( s \right),l} k_{\operatorname{pa} \left( s \right), l}^{-1}, \\
&= \exp\left( -\frac{d_{\mathcal T} \left( \boldsymbol t_s, \boldsymbol t_{\operatorname{pa} \left( s \right)} \right) }{\ell}\right), \\
\eta_{sl}^2 &= k_{s, l} - \nu_{s, l} \, k_{\operatorname{pa} \left( s \right), s, l}, \\
&= h_l^2 \left( 1 - \exp \left( - 2 \frac{d_{\mathcal T} \left( \boldsymbol t_s, \boldsymbol t_{\operatorname{pa} \left( s \right)} \right) }{\ell} \right) \right) + \left( 1 - h_l^2\right) \delta_{\boldsymbol t_s}
\end{aligned}
")

for
![s = 1, \\dots 2S - 1](https://latex.codecogs.com/png.latex?s%20%3D%201%2C%20%5Cdots%202S%20-%201 "s = 1, \dots 2S - 1").
We also define
![\\boldsymbol \\nu\_s = \\left(\\nu\_{s1}, \\dots, \\nu\_{sL} \\right)^\\top](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20%5Cnu_s%20%3D%20%5Cleft%28%5Cnu_%7Bs1%7D%2C%20%5Cdots%2C%20%5Cnu_%7BsL%7D%20%5Cright%29%5E%5Ctop "\boldsymbol \nu_s = \left(\nu_{s1}, \dots, \nu_{sL} \right)^\top")
and
![\\boldsymbol \\eta\_s = \\left(\\eta\_{s1}, \\dots, \\eta\_{sL} \\right)^\\top](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20%5Ceta_s%20%3D%20%5Cleft%28%5Ceta_%7Bs1%7D%2C%20%5Cdots%2C%20%5Ceta_%7BsL%7D%20%5Cright%29%5E%5Ctop "\boldsymbol \eta_s = \left(\eta_{s1}, \dots, \eta_{sL} \right)^\top")
at this point.

The PLVM admits conjugate priors for the loading
![\\boldsymbol W](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20W "\boldsymbol W")
and free auxiliary trait precision parameters
![\\boldsymbol \\lambda](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20%5Clambda "\boldsymbol \lambda")
when we have

![
\\begin{aligned}
\\boldsymbol w\_l \\mid \\alpha\_l &\\sim \\mathcal N \\left( \\boldsymbol 0, \\alpha\_l^{-1} \\boldsymbol C\_w \\right), \\\\
\\lambda\_p &\\sim
\\begin{cases}
\\delta \\left( 1 \\right), & \\text{if } p \\text{ is a discrete trait}, \\\\
\\operatorname{Gamma} \\left( a\_\\lambda, b\_\\lambda \\right), &\\text{otherwise},
\\end{cases}
\\end{aligned}
](https://latex.codecogs.com/png.latex?%0A%5Cbegin%7Baligned%7D%0A%5Cboldsymbol%20w_l%20%5Cmid%20%5Calpha_l%20%26%5Csim%20%5Cmathcal%20N%20%5Cleft%28%20%5Cboldsymbol%200%2C%20%5Calpha_l%5E%7B-1%7D%20%5Cboldsymbol%20C_w%20%5Cright%29%2C%20%5C%5C%0A%5Clambda_p%20%26%5Csim%0A%5Cbegin%7Bcases%7D%0A%5Cdelta%20%5Cleft%28%201%20%5Cright%29%2C%20%26%20%5Ctext%7Bif%20%7D%20p%20%5Ctext%7B%20is%20a%20discrete%20trait%7D%2C%20%5C%5C%0A%5Coperatorname%7BGamma%7D%20%5Cleft%28%20a_%5Clambda%2C%20b_%5Clambda%20%5Cright%29%2C%20%26%5Ctext%7Botherwise%7D%2C%0A%5Cend%7Bcases%7D%0A%5Cend%7Baligned%7D%0A "
\begin{aligned}
\boldsymbol w_l \mid \alpha_l &\sim \mathcal N \left( \boldsymbol 0, \alpha_l^{-1} \boldsymbol C_w \right), \\
\lambda_p &\sim
\begin{cases}
\delta \left( 1 \right), & \text{if } p \text{ is a discrete trait}, \\
\operatorname{Gamma} \left( a_\lambda, b_\lambda \right), &\text{otherwise},
\end{cases}
\end{aligned}
")

where
![\\alpha\_l \\in \\mathbb R\_+](https://latex.codecogs.com/png.latex?%5Calpha_l%20%5Cin%20%5Cmathbb%20R_%2B "\alpha_l \in \mathbb R_+")
is the Automatic Relevance Determination (ARD) precision hyperparameter
such that
![\\boldsymbol \\alpha = \\left( \\alpha\_1, \\dots, \\alpha\_L \\right)^\\top](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20%5Calpha%20%3D%20%5Cleft%28%20%5Calpha_1%2C%20%5Cdots%2C%20%5Calpha_L%20%5Cright%29%5E%5Ctop "\boldsymbol \alpha = \left( \alpha_1, \dots, \alpha_L \right)^\top"),
![\\boldsymbol C\_w](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20C_w "\boldsymbol C_w")
is a fixed covariance matrix enforcing prior smoothness constraints on
![\\boldsymbol w\_q](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20w_q "\boldsymbol w_q"),
and
![a\_\\lambda](https://latex.codecogs.com/png.latex?a_%5Clambda "a_\lambda")
and
![b\_\\lambda](https://latex.codecogs.com/png.latex?b_%5Clambda "b_\lambda")
are shape and rate hyperparameters for trait
![p](https://latex.codecogs.com/png.latex?p "p") taking values along a
continuum, and we let
![\\delta \\left( 1 \\right)](https://latex.codecogs.com/png.latex?%5Cdelta%20%5Cleft%28%201%20%5Cright%29 "\delta \left( 1 \right)")
denote the one point distribution with all its mass at
![1](https://latex.codecogs.com/png.latex?1 "1").

The PLVM also provides a direct approach to approach to ancestral
reconstruction. Firstly, note that this model for trait evolution allows
for within-taxon variation on the latent traits. Stripping away this
stochastic effect on trait evolution allows us define taxon-specific
manifest and auxiliary traits,
![\\boldsymbol y\_s](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20y_s "\boldsymbol y_s")
and
![\\boldsymbol x\_s](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20x_s "\boldsymbol x_s"),
such that

![
\\boldsymbol x\_s \\sim \\mathcal N \\left( \\boldsymbol W \\left( \\boldsymbol I\_L - \\Sigma\_{\\tau} \\right)^{1/2} \\boldsymbol f\_s, \\Lambda^{-1} \\right)
](https://latex.codecogs.com/png.latex?%0A%5Cboldsymbol%20x_s%20%5Csim%20%5Cmathcal%20N%20%5Cleft%28%20%5Cboldsymbol%20W%20%5Cleft%28%20%5Cboldsymbol%20I_L%20-%20%5CSigma_%7B%5Ctau%7D%20%5Cright%29%5E%7B1%2F2%7D%20%5Cboldsymbol%20f_s%2C%20%5CLambda%5E%7B-1%7D%20%5Cright%29%0A "
\boldsymbol x_s \sim \mathcal N \left( \boldsymbol W \left( \boldsymbol I_L - \Sigma_{\tau} \right)^{1/2} \boldsymbol f_s, \Lambda^{-1} \right)
")

for
![s = 1, \\dots, 2S - 1](https://latex.codecogs.com/png.latex?s%20%3D%201%2C%20%5Cdots%2C%202S%20-%201 "s = 1, \dots, 2S - 1"),
![g : \\boldsymbol x\_s \\to \\boldsymbol y\_s](https://latex.codecogs.com/png.latex?g%20%3A%20%5Cboldsymbol%20x_s%20%5Cto%20%5Cboldsymbol%20y_s "g : \boldsymbol x_s \to \boldsymbol y_s")
as before, and
![\\Sigma\_{\\tau}](https://latex.codecogs.com/png.latex?%5CSigma_%7B%5Ctau%7D "\Sigma_{\tau}")
is the diagonal covariance matrix parametrised by
![\\left( \\tau\_1^2, \\dots, \\tau\_L^2\\right)^\\top](https://latex.codecogs.com/png.latex?%5Cleft%28%20%5Ctau_1%5E2%2C%20%5Cdots%2C%20%5Ctau_L%5E2%5Cright%29%5E%5Ctop "\left( \tau_1^2, \dots, \tau_L^2\right)^\top").

## Sampling the Generative Model

Before outlining our approach inference for the PLVM, we sample from the
model. This provides the synthetic dataset we use to validate the
inference scheme. We require the following packages.

``` r
library(vbar)
library(ape)
library(ggplot2)
library(magrittr)
library(reshape2)
library(dplyr)
library(mvnfast)
library(HDInterval)
library(ggpubr)
```

Our analysis considers
![S = 128](https://latex.codecogs.com/png.latex?S%20%3D%20128 "S = 128")
extant taxa with a shared evolutionary history defined by the phylogeny
below. We label the root node and its children. Not that the `ape`
package labels the root node as
![S+1](https://latex.codecogs.com/png.latex?S%2B1 "S+1") by default.

``` r
S <- 2^7
set.seed(98)
phy <- rcoal(S)
phy <- scale_phylo(phy, max_dist = 1)
#> Loading required namespace: testthat
```

<img src="man/figures/README-plot_phylogeny-1.png" width="100%" />

![P = 4](https://latex.codecogs.com/png.latex?P%20%3D%204 "P = 4")
traits are recorded for
![N\_s = 3](https://latex.codecogs.com/png.latex?N_s%20%3D%203 "N_s = 3")
individuals in the extant taxon
![s = 1, \\dots, S](https://latex.codecogs.com/png.latex?s%20%3D%201%2C%20%5Cdots%2C%20S "s = 1, \dots, S").
These traits consist of an ordinal trait with
![K = 4](https://latex.codecogs.com/png.latex?K%20%3D%204 "K = 4")
ordered categories, a nominal trait with
![K = 3](https://latex.codecogs.com/png.latex?K%20%3D%203 "K = 3")
unordered categories, a scalar-valued continuous trait, and a strictly
positive function-valued trait registered to 32 equidistant points along
the unit interval. Thus,
![\\boldsymbol y\_{sn}](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20y_%7Bsn%7D "\boldsymbol y_{sn}")
and
![\\boldsymbol x\_{sn}](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20x_%7Bsn%7D "\boldsymbol x_{sn}")
are
![D = 35](https://latex.codecogs.com/png.latex?D%20%3D%2035 "D = 35")
and
![D = 37](https://latex.codecogs.com/png.latex?D%20%3D%2037 "D = 37")
dimensional vectors respectively.

``` r
P <- 4
N_s <- rep(3, S)
N <- sum(N_s)
K_nom <- 3
K_ord <- 4
D <- 35
D_prime <- 37
```

Auxiliary traits are generated by a PLVM with
![L = 3](https://latex.codecogs.com/png.latex?L%20%3D%203 "L = 3")
latent traits, parameterised by
![\\boldsymbol W](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20W "\boldsymbol W")
and
![\\boldsymbol \\lambda = \\left(1, 1, 10, 10000 \\right)^\\top](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20%5Clambda%20%3D%20%5Cleft%281%2C%201%2C%2010%2C%2010000%20%5Cright%29%5E%5Ctop "\boldsymbol \lambda = \left(1, 1, 10, 10000 \right)^\top")
such that auxiliary traits associated with discrete traits have
precision 1, the scalar-valued continuous trait has precision 10, and
the function-valued trait has precision of 10000. The loading matrix
![\\boldsymbol W](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20W "\boldsymbol W")
is presented below.

``` r
auxiliary_traits_precision <- c(1, 1, 10, 10000)
L <- 3 
# FVT Loading
K_fvt <- 2^5 
x_fvt <- seq(0, 1, length.out = K_fvt) 
loading_fvt_means <- c(.2, .5, .75)
loading_fvt_sds <- c(.06, .12, .07)
loading_fvt <- sapply(
  seq_along(loading_fvt_means), 
  function(i) dnorm(x_fvt, mean = loading_fvt_means[i], sd = loading_fvt_sds[i])
  )
loading_fvt <- loading_fvt %*% diag(1 / c(5, 4, 6))
# Continuous Trait Loading
loading_con <- c(0.5, 0.25, -0.25)
# Nominal Trait loading
loading_nom <- matrix(c(1, 0, 0, 0, -0.5, 0, 0, 0, 1), byrow = TRUE, nrow = 3, ncol = L)
# Ordinal Trait Loading
loading_ord <- c(0, 1, 0.5)
# Full Loading
loading <- rbind(loading_ord, loading_nom, loading_con, loading_fvt)
rownames(loading) <- c(
  "ord", 
  paste("nom", 1:K_nom, sep = "."),
  "con", 
  paste("fvt", 1:K_fvt, sep = ".")
)
colnames(loading) <- paste("loading", 1:L, sep = ".")
```

<img src="man/figures/README-plot_loading-1.png" width="100%" />

The distribution of taxon- and individual-specific latent traits is
defined by the phylogenetic hyperparameters
![\\boldsymbol \\tau = \\left( \\sqrt{0.001}, \\sqrt{0.05}, \\sqrt{0.1} \\right)^\\top](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20%5Ctau%20%3D%20%5Cleft%28%20%5Csqrt%7B0.001%7D%2C%20%5Csqrt%7B0.05%7D%2C%20%5Csqrt%7B0.1%7D%20%5Cright%29%5E%5Ctop "\boldsymbol \tau = \left( \sqrt{0.001}, \sqrt{0.05}, \sqrt{0.1} \right)^\top")
and
![\\boldsymbol h = \\left( \\sqrt{0.95}, \\sqrt{0.66}, \\sqrt{0.25} \\right)^\\top](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20h%20%3D%20%5Cleft%28%20%5Csqrt%7B0.95%7D%2C%20%5Csqrt%7B0.66%7D%2C%20%5Csqrt%7B0.25%7D%20%5Cright%29%5E%5Ctop "\boldsymbol h = \left( \sqrt{0.95}, \sqrt{0.66}, \sqrt{0.25} \right)^\top")
and so we sample
![\\boldsymbol z\_{sn}](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20z_%7Bsn%7D "\boldsymbol z_{sn}").

``` r
tau <- sqrt(c(0.001, 0.05, 0.1))
h <- sqrt(c(0.95, 0.66, 0.25))

set.seed(99)
latent_taxon_traits <- sapply(
  h, function(x) {
    simulate_phylogenetic_ou(
      phy = phy, 
      heritable_amplitude = x, length_scale = 2,
      environmental_amplitude = sqrt(1 - x^2),
      internal = TRUE
    )
  }
)

taxon_id <- sapply(1:S, function(i) rep(phy$tip.label[i], N_s[i])) %>% 
  c()

latent_individual_traits <- sapply(
  1:L, function(l){
    sqrt(1 - tau[l]^2) * latent_taxon_traits[taxon_id, l] + 
      tau[l] * rnorm(N)
  })
```

Given
![\\boldsymbol W](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20W "\boldsymbol W"),
![\\boldsymbol z\_{sn}](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20z_%7Bsn%7D "\boldsymbol z_{sn}"),
and
![\\boldsymbol \\lambda](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20%5Clambda "\boldsymbol \lambda"),
we can define
![\\boldsymbol \\mu\_{sn}](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20%5Cmu_%7Bsn%7D "\boldsymbol \mu_{sn}")
and sample
![\\boldsymbol x\_{sn}](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20x_%7Bsn%7D "\boldsymbol x_{sn}").

``` r
auxiliary_traits_expectation <- latent_individual_traits %*% t(loading)
auxiliary_traits_precision_vec <- c(
  auxiliary_traits_precision[1], 
  rep(auxiliary_traits_precision[2], K_nom), 
  auxiliary_traits_precision[3], 
  rep(auxiliary_traits_precision[4], K_fvt)
)

set.seed(100)
auxiliary_traits <- auxiliary_traits_expectation + 
  t(
    sapply(1:N, function(i){
      rnorm(D_prime, sd = sqrt(1 / auxiliary_traits_precision_vec))
      }
      )
    )
```

We map
![\\boldsymbol x\_{sn}](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20x_%7Bsn%7D "\boldsymbol x_{sn}")
to
![\\boldsymbol y\_{sn}](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20y_%7Bsn%7D "\boldsymbol y_{sn}")
by ![g](https://latex.codecogs.com/png.latex?g "g") given the ordinal
trait cut-off points
![\\boldsymbol \\gamma = \\left( -\\infty, 0, 0.5, 1, \\infty \\right)^\\top](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20%5Cgamma%20%3D%20%5Cleft%28%20-%5Cinfty%2C%200%2C%200.5%2C%201%2C%20%5Cinfty%20%5Cright%29%5E%5Ctop "\boldsymbol \gamma = \left( -\infty, 0, 0.5, 1, \infty \right)^\top").
Furthermore, we set
![g \\left( x\_{snd} \\right) = \\exp \\left( x\_{snd} \\right) = y\_{snd}](https://latex.codecogs.com/png.latex?g%20%5Cleft%28%20x_%7Bsnd%7D%20%5Cright%29%20%3D%20%5Cexp%20%5Cleft%28%20x_%7Bsnd%7D%20%5Cright%29%20%3D%20y_%7Bsnd%7D "g \left( x_{snd} \right) = \exp \left( x_{snd} \right) = y_{snd}")
when the index ![d](https://latex.codecogs.com/png.latex?d "d")
corresponds to a function-valued trait.

``` r
cut_off_points <- c(-Inf, 0, 0.5, 1, Inf)

auxiliary_index_ord <- 1
auxiliary_index_nom <- max(auxiliary_index_ord) + (1:K_nom)
auxiliary_index_con <- max(auxiliary_index_nom) + 1
auxiliary_index_fvt <- max(auxiliary_index_con) + 1:K_fvt

manifest_traits <- data.frame(
  taxon_id = taxon_id,
  ord = sapply(auxiliary_traits[, auxiliary_index_ord], function(x){
    x > cut_off_points
  }) %>% 
    colSums() %>% 
    factor(ordered = TRUE),
  nom = apply(auxiliary_traits[, auxiliary_index_nom], 1, which.max) %>% 
    factor(ordered = FALSE),
  con = auxiliary_traits[, auxiliary_index_con],
  fvt = exp(unname(auxiliary_traits[, auxiliary_index_fvt]))
)
```

<img src="man/figures/README-plot_manifest_traits-1.png" width="100%" />

Given the mapping ![g](https://latex.codecogs.com/png.latex?g "g"),
loading
![\\boldsymbol W](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20W "\boldsymbol W"),
auxiliary trait precision parameters
![\\boldsymbol \\lambda](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20%5Clambda "\boldsymbol \lambda"),
individual-specific latent traits, and phylogenetic hyperparameters
![\\boldsymbol \\tau](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20%5Ctau "\boldsymbol \tau")
and
![\\boldsymbol h](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20h "\boldsymbol h"),
we can explore the distribution of taxon-specific manifest traits
![\\boldsymbol y\_s](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20y_s "\boldsymbol y_s")
for any
![s = 1, \\dots, 2S - 1](https://latex.codecogs.com/png.latex?s%20%3D%201%2C%20%5Cdots%2C%202S%20-%201 "s = 1, \dots, 2S - 1").

``` r
covariance_ou <- ou_kernel(
  dist.nodes(phy),
  amplitude = 1, length_scale = 2
)
covariance_f <- lapply(h, function(x){
  tmp <- x^2 * covariance_ou
  tmp[1:S, 1:S] <- tmp[1:S, 1:S] + (1 - x^2)* diag(S)
  colnames(tmp) <- rownames(tmp) <- c(phy$tip.label, S + (1:(S-1)))
  tmp
})
covariance_zf <- lapply(1:L, function(i){
  tmp <- (1 - tau[i]^2) * covariance_f[[i]]
  tmp <- tmp[
    c(taxon_id, phy$tip.label, S + (1:(S-1))), 
    c(taxon_id, phy$tip.label, S + (1:(S-1)))
    ]
  tmp[1:N, 1:N] <- tmp[1:N, 1:N] + tau[i]^2 * diag(N)
  tmp
})

ancestral_mean <- sapply(1:L, function(i){
  (covariance_zf[[i]][-(1:N), (1:N)] %*% 
     chol2inv(chol(covariance_zf[[i]][(1:N), (1:N)])) %*% 
     latent_individual_traits[, i])[, 1]
})
ancestral_cov <- lapply(1:L, function(i){
  covariance_zf[[i]][-(1:N), -(1:N)] - 
    covariance_zf[[i]][-(1:N), (1:N)] %*% 
    chol2inv(chol(covariance_zf[[i]][(1:N), (1:N)])) %*% 
    covariance_zf[[i]][(1:N), -(1:N)]
})
ancestral_sd <- sapply(ancestral_cov, function(X){
  sqrt(diag(X))
})
```

``` r
a <- "130" # ancestor
M <- 1000 # samples for ancestral reconstruction

scaled_taxon_latent_trait_ar <- rmvn(
  M, mu = ancestral_mean[a,], sigma = diag(ancestral_sd[a, ]), isChol = TRUE
  )
auxiliary_trait_ar <- scaled_taxon_latent_trait_ar %*% t(loading) + 
  t(
    sapply(1:M, function(i){
      rnorm(D_prime, sd = sqrt(1 / auxiliary_traits_precision_vec))
      }
      )
    )

manifest_trait_ar <- data.frame(
  ord = sapply(auxiliary_trait_ar[, auxiliary_index_ord], function(x){
    x > cut_off_points
  }) %>% 
    colSums() %>% 
    factor(ordered = TRUE),
  nom = apply(auxiliary_trait_ar[, auxiliary_index_nom], 1, which.max) %>% 
    factor(ordered = FALSE),
  con = auxiliary_trait_ar[, auxiliary_index_con],
  fvt = exp(unname(auxiliary_trait_ar[, auxiliary_index_fvt]))
)
```

<img src="man/figures/README-plot_ancestral_reconstruction-1.png" width="100%" />

## Approximate Bayesian Inference

Our objective is to infer a posterior distribution over
individual-specific auxiliary traits, the PLVM loading, individual- and
taxon-specific latent traits, and free auxiliary trait precision
parameters given manifest traits, a fixed phylogeny, and priors on the
model parameters. That is

![
\\begin{aligned}
p \\left( \\boldsymbol X, \\boldsymbol W, \\boldsymbol Z, \\boldsymbol F \\mid \\boldsymbol Y, \\mathcal T, \\boldsymbol \\psi \\right) &\\propto p \\left(\\boldsymbol Y \\mid \\boldsymbol X, \\boldsymbol \\gamma \\right) p \\left( \\boldsymbol X \\mid \\boldsymbol W, \\boldsymbol Z, \\boldsymbol \\lambda \\right) p \\left( \\boldsymbol Z \\mid \\boldsymbol F, \\boldsymbol \\tau \\right) p \\left( \\boldsymbol F \\mid \\boldsymbol h \\right) p \\left( \\boldsymbol W \\mid \\boldsymbol \\alpha \\right) p \\left( \\boldsymbol \\lambda \\mid a\_\\lambda, b\_\\lambda \\right), \\\\
&= \\left( \\prod\_{n = 1}^N p \\left(\\boldsymbol y\_n \\mid \\boldsymbol x\_n \\right) p \\left( \\boldsymbol x\_n \\mid \\boldsymbol W, \\boldsymbol z\_n, \\boldsymbol \\lambda \\right) \\right) 
\\left( \\prod\_{s = 1}^S \\prod\_{n = 1}^{N\_s} p \\left( \\boldsymbol z\_{sn} \\mid \\boldsymbol f\_s, \\boldsymbol \\tau \\right) \\right) \\\\
&\\left( \\prod\_{s = 1}^{2S - 2} p \\left( \\boldsymbol f\_s \\mid \\boldsymbol f\_{\\operatorname{pa} \\left( s \\right)}, \\boldsymbol h \\right) p \\left( \\boldsymbol f\_{2S - 1} \\mid \\boldsymbol h \\right) \\right) 
 \\prod\_{l = 1}^L p \\left( \\boldsymbol w\_l \\mid \\alpha\_l, \\boldsymbol C\_w \\right) \\prod\_{p = 1}^{P} p \\left( \\lambda\_p \\mid a\_\\lambda, b\_\\lambda \\right), 
\\end{aligned}
](https://latex.codecogs.com/png.latex?%0A%5Cbegin%7Baligned%7D%0Ap%20%5Cleft%28%20%5Cboldsymbol%20X%2C%20%5Cboldsymbol%20W%2C%20%5Cboldsymbol%20Z%2C%20%5Cboldsymbol%20F%20%5Cmid%20%5Cboldsymbol%20Y%2C%20%5Cmathcal%20T%2C%20%5Cboldsymbol%20%5Cpsi%20%5Cright%29%20%26%5Cpropto%20p%20%5Cleft%28%5Cboldsymbol%20Y%20%5Cmid%20%5Cboldsymbol%20X%2C%20%5Cboldsymbol%20%5Cgamma%20%5Cright%29%20p%20%5Cleft%28%20%5Cboldsymbol%20X%20%5Cmid%20%5Cboldsymbol%20W%2C%20%5Cboldsymbol%20Z%2C%20%5Cboldsymbol%20%5Clambda%20%5Cright%29%20p%20%5Cleft%28%20%5Cboldsymbol%20Z%20%5Cmid%20%5Cboldsymbol%20F%2C%20%5Cboldsymbol%20%5Ctau%20%5Cright%29%20p%20%5Cleft%28%20%5Cboldsymbol%20F%20%5Cmid%20%5Cboldsymbol%20h%20%5Cright%29%20p%20%5Cleft%28%20%5Cboldsymbol%20W%20%5Cmid%20%5Cboldsymbol%20%5Calpha%20%5Cright%29%20p%20%5Cleft%28%20%5Cboldsymbol%20%5Clambda%20%5Cmid%20a_%5Clambda%2C%20b_%5Clambda%20%5Cright%29%2C%20%5C%5C%0A%26%3D%20%5Cleft%28%20%5Cprod_%7Bn%20%3D%201%7D%5EN%20p%20%5Cleft%28%5Cboldsymbol%20y_n%20%5Cmid%20%5Cboldsymbol%20x_n%20%5Cright%29%20p%20%5Cleft%28%20%5Cboldsymbol%20x_n%20%5Cmid%20%5Cboldsymbol%20W%2C%20%5Cboldsymbol%20z_n%2C%20%5Cboldsymbol%20%5Clambda%20%5Cright%29%20%5Cright%29%20%0A%5Cleft%28%20%5Cprod_%7Bs%20%3D%201%7D%5ES%20%5Cprod_%7Bn%20%3D%201%7D%5E%7BN_s%7D%20p%20%5Cleft%28%20%5Cboldsymbol%20z_%7Bsn%7D%20%5Cmid%20%5Cboldsymbol%20f_s%2C%20%5Cboldsymbol%20%5Ctau%20%5Cright%29%20%5Cright%29%20%5C%5C%0A%26%5Cleft%28%20%5Cprod_%7Bs%20%3D%201%7D%5E%7B2S%20-%202%7D%20p%20%5Cleft%28%20%5Cboldsymbol%20f_s%20%5Cmid%20%5Cboldsymbol%20f_%7B%5Coperatorname%7Bpa%7D%20%5Cleft%28%20s%20%5Cright%29%7D%2C%20%5Cboldsymbol%20h%20%5Cright%29%20p%20%5Cleft%28%20%5Cboldsymbol%20f_%7B2S%20-%201%7D%20%5Cmid%20%5Cboldsymbol%20h%20%5Cright%29%20%5Cright%29%20%0A%20%5Cprod_%7Bl%20%3D%201%7D%5EL%20p%20%5Cleft%28%20%5Cboldsymbol%20w_l%20%5Cmid%20%5Calpha_l%2C%20%5Cboldsymbol%20C_w%20%5Cright%29%20%5Cprod_%7Bp%20%3D%201%7D%5E%7BP%7D%20p%20%5Cleft%28%20%5Clambda_p%20%5Cmid%20a_%5Clambda%2C%20b_%5Clambda%20%5Cright%29%2C%20%0A%5Cend%7Baligned%7D%0A "
\begin{aligned}
p \left( \boldsymbol X, \boldsymbol W, \boldsymbol Z, \boldsymbol F \mid \boldsymbol Y, \mathcal T, \boldsymbol \psi \right) &\propto p \left(\boldsymbol Y \mid \boldsymbol X, \boldsymbol \gamma \right) p \left( \boldsymbol X \mid \boldsymbol W, \boldsymbol Z, \boldsymbol \lambda \right) p \left( \boldsymbol Z \mid \boldsymbol F, \boldsymbol \tau \right) p \left( \boldsymbol F \mid \boldsymbol h \right) p \left( \boldsymbol W \mid \boldsymbol \alpha \right) p \left( \boldsymbol \lambda \mid a_\lambda, b_\lambda \right), \\
&= \left( \prod_{n = 1}^N p \left(\boldsymbol y_n \mid \boldsymbol x_n \right) p \left( \boldsymbol x_n \mid \boldsymbol W, \boldsymbol z_n, \boldsymbol \lambda \right) \right) 
\left( \prod_{s = 1}^S \prod_{n = 1}^{N_s} p \left( \boldsymbol z_{sn} \mid \boldsymbol f_s, \boldsymbol \tau \right) \right) \\
&\left( \prod_{s = 1}^{2S - 2} p \left( \boldsymbol f_s \mid \boldsymbol f_{\operatorname{pa} \left( s \right)}, \boldsymbol h \right) p \left( \boldsymbol f_{2S - 1} \mid \boldsymbol h \right) \right) 
 \prod_{l = 1}^L p \left( \boldsymbol w_l \mid \alpha_l, \boldsymbol C_w \right) \prod_{p = 1}^{P} p \left( \lambda_p \mid a_\lambda, b_\lambda \right), 
\end{aligned}
")

where we have introduced the notation
![\\boldsymbol Y = \\left( \\boldsymbol y\_{11}, \\dots, \\boldsymbol y\_{1N\_1}, \\dots, \\boldsymbol y\_{S1}, \\dots, \\boldsymbol y\_{SN\_S} \\right)^\\top = \\left( \\boldsymbol y\_1, \\dots, \\boldsymbol y\_N \\right)^\\top](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20Y%20%3D%20%5Cleft%28%20%5Cboldsymbol%20y_%7B11%7D%2C%20%5Cdots%2C%20%5Cboldsymbol%20y_%7B1N_1%7D%2C%20%5Cdots%2C%20%5Cboldsymbol%20y_%7BS1%7D%2C%20%5Cdots%2C%20%5Cboldsymbol%20y_%7BSN_S%7D%20%5Cright%29%5E%5Ctop%20%3D%20%5Cleft%28%20%5Cboldsymbol%20y_1%2C%20%5Cdots%2C%20%5Cboldsymbol%20y_N%20%5Cright%29%5E%5Ctop "\boldsymbol Y = \left( \boldsymbol y_{11}, \dots, \boldsymbol y_{1N_1}, \dots, \boldsymbol y_{S1}, \dots, \boldsymbol y_{SN_S} \right)^\top = \left( \boldsymbol y_1, \dots, \boldsymbol y_N \right)^\top")
for
![N = \\sum\_{s = 1}^S N\_s](https://latex.codecogs.com/png.latex?N%20%3D%20%5Csum_%7Bs%20%3D%201%7D%5ES%20N_s "N = \sum_{s = 1}^S N_s")
which suppresses the index
![s](https://latex.codecogs.com/png.latex?s "s"),
![\\boldsymbol X](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20X "\boldsymbol X")
and
![\\boldsymbol Z](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20Z "\boldsymbol Z")
have analogous definitions, and
![\\boldsymbol \\psi = \\left\\{\\boldsymbol \\gamma, \\boldsymbol \\tau, \\boldsymbol h, \\boldsymbol \\alpha, \\boldsymbol C\_w, a\_\\lambda, b\_\\lambda \\right\\}](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20%5Cpsi%20%3D%20%5Cleft%5C%7B%5Cboldsymbol%20%5Cgamma%2C%20%5Cboldsymbol%20%5Ctau%2C%20%5Cboldsymbol%20h%2C%20%5Cboldsymbol%20%5Calpha%2C%20%5Cboldsymbol%20C_w%2C%20a_%5Clambda%2C%20b_%5Clambda%20%5Cright%5C%7D "\boldsymbol \psi = \left\{\boldsymbol \gamma, \boldsymbol \tau, \boldsymbol h, \boldsymbol \alpha, \boldsymbol C_w, a_\lambda, b_\lambda \right\}")
model hyperparameters.
![\\boldsymbol z\_{sn}](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20z_%7Bsn%7D "\boldsymbol z_{sn}")
and
![\\boldsymbol z\_n](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20z_n "\boldsymbol z_n")
are used interchangeably as is most convenient. This hierarchical model
is defined by

![
\\begin{aligned}
\\boldsymbol y\_n \\mid \\boldsymbol x\_n, \\boldsymbol \\gamma  &= g \\left( \\boldsymbol x\_n \\right), \\\\
\\boldsymbol x\_n \\mid \\boldsymbol W, \\boldsymbol z\_n, \\boldsymbol \\lambda &\\sim \\mathcal N \\left( \\boldsymbol W \\boldsymbol z\_n, \\Lambda^{-1}\\right), \\\\
\\boldsymbol z\_{sn} \\mid \\boldsymbol f\_s, \\boldsymbol \\tau  &\\sim \\mathcal N \\left(\\left( \\boldsymbol I\_L - \\Sigma\_\\tau \\right)^{1/2} \\boldsymbol f\_s, \\Sigma\_\\tau\\right), \\\\
\\boldsymbol f\_s \\mid \\boldsymbol f\_{\\operatorname{pa} \\left( s \\right)}, \\boldsymbol h &\\sim \\mathcal N \\left( \\operatorname{diag} \\left( \\boldsymbol \\nu\_s\\right) \\boldsymbol f\_{\\operatorname{pa} \\left( s \\right)}, \\Sigma\_{\\eta, s} \\right), \\\\
\\boldsymbol f\_{2S - 1} \\mid \\boldsymbol h &\\sim \\mathcal N \\left( \\boldsymbol 0, \\Sigma\_h \\right), \\\\
\\boldsymbol w\_q \\mid \\boldsymbol \\alpha, \\boldsymbol C\_w &\\sim \\mathcal N \\left( \\boldsymbol 0, \\alpha\_l^{-1} \\boldsymbol C\_w \\right), \\\\
\\lambda\_p &\\sim
\\begin{cases}
\\delta \\left( 1 \\right), & \\text{if } p \\text{ is a discrete trait}, \\\\
\\operatorname{Gamma} \\left( a\_\\lambda, b\_\\lambda \\right), &\\text{otherwise},
\\end{cases}
\\end{aligned}
](https://latex.codecogs.com/png.latex?%0A%5Cbegin%7Baligned%7D%0A%5Cboldsymbol%20y_n%20%5Cmid%20%5Cboldsymbol%20x_n%2C%20%5Cboldsymbol%20%5Cgamma%20%20%26%3D%20g%20%5Cleft%28%20%5Cboldsymbol%20x_n%20%5Cright%29%2C%20%5C%5C%0A%5Cboldsymbol%20x_n%20%5Cmid%20%5Cboldsymbol%20W%2C%20%5Cboldsymbol%20z_n%2C%20%5Cboldsymbol%20%5Clambda%20%26%5Csim%20%5Cmathcal%20N%20%5Cleft%28%20%5Cboldsymbol%20W%20%5Cboldsymbol%20z_n%2C%20%5CLambda%5E%7B-1%7D%5Cright%29%2C%20%5C%5C%0A%5Cboldsymbol%20z_%7Bsn%7D%20%5Cmid%20%5Cboldsymbol%20f_s%2C%20%5Cboldsymbol%20%5Ctau%20%20%26%5Csim%20%5Cmathcal%20N%20%5Cleft%28%5Cleft%28%20%5Cboldsymbol%20I_L%20-%20%5CSigma_%5Ctau%20%5Cright%29%5E%7B1%2F2%7D%20%5Cboldsymbol%20f_s%2C%20%5CSigma_%5Ctau%5Cright%29%2C%20%5C%5C%0A%5Cboldsymbol%20f_s%20%5Cmid%20%5Cboldsymbol%20f_%7B%5Coperatorname%7Bpa%7D%20%5Cleft%28%20s%20%5Cright%29%7D%2C%20%5Cboldsymbol%20h%20%26%5Csim%20%5Cmathcal%20N%20%5Cleft%28%20%5Coperatorname%7Bdiag%7D%20%5Cleft%28%20%5Cboldsymbol%20%5Cnu_s%5Cright%29%20%5Cboldsymbol%20f_%7B%5Coperatorname%7Bpa%7D%20%5Cleft%28%20s%20%5Cright%29%7D%2C%20%5CSigma_%7B%5Ceta%2C%20s%7D%20%5Cright%29%2C%20%5C%5C%0A%5Cboldsymbol%20f_%7B2S%20-%201%7D%20%5Cmid%20%5Cboldsymbol%20h%20%26%5Csim%20%5Cmathcal%20N%20%5Cleft%28%20%5Cboldsymbol%200%2C%20%5CSigma_h%20%5Cright%29%2C%20%5C%5C%0A%5Cboldsymbol%20w_q%20%5Cmid%20%5Cboldsymbol%20%5Calpha%2C%20%5Cboldsymbol%20C_w%20%26%5Csim%20%5Cmathcal%20N%20%5Cleft%28%20%5Cboldsymbol%200%2C%20%5Calpha_l%5E%7B-1%7D%20%5Cboldsymbol%20C_w%20%5Cright%29%2C%20%5C%5C%0A%5Clambda_p%20%26%5Csim%0A%5Cbegin%7Bcases%7D%0A%5Cdelta%20%5Cleft%28%201%20%5Cright%29%2C%20%26%20%5Ctext%7Bif%20%7D%20p%20%5Ctext%7B%20is%20a%20discrete%20trait%7D%2C%20%5C%5C%0A%5Coperatorname%7BGamma%7D%20%5Cleft%28%20a_%5Clambda%2C%20b_%5Clambda%20%5Cright%29%2C%20%26%5Ctext%7Botherwise%7D%2C%0A%5Cend%7Bcases%7D%0A%5Cend%7Baligned%7D%0A "
\begin{aligned}
\boldsymbol y_n \mid \boldsymbol x_n, \boldsymbol \gamma  &= g \left( \boldsymbol x_n \right), \\
\boldsymbol x_n \mid \boldsymbol W, \boldsymbol z_n, \boldsymbol \lambda &\sim \mathcal N \left( \boldsymbol W \boldsymbol z_n, \Lambda^{-1}\right), \\
\boldsymbol z_{sn} \mid \boldsymbol f_s, \boldsymbol \tau  &\sim \mathcal N \left(\left( \boldsymbol I_L - \Sigma_\tau \right)^{1/2} \boldsymbol f_s, \Sigma_\tau\right), \\
\boldsymbol f_s \mid \boldsymbol f_{\operatorname{pa} \left( s \right)}, \boldsymbol h &\sim \mathcal N \left( \operatorname{diag} \left( \boldsymbol \nu_s\right) \boldsymbol f_{\operatorname{pa} \left( s \right)}, \Sigma_{\eta, s} \right), \\
\boldsymbol f_{2S - 1} \mid \boldsymbol h &\sim \mathcal N \left( \boldsymbol 0, \Sigma_h \right), \\
\boldsymbol w_q \mid \boldsymbol \alpha, \boldsymbol C_w &\sim \mathcal N \left( \boldsymbol 0, \alpha_l^{-1} \boldsymbol C_w \right), \\
\lambda_p &\sim
\begin{cases}
\delta \left( 1 \right), & \text{if } p \text{ is a discrete trait}, \\
\operatorname{Gamma} \left( a_\lambda, b_\lambda \right), &\text{otherwise},
\end{cases}
\end{aligned}
")

where
![\\Sigma\_{\\tau}](https://latex.codecogs.com/png.latex?%5CSigma_%7B%5Ctau%7D "\Sigma_{\tau}")
is the diagonal covariance matrix parametrised by
![\\left( \\tau\_1^2, \\dots, \\tau\_L^2\\right)^\\top](https://latex.codecogs.com/png.latex?%5Cleft%28%20%5Ctau_1%5E2%2C%20%5Cdots%2C%20%5Ctau_L%5E2%5Cright%29%5E%5Ctop "\left( \tau_1^2, \dots, \tau_L^2\right)^\top"),
![\\Sigma\_{\\eta, }](https://latex.codecogs.com/png.latex?%5CSigma_%7B%5Ceta%2C%20%7D "\Sigma_{\eta, }")
is the diagonal covariance matrix parametrised by
![\\left( \\tau\_1^2, \\dots, \\tau\_L^2\\right)^\\top](https://latex.codecogs.com/png.latex?%5Cleft%28%20%5Ctau_1%5E2%2C%20%5Cdots%2C%20%5Ctau_L%5E2%5Cright%29%5E%5Ctop "\left( \tau_1^2, \dots, \tau_L^2\right)^\top")..
