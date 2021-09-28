[![issues](https://img.shields.io/github/issues/sparktseung/LRMoE)](https://github.com/sparktseung/LRMoE/issues)
[![forks](https://img.shields.io/github/forks/sparktseung/LRMoE)](https://github.com/sparktseung/LRMoE/network/members)
[![stars](https://img.shields.io/github/stars/sparktseung/LRMoE)](https://github.com/sparktseung/LRMoE/stargazers)

# LRMoE Package

## Introduction

**LRMoE** is a package tailor-made for actuarial applications which allows actuarial researchers and practitioners to model and analyze insurance loss frequencies and severities using the Logit-weighted Reduced Mixture-of-Experts (LRMoE) model. The flexibility of LRMoE models is theoretically justified in [Fung et al. (2019)](https://www.sciencedirect.com/science/article/pii/S0167668719303956), and an application of LRMoE for modelling correlated insurance claim frequencies is in [Fung et al. (2019)](https://www.cambridge.org/core/journals/astin-bulletin-journal-of-the-iaa/article/class-of-mixture-of-experts-models-for-general-insurance-application-to-correlated-claim-frequencies/E9FCCAD03E68C3908008448B806BAF8E).

The package **LRMoE** offers several new distinctive features which are motivated by various actuarial applications and mostly cannot be achieved using existing packages for mixture models. Key features include:
* A wider coverage on frequency and severity distributions and their zero inflation;
* The flexibility to vary classes of distributions across components;
* Parameter estimation under data censoring and truncation;
* A collection of insurance rate making and reserving functions; and
* Model selection and visualization tools.

While **LRMoE** was initially developed for actuarial application, this package also allows for customized expert functions for various modelling problems outside of the insurance context.

## Package Installation and Usage

While we prepare for submitting the package to Comprehensive R Archive Network ([CRAN](https://cran.r-project.org/)), the package is available in this github repository. The package manual documenting all functions in the package can also be found in the current repository. 

In R, the **LRMoE** package can be downloaded and installed by running the following code.

```R
library(devtools)
install_github("sparktseung/LRMoE")
```

To get started, see [here](./articles/LRMoE.html).

## Issue Report and Suggestions

Issues and suggestions can be posted on [https://github.com/sparktseung/LRMoE/issues](https://github.com/sparktseung/LRMoE/issues).

## Development Planned and In Progress

See [https://github.com/sparktseung/LRMoE/projects](https://github.com/sparktseung/LRMoE/projects).
