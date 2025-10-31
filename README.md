# Hyperevent network modelling of partially observed gossip data

This repository contains the code used for the simulation study and empirical application presented in the paper *Hyperevent network modelling of partially observed
gossip data*.

Traditional inferential techniques for relational event models typically rely on fully observed event sequences. However, in many real-world networks, events are only partially observed, posing challenges for standard likelihood-based methods. By extending traditional inferential approaches for RHEM to right-censored interval-time data, we show how flexible and efficient generalized additive models can be used to estimate effects of interest. The proposed method is validated through a simulation study. Finally, our analysis of longitudinal survey data from over 40 secondary school communities in Hungary illustrates how a model that accounts for linear, smooth and random effects can identify the social drivers of gossiping, while revealing complex temporal dynamics.



## Repository structure:

- **`gossipdata_analysis/`** : contains survey data and scripts for data pre-processing, model fitting, and plotting related to the analysis of gossiping among secondary school students in Hungary, as presented in the final section of the paper.
- **`gossipdata_simulation/`**: includes scripts for data simulation, model fitting, and plotting  for the simulation study.
- **`utils/`** : provides the definitions and implementations of the covariates used in both the data analysis and the simulation study.
