# Hyperevent network modelling of partially observed gossip data

This repository contains the code used for the simulation study and empirical application presented in the paper *Hyperevent network modelling of partially observed
gossip data*.

Traditional inferential techniques for relational event models typically rely on fully observed event sequences. However, in many real-world networks, events are only partially observed, posing challenges for standard likelihood-based methods. In this work, we build on the Relational Hyperevent Model (RHEM) to analyze higher-order interactions, such as gossip involving multiple individuals, within a framework that accommodates partially observed data.

By introducing a *right-censored full likelihood* approach, we estimate relevant effects in a gossip case study, capturing both linear and smooth components as well as unobserved heterogeneity. The proposed method is validated through a simulation study.

Finally, we we apply the model to longitudinal survey data from over 40 secondary school communities in Hungary, uncovering the temporal dynamics and social drivers of gossiping behavior.

Repository structure:
- gossipdata_analysis: contains survey data and scripts for data pre-processing, model fitting, and plotting related to the analysis of gossiping among secondary school students in Hungary, as presented in the final section of the paper.
- gossipdata_simulation: includes scripts for data simulation, model fitting, and plotting  for the simulation study.
- utils: provides the definitions and implementations of the covariates used in both the data analysis and the simulation study.
