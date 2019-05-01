# GcompBART
Accelerated Multinomial Probit BART and an Integrated Tool  for Bayesian Prediction and Causal Comparison of Dynamic Policies

This package combines two work:

1. Bayesian Framework for Predictive and Causal Modeling with Application to HIV care Cascade 
This work studies the Bayesian formulation of G computation algorithm, and incorporates Bayesian additive regression trees (BART) into the generative components of the causal framework for posterior sampling of counterfactual longitudinal paths over time under certain dynamic policies of interest.
This integrated tool can be used when the time-varying confounders and the outcomes are continuous, binary, or multinomial. 

2. Accelerated Multinomial Probit Bayesian Additive Regression Trees (MPBART)
The model for multinomial response uses the multinomial probit (MNP) regression framework (Imai and van Dyk 2005). In this package I improved the original MPBART (Kindo et al 2016) so the MCMC convergence is much improved.
