# Particle-filters-for-multi-scale-CRNs-v.1
README

This project provides sequential importance resampling particle filters for two specific gene reaction networks. The filter is established based on hybrid approximations of those models and can accurately approximate the actual filter when the particle population is large. Please see the literature [1] for the derivation and analysis of these filters.

In each folder, the file "main.m" is the main function, which solves the filtering problem for the considered gene network and compares the performance of different filters. Also, we enclosed three sub-folders showing the performance of the filters when the assumption of Gaussian observation noise is violated.     

Reference:
[1] Z Fang, A Gupta, M Khammash, "Stochastic filtering for multiscale stochastic reaction networks based on hybrid approximations", Journal of Computational Physics 467 (2022): 111441.
