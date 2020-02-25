# Synthetic_likelihood_MRF
This project aims to tackle intractable normlizing contant in the Potts model. We propose a novel method to achieve this aim. 
We take advantage of conditional independence of the Potts model to simplify the likelihood of Potts model to products of some small terms. These small terms can be approximated by Monte Carlo methods. The procedure is refered to as SLPCD (synthetic likelihood for partial conditional denpendence).

# Instructions


## Example of grass image
Run SLPCDforgrass.m to obtain the posterior samples of all the parameter.
In this example, the Potts model is used as hidden variables. And the neighouhood structure is the first neighbourhood structure.
Run SLPCDforgrasssecond.m in folder secondorder to obtain the posterior samples of all the parameter.
In this example, the Potts model is used as hidden variables. And the neighouhood structure is the second neighbourhood structure.


