# Synthetic_likelihood_MRF
This project aims to tackle intractable normlizing contant in the Potts model. We propose a novel method to achieve this aim. 
We take advantage of conditional independence of the Potts model to simplify the likelihood of Potts model to products of some small terms. These small terms can be approximated by Monte Carlo methods. The procedure is refered to as SLPCD (synthetic likelihood for partial conditional denpendence).

# Instructions


## Example of grass image
Run SLPCDforgrass.m to obtain the posterior samples of all the parameter.
In this example, the Potts model is used as hidden variables. And the neighouhood structure is the first neighbourhood structure.
Run SLPCDforgrasssecond.m in folder secondorder to obtain the posterior samples of all the parameter.
In this example, the Potts model is used as hidden variables. And the neighouhood structure is the second neighbourhood structure.

## The explaination of the TXT files
p1tableCT=1q=2.txt is the results obtained from Monte Carlo simulations. 'q=2' indicates that this is for Ising model. If 'q=3', the file is for Potts model with q=3. 'CT=1' indicates that the two conditional items have the same label. While "CT=2" indicates that the two conditional items have different labels. 'p1table' indicates that there is only one conditional item. While 'p2table' indicates that there are two conditional items in the conditional density.

The other TXT files can be interpreted similarly.

# Reference

