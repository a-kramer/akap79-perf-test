# AKAP7979 MCMC Performance Tests

This Respository holds several R scripts that sample the AKAP79 model
using several different approaches and method combinations. 

We tests which has the highest effective sampling speed:

$$ v = \frac{N}{2 \tau_{l,\text{int.}}}\,,$$

where $\tau_{l,\text{int.}}$ is the integrated auto-correlation length
(Markov chain time), and $N$ the sample size.
