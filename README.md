# Discretisation-of-Stochastic-Differential-Equations
Main code in C++ and graphs in Pylab

## Objectives

The goal of this project is to implement and compare the Euler schema and the Milshtein schema on the Black-Scholes model. We denote <img src="https://latex.codecogs.com/svg.image?T" title="T" /> the time horizon, we have the following :

The discretisation step is <img src="https://latex.codecogs.com/svg.image?\Delta&space;t&space;=&space;\frac{T}{N}" title="\Delta t = \frac{T}{N}" /> with <img src="https://latex.codecogs.com/svg.image?N&space;\in&space;\mathbb{N}^*" title="N \in \mathbb{N}^*" />. 

The k-th discretisation step is <img src="https://latex.codecogs.com/svg.image?t_k&space;=&space;k\Delta&space;t" title="t_k = k\Delta t" /> with <img src="https://latex.codecogs.com/svg.image?k&space;\in&space;\left\{&space;0,...,N&space;\right\}" title="k \in \left\{ 0,...,N \right\}" />.

## Strong convergence

We consider the Black-Scholes model :

<img src="https://latex.codecogs.com/svg.image?dS_t&space;=&space;&space;\sigma&space;S_t&space;dW_t&space;&plus;&space;rS_t&space;dt" title="dS_t = \sigma S_t dW_t + rS_t dt" />

We can do the exact simulation of the underlying using the same brownian increments than the ones used for the generation of the discretised processes.

The values of the underlying are as follows :
- Exact solution : <img src="https://latex.codecogs.com/svg.image?S_{t_{k&plus;1}}&space;=&space;S_{t_k}&space;e^{(r&space;-&space;\frac{\sigma^2}{2})(t_{k&plus;1}&space;-&space;t_k)&space;&plus;&space;\sigma(W_{t_{k&plus;1}}&space;-&space;W_{t_k})}" title="S_{t_{k+1}} = S_{t_k} e^{(r - \frac{\sigma^2}{2})(t_{k+1} - t_k) + \sigma(W_{t_{k+1}} - W_{t_k})}" />
- Euler schema : <img src="https://latex.codecogs.com/svg.image?S^{e}_{t_{k&plus;1}}&space;=&space;S^{e}_{t_k}(1&plus;\sigma(W_{t_{k&plus;1}}&space;-&space;W_{t_k})&plus;r(t_{k&plus;1}&space;-&space;t_k)))" title="S^{e}_{t_{k+1}} = S^{e}_{t_k}(1+\sigma(W_{t_{k+1}} - W_{t_k})+r(t_{k+1} - t_k)))" />
- Milshtein schema : <img src="https://latex.codecogs.com/svg.image?S^{m}_{t_{k&plus;1}}&space;=&space;S^{e}_{t_k}(1&plus;\sigma(W_{t_{k&plus;1}}&space;-&space;W_{t_k})&plus;&space;\frac{1}{2}\sigma^2(W_{t_{k&plus;1}}&space;-&space;W_{t_k})^2&space;&plus;&space;(r-\frac{\sigma^2}{2})(t_{k&plus;1}&space;-&space;t_k))" title="S^{m}_{t_{k+1}} = S^{m}_{t_k}(1+\sigma(W_{t_{k+1}} - W_{t_k})+ \frac{1}{2}\sigma^2(W_{t_{k+1}} - W_{t_k})^2 + (r-\frac{\sigma^2}{2})(t_{k+1} - t_k))" />

The _vitfort.cpp_ file computes the simulation using those equations and write the results in the _vitfort.csv_ file. Then, the file _plot_vitfort.py_ loads this result and plot it using Pylab.

Theoretically, the quantity

<img src="https://latex.codecogs.com/svg.image?\mathbb{E}(\max_{1\leqslant&space;k\leqslant&space;&space;N}&space;(S_{\frac{kT}{N}}&space;-&space;S^{e}_{\frac{kT}{N}})^2)" title="\mathbb{E}(\max_{1\leqslant k\leqslant N} (S_{\frac{kT}{N}} - S^{e}_{\frac{kT}{N}})^2)" />

depends on <img src="https://latex.codecogs.com/svg.image?\frac{1}{N}" title="\frac{1}{N}" /> and the quantity 

<img src="https://latex.codecogs.com/svg.image?\mathbb{E}(\max_{1\leqslant&space;k\leqslant&space;&space;N}&space;(S_{\frac{kT}{N}}&space;-&space;S^{m}_{\frac{kT}{N}})^2)" title="\mathbb{E}(\max_{1\leqslant k\leqslant N} (S_{\frac{kT}{N}} - S^{m}_{\frac{kT}{N}})^2)" />

depends on <img src="https://latex.codecogs.com/svg.image?\frac{1}{N^2}" title="\frac{1}{N^2}" />.

## Weak convergence

### Weak speed

Now, we focus our study on the weak speed of Euler and Milshtein schemas to compute a European Put of maturity <img src="https://latex.codecogs.com/svg.image?T" title="T" /> in the Black-Scholes model.

We want to identify the dependance in <img src="https://latex.codecogs.com/svg.image?N" title="N" /> of the quantities <img src="https://latex.codecogs.com/svg.image?\mathbb{E}(e^{-rT}&space;(K&space;-&space;S^{e}_{T})^&plus;)&space;-&space;\mathbb{E}(e^{-rT}&space;(K&space;-&space;S_{T})^&plus;)" title="\mathbb{E}(e^{-rT} (K - S^{e}_{T})^+) - \mathbb{E}(e^{-rT} (K - S_{T})^+)" /> and <img src="https://latex.codecogs.com/svg.image?\mathbb{E}(e^{-rT}&space;(K&space;-&space;S^{m}_{T})^&plus;)&space;-&space;\mathbb{E}(e^{-rT}&space;(K&space;-&space;S_{T})^&plus;)" title="\mathbb{E}(e^{-rT} (K - S^{m}_{T})^+) - \mathbb{E}(e^{-rT} (K - S_{T})^+)" />.

We first compute <img src="https://latex.codecogs.com/svg.image?\mathbb{E}(e^{-rT}&space;(K&space;-&space;S_{T})^&plus;)" title="\mathbb{E}(e^{-rT} (K - S_{T})^+)" /> with Black-Scholes formula and then we can approach this expectation with a Monte-Carlo simulation, using the same brownian increments to compute <img src="https://latex.codecogs.com/svg.image?\mathbb{E}(e^{-rT}&space;(K&space;-&space;S^{e}_{T})^&plus;)" title="\mathbb{E}(e^{-rT} (K - S^{e}_{T})^+)" /> and <img src="https://latex.codecogs.com/svg.image?\mathbb{E}(e^{-rT}&space;(K&space;-&space;S^{m}_{T})^&plus;)" title="\mathbb{E}(e^{-rT} (K - S^{m}_{T})^+)" />.

This is a variance reduction technique by control variates method.

The _vitfaible.cpp_ file computes the simulation using those equations and write the results in the _vitfaible.csv_ file. Then, the file _plot_vitfaible.py_ loads this result and plot it using Pylab.

The theoretical behavior of <img src="https://latex.codecogs.com/svg.image?\mathbb{E}(e^{-rT}&space;(K&space;-&space;S^{e}_{T})^&plus;)&space;-&space;\mathbb{E}(e^{-rT}&space;(K&space;-&space;S_{T})^&plus;)" title="\mathbb{E}(e^{-rT} (K - S^{e}_{T})^+) - \mathbb{E}(e^{-rT} (K - S_{T})^+)" /> is to be dependent of <img src="https://latex.codecogs.com/svg.image?\frac{1}{N}" title="\frac{1}{N}" />.

As we conduct our study in the context of the Black-Scholes model, the coefficients <img src="https://latex.codecogs.com/svg.image?\sigma" title="\sigma" /> and <img src="https://latex.codecogs.com/svg.image?b" title="b" /> are constants about the time, which leads to the hypothesis for both weak and strong speed to be always verified. Therefore, the model is too general to distinguish between both methods. 

### Romberg's extrapolation

This section focus on the study of the acceleration of weak convergence by Romberg's extrapolation method. We denote respectively <img src="https://latex.codecogs.com/svg.image?S^{e,N}" title="S^{e,N}" /> and <img src="https://latex.codecogs.com/svg.image?S^{e,2N}" title="S^{e,2N}" /> for Euler schemas at time steps <img src="https://latex.codecogs.com/svg.image?N" title="N" /> and <img src="https://latex.codecogs.com/svg.image?2N" title="2N" /> and <img src="https://latex.codecogs.com/svg.image?S^{m,N}" title="S^{m,N}" /> and <img src="https://latex.codecogs.com/svg.image?S^{m,2N}" title="S^{m,2N}" /> for Milshtein schemas at time steps <img src="https://latex.codecogs.com/svg.image?N" title="N" /> and <img src="https://latex.codecogs.com/svg.image?2N" title="2N" />.

We evaluate the following quantities, respectively for Euler and Milshtein schemas :
- <img src="https://latex.codecogs.com/svg.image?\mathbb{E}(e^{-rT}[2(K&space;-&space;S^{e,2N}_{T})^&plus;&space;-&space;(K&space;-&space;S^{e,N}_{T})^&plus;&space;-&space;(K&space;-&space;S_{T})^&plus;])" title="\mathbb{E}(e^{-rT}[2(K - S^{e,2N}_{T})^+ - (K - S^{e,N}_{T})^+ - (K - S_{T})^+])" />
- <img src="https://latex.codecogs.com/svg.image?\mathbb{E}(e^{-rT}[2(K&space;-&space;S^{m,2N}_{T})^&plus;&space;-&space;(K&space;-&space;S^{m,N}_{T})^&plus;&space;-&space;(K&space;-&space;S_{T})^&plus;])" title="\mathbb{E}(e^{-rT}[2(K - S^{m,2N}_{T})^+ - (K - S^{m,N}_{T})^+ - (K - S_{T})^+])" />

As the previous sections, we keep the same brownian increments for all terms in the expectation to reduce the variance.

The evolution equations are the following :
- <img src="https://latex.codecogs.com/svg.image?S_{\frac{(k&plus;1)T}{N}}&space;=&space;S_{\frac{kT}{N}}e^{\sigma(W_{\frac{(k&plus;1)T}{N}}&space;-&space;W_{\frac{kT}{N}})&space;&plus;&space;(r-\frac{\sigma^2}{2})(\frac{(k&plus;1)T}{N}&space;-&space;\frac{kT}{N})}" title="S_{\frac{(k+1)T}{N}} = S_{\frac{kT}{N}}e^{\sigma(W_{\frac{(k+1)T}{N}} - W_{\frac{kT}{N}}) + (r-\frac{\sigma^2}{2})(\frac{(k+1)T}{N} - \frac{kT}{N})}" />
- <img src="https://latex.codecogs.com/svg.image?S^{e,N}_{\frac{(k&plus;1)T}{N}}&space;=&space;S^{e,N}_{\frac{kT}{N}}(1&plus;\sigma(W_{\frac{(k&plus;1)T}{N}}&space;-&space;W_{\frac{kT}{N}})&space;&plus;&space;r(\frac{(k&plus;1)T}{N}&space;-&space;\frac{kT}{N})" title="S^{e,N}_{\frac{(k+1)T}{N}} = S^{e,N}_{\frac{kT}{N}}(1+\sigma(W_{\frac{(k+1)T}{N}} - W_{\frac{kT}{N}}) + r(\frac{(k+1)T}{N} - \frac{kT}{N})" />
- <img src="https://latex.codecogs.com/svg.image?S^{e,2N}_{\frac{(k&plus;1)T}{N}}&space;=&space;S^{e,2N}_{\frac{kT}{N}}(1&plus;\sigma(W_{\frac{(k&plus;1)T}{N}}&space;-&space;W_{\frac{(2k&plus;1)T}{2N}})&space;&plus;&space;\frac{r}{2}(\frac{(k&plus;1)T}{N}&space;-&space;\frac{kT}{N}))&space;(1&space;&plus;&space;\sigma(W_{\frac{(2k&plus;1)T}{2N}}&space;-&space;W_{\frac{kT}{N}})&space;&plus;&space;\frac{r}{2}(\frac{(k&plus;1)T}{N}&space;-&space;\frac{kT}{N}))" title="S^{e,2N}_{\frac{(k+1)T}{N}} = S^{e,2N}_{\frac{kT}{N}}(1+\sigma(W_{\frac{(k+1)T}{N}} - W_{\frac{(2k+1)T}{2N}}) + \frac{r}{2}(\frac{(k+1)T}{N} - \frac{kT}{N})) (1 + \sigma(W_{\frac{(2k+1)T}{2N}} - W_{\frac{kT}{N}}) + \frac{r}{2}(\frac{(k+1)T}{N} - \frac{kT}{N}))" />
- <img src="https://latex.codecogs.com/svg.image?S^{m,N}_{\frac{(k&plus;1)T}{N}}&space;=&space;S^{m,N}_{\frac{kT}{N}}(1&plus;\sigma(W_{\frac{(k&plus;1)T}{N}}&space;-&space;W_{\frac{kT}{N}})&plus;&space;\frac{1}{2}\sigma^2(W_{\frac{(k&plus;1)T}{N}}&space;-&space;W_{\frac{kT}{N}})^2&space;&plus;&space;(r-\frac{\sigma^2}{2})(\frac{(k&plus;1)T}{N}&space;-&space;\frac{kT}{N}))" title="S^{m,N}_{\frac{(k+1)T}{N}} = S^{m,N}_{\frac{kT}{N}}(1+\sigma(W_{\frac{(k+1)T}{N}} - W_{\frac{kT}{N}})+ \frac{1}{2}\sigma^2(W_{\frac{(k+1)T}{N}} - W_{\frac{kT}{N}})^2 + (r-\frac{\sigma^2}{2})(\frac{(k+1)T}{N} - \frac{kT}{N}))" />
- <img src="https://latex.codecogs.com/svg.image?S^{m,2N}_{\frac{(k&plus;1)T}{N}}&space;=&space;S^{m,2N}_{\frac{kT}{N}}(1&plus;\sigma(W_{\frac{(k&plus;1)T}{N}}&space;-&space;W_{\frac{(2k&plus;1)T}{2N}})&plus;&space;\frac{1}{2}\sigma^2(W_{\frac{(k&plus;1)T}{N}}&space;-&space;W_{\frac{(2k&plus;1)T}{2N}})^2&space;&plus;&space;\frac{1}{2}(r-\frac{\sigma^2}{2})(\frac{(k&plus;1)T}{N}&space;-&space;\frac{(2k&plus;1)T}{2N}))(1&plus;\sigma(W_{\frac{(2k&plus;1)T}{2N}}&space;-&space;W_{\frac{kT}{N}})&plus;&space;\frac{1}{2}\sigma^2(W_{\frac{(2k&plus;1)T}{2N}}&space;-&space;W_{\frac{kT}{N}})^2&space;&plus;&space;\frac{1}{2}(r-\frac{\sigma^2}{2})(\frac{(2k&plus;1)T}{2N}&space;-&space;\frac{kT}{N})" title="S^{m,2N}_{\frac{(k+1)T}{N}} = S^{m,2N}_{\frac{kT}{N}}(1+\sigma(W_{\frac{(k+1)T}{N}} - W_{\frac{(2k+1)T}{2N}})+ \frac{1}{2}\sigma^2(W_{\frac{(k+1)T}{N}} - W_{\frac{(2k+1)T}{2N}})^2 + \frac{1}{2}(r-\frac{\sigma^2}{2})(\frac{(k+1)T}{N} - \frac{(2k+1)T}{2N}))(1+\sigma(W_{\frac{(2k+1)T}{2N}} - W_{\frac{kT}{N}})+ \frac{1}{2}\sigma^2(W_{\frac{(2k+1)T}{2N}} - W_{\frac{kT}{N}})^2 + \frac{1}{2}(r-\frac{\sigma^2}{2})(\frac{(2k+1)T}{2N} - \frac{kT}{N})" />

The _romberg.cpp_ file computes the simulation using those equations and write the results in the _romberg.csv_ file. Then, the file _plot_romberg.py_ loads this result and plot it using Pylab.

As there are more quantities to compute regarding the previous methods, the execution time of this code is higher.