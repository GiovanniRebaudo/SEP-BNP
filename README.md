# SEP-BNP

R codes for separately exchangeable Bayesian non-parametric model.

**Authors**: [Giovanni Rebaudo](https://giovannirebaudo.github.io), [Qiaohui Lin](https://qiaohuilin.github.io) and [Peter Müller](https://web.ma.utexas.edu/users/pmueller).

#### Overview 
This repository is associated with the article [Rebaudo, G., Lin Q. and Müller, P. (2024). Separate exchangeability as modeling principle in Bayesian nonparametrics.]()
The key contribution of the paper is outlined below.
 
> [...] We argue for the use of separate exchangeability as a modeling principle in Bayesian nonparametric (BNP) inference. 
Separate exchangeability is \emph{de facto} widely applied in the Bayesian parametric case, e.g., it naturally arises in simple mixed models.
However, while in some areas, such as random graphs, separate and (closely related) joint exchangeability are widely used, it is curiously underused for several other applications in BNP.
We briefly review the definition of separate exchangeability focusing on the implications of such a definition in Bayesian modeling.
We then discuss two tractable classes of models that implement separate exchangeability that are the natural counterparts of familiar partially exchangeable BNP models.

This repository provides codes to replicate the results in Rebaudo, G., Lin Q., and Müller, P. (2024). Separate exchangeability as modeling principle in Bayesian nonparametrics.

In particular, we provide the `R` code to implement the MCMC to perform posterior inference under the GARP model.

The repository contains the following:

1. `SEP_RPM.R` code to reproduce the main results in the article obtained with the separately exchangeable random partition model;
2. `SEP_Reg.R` code to reproduce the main results in the article obtained with the separately exchangeable regression BNP model;
3. `SEP_fcts.R` functions needed to run the main code;
4. `SEP_Splines.R` code to reproduce the spline plot;
5. `Data-and-Results` folder with data and some results of the analyses.

#### Citation
Please cite the following publication if you use this repository in your research: [Rebaudo, G., Lin Q. and Müller, P. (2024). Separate exchangeability as modeling principle in Bayesian nonparametrics.]()
