# Introduction to Bayesian statistics in R & brms

4 day course: 24-27 February 2025  
German Centre for Integrative Biodiversity Research (iDiv) Halle-Jena-Leipzig  
benjamin.rosenbaum@idiv.de  

### Outline

The course offers a straightforward and practical approach to applied statistics using Bayesian inference for ecologists. It starts with a general introduction to statistical modeling and the concepts of Bayesian statistics (likelihood, priors, posterior distribution, MCMC sampling). We will move step-by-step from classical ANOVA and linear regression to generalized, nonlinear, or mixed-effects models, with a strong conceptual focus on the building blocks of statistical models.

While previous software required users to code in specific modeling languages (JAGS, NIMBLE, Stan), we are focusing on the user-friendly and flexible R-package ‘brms’, which makes the transition easy for people familiar with ‘lm’ or ‘lme4’. An additional introduction to coding in Stan will be provided for interested participants.

Participants learn how to practically conceptualize their research questions into statistical models. They learn how to specify and critically interpret models of varying complexity in R. The course prepares participants to analyze their own data with state-of-the-art methodology.

### Curriculum

|   | Lecture | Practical |
| ------------- | ------------- | ------------- |
| (1) Statistical modeling      | [pdf](https://benjamin-rosenbaum.github.io/bayesian-intro/Lecture_01.pdf) | [pdf](https://benjamin-rosenbaum.github.io/bayesian-intro/Practical_01.pdf) &nbsp; [html](https://benjamin-rosenbaum.github.io/bayesian-intro/Practical_01.html) &nbsp; [Rcode](https://benjamin-rosenbaum.github.io/bayesian-intro/Practical_01.R)  |
| (2) Bayesian principles       | [pdf](https://benjamin-rosenbaum.github.io/bayesian-intro/Lecture_02.pdf) | [pdf](https://benjamin-rosenbaum.github.io/bayesian-intro/Practical_02.pdf) &nbsp; [html](https://benjamin-rosenbaum.github.io/bayesian-intro/Practical_02.html) &nbsp; [Rcode](https://benjamin-rosenbaum.github.io/bayesian-intro/Practical_02.R)  |
| (3) Priors and posteriors     | [pdf](https://benjamin-rosenbaum.github.io/bayesian-intro/Lecture_03.pdf) | [pdf](https://benjamin-rosenbaum.github.io/bayesian-intro/Practical_03.pdf) &nbsp; [html](https://benjamin-rosenbaum.github.io/bayesian-intro/Practical_03.html) &nbsp; [Rcode](https://benjamin-rosenbaum.github.io/bayesian-intro/Practical_03.R)  |
| (4) Linear models             | [pdf](https://benjamin-rosenbaum.github.io/bayesian-intro/Lecture_04.pdf) | [pdf](https://benjamin-rosenbaum.github.io/bayesian-intro/Practical_04.pdf) &nbsp; [html](https://benjamin-rosenbaum.github.io/bayesian-intro/Practical_04.html) &nbsp; [Rcode](https://benjamin-rosenbaum.github.io/bayesian-intro/Practical_04.R)  |
| (5) Generalized linear models | [pdf](https://benjamin-rosenbaum.github.io/bayesian-intro/Lecture_05.pdf) | [pdf](https://benjamin-rosenbaum.github.io/bayesian-intro/Practical_05.pdf) &nbsp; [html](https://benjamin-rosenbaum.github.io/bayesian-intro/Practical_05.html) &nbsp; [Rcode](https://benjamin-rosenbaum.github.io/bayesian-intro/Practical_05.R)  |
| (6) Mixed effects models      | [pdf](https://benjamin-rosenbaum.github.io/bayesian-intro/Lecture_06.pdf) | [pdf](https://benjamin-rosenbaum.github.io/bayesian-intro/Practical_06.pdf) &nbsp; [html](https://benjamin-rosenbaum.github.io/bayesian-intro/Practical_06.html) &nbsp; [Rcode](https://benjamin-rosenbaum.github.io/bayesian-intro/Practical_06.R)  |
| (7) Stan introduction         | [pdf](https://benjamin-rosenbaum.github.io/bayesian-intro/Lecture_07.pdf) | [pdf](https://benjamin-rosenbaum.github.io/bayesian-intro/Practical_07.pdf) &nbsp; [html](https://benjamin-rosenbaum.github.io/bayesian-intro/Practical_07.html) &nbsp; [Rcode](https://benjamin-rosenbaum.github.io/bayesian-intro/Practical_07.R)  |
| (8) Conclusions               | [pdf](https://benjamin-rosenbaum.github.io/bayesian-intro/Lecture_08.pdf) | pdf &nbsp; html &nbsp; [Rcode](https://benjamin-rosenbaum.github.io/bayesian-intro/Practical_08.R)  |

### Software requirements

- Rstudio: [link](https://posit.co/download/rstudio-desktop/)
- C-toolchain: [link](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started#configuring-c-toolchain),
- which involves RTools (for Windows): [link](https://cran.r-project.org/bin/windows/Rtools/)
- Some R-packages:

```r
update.packages()
install.packages("devtools")    # install packages from github
install.packages("brms")        # our main software package
install.packages("ggplot2")     # plotting
install.packages("bayesplot")   # additional plotting tools
install.packages("sfsmisc")     # mathematical integration through data points
install.packages("performance") # model evaluation
install.packages("arm")         # model evaluation
install.packages("GGally")      # pairs plots
install.packages("emmeans")     # post-hoc analysis
install.packages("ecostats")    # some datasets
devtools::install_github("jfieberg/Data4Ecologists") # more datasets
```
