# MCknownCured (A Generalized Mixture Cure Model Incorporating Known Cured Individuals)

## Overview

The Mixture Cure (MC) models constitute an appropriate and easily interpretable method when studying a time-to-event variable in a population comprised of both susceptible and cured individuals.
In literature, those models usually assume that the latter are unobservable. However, there are cases in which a cured individual may be identified.
For example, when studying the distant metastasis during the lifetime or the miscarriage during pregnancy, individuals that have died without a metastasis or have given birth are certainly non-susceptible.
The same also holds when studying the x-year overall survival or the death during hospital stay. Common MC models ignore this information and consider them all censored, thus yielding in risk of assigning low immune probabilities to cured individuals.
In this study, we consider a MC model that incorporates known information on cured individuals, with the time to cure identification being deterministic, stochastic or even irrelevant.
Furthermore, we compare different strategies that account for cure information such as:
* assigning infinite times to event for known cured cases and adjusting the traditional model
* considering only the probability of cure identification but ignoring the time until that happens.

## Project Structure

- ***`Rscripts`***: R source code/scripts for the data generation (simulations) and building the different models.

