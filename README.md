# pycopanbehave
Python implementation of COPAN:BEHAVE model

Provides code for reproducing the analysis reported in:

C.-F. Schleussner (#), J.F. Donges (#), D.A. Engemann (#), and A. Levermann,
Co-evolutionary behaviour selection in adaptive social networks predicts clustered marginalization of minorities,
Preprint: arxiv.org:1512.05013 [physics.soc-ph] (2015), in review.
(#) The first three authors share the lead authorship.

Please see the above paper for a detailed mathematical description of the model (Methods section) and references to the relevant scientific literature.

## Dependencies
pyunicorn (Unified Complex Network and RecurreNce analysis toolbox)
http://www.pik-potsdam.de/~donges/pyunicorn/
get it here: https://github.com/pik-copan/pyunicorn
or via pip: pip install pyunicorn

## Code structure

The modelling approach consists of two parts:

- The dynamic social network model (pysoc)
    *  Contains the main functionalities of the underlying adaptive network
    *  provides interface to specific applications e.g. for the smoker evolution case

- The smoker transition simulations (bin)
    *  Specifically adapted to test case
    * Add simulation specific functionalities e.g. the smoking disposition and the respective transition
    ** Optimised to be run in parallel setup

## Run the Model
Being rather computationally intensive, the modelling setup is optimised to run in a parallel MPI setup and individual ensemble members will be saved as separate outputs (requires post-processing, see below)

- A run_ens_smoker_*.py contains a dictionary including all relevant configuration parameters
- The setup in manuscript_runs/ resembles the full 1000 Ensemble member setups underlying the manuscript
- The setup under tests/ resembles a computationally less expensive test case of the model setup

## Postprocessing and analysis
- The parallelisation procedure saves individual ensemble member representations
- bin/derive_percentiles_from_output.py derives percentiles from the respective output that can be used for further analysis
