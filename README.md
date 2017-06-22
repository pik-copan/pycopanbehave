# pycopanbehave
Python implementation of COPAN:BEHAVE model

Provides code for reproducing the analysis reported in:

C.-F. Schleussner*, J.F. Donges*, D.A. Engemann*, and A. Levermann,
Clustered marginalization of minorities during social transitions induced by co-evolution of behaviour and network structure,
Nature Scientific Reports 6, 30790 (2016),
DOI: 10.1038/srep30790,
* The first three authors share the lead authorship.

Please see the above paper for a detailed mathematical description of the model (Methods section) and references to the relevant scientific literature.

## Dependencies

Main:

- numpy > 1.11.1 (http://www.numpy.org/)
- scipy > 0.17.1 (http://www.scipy.org)
- igraph > 0.7.0 (http://igraph.org/python/)
- pyunicorn > 0.5.1 (Unified Complex Network and RecurreNce analysis toolbox),
  http://www.pik-potsdam.de/~donges/pyunicorn/,
  get it here: https://github.com/pik-copan/pyunicorn,
  or via pip: pip install pyunicorn

For plotting:

- matplotlib > 1.5.1 (http://www.matplotlib.org/)
- seaborn > 0.7.1 (https://web.stanford.edu/~mwaskom/software/seaborn/)

## Code structure

The modelling approach consists of two parts:

- The dynamic social network model (pycopanbehave/, class CopanBehaveModel)
    *  Contains the main functionalities of the underlying adaptive network model
    *  provides interfaces to specific applications, e.g. for the smoker evolution case

- The smoker transition simulations (bin/)
    *  Specifically adapted to test case
    *  Add simulation specific functionalities, e.g. the smoking disposition and the respective transition
    * Optimised to be run in parallel setup

## Run the Model
Being rather computationally intensive, the modelling setup is optimised to run in a parallel MPI setup and individual ensemble members will be saved as separate outputs (requires post-processing, see below)

- A run_ens_smoker_*.py contains a dictionary including all relevant configuration parameters
- The setup in manuscript_runs/ resembles the full 1000 Ensemble member setups underlying the manuscript
- The setup under tests/ resembles a computationally less expensive test case of the model setup

## Postprocessing and analysis
- The parallelisation procedure saves individual ensemble member representations
- bin/derive_percentiles_from_output.py derives percentiles from the respective output that can be used for further analysis

## Plotting results
- manuscript_plot_scripts/plot_manuscript_figures.py allows to reproduce the manuscript figures from model runs
