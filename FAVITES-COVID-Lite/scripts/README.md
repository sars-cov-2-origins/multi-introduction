- `FAVITES-COVID-Lite_noSeqgen.py` - an updated version FAVITES-COVID-Lite that does not simulate sequences. The commands for the simulations are in the `commands` directory. 

- `FAVITES-COVID-Lite-followup.py` - a script that subsamples the transmission network to the first 50k infected people, and from those 50k, subsamples the index case and ascertained individuals after the date of the first hospitalization to subsequently create a coalescent tree

- `GEMF_firsts.py` - a script to pull the first ascertained, unascertained, and hospitalized individual from the simulation

- `stableCoalescence_cladeAnalysis.py` - a script to infer the time of stable coalescence from the coalescent trees, evolve individual mutations down the coalescent tree, and examine the phylogenetic structure of the resulting "mutation" tree

- `collect_FAVITES_results.py` - a script to collect the results of `GEMF_firsts.py` and the stable coalescence for each simulation, putting them in one text for the rejection sampling

- `create_gemf_dict.py` - a script to pool the GEMF results from each simulation into one pickle file for post-processing
