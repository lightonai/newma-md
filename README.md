# Conformational exploration SARS-CoV-2 (coronavirus responsible for COVID-19)

This repository contains the script to reproduce the results of Figure 4 of the blog post [**Coronavirus “SARS-CoV-2” conformational change-points detected with optical random features**](https://medium.com/@LightOnIO/accelerating-sars-cov2-molecular-dynamics-studies-with-optical-random-features-b8cffdb99b01). 

## How to install 

We advise creating a `virtualenv` before running these commands. You can create one with `python3 -m venv <venv_name>`. Activate it with source `<path_to_venv>/bin/activate`  before proceeding. We used `python 3.7`for all the simulations.

- Clone the repository and then do `pip install <path_to_repo>`.
- Download the dataset from [this page](https://figshare.com/articles/6_molecular_dynamics_simulations_of_coronavirus_2019-nCoV_protease_model_in_complex_with_different_conformations_of_lopinavir_/11764158). You should put the dataset in the  same folder as the repository. Otherwise, make sure to update the variable ```path traj```.



## Replicate our results

Run 

```
python newma_coronavirus.py
```

This outputs a numpy zipped archive `newmaopu_corona_natoms_{n_atoms}.npz`. To analyse the results and create your own plot, use the notebook `corona_exploration.ipynb`.

## Hardware Specifications

All our results were obtains on a Intel(R) Xeon(R) Gold 6128 CPU @3.40GHz with 12 cores. 



