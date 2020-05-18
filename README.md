# <img src="_static/lighton_small.png" width=60/> Conformational exploration SARS-CoV-2 (coronavirus responsible for COVID-19)

[![GitHub license](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)  [![Twitter](https://img.shields.io/twitter/follow/LightOnIO?style=social)](https://twitter.com/LightOnIO)

This repository contains the script to reproduce the results of the blog post [**Coronavirus “SARS-CoV-2” conformational change-points detected with optical random features**](https://medium.com/@LightOnIO/accelerating-sars-cov2-molecular-dynamics-studies-with-optical-random-features-b8cffdb99b01) as well as the results of [**Optical Random Features versus SARS-CoV-2 Glycoprotein Trajectory: Round#2**](https://medium.com/@LightOnIO/optical-random-features-versus-sars-cov-2-glycoprotein-trajectory-round-2-adcf04d6036d).

## 

## Requirements

We advise creating a `virtualenv` before running these commands. You can create one with `python3 -m venv <venv_name>`. Activate it with source `<path_to_venv>/bin/activate`  before proceeding. We used `python 3.7`for all the simulations.

- Clone the repository and then do `pip install -r requirements.txt`.
- Download the dataset from [this page](https://figshare.com/articles/6_molecular_dynamics_simulations_of_coronavirus_2019-nCoV_protease_model_in_complex_with_different_conformations_of_lopinavir_/11764158). You should put the dataset in the  same folder as the repository. Otherwise, make sure to update the variable `path traj`.

## 

## Replicate our results

For replicating the results of [**Coronavirus “SARS-CoV-2” conformational change-points detected with optical random features**](https://medium.com/@LightOnIO/accelerating-sars-cov2-molecular-dynamics-studies-with-optical-random-features-b8cffdb99b01): start by running

```
python newma_coronavirus.py
```

This outputs a numpy zipped archive. To analyse the results and create your own plot, use the notebook `notebooks/corona_exploration.ipynb`.

For replicating the results of [**Optical Random Features versus SARS-CoV-2 Glycoprotein Trajectory: Round#2**](https://medium.com/@LightOnIO/optical-random-features-versus-sars-cov-2-glycoprotein-trajectory-round-2-adcf04d6036d): start by running

```
python newma_coronavirus_DESRES.py
```

This outputs a numpy zipped archive. To analyse the results and create your own plot, use the notebook `notebooks/corona_exploration_DESRES.ipynb`.

## <img src="_static/lighton_cloud_small.png" width=120/> Access to Optical Processing Units

To request access to LightOn Cloud and try our photonic co-processor, please visit: <https://cloud.lighton.ai/>

For researchers, we also have a LightOn Cloud for Research program, please visit <https://cloud.lighton.ai/lighton-research/> for more information.

## 

## Hardware Specifications

All our results were obtains on a Intel(R) Xeon(R) Gold 6128 CPU @3.40GHz with 12 cores.