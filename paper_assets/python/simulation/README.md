# Simulation code for COMB

This code allows us to run simulations for the COMB paper, and keep track of the results in a systematic way.

The main idea is that we have many different combinations of parameters to try, but these large number of combinations can end up becoming unwieldy to deal with. This library is designed to streamline the process of running these simulations and storing the outputs for further analysis.

## Running simulations

### Prerequisites

JAGS must be installed, along with [uv](https://github.com/astral-sh/uv).

### Configs

The configs are located in `sim_lib/sim_configs.py`. For each parameter, add the different values you would like that parameter to take on. Be mindful that adding too many values quickly leads to a large number of combinations.

### Running

To run, use the following command

```bash
make run CONFIG=<config-name> NUM_SIMS=<nsims> NUM_PROCESSES=<nproc>
```

`NUM_SIMS` corresponds to the number of simulations you want to run for a single combination. `NUM_PROCESSES` corresponds to the number of processes you want to use when running (they will run concurrently).

## Outputs

The outputs will be located in `data/<config-name>`. Each combination gets assigned a hash, which allows you to rerun the same config and just add more simulations without removing old simulations.
