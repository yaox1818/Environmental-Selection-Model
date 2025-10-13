# Environmental-Selection-Model
This is a model about the environment having a preference for microorganisms.

###Parameters

Under our neutral model, several parameters are adjustable:
  1. pct_evn: Percentage of environmental acquisition, 1-pct_evn is the proportion of parental inheritance.
  2. pct_pool: Percentage of pooled component in the environment.
  3. population size: the population size of hosts.
  4. microbe size: number of microbes associated with one host
  5. number of species: the total number of species in the environment
  6. number of generations: the number of host generations that will be simulated
  7. number generation for observation: every this many generations the diversities and other summary statistics are calculated
  8. replication: the number of simulation with the same parameters you want to repeat  



### Requirement:
   * [**Java 8**](https://www.java.com/)

### Install

Clone git repository and the jar file can be found under the `release` folder.
Alternative, you can recompile it with [ant](http://ant.apache.org/).
```bash
cd microbiosima/java
ant release
cd release/microbiosima*
./bin/microbiosima
```


### Usage
##### Help Menu
```bash
./bin/microbiosima -h
```

##### Case 1

Two arguments provided for percentage of environmental acquisition, and percentage of pooled component in the environment.
```bash
./bin/microbiosima 0.2 0.5
#./bin/microbiosima  pct_env pct_pool
#Effectively equals
#./bin/microbiosima  0.2 0.5 500 1000 150 10000
```
the default settings for other parameters are following:
  - population size=500
  - microbe size=1000
  - number of species=150
  - number of generations=10000

##### Case 2

```bash
./bin/microbiosima 0.2 0.5 -c 50 200 20 50
#./bin/microbiosima pctEnv pctPool -c Pop Micro Spec Gen
```
To run the simulation from terminal with six arguments taken.
- pctEnv: percentage of environmental acquisition
- pctPool: percentage of pooled component in the environment
- Pop: population size
- Micro: microbe size
- Spce: number of species
- Gen: number of generations



##### Additional parameters
  - `--obs` Number generation for observation [default: 200]
  - `--rep`Number of replication [default: 5]




##Development

Our selection and horizontal gene transfer (HGT) model are still under developing process.
