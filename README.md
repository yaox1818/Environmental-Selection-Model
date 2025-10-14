# Environmental-Selection-Model
This is a model about the environment having a preference for microorganisms.

###Parameters

Under our neutral model, several parameters are adjustable:
  1. pct_evn: Percentage of environmental acquisition, 1-pct_evn is the proportion of parental inheritance.
  2. pct_pool: Percentage of pooled component in the environment.
  3. msCoeffInHost: parameter related to microbe selection strength in host.
  4. msCoeffInEnv: parameter related to microbe selection strength in environment.
  5. hsCoeff: parameter related to host selection strength.
  6. HMS_or_TMS: string HMS or TMS to specify host-mediated or trait-mediated microbe selection.
  7. population size: the population size of hosts.
  8. microbe size: number of microbes associated with one host.
  9. number of species: the total number of species in the environment.
  10. number of generations: the number of host generations that will be simulated.
  11. number generation for observation: every this many generations the diversities and other summary statistics are calculated.
  12. number of total traits: the total number of available microbial traits.
  13. number of traits per microbe: the number of traits for each microbe.
  14. replication: the number of simulation with the same parameters you want to repeat.
  15. resource provisioning process: hosts incur fitness costs to provide resources to recruit beneficial microbes and suppress harmful microbes(Reference: [Models of microbiome evolution incorporating host resource provisioning](https://academic.oup.com/ismecommun/article/5/1/ycaf059/8106546)).



### Requirement:
   * [**Java 8**](https://www.java.com/)

### Install

Clone git repository and compile it with [ant](http://ant.apache.org/).
```bash
cd Environmental-Selection-Model/modelE
ant release
cd release/environmental_microbiosima_Java_v1.0
./bin/environmentalMicrobiosima 
```


### Usage
##### Help Menu
```bash
./bin/environmentalMicrobiosima -h
```

##### Case 1

Six arguments provided for percentage of environmental acquisition, percentage of pooled component in the environment, parameter related to microbe selection strength in host, parameter related to microbe selection strength in environment, parameter related to host selection strength, and string HMS or TMS to specify host-mediated or trait-mediated microbe selection.
```bash
./bin/microbiosima 0.9 0.9 1 1 10 TMS
#./bin/microbiosima  pct_env pct_pool msCoeffInHost msCoeffInEnv hsCoeff HMS_or_TMS
#Effectively equals
#./bin/microbiosima  0.9 0.9 1 1 10 TMS 5000 1000000000 150 20000 25 5
```
the default settings for other parameters are following:
  - population size=5000
  - microbe size=1000000000
  - number of species=150
  - number of generations=20000
  - number of total traits=25
  - number of traits per microbe=5

##### Case 2

```bash
./bin/microbiosima 0.9 0.9 1 1 10 TMS -c 50 200 20 50 10 2
#./bin/microbiosima pctEnv pctPool -c Pop Micro Spec Gen
```
To run the simulation from terminal with twelve arguments taken.
- pctEnv: percentage of environmental acquisition
- pctPool: percentage of pooled component in the environment
- msCoeffInHost: parameter related to microbe selection strength in host
- msCoeffInEnv: parameter related to microbe selection strength in environment
- hsCoeff: parameter related to host selection strength
- HMS_or_TMS: string HMS or TMS to specify host-mediated or trait-mediated microbe selection
- Pop: population size
- Micro: microbe size
- Spce: number of species
- Gen: number of generations
- Ngene: number of total traits
- Ngenepm: number of traits per microbe



##### Additional parameters
  - `--obs` Number generation for observation [default: 200]
  - `--rep`Number of replication [default: 5]
  - `--rpp` Use resource provisioning process

