# Fit-Seq

## What is Fit-Seq?

Fit-Seq is a MATLAB-based fitness estimation tool for pooled amplicon sequencing studies. Fit-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

It currently has three main functions:
1. evolution_0mut_simplified.m performs simplified simulations of competitve pooled growth of a population of genotypes.
2. evolution_0mut_complex.m performs complex simulations of competitve pooled growth of a population of genotypes.
3. FitSeq.m calculates the fitness of each genotype from read-count time-series data.


## Installing

You must have MATLAB (version 2018a or newer) installed that includes the Optimization Toolbox. 


## Getting Started

A walk-through of how to perform each simulation and FitSeq fitness estimates on these simulations is provided in walk_through.m

### Simplified Simulations 
Models competative pooled growth of a population of genotypes with different fitnesses. This simulation include only growth noise, but not noise from cell transfers, DNA extraction, PCR, or sequencing.

#### INPUTS

-- lineage: number of genotypes of the population

-- t_evo: total number of generations grown

 -- cell_num_ini: a vector of the initial cell number of each genotype at 
                  generation 0,
                  size = lineage * 1

 -- x_ini: a vector of the fitness of each genotype,
           size = lineage * 1

 -- noise_option: options of whether cell growth noise is simulated, 
                  logical (0-1) scaler value, 
                  1 means that the cell growth noise is included
                  0 means that the cell growth noise is not included

 -- 'format': optional, file format of the output file, 'csv'(default) or 'mat'


#### OUTPUTS

-- file_name: the name of the file(s) written by the function.
               When 'format' is set to 'mat', output will be:  
                  'data_evo_simu_0mut_simplified_*Time*.mat'   
               When 'format' is set to 'csv', output will be:  
                   'data_evo_simu_0mut_simplified_*Time*_MeanFitness.mat'   
                   'data_evo_simu_0mut_simplified_*Time*_CellNumber.mat'   
                   'data_evo_simu_0mut_simplified_*Time*_EffectiveCellDepth.mat'   
                   'data_evo_simu_0mut_simplified_*Time*_Parameters.mat'   

#### Example
```
lineage = 1e4;  % number of genotypes of the population
t_evo = 24;  % total number of generations grown
cell_num_ini = 1e2*ones(lineage,1);  % a vector of initial cell number of each genotype at generation 0
x_ini = random('Normal',0,0.2, [lineage,1]);  % a vector of the fitness of each genotype. Here is a Gaussian distribution with mean = 0 and standard devistion = 0.2
noise_option = 1;  %cell growth noise is simulated

% Execute
[ test.simple.csv ] = evolution_0mut_simplified(lineage, t_evo, cell_num_ini, x_ini, noise_option);
```

### Complex Simulations 
Models competative pooled growth of a population of genotypes with different fitnesses. This simulation may include many sources of noise, including growth noise, noise from cell transfers, DNA extraction, PCR, and sequencing.

#### INPUTS
 
 -- lineage: number of genotypes of the population

 -- t_evo: total number of generations grown 

 -- cell_num_ini: a vector of the initial cell number of each genotype at 
                  generation 0,
                  size = lineage * 1

 -- x_ini: a vector of the fitness of each genotype,
           size = lineage * 1

 -- deltat: number of generations between two successive cell transfers

 -- read_depth_average: average number of reads per genotype per
                        sequencing time point

-- cell_num_ini: a vector of the initial cell number of every genotype at the 
                  generation 0, 
                  size = lineage * 1

 -- x_ini: a vector of the fitness of every genotype, 
           size = lineage * 1

 -- noise_option: a vector of whether five types of noise are simulated
                  (cell growth, bottleneck transfer, DNA extraction, PCR, sequencing)
                  size = 1*5 logical (0-1) vector, 
                  1 for the first element means that the cell growth noise is included
                  0 for the first element means that the cell growth noise is not included
                  1 for the second position means that the bottleneck transfer noise is included
                  0 for the second position means that the bottleneck transfer noise is not included
                  1 for the third element means that the DNA extraction noise is included
                  0 for the third element means that the DNA extraction noise is not included
                  1 for the fourth position means that the PCR noise is included
                  0 for the fourth position means that the PCR is not included
                  1 for the fifth element means that the sequencing noise is included
                  0 for the fifth element means that the sequencing noise is not included

 -- 'format': optional, file format of the output file, 'csv'(default) or 'mat'


#### OUTPUTS

 -- file_name: the name of the file(s) written by the function.
               When 'format' is set to 'mat', output will be:
                  'data_evo_simu_0mut_complex_*Time*.mat' 
               When 'format' is set to 'csv', output will be:
                   'data_evo_simu_0mut_complex_*Time*_MeanFitness.mat' 
                   'data_evo_simu_0mut_complex_*Time*_CellNumber.mat' 
                   'data_evo_simu_0mut_complex_*Time*_EffectiveCellDepth.mat' 
                   'data_evo_simu_0mut_complex_*Time*_Parameters.mat'
                   'data_evo_simu_0mut_complex_*Time*_SequencedTimepoints.mat'
                   'data_evo_simu_0mut_complex_*Time*_Reads.mat' 
#### Example
```
lineage = 1e4;  % number of genotypes of the population
t_evo = 24;  % total number of generations grown
cell_num_ini = 1e2*ones(lineage,1);  % a vector of initial cell number of each genotype at generation 0
x_ini = random('Normal',0,0.15, [lineage,1]);  % a vector of the fitness of each genotype. Here, this is a Gaussian distribution with mean = 0 and standard devistion = 0.15
read_depth_average = 100;  % average number of reads per genotype per sequencing
deltat = 8;  % number of generations between two successive cell transfers
noise_option = [1,1,1,1,1];  % a vector of whether to simulate each of five types of noise (cell growth, bottleneck transfer, DNA extraction, PCR, sequencing). Here, all five types of noise are simulated.
                             
% Execute
[ test.complex.csv] = evolution_0mut_complex(lineage, t_evo, cell_num_ini, x_ini, deltat, read_depth_average, noise_option);
```


### Fitness Estimation
Estimates the fitness of each genotype from read-count time-series data.

 INPUTS
 -- t_seq_vec: a vector of all sequencing time points

 -- BC_num_mat_original: a matrix of read number of each genotype at each
                         sequencing time point,
                         size = genotypes * length(t_seq_vec) 

 -- effective_cell_depth: a vector of the effective cell number (number of 
                          cells transferred at the bottleneck) of the population 
                          at each sequencing time point,
                          size = 1 * length(t_seq_vec)

 --  deltat: number of generations between two succesive cell transfers
             This is required in addition to t_seq_vec because not every 
             cell transfer is necessarely sequenced (e.g. t_seq_vec = [0, 3, 6,
             9, 15])

 -- 'format': optional, file format of the output file, 'csv'(default) or 'mat'

 -- 'kappa': optional, a noise parameter that characterizes the total noise 
             introduced by growth, cell transfer, DNA extraction, PCR, and 
             sequencing, default value is 2.5, from Levy et al. Nature 2015 519, 
             181-186. To measure kappa empirically, see that reference. 

 -- 'opt_cycle': optional, the number of cycles used when using likelihood 
                 optimization method to estimate fitness, default value is 2

 -- file_name: the name of the file(s) written by the function.
               When 'format' is set to 'mat', output will be:
                   Fit-Seq_result_*Time*.mat' 
               When 'format' is set to 'csv', output will be:
                   Fit-Seq_result_EstimatedFitness_*Time*.csv'
                   Fit-Seq_result_EstimatedReads_*Time*.csv'

 -- x_estimate_result: a vector of the estimated fitness of each genotype

 -- r_estimate_result: a matrix of the estimated read number of each genotype 
                       at each sequencing time point,
                       size = genotypes * length(t_seq_vec)

 -- x_mean_est: a vector of the estimated mean fitness of the population
                at each sequencing time point,
                size = 1 * length(t_seq_vec)

Inputting data from simulation

```
t_seq_vec = csvread([test.complex.csv, '_SeuqencedTimepoints.csv']); % a vector of all sequencing time points
BC_num_mat_original = csvread([test.complex.csv, '_Reads.csv']); % a matrix of the read number of each genotype at each sequencing time point
effective_cell_depth = csvread([test.complex.csv, '_EffectiveCellDepth.csv'); % a vector of the effective cell number (number of cells transferred at the bottleneck) of population at each sequencing time point
deltat_temp = textscan(fopen([test.complex.csv, '_Paramaters.csv']),'%*f %*f %*f %*f %f %*f %*s','Delimiter',',','headerLines', 1);
deltat = deltat_temp{1}(1); % number of generations between two succesive cell transfers
```


Inputting data from a file (Simulated-Pooled_Growth_Reads.csv) 



### Pooled growth simulation
A numerical simulation framework is developed to simulate the competitive pooled cell growth of a population of genotypes with different fitnesses, and includes all sources of experimental noise, such as growth noise, sampling during bottlenecks, DNA extraction, PCR, and sampling on the sequencer. The simulation is performed in MATLAB using “evolution_0mut_complex.m”.


#### Inputs
-- lineage: number of genotypes of the population

-- t_evo: total number of growth generations

-- cell_num_ini: a vector of initial cell number of each genotype at the 0-th generation, size = lineage * 1

-- x_ini: a vector of the fitness of each genotype, size = lineage * 1

-- deltat: number of generations between successive cell transfers

-- read_depth_average: average number of reads per genotype per sequencing time point

-- noise_option: a vector of options of whether five types of noise, cell growth, cell transfer at the bottleneck, DNA extraction, PCR, and sequencing (in that order in the vector) are simulated, size = 1 * 5, 1 means the noise is included and 0 means the noise is not included.

-- 'format': optional, file format of the output file, 'csv'(default) or 'mat'

#### Outputs
-- file_name: name of the file generated, the standard file generated 'data_evo_simu_0mut_complex_********-*********.mat' when 'format' is set to be 'mat', and 'data_evo_simu_0mut_complex_********-*********.csv' when 'format' is set to be 'csv'

```
lineage = 1e4;
t_evo = 24;
cell_num_ini = 1e2*ones(lineage,1);  
x_ini = random('Normal',0,0.15, [lineage,1]);  % Gaussian with mean = 0 and standard devistion = 0.15
read_depth_average = 100;  
deltat = 8;  
noise_option = [1,1,1,1,1];  % All five types of noise are simulated

[ file_name ] = evolution_0mut_simplified(lineage, t_evo, cell_num_ini, x_ini, noise_option);
```

#### Complex version

### Fitness estimation using Fit-Seq

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone who's code was used
* Inspiration
* etc

