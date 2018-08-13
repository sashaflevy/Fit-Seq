# Fit-Seq

## What is Fit-Seq?

Fit-Seq is a MATLAB-based fitness estimation tool for pooled amplicon sequencing studies. Fit-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. It is GNU Octove compatible (https://www.gnu.org/software/octave/). Five packages (io, nan, struct, statistics, and optim) are required to run Fit-Seq on GNU Octave platform. Computaitonal performace is severely compromised with GNU Octave relative MATLAB. We therefore recommend limiting the use of GNU Octave to small genotype libraries (<1000).  

It currently has three main functions:
1. evolution_0mut_simplified.m performs simplified simulations of competitve pooled growth of a population of genotypes.
2. evolution_0mut_complex.m performs complex simulations of competitve pooled growth of a population of genotypes.
3. FitSeq.m calculates the fitness of each genotype from read-count time-series data.


## Installing

You must have MATLAB (version 2018a or newer) installed that includes the Optimization Toolbox. Download all files and make sure they are in your MATLAB path.


## Getting Started

A walk-through of how to perform each simulation and FitSeq fitness estimates on these simulations is provided in walk_through.m

### Simplified Simulations

Models competative pooled growth of a population of genotypes with different fitnesses. This simulation include only growth noise, but not noise from cell transfers, DNA extraction, PCR, or sequencing.

#### INPUTS

+ lineage: number of genotypes of the population

+ t_evo: total number of generations grown

+ cell_num_ini: a vector of the initial cell number of each genotype at generation 0, size = lineage * 1

+ x_ini: a vector of the fitness of each genotype, size = lineage * 1

+ noise_option: options of whether cell growth noise is simulated, logical (0-1) scaler value, 1 means that the cell growth noise is included and 0 means that the cell growth noise is not included  

+ 'format': optional, file format of the output file, 'csv'(default) or 'mat'

+ 'platform': optional, the platform that run Fit-Seq, 'MATLAB (default)' or 'OCTAVE'


#### OUTPUTS

+ file_name: the name of the file(s) written by the function<br/>
    When 'format' is set to 'mat', output will be:<br/>
    - 'data_evo_simu_0mut_simplified_*Time*.mat'<br/>
    
    When 'format' is set to 'csv', output will be:<br/>
    - 'data_evo_simu_0mut_simplified_*Time*_MeanFitness.csv'<br/>
    - 'data_evo_simu_0mut_simplified_*Time*_CellNumber.csv'<br/>
    - 'data_evo_simu_0mut_simplified_*Time*_EffectiveCellDepth.csv'<br/>
    - 'data_evo_simu_0mut_simplified_*Time*_Parameters.csv' 

#### Example

```
lineage = 1e4;   % number of genotypes of the population
t_evo = 24;   % total number of generations grown
cell_num_ini = 1e2*ones(lineage,1);   % a vector of initial cell number of each genotype at generation 0
x_ini = random('Normal',0,0.2, [lineage,1]);   % a vector of the fitness of each genotype. Here is a Gaussian distribution with mean = 0 and standard devistion = 0.2
noise_option = 1;   % cell growth noise is simulated

% Execute
[ TestSimplified_csv ] = evolution_0mut_simplified(lineage, t_evo, cell_num_ini, x_ini, noise_option);
```

### Complex Simulations 

Models competative pooled growth of a population of genotypes with different fitnesses. This simulation may include many sources of noise, including growth noise, noise from cell transfers, DNA extraction, PCR, and sequencing.

#### INPUTS

+ lineage: number of genotypes of the population

+ t_evo: total number of generations grown

+ cell_num_ini: a vector of the initial cell number of each genotype at generation 0, size = lineage * 1

+ x_ini: a vector of the fitness of each genotype, size = lineage * 1

+ deltat: number of generations between two successive cell transfers

+ read_depth_average: average number of reads per genotype per sequencing time point

+ noise_option: a vector of whether five types of noise are simulated (cell growth, bottleneck transfer, DNA extraction, PCR, sequencing), size = 1*5 logical (0-1) vector,
    - 1 or 0 at the 1st position determines if the cell growth noise is included or not
    - 1 or 0 at the 2nd position determines if the bottleneck cell transfer noise is included or not
    - 1 or 0 at the 3rd position determines if the DNA extraction noise is included or not
    - 1 or 0 at the 4th position determines if the PCR noise is included or not
    - 1 or 0 at the 5th position determines if the sequencing noise is included or not
    
+  'format': optional, file format of the output file, 'csv'(default) or 'mat'

+ 'gDNA_copy': optional, average copy number of genome DNA per genotype as template in PCR, default value is 500

+ 'PCR_cycle': optional, number of cycles in PCR, default value is 25

+ 'platform': optional, the platform that run Fit-Seq, 'MATLAB (default)' or 'OCTAVE'


#### OUTPUTS

+  file_name: the name of the file(s) written by the function<br/>
    When 'format' is set to 'mat', output will be:<br/>
    - 'data_evo_simu_0mut_complex_*Time*.mat'<br/>
    
    When 'format' is set to 'csv', output will be:<br/> 
    - 'data_evo_simu_0mut_complex_*Time*_MeanFitness.csv'<br/>
    - 'data_evo_simu_0mut_complex_*Time*_CellNumber.csv'<br/>
    - 'data_evo_simu_0mut_complex_*Time*_EffectiveCellDepth.csv'<br/>
    - 'data_evo_simu_0mut_complex_*Time*_Parameters.csv'<br/>
    - 'data_evo_simu_0mut_complex_*Time*_SequencedTimepoints.csv'<br/>
    - 'data_evo_simu_0mut_complex_*Time*_Reads.csv'
                   
                   
#### Example

```
lineage = 1e4;   % number of genotypes of the population
t_evo = 24;   % total number of generations grown
cell_num_ini = 1e2*ones(lineage,1);  % a vector of initial cell number of each genotype at generation 0
x_ini = random('Normal',0,0.15, [lineage,1]);   % a vector of the fitness of each genotype. Here, this is a Gaussian distribution with mean = 0 and standard devistion = 0.15
read_depth_average = 100;   % average number of reads per genotype per sequencing
deltat = 8;   % number of generations between two successive cell transfers
noise_option = [1,1,1,1,1];   % a vector of whether to simulate each of five types of noise (cell growth, bottleneck transfer, DNA extraction, PCR, sequencing). Here, all five types of noise are simulated.
                             
% Execute
[ TestComplex_csv ] = evolution_0mut_complex(lineage, t_evo, cell_num_ini, x_ini, deltat, read_depth_average, noise_option);
```


### Fitness Estimation
Estimates the fitness of each genotype from read-count time-series data.


#### INPUTS
+ BC_num_mat_original: a matrix of read number of each genotype at each sequencing time point, size = number of genotypes * length(t_seq_vec). This is a required input.

+ t_seq_vec: a vector of all sequencing time points. Either t_seq_vec or cell_depth is a required input. Use [] as input if it is unknown. 

+ cell_depth: a matrix of the cell number across growth cycles. The first row of this matrix is the cell number after bottleneck (before growth) for each time point. The second row of this matrix is the cell number before bottleneck (after growth). t_seq_vec is calculated from this matrix, and length(t_seq_vec) = size(cell_depth,2) + 1 . If t_seq_vec is also supplied by the user, the t_seq_vec calculated from this matrix will be used instead. Use [] as input if it is unknown.

+ file_name: the name of the file(s) written by the function. Use [] as input if it is not given.<br/>
    When 'format' is set to 'mat', output will be:<br/>
    - *file_name*_Fit-Seq_result_*Time*.mat'<br/>

    When 'format' is set to 'csv', output will be:<br/>
    - *file_name*_Fit-Seq_result_EstimatedFitness_*Time*.csv'<br/>
    - *file_name*_Fit-Seq_result_EstimatedReads_*Time*.csv'<br/>

+ 'format': optional, file format of the output file, 'csv'(default) or 'mat'

+ 'kappa': optional, a noise parameter that characterizes the total noise introduced by growth, cell transfer, DNA extraction, PCR, and sequencing, default value is 2.5, from Levy et al. Nature 2015 519, 181-186. To measure kappa empirically, see that reference. 

+ 'opt_cycle_max': optional, maximum number of cycles used when using likelihood optimization method to estimate fitness, default value is 10. Increse this cycle number might increase the fitness estimation accuracy, but extend the compute time.

+ 'platform': optional, the platform that run Fit-Seq, 'MATLAB (default)' or 'OCTAVE'


#### OUTPUTS
+ x_estimate_result: a vector of the estimated fitness of each genotype

+ r_estimate_result: a matrix of the estimated read number of each genotype at each sequencing time point, size = number of genotypes * length(t_seq_vec)

+ x_mean_estimate_result: a vector of the estimated mean fitness of the population at each sequencing time point, size = 1 * length(t_seq_vec)


#### Inputting data from simulation

```
t_seq_vec = csvread([TestComplex_csv(1:end-4), '_SeuqencedTimepoints.csv']);   % a vector of all sequencing time points
BC_num_mat_original = csvread([TestComplex_csv(1:end-4), '_Reads.csv']);   % a matrix of the read number of each genotype at each sequencing time point
cell_depth = [];   % a matrix of the cell number after bottleneck (before growth, first row), and the cell number before bottleneck (after growth, second row), size = 2 * (length(t_seq_vec)-1). Use [] as input if it is unknown
file_name = 'Test';
```


#### Inputting data from a file (Simulated-Pooled-Growth_Reads_10000Genotypes.csv, on MATLAB) 

```
t_seq_vec = 0:8:24;   % a vector of all sequencing time points
BC_num_mat_original = csvread('Simulated-Pooled-Growth_Reads_10000Genotypes.csv');   % a matrix of read number of each genotype at each sequencing time point
cell_depth = [];   % a matrix of the cell number after bottleneck (before growth, first row), and the cell number before bottleneck (after growth, second row), size = 2 * (length(t_seq_vec)-1). Use [] as input if it is unknown.
file_name = 'Test';
```

#### Inputting data from a file (Simulated-Pooled-Growth_Reads_1000Genotypes.csv, GNU Octave) 

```
t_seq_vec = 0:8:24;   % a vector of all sequencing time points
BC_num_mat_original = csvread('Simulated-Pooled-Growth_Reads_10000enotypes.csv');   % a matrix of read number of each genotype at each sequencing time point
cell_depth = [];   % a matrix of the cell number after bottleneck (before growth, first row), and the cell number before bottleneck (after growth, second row), size = 2 * (length(t_seq_vec)-1). Use [] as input if it is unknown.
file_name = 'Test';
```

#### Running Fit-Seq (on MATLAB)
```
[x_estimate_result, r_estimate_result, x_mean_estimate_result] =  FitSeq(BC_num_mat_original, t_seq_vec, cell_depth, file_name);
```

#### Running Fit-Seq (on GNU Octave)
```
[x_estimate_result, r_estimate_result, x_mean_estimate_result] =  FitSeq(BC_num_mat_original, t_seq_vec, cell_depth, file_name, 'platform','OCTAVE');
```
