%% An expample of evolution_0mut_simplified.m (run on MATLAB)

% INPUTS
lineage = 1e4;    % number of genotypes of the population
t_evo = 24;   % total number of generations grown
cell_num_ini = 1e2*ones(lineage,1);   % a vector of initial cell number of each genotype at generation 0
x_ini = random('Normal',0,0.2, [lineage,1]);   % a vector of the fitness of each genotype. Here is a Gaussian distribution with mean = 0 and standard devistion = 0.2
noise_option = 1;   % cell growth noise is simulated

% Execute
[ TestSimplified_csv ] = evolution_0mut_simplified(lineage, t_evo, cell_num_ini, x_ini, noise_option);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% An expample of evolution_0mut_complex.m (run on MATLAB)

% INPUTS
lineage = 1e4;   % number of genotypes of the population
t_evo = 24;   % total number of generations grown
cell_num_ini = 1e2*ones(lineage,1);  % a vector of initial cell number of each genotype at generation 0
x_ini = random('Normal',0,0.15, [lineage,1]);   % a vector of the fitness of each genotype. Here, this is a Gaussian distribution with mean = 0 and standard devistion = 0.15
read_depth_average = 100;   % average number of reads per genotype per sequencing
deltat = 8;   % number of generations between successive cell transfers
noise_option = [1,1,1,1,1];   % a vector of whether to simulate each of five types of noise (cell growth, bottleneck transfer, DNA extraction, PCR, sequencing). Here, all five types of noise are simulated.
                             
% Execute
[ TestComplex_csv ] = evolution_0mut_complex(lineage, t_evo, cell_num_ini, x_ini, deltat, read_depth_average, noise_option);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% An expample of Fit-Seq.m (run on MATLAB)

% There are two way of inputting
%% The first is by loading data from simulation at previouse step,
t_seq_vec = csvread([TestComplex_csv(1:end-4), '_SeuqencedTimepoints.csv']);   % a vector of all sequencing time points
BC_num_mat_original = csvread([TestComplex_csv(1:end-4), '_Reads.csv']);   % a matrix of the read number of each genotype at each sequencing time point
cell_depth = [];   % a matrix of the cell number after bottleneck (before growth, first row), and the cell number before bottleneck (after growth, second row). size = 2 * (length(t_seq_vec)-1). Use [] as input if it is unknown.
file_name = 'Test';

%% The second is by loading'Simulated-Pooled_Growth_Reads_10000Genotypes.csv' provided
t_seq_vec = 0:8:24;   % a vector of all sequencing time points
BC_num_mat_original = csvread('Simulated-Pooled-Growth_Reads_10000Genotypes.csv');   % a matrix of read number of each genotype at each sequencing time point,
cell_depth = [];  % a matrix of the cell number after bottleneck (before growth, first row), and the cell number before bottleneck (after growth, second row). size = 2 * (length(t_seq_vec)-1). Use [] as input if it is unknown.
file_name = 'Test';

%%
% Execute
[x_estimate_result, r_estimate_result, x_mean_estimate_result] =  FitSeq(BC_num_mat_original, t_seq_vec, cell_depth, file_name);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%% An expample of evolution_0mut_simplified.m (run on MATLAB)

% INPUTS
lineage = 1e4;    % number of genotypes of the population
t_evo = 24;   % total number of generations grown
cell_num_ini = 1e2*ones(lineage,1);   % a vector of initial cell number of each genotype at generation 0
x_ini = random('Normal',0,0.2, [lineage,1]);   % a vector of the fitness of each genotype. Here is a Gaussian distribution with mean = 0 and standard devistion = 0.2
noise_option = 1;   % cell growth noise is simulated

% Execute
[ TestSimplified_csv ] = evolution_0mut_simplified(lineage, t_evo, cell_num_ini, x_ini, noise_option);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%
% Need to install and load the following packages before running GNU Octave
pkg load optim
pkg load struct
pkg load statistics
pkg load io
pkg load nan

%% An expample of evolution_0mut_simplified.m (run on GNU Octave)

% INPUTS
lineage = 1e3;    % number of genotypes of the population
t_evo = 24;   % total number of generations grown
cell_num_ini = 1e2*ones(lineage,1);   % a vector of initial cell number of each genotype at generation 0
x_ini = random('Normal',0,0.2, [lineage,1]);   % a vector of the fitness of each genotype. Here is a Gaussian distribution with mean = 0 and standard devistion = 0.2
noise_option = 1;   % cell growth noise is simulated

% Execute
[ TestSimplified_csv ] = evolution_0mut_simplified(lineage, t_evo, cell_num_ini, x_ini, noise_option,'platform','OCTAVE');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% An expample of evolution_0mut_complex.m (run on GNU Octave)

% INPUTS
lineage = 1e3;   % number of genotypes of the population
t_evo = 24;   % total number of generations grown
cell_num_ini = 1e2*ones(lineage,1);  % a vector of initial cell number of each genotype at generation 0
x_ini = random('Normal',0,0.15, [lineage,1]);   % a vector of the fitness of each genotype. Here, this is a Gaussian distribution with mean = 0 and standard devistion = 0.15
read_depth_average = 100;   % average number of reads per genotype per sequencing
deltat = 8;   % number of generations between successive cell transfers
noise_option = [1,1,1,1,1];   % a vector of whether to simulate each of five types of noise (cell growth, bottleneck transfer, DNA extraction, PCR, sequencing). Here, all five types of noise are simulated.
                             
% Execute
[ TestComplex_csv ] = evolution_0mut_complex(lineage, t_evo, cell_num_ini, x_ini, deltat, read_depth_average, noise_option,'platform','OCTAVE');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% An expample of Fit-Seq.m (run on GNU Octave)

% There are two way of inputting
%% The first is by loading data from simulation at previouse step,
t_seq_vec = csvread([TestComplex_csv(1:end-4), '_SeuqencedTimepoints.csv']);   % a vector of all sequencing time points
BC_num_mat_original = csvread([TestComplex_csv(1:end-4), '_Reads.csv']);   % a matrix of the read number of each genotype at each sequencing time point
cell_depth = [];   % a matrix of the cell number after bottleneck (before growth, first row), and the cell number before bottleneck (after growth, second row). size = 2 * (length(t_seq_vec)-1). Use [] as input if it is unknown.
file_name = 'Test';

%% The second is by loading'Simulated-Pooled_Growth_Reads_1000Genotypes.csv' provided
t_seq_vec = 0:8:24;   % a vector of all sequencing time points
BC_num_mat_original = csvread('Simulated-Pooled-Growth_Reads_1000Genotypes.csv');   % a matrix of read number of each genotype at each sequencing time point,
cell_depth = [];  % a matrix of the cell number after bottleneck (before growth, first row), and the cell number before bottleneck (after growth, second row). size = 2 * (length(t_seq_vec)-1). Use [] as input if it is unknown.
file_name = 'Test';

%%
% Execute
[x_estimate_result, r_estimate_result, x_mean_estimate_result] =  FitSeq(BC_num_mat_original, t_seq_vec, cell_depth, file_name,'platform','OCTAVE');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%