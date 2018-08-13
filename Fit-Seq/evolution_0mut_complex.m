function [ file_name ] = evolution_0mut_complex(lineage, t_evo, ...
    cell_num_ini, x_ini, deltat, read_depth_average, noise_option, varargin)
% -------------------------------------------------------------------------
% evolution_0mut_complex
% COMPLEX VERSION OF SIMULATED COMPETETIVE POOLED GROWTH OF A POPULATION OF
% GENOTYPES WITH DIFFERENT FITNESSES. THESE SIMULATIONS INCLUDE EXPERIMENTAL 
% NOISE SUCH AS GROWTH NOISE, SAMPLING DURING BOTTLENECKS, DNA EXTRACTION, 
% PCR, AND SAMPLING ON SEQUENCER.
%
%
% INPUTS
% -- lineage: number of genotypes of the population
%
% -- t_evo: total number of generations grown
%
% -- cell_num_ini: a vector of the initial cell number of each genotype at
%                  generation 0,
%                  size = lineage * 1
%
% -- x_ini: a vector of the fitness of each genotype,
%           size = lineage * 1
%
% -- deltat: number of generations between two successive cell transfers
%
% -- read_depth_average: average number of reads per genotype per
%                        sequencing time point
%
% -- noise_option: a vector of whether five types of noise are simulated
%                  (cell growth, bottleneck transfer, DNA extraction, PCR, sequencing)
%                  size = 1*5 logical (0-1) vector,
%                  1 or 0 at the 1st position determines if the cell growth noise is included or not
%                  1 or 0 at the 2nd position determines if the bottleneck cell transfer noise is included or not
%                  1 or 0 at the 3rd position determines if the DNA extraction noise is included or not
%                  1 or 0 at the 4th position determines if the PCR noise is included or not
%                  1 or 0 at the 5th position determines if the sequencing noise is included or not

% -- 'format': optional, file format of the output file, 'csv'(default) or 'mat'
%
% -- 'gDNA_copy': optional, average copy number of genome DNA per genotype
%                 as template in PCR, default value is 500
%
% -- 'PCR_cycle': optional, number of cycles in PCR, default value is 25
% 
% -- 'platform': optional, the platform that run Fit-Seq, 'MATLAB (default)' or 'OCTAVE'
%
%
% OUTPUTS
% -- file_name: the name of the file(s) written by the function.
%               When 'format' is set to 'mat', output will be:
%                  'data_evo_simu_0mut_complex_*Time*.mat'
%               When 'format' is set to 'csv', output will be:
%                   'data_evo_simu_0mut_complex_*Time*_MeanFitness.mat'
%                   'data_evo_simu_0mut_complex_*Time*_CellNumber.mat'
%                   'data_evo_simu_0mut_complex_*Time*_EffectiveCellDepth.mat'
%                   'data_evo_simu_0mut_complex_*Time*_Parameters.mat'
%                   'data_evo_simu_0mut_complex_*Time*_SequencedTimepoints.mat'
%                   'data_evo_simu_0mut_complex_*Time*_Reads.mat'
%
% -------------------------------------------------------------------------
% Parse inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numvarargs = length(varargin);
if numvarargs > 8
    error('evolution_0mut_complex:TooManyInputs', 'requires at most 4 optional inputs');
end

vararg_update = {'format','csv', 'gDNA_copy',500, 'PCR_cycle',25, 'platform','MATLAB'};
pos1 = find(strcmpi(varargin,'format'))+1;
pos2 = find(strcmpi(varargin,'gDNA_copy'))+1;
pos3 = find(strcmpi(varargin,'PCR_cycle'))+1;
pos4 = find(strcmpi(varargin,'platform'))+1;
if numvarargs/2 == 4
    vararg_update([2,4,6,8]) = varargin([pos1,pos2,pos3,pos4]);
elseif numvarargs/2 == 3
    if isempty(pos1)
    vararg_update([4,6,8]) = varargin([pos2,pos3,pos4]);
    elseif isempty(pos2)
        vararg_update([2,6,8]) = varargin([pos1,pos3,pos4]);
    elseif isempty(pos3)
        vararg_update([2,4,8]) = varargin([pos1,pos2,pos4]);
    elseif isempty(pos4)    
        vararg_update([2,4,6]) = varargin([pos1,pos2,pos3]);
    end  
elseif numvarargs/2 == 2
    if isempty(pos1) && isempty(pos2)
        vararg_update([6,8]) = varargin([pos3,pos4]);
    elseif isempty(pos1) && isempty(pos3)
        vararg_update([4,8]) = varargin([pos2,pos4]);
    elseif isempty(pos1) && isempty(pos4)
        vararg_update([4,6]) = varargin([pos2,pos3]);
    elseif isempty(pos2) && isempty(pos3)
        vararg_update([2,8]) = varargin([pos1,pos4]);
    elseif isempty(pos2) && isempty(pos4)
        vararg_update([2,6]) = varargin([pos1,pos3]);
    elseif isempty(pos3) && isempty(pos4)
        vararg_update([2,4]) = varargin([pos1,pos2]);
    end
elseif numvarargs/2 == 1
    if ~isempty(pos1)
        vararg_update(2) = varargin(pos1);
    elseif ~isempty(pos2)
        vararg_update(4) = varargin(pos2);
    elseif ~isempty(pos3)
        vararg_update(6) = varargin(pos3);
    elseif ~isempty(pos4)
        vararg_update(8) = varargin(pos4);
    end
end

ouptput_format = vararg_update{2};
gDNA_copy_number = vararg_update{4};
PCR_cycle_number = vararg_update{6};
platform_choice = vararg_update{8};

noise_growth = noise_option(1);
noise_bottleneck_transfer = noise_option(2);
noise_genome_DNA_extraction = noise_option(3);
noise_barcode_PCR = noise_option(4);
noise_sequencing = noise_option(5);

cell_num_evo = zeros(lineage, t_evo+1);
cell_num_evo(:,1) = cell_num_ini;
cell_num_evo_flask = cell_num_evo; % number of cells transferred at bottleneck
x_mean = nan(1,t_evo+1);
x_mean(1) = cell_num_ini'*x_ini/sum(cell_num_ini);

t_seq_vec = 0:deltat:t_evo;


% Simulate pooled growth %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Calulate the minimum cell number that is expected for a genotype
%         that is not growing and being lost to dilution
cell_min = zeros(lineage, t_evo+1);
cell_min(:,1) = cell_num_ini;
if noise_bottleneck_transfer == 0
    for j1 = 2:(t_evo+1)
        cell_min(:,j1) = round(cell_num_ini/(2^(floor((j1-2)/deltat)*deltat)));
    end
elseif noise_bottleneck_transfer == 1
    for j1 = 2:(t_evo+1)
        cell_min(:,j1) = poissrnd(cell_num_ini/(2^(floor((j1-2)/deltat)*deltat)));
    end
end

% Step 2: Growth under four possible conditions:
%         1. Without growth noise and without sampling noise during bottlenecks
%         2. With growth noise but without sampling noise during bottlenecks
%         3. Without growth noise but with sampling noise during bottlenecks
%         4. With growth noise and with sampling noise during bottlenecks

tstart0 = tic;
if noise_growth==0 && noise_bottleneck_transfer==0
    if deltat>=2
        for j = 2:deltat
            if mod(j,20) == 1
                fprintf('Current generation: %i\n', j-1)
            end
            pos = cell_num_evo(:,j-1)~=0;
            n_prev = cell_num_evo(pos,j-1);
            x_ini_rela = max((1+x_ini(pos))/(1+x_mean(j-1)),0);
            cell_num_evo(pos,j) = max(round(2*n_prev.*x_ini_rela),cell_min(pos,j));
            cell_num_evo_flask(:,j) = cell_num_evo(:,j);
            x_mean(j) = cell_num_evo(:,j)'*x_ini/sum(cell_num_evo(:,j));
        end
        for j = (deltat+1):(t_evo+1)
            if mod(j,20) == 1
                fprintf('Current generation: %i\n', j-1)
            end
            if mod(j-1,deltat) == 0
                pos = cell_num_evo(:,j-1)~=0;
                n_prev = cell_num_evo(pos,j-1);
                x_ini_rela = max((1+x_ini(pos))/(1+x_mean(j-1)),0);
                cell_num_evo(pos,j) = max(round(2*n_prev.*x_ini_rela),cell_min(pos,j));
                cell_num_evo_flask(:,j) = round(cell_num_evo(:,j)/2^deltat);
            elseif mod(j-1,deltat) ~= 0
                pos = cell_num_evo_flask(:,j-1)~=0;
                n_prev = cell_num_evo_flask(pos,j-1);
                x_ini_rela = max((1+x_ini(pos))/(1+x_mean(j-1)),0);
                cell_num_evo(pos,j) = max(round(2*n_prev.*x_ini_rela),cell_min(pos,j));
                cell_num_evo_flask(:,j) = cell_num_evo(:,j);
            end
            x_mean(j) = cell_num_evo(:,j)'*x_ini/sum(cell_num_evo(:,j));
        end
    elseif deltat == 1
        for j = 2:(t_evo+1)
            if mod(j,20) == 1
                fprintf('Current generation: %i\n', j-1)
            end
            pos = cell_num_evo_flask(:,j-1)~=0;
            n_prev = cell_num_evo_flask(pos,j-1);
            x_ini_rela = max((1+x_ini(pos))/(1+x_mean(j-1)),0);
            cell_num_evo(pos,j) = max(round(2*n_prev.*x_ini_rela),cell_min(pos,j));
            cell_num_evo_flask(:,j) = round(cell_num_evo(:,j)/2);
            x_mean(j) = cell_num_evo(:,j)'*x_ini/sum(cell_num_evo(:,j));
        end
    end
    
    
elseif noise_growth==1 && noise_bottleneck_transfer==0
    if deltat>=2
        for j = 2:deltat
            if mod(j,20) == 1
                fprintf('Current generation: %i\n', j-1)
            end
            pos = cell_num_evo(:,j-1)~=0;
            n_prev = cell_num_evo(pos,j-1);
            x_ini_rela = max((1+x_ini(pos))/(1+x_mean(j-1)),0);
            cell_num_evo(pos,j) = max(poissrnd(2*n_prev.*x_ini_rela),cell_min(pos,j));
            cell_num_evo_flask(:,j) = cell_num_evo(:,j);
            x_mean(j) = cell_num_evo(:,j)'*x_ini/sum(cell_num_evo(:,j));
        end
        for j = (deltat+1):(t_evo+1)
            if mod(j,20) == 1
                fprintf('Current generation: %i\n', j-1)
            end
            if mod(j-1,deltat) == 0
                pos = cell_num_evo(:,j-1)~=0;
                n_prev = cell_num_evo(pos,j-1);
                x_ini_rela = max((1+x_ini(pos))/(1+x_mean(j-1)),0);
                cell_num_evo(pos,j) = max(poissrnd(2*n_prev.*x_ini_rela),cell_min(pos,j));
                cell_num_evo_flask(:,j) = round(cell_num_evo(:,j)/2^deltat);
            elseif mod(j-1,deltat) ~= 0
                pos = cell_num_evo_flask(:,j-1)~=0;
                n_prev = cell_num_evo_flask(pos,j-1);
                x_ini_rela = max((1+x_ini(pos))/(1+x_mean(j-1)),0);
                cell_num_evo(pos,j) = max(poissrnd(2*n_prev.*x_ini_rela),cell_min(pos,j));
                cell_num_evo_flask(:,j) = cell_num_evo(:,j);
            end
            x_mean(j) = cell_num_evo(:,j)'*x_ini/sum(cell_num_evo(:,j));
        end
    elseif deltat == 1
        for j = 2:(t_evo+1)
            if mod(j,20) == 1
                fprintf('Current generation: %i\n', j-1)
            end
            pos = cell_num_evo_flask(:,j-1)~=0;
            n_prev = cell_num_evo_flask(pos,j-1);
            x_ini_rela = max((1+x_ini(pos))/(1+x_mean(j-1)),0);
            cell_num_evo(pos,j) = max(poissrnd(2*n_prev.*x_ini_rela),cell_min(pos,j));
            cell_num_evo_flask(:,j) = round(cell_num_evo(:,j)/2);
            x_mean(j) = cell_num_evo(:,j)'*x_ini/sum(cell_num_evo(:,j));
        end
    end
    
    
elseif noise_growth==0 && noise_bottleneck_transfer==1
    if deltat >=2
        for j = 2:deltat
            if mod(j,20) == 1
                fprintf('Current generation: %i\n', j-1)
            end
            pos = cell_num_evo(:,j-1)~=0;
            n_prev = cell_num_evo(pos,j-1);
            x_ini_rela = max((1+x_ini(pos))/(1+x_mean(j-1)),0);
            cell_num_evo(pos,j) = max(round(2*n_prev.*x_ini_rela),cell_min(pos,j));
            cell_num_evo_flask(:,j) = cell_num_evo(:,j);
            x_mean(j) = cell_num_evo(:,j)'*x_ini/sum(cell_num_evo(:,j));
        end
        for j = (deltat+1):(t_evo+1)
            if mod(j,20) == 1
                fprintf('Current generation: %i\n', j-1)
            end
            if mod(j-1,deltat) == 0
                pos = cell_num_evo(:,j-1)~=0;
                n_prev = cell_num_evo(pos,j-1);
                x_ini_rela = max((1+x_ini(pos))/(1+x_mean(j-1)),0);
                cell_num_evo(pos,j) = max(round(2*n_prev.*x_ini_rela),cell_min(pos,j));
                cell_num_evo_flask(:,j) = poissrnd(cell_num_evo(:,j)/2^deltat);
            elseif mod(j-1,deltat) ~= 0
                pos = cell_num_evo_flask(:,j-1)~=0;
                n_prev = cell_num_evo_flask(pos,j-1);
                x_ini_rela = max((1+x_ini(pos))/(1+x_mean(j-1)),0);
                cell_num_evo(pos,j) = max(round(2*n_prev.*x_ini_rela),cell_min(pos,j));
                cell_num_evo_flask(:,j) = cell_num_evo(:,j);
            end
            x_mean(j) = cell_num_evo(:,j)'*x_ini/sum(cell_num_evo(:,j));
        end
    elseif deltat == 1
        for j = 2:(t_evo+1)
            if mod(j,20) == 1
                fprintf('Current generation: %i\n', j-1)
            end
            pos = cell_num_evo_flask(:,j-1)~=0;
            n_prev = cell_num_evo_flask(pos,j-1);
            x_ini_rela = max((1+x_ini(pos))/(1+x_mean(j-1)),0);
            cell_num_evo(pos,j) = max(round(2*n_prev.*x_ini_rela),cell_min(pos,j));
            cell_num_evo_flask(:,j) = poissrnd(cell_num_evo(:,j)/2);
            x_mean(j) = cell_num_evo(:,j)'*x_ini/sum(cell_num_evo(:,j));
        end
    end
    
    
elseif noise_growth==1 && noise_bottleneck_transfer==1
    if deltat >=2
        for j = 2:deltat
            if mod(j,20) == 1
                fprintf('Current generation: %i\n', j-1)
            end
            pos = cell_num_evo(:,j-1)~=0;
            n_prev = cell_num_evo(pos,j-1);
            x_ini_rela = max((1+x_ini(pos))/(1+x_mean(j-1)),0);
            cell_num_evo(pos,j) = max(poissrnd(2*n_prev.*x_ini_rela),cell_min(pos,j));
            cell_num_evo_flask(:,j) = cell_num_evo(:,j);
            x_mean(j) = cell_num_evo(:,j)'*x_ini/sum(cell_num_evo(:,j));
        end
        for j = (deltat+1):(t_evo+1)
            if mod(j,20) == 1
                fprintf('Current generation: %i\n', j-1)
            end
            if mod(j-1,deltat) == 0
                pos = cell_num_evo(:,j-1)~=0;
                n_prev = cell_num_evo(pos,j-1);
                x_ini_rela = max((1+x_ini(pos))/(1+x_mean(j-1)),0);
                cell_num_evo(pos,j) = max(poissrnd(2*n_prev.*x_ini_rela),cell_min(pos,j));
                cell_num_evo_flask(:,j) = poissrnd(cell_num_evo(:,j)/2^deltat);
            elseif mod(j-1,deltat) ~= 0
                pos = cell_num_evo_flask(:,j-1)~=0;
                n_prev = cell_num_evo_flask(pos,j-1);
                x_ini_rela = max((1+x_ini(pos))/(1+x_mean(j-1)),0);
                cell_num_evo(pos,j) = max(poissrnd(2*n_prev.*x_ini_rela),cell_min(pos,j));
                cell_num_evo_flask(:,j) = cell_num_evo(:,j);
            end
            x_mean(j) = cell_num_evo(:,j)'*x_ini/sum(cell_num_evo(:,j));
        end
    elseif deltat == 1
        for j = 2:(t_evo+1)
            if mod(j,20) == 1
                fprintf('Current generation: %i\n', j-1)
            end
            pos = cell_num_evo_flask(:,j-1)~=0;
            n_prev = cell_num_evo_flask(pos,j-1);
            x_ini_rela = max((1+x_ini(pos))/(1+x_mean(j-1)),0);
            cell_num_evo(pos,j) = max(poissrnd(2*n_prev.*x_ini_rela),cell_min(pos,j));
            cell_num_evo_flask(:,j) = poissrnd(cell_num_evo(:,j)/2);
            x_mean(j) = cell_num_evo(:,j)'*x_ini/sum(cell_num_evo(:,j));
        end
    end
end

effective_cell_depth = sum(cell_num_evo_flask(:,t_seq_vec+1));

telaps0=toc(tstart0);
fprintf('Computing time for %i generations: %f seconds.\n', t_evo, telaps0)


% Step 3: Sample prep under eight possible conditions:
%         1. Without DNA extraction noise, without PCR noise, and without sequencing noise
%         2. With DNA extraction noise, without PCR noise, and without sequencing noise
%         3. Without DNA extraction noise, with PCR noise, and without sequencing noise
%         4. Without DNA extraction noise, without PCR noise, and with sequencing noise
%         5. With DNA extraction noise, with PCR noise, and without sequencing noise
%         6. With DNA extraction noise, without PCR noise, and with sequencing noise
%         7. Without DNA extraction noise, with PCR noise, and with sequencing noise
%         8. With DNA extraction noise, with PCR noise, and with sequencing noise
cell_num_mat_data = cell_num_evo(:,t_seq_vec+1);
genomeDNA_percentage = gDNA_copy_number/(mean(cell_num_ini)*2^deltat);
PCR_percentage = 0.01; % PCR_percentage: percantage of total PCR product
% sent for sequencing

if noise_genome_DNA_extraction==0 && noise_barcode_PCR==0 && ...
        noise_sequencing==0
    num_mat_DNA = floor(cell_num_mat_data.*genomeDNA_percentage);
    num_mat_DNA(:,1) = floor(cell_num_mat_data(:,1)*2^deltat.*...
        genomeDNA_percentage);
    num_mat_PCR_tempt = num_mat_DNA;
    for i1 = 1:PCR_cycle_number
        num_mat_PCR_tempt = floor(2*num_mat_PCR_tempt);
    end
    num_mat_PCR = floor(num_mat_PCR_tempt.*PCR_percentage);
    num_mat_sequencing = floor(num_mat_PCR*read_depth_average*lineage./...
        sum(num_mat_PCR));
    
elseif noise_genome_DNA_extraction==1 && noise_barcode_PCR==0 && ...
        noise_sequencing==0
    num_mat_DNA = poissrnd(cell_num_mat_data.*genomeDNA_percentage);
    num_mat_DNA(:,1) = poissrnd(cell_num_mat_data(:,1)*2^deltat.*...
        genomeDNA_percentage);
    num_mat_PCR_tempt = num_mat_DNA;
    for i1 = 1:PCR_cycle_number
        num_mat_PCR_tempt = floor(2*num_mat_PCR_tempt);
    end
    num_mat_PCR = floor(num_mat_PCR_tempt.*PCR_percentage);
    num_mat_sequencing = floor(num_mat_PCR*read_depth_average*lineage./...
        sum(num_mat_PCR));
    
elseif noise_genome_DNA_extraction==0 && noise_barcode_PCR==1 && ...
        noise_sequencing==0
    num_mat_DNA = floor(cell_num_mat_data.*genomeDNA_percentage);
    num_mat_DNA(:,1) = floor(cell_num_mat_data(:,1)*2^deltat.*...
        genomeDNA_percentage);
    num_mat_PCR_tempt = num_mat_DNA;
    for i1 = 1:PCR_cycle_number
        num_mat_PCR_tempt = poissrnd(2*num_mat_PCR_tempt);
    end
    num_mat_PCR = poissrnd(num_mat_PCR_tempt.*PCR_percentage);
    num_mat_sequencing = floor(num_mat_PCR*read_depth_average*lineage./...
        sum(num_mat_PCR));
    
elseif noise_genome_DNA_extraction==0 && noise_barcode_PCR==0 && ...
        noise_sequencing==1
    num_mat_DNA = floor(cell_num_mat_data.*genomeDNA_percentage);
    num_mat_DNA(:,1) = floor(cell_num_mat_data(:,1)*2^deltat.*...
        genomeDNA_percentage);
    num_mat_PCR_tempt = num_mat_DNA;
    for i1 = 1:PCR_cycle_number
        num_mat_PCR_tempt = floor(2*num_mat_PCR_tempt);
    end
    num_mat_PCR = floor(num_mat_PCR_tempt.*PCR_percentage);
    num_mat_sequencing = poissrnd(num_mat_PCR*read_depth_average*lineage./...
        sum(num_mat_PCR));
    
elseif noise_genome_DNA_extraction==1 && noise_barcode_PCR==1 && ...
        noise_sequencing==0
    num_mat_DNA = poissrnd(cell_num_mat_data.*genomeDNA_percentage);
    num_mat_DNA(:,1) = poissrnd(cell_num_mat_data(:,1)*2^deltat.*...
        genomeDNA_percentage);
    num_mat_PCR_tempt = num_mat_DNA;
    for i1 = 1:PCR_cycle_number
        num_mat_PCR_tempt = poissrnd(2*num_mat_PCR_tempt);
    end
    num_mat_PCR = poissrnd(num_mat_PCR_tempt.*PCR_percentage);
    num_mat_sequencing = floor(num_mat_PCR*read_depth_average*lineage./...
        sum(num_mat_PCR));
    
elseif noise_genome_DNA_extraction==1 && noise_barcode_PCR==0 && ...
        noise_sequencing==1
    num_mat_DNA = poissrnd(cell_num_mat_data.*genomeDNA_percentage);
    num_mat_DNA(:,1) = poissrnd(cell_num_mat_data(:,1)*2^deltat.*...
        genomeDNA_percentage);
    num_mat_PCR_tempt = num_mat_DNA;
    for i1 = 1:PCR_cycle_number
        num_mat_PCR_tempt = floor(2*num_mat_PCR_tempt);
    end
    num_mat_PCR = floor(num_mat_PCR_tempt.*PCR_percentage);
    num_mat_sequencing = poissrnd(num_mat_PCR*read_depth_average*lineage./...
        sum(num_mat_PCR));
    
elseif noise_genome_DNA_extraction==0 && noise_barcode_PCR==1 && ...
        noise_sequencing==1
    num_mat_DNA = floor(cell_num_mat_data.*genomeDNA_percentage);
    num_mat_DNA(:,1) = floor(cell_num_mat_data(:,1)*2^deltat.*...
        genomeDNA_percentage);
    num_mat_PCR_tempt = num_mat_DNA;
    for i1 = 1:PCR_cycle_number
        num_mat_PCR_tempt = poissrnd(2*num_mat_PCR_tempt);
    end
    num_mat_PCR = poissrnd(num_mat_PCR_tempt.*PCR_percentage);
    num_mat_sequencing = poissrnd(num_mat_PCR*read_depth_average*lineage./...
        sum(num_mat_PCR));
    
elseif noise_genome_DNA_extraction==1 && noise_barcode_PCR==1 && ...
        noise_sequencing==1
    num_mat_DNA = poissrnd(cell_num_mat_data.*genomeDNA_percentage);
    num_mat_DNA(:,1) = poissrnd(cell_num_mat_data(:,1)*2^deltat.*...
        genomeDNA_percentage);
    num_mat_PCR_tempt = num_mat_DNA;
    for i1 = 1:PCR_cycle_number
        num_mat_PCR_tempt = poissrnd(2*num_mat_PCR_tempt);
    end
    num_mat_PCR = poissrnd(num_mat_PCR_tempt.*PCR_percentage);
    num_mat_sequencing = poissrnd(num_mat_PCR*read_depth_average*lineage./...
        sum(num_mat_PCR));
end
% End of the main function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Report Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = datestr(now,'yyyymmdd-HHMMSSFFF');
switch lower(ouptput_format)
    case {'mat'}
        file_name = ['data_evo_simu_0mut_complex_' dt '.mat'];
        switch lower(platform_choice)
            case{'matlab'}
                save(file_name, 'num_mat_sequencing','t_seq_vec',...
                    'effective_cell_depth','cell_num_evo','x_mean',...
                    'x_ini','cell_num_ini','lineage','t_evo','deltat',...
                    'read_depth_average','noise_option', ...
                    'gDNA_copy_number', 'PCR_cycle_number')
            case{'octave'}
                save('-mat7-binary', file_name, 'num_mat_sequencing','t_seq_vec',...
                    'effective_cell_depth','cell_num_evo','x_mean',...
                    'x_ini','cell_num_ini','lineage','t_evo','deltat',...
                    'read_depth_average','noise_option', ...
                    'gDNA_copy_number', 'PCR_cycle_number')
        end

    case {'csv'}
        file_name_1 = ['data_evo_simu_0mut_complex_' dt '_Reads.csv'];
        file_name_2 = ['data_evo_simu_0mut_complex_' dt '_SeuqencedTimePoints.csv'];
        file_name_3 = ['data_evo_simu_0mut_complex_' dt '_EffectiveCellDepth.csv'];
        file_name_4 = ['data_evo_simu_0mut_complex_' dt '_CellNumber.csv'];
        file_name_5 = ['data_evo_simu_0mut_complex_' dt '_MeanFitness.csv'];
        file_name_6 = ['data_evo_simu_0mut_complex_' dt '_Paramaters.csv'];
             
        csvwrite(file_name_1,num_mat_sequencing)
        csvwrite(file_name_2,t_seq_vec)
        csvwrite(file_name_3,effective_cell_depth)
        csvwrite(file_name_4,cell_num_evo)
        csvwrite(file_name_5,x_mean)
        output_parameters = cell(lineage+1,9);
        output_parameters{1,1} = 'Fitness of each genotype (x_i)';
        output_parameters{1,2} = ...
            'Initial cell number of each genotype (n0_i)';
        output_parameters{1,3} = 'Number of genotypes (L)';
        output_parameters{1,4} = 'Total number of generations grown (T)';
        output_parameters{1,5} = ...
            'Number of generations between successive cell transfers (Delta t)';
        output_parameters{1,6} = ...
            'Average read number per genotype per time point (R/L)';
        output_parameters{1,7} = 'Average copy number of genome DNA per genotype';
        output_parameters{1,8} = 'Number of cycles in PCR';
        output_parameters{1,9} = 'Optional noise (noise_option)';
        
        output_parameters(2:end,1) = num2cell(x_ini);
        output_parameters(2:end,2) = num2cell(cell_num_ini);
        output_parameters{2,3} = lineage;
        output_parameters{2,4} = t_evo;
        output_parameters{2,5} = deltat;
        output_parameters{2,6} = read_depth_average;
        output_parameters{2,7} = gDNA_copy_number;
        output_parameters{2,8} = PCR_cycle_number;
        output_parameters{2,9} = ['[' num2str(noise_option) ']'];
        fileID = fopen(file_name_6,'wt');
        for k1 = 1:(size(output_parameters,2)-1)
            fprintf(fileID,'%s,', output_parameters{1,k1});
        end
        fprintf(fileID,'%s\n', output_parameters{1,size(output_parameters,2)});
        for k2 = 1:(size(output_parameters,2)-1)
            fprintf(fileID,'%f,', output_parameters{2,k2});
        end
        fprintf(fileID,'%s\n', output_parameters{2,size(output_parameters,2)});
        for k3 = 2:lineage
            fprintf(fileID,'%f,', output_parameters{k3+1,1});
            fprintf(fileID,'%f\n', output_parameters{k3+1,2});
        end
        fclose(fileID);     
        file_name = ['data_evo_simu_0mut_complex_' dt '.csv'];   
end
end
