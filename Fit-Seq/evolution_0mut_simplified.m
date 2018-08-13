function [ file_name ] = evolution_0mut_simplified(lineage, t_evo, ...
    cell_num_ini, x_ini, noise_option, varargin)
% -------------------------------------------------------------------------
% evolution_0mut_simplified
% SIMPLIFIED VERSION OF SIMULATED COMPETETIVE POOLED GROWTH OF A POPULATION
% OF GENOTYPES WITH DIFFERENT FITNESSES. THESE SIMULATIONS INCLUDE ONLY GROWTH
% NOISE, BUT NOT NOISE FROM CELL TRANSFERS, DNA EXTRACTION, PCR OR SEQUENCING.
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
% -- noise_option: options of whether cell growth noise is simulated,
%                  logical (0-1) scaler value,
%                  1 means that the cell growth noise is included
%                  0 means that the cell growth noise is not included
%
% -- 'format': optional, file format of the output file, 'csv'(default) or 'mat'
%
% -- 'platform': optional, the platform that run Fit-Seq, 'MATLAB (default)' or 'OCTAVE'
%
%
% OUTPUTS
% -- file_name: the name of the file(s) written by the function.
%               When 'format' is set to 'mat', output will be:
%                  'data_evo_simu_0mut_simplified_*Time*.mat'
%               When 'format' is set to 'csv', output will be:
%                   'data_evo_simu_0mut_simplified_*Time*_MeanFitness.mat'
%                   'data_evo_simu_0mut_simplified_*Time*_CellNumber.mat'
%                   'data_evo_simu_0mut_simplified_*Time*_EffectiveCellDepth.mat'
%                   'data_evo_simu_0mut_simplified_*Time*_Parameters.mat'
%% -------------------------------------------------------------------------
% Parse inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numvarargs = length(varargin);
if numvarargs > 4
    error('evolution_0mut_simplified:TooManyInputs', ...
        'requires at most 2 optional inputs');
end
vararg_update = {'format', 'csv', 'platform','MATLAB'};
pos1 = find(strcmpi(varargin,'format'))+1;
pos2 = find(strcmpi(varargin,'platform'))+1;
if numvarargs/2 == 2
    vararg_update([2,4]) = varargin([pos1,pos2]);
elseif numvarargs/2 == 1
    if ~isempty(pos1)
        vararg_update(2) = varargin(pos1);
    elseif ~isempty(pos2)
        vararg_update(4) = varargin(pos2);
    end
end
ouptput_format = vararg_update{2};
platform_choice = vararg_update{4};

noise_growth = noise_option;

cell_num_evo = zeros(lineage, t_evo+1);
cell_num_evo(:,1) = cell_num_ini;
x_mean = nan(1,t_evo+1);
x_mean(1) = cell_num_ini'*x_ini/sum(cell_num_ini);


% Simulate pooled growth %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Grow under two possible conditions:
%         1. Without growth noise
%         2. With growth noise


tstart0 = tic;
if noise_growth == 0
    for j = 2:(t_evo+1)
        if mod(j,20) == 1
            fprintf('Current generation: %i\n', j-1)
        end
        pos = cell_num_evo(:,j-1)~=0;
        n_prev = cell_num_evo(pos,j-1);
        x_ini_rela = max((1+x_ini(pos))/(1+x_mean(j-1)),0);
        cell_num_evo(pos,j) = round(n_prev.*x_ini_rela);
        x_mean(j) = cell_num_evo(:,j)'*x_ini/sum(cell_num_evo(:,j));
    end
elseif noise_growth == 1
    for j = 2:(t_evo+1)
        if mod(j,20) == 1
            fprintf('Current generation: %i\n', j-1)
        end
        pos = cell_num_evo(:,j-1)~=0;
        n_prev = cell_num_evo(pos,j-1);
        x_ini_rela = max((1+x_ini(pos))/(1+x_mean(j-1)),0);
        cell_num_evo(pos,j) = poissrnd(n_prev.*x_ini_rela);
        x_mean(j) = cell_num_evo(:,j)'*x_ini/sum(cell_num_evo(:,j));
    end
end

effective_cell_depth = sum(cell_num_evo);
telaps0=toc(tstart0);
fprintf('Computing time for %i generations: %f seconds.\n', t_evo, telaps0)
% End of the main function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Report Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = datestr(now,'yyyymmdd-HHMMSSFFF');
switch lower(ouptput_format)
    case {'mat'}
        file_name = ['data_evo_simu_0mut_simplified_' dt '.mat'];
        switch lower(platform_choice)
            case{'matlab'}
                save(file_name, 'effective_cell_depth','cell_num_evo','x_mean',...
                    'x_ini','cell_num_ini','lineage','t_evo','noise_option')
            case{'octave'}
                save('-mat7-binary', file_name, 'effective_cell_depth',...
                    'cell_num_evo','x_mean','x_ini','cell_num_ini',...
                    'lineage','t_evo','noise_option')
        end
        
    case {'csv'}
        file_name_3 = ['data_evo_simu_0mut_simplified_' dt '_EffectiveCellDepth.csv'];
        file_name_4 = ['data_evo_simu_0mut_simplified_' dt '_CellNumber.csv'];
        file_name_5 = ['data_evo_simu_0mut_simplified_' dt '_MeanFitness.csv'];
        file_name_6 = ['data_evo_simu_0mut_simplified_' dt '_Paramaters.csv'];
        
        csvwrite(file_name_3,effective_cell_depth)
        csvwrite(file_name_4,cell_num_evo)
        csvwrite(file_name_5,x_mean)
        output_parameters = cell(lineage+1,5);
        output_parameters{1,1} = 'Fitness of each genotype (x_i)';
        output_parameters{1,2} = ...
            'Initial cell number of each genotype (n0_i)';
        output_parameters{1,3} = 'Number of genotypes (L)';
        output_parameters{1,4} = 'Total number of generations grown (T)';
        output_parameters{1,5} = 'Optional noise (noise_option)';
        output_parameters(2:end,1) = num2cell(x_ini);
        output_parameters(2:end,2) = num2cell(cell_num_ini);
        output_parameters{2,3} = lineage;
        output_parameters{2,4} = t_evo;
        output_parameters{2,5} = [num2str(noise_option)];
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
        file_name = ['data_evo_simu_0mut_simplified_' dt '.csv'];
end
end
