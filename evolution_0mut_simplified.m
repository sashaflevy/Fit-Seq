function [ file_name ] = evolution_0mut_simplified(lineage, t_evo, ...
    cell_num_ini, x_ini, noise_option, varargin)
% -------------------------------------------------------------------------
% evolution_0mut_simplified
% SIMPLIFIED VERSION OF SIMULATED COMPETETIVE POOLED GROWTH OF A POPULATION OF 
% GENOTYPES WITH DIFFERENT FITNESSES. [DEFINE WHAT MAKES THIS SIMPLIFIED] 
%
% INPUTS
% -- lineage: number of genotypes of the population
% -- t_evo: total number of growth generations
% -- cell_num_ini: vector of the initial cell number of each genotype at the 0-th
%                  generation
% -- x_ini: vector of the fitness of each genotype
% -- noise_option: options of whether cell growth noise is simulated, 
%                  logical (0-1) scaler value
% -- 'format': optional, file format of outputs, 'csv'(default) or 'mat'
%
% OUTPUTS
% -- file_name: name of the file generated [WHAT IS OUTPUT TO THIS FILE?]
% -------------------------------------------------------------------------

% Parse inputs
numvarargs = length(varargin);
if numvarargs > 2
    error('evolution_0mut_complex:TooManyInputs', ...
        'requires at most 1 optional inputs');
end
optargs = {'format', 'csv'};
optargs(1:numvarargs) = varargin;
ouptput_format = optargs{2};

noise_growth = noise_option;

cell_num_evo = zeros(lineage, t_evo+1);
cell_num_evo(:,1) = cell_num_ini;
x_mean = nan(1,t_evo+1);
x_mean(1) = cell_num_ini'*x_ini/sum(cell_num_ini);


% Simulate Pooled growth [CAN YOU GIVE A FEW MORE DESCRIPTIONS THROUGHOUT?]
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

dt = datestr(now,'yyyymmdd-HHMMSSFFF');
switch lower(ouptput_format)
    case {'mat'}
        file_name = ['data_evo_simu_0mut_simplified_' dt '.mat'];
        save(file_name, 'effective_cell_depth','cell_num_evo','x_mean',...
            'x_ini','cell_num_ini','lineage','t_evo','noise_option')       
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
        output_parameters{2,5} = ['[' num2str(noise_option) ']'];
        fileID = fopen(file_name_6,'wt');
        for k1 = 1:4
            fprintf(fileID,'%s,', output_parameters{1,k1});
        end
        fprintf(fileID,'%s\n', output_parameters{1,5});
        for k2 = 1:4
            fprintf(fileID,'%f,', output_parameters{2,k2});
        end
        fprintf(fileID,'%s\n', output_parameters{2,5});
        for k3 = 2:lineage
            fprintf(fileID,'%f,', output_parameters{k3+1,1});
            fprintf(fileID,'%f\n', output_parameters{k3+1,2});
        end
        fclose(fileID);
        file_name = ['data_evo_simu_0mut_simplified_' dt];
end
return