function [ x_ini_est_opt ] = xini_opt_command(t_seq_vec, BC_num_mat, ...
    read_depth, cell_depth, x_ini0, x_mean_est, kappa_vec, deltat)
%--------------------------------------------------------------------------
% xini_opt_command
% A SUB-FUNCTION CALLED BY Fit-Seq.m TO ESTIMATE FITNESSES OF ALL GENOTYPES
% USING LIKELIHOOD OPTIMIZATION METHOD
%
%
% INPUTS
% -- t_seq_vec: a vector of all sequencing time points
%
% -- BC_num_mat: a matrix of the read number of each genotype at each 
%                sequencing time point,
%                size = genotypes * length(t_seq_vec)
%
% -- read_depth: a vector of the total read number of the population at each 
%                sequencing time point, 
%                size = 1 * length(t_seq_vec)
%
% -- cell_depth: a vector of the effective cell number of the population 
%                at each sequencing time point, 
%                size = 1 * length(t_seq_vec)
%
% -- x_ini0: a vector of initial value of fitness of each genotype used in the
%            optimization,
%            size = genotypes * 1 
%
% -- x_mean_est: a vector of the mean fitness of the population at each sequencing 
%                time point, 
%                size = 1 * length(t_seq_vec)
%
% -- kappa_vec: a vector of the kappa value at each sequencing time point, 
%               kappa is a noise parameter that characterizes the total noise 
%               introduced by growth, cell transfer, DNA extraction, PCR, and 
%               sequencing, 
%               size = 1 * length(t_seq_vec)
%
% -- deltat: number of generations between two succesive cell transfers
%            (This is not redundant with t_seq_vec when t_seq_vec contains 
%            arbitrary sampling times)
%
%
% OUTPUTS
% -- x_ini_est_opt: a vector of the estimated fitness of each genotype, 
%                   size = genotypes * 1
%
%--------------------------------------------------------------------------
global t_seq_vec_tempt read_depth_tempt cell_depth_tempt x_mean_est_tempt ...
    r_exp_l kappa_vec_tempt r_est_l_min

t_seq_vec_tempt = t_seq_vec;
read_depth_tempt = read_depth;
cell_depth_tempt = cell_depth;
x_mean_est_tempt = x_mean_est;
kappa_vec_tempt = kappa_vec;
lineage = size(BC_num_mat,1);
x_ini_est_opt = nan(lineage,1);

vec_length = length(t_seq_vec_tempt);
r_est_l_min = zeros(1, vec_length);

for i = 1:lineage
    if mod(i,5e3) == 0
        fprintf('The %i-th lineage of total %i lineages.\n', i, lineage)
    end
    if ~isnan(x_ini0(i))
        % Calculate the minimum read number for a genotype at each sequencing time point
        % This minimum represents the number of reads that are expected for a genotype 
        % that is not growing and being lost to dilution  
        r_exp_l = BC_num_mat(i,:);  % r_exp_l: true read number of a genotype across all sequencing time points      
        if t_seq_vec_tempt(1) == 0
            r_est_l_min(1) = r_exp_l(1)/2^deltat;
            r_est_l_min(2) = r_est_l_min(1)/2^(t_seq_vec_tempt(2)-t_seq_vec_tempt(1)-deltat);
            for j1 = 3:vec_length
                r_est_l_min(j1) = r_est_l_min(j1-1)/...
                    2^(t_seq_vec_tempt(j1)-t_seq_vec_tempt(j1-1)-t_seq_vec_tempt(1))*...
                    (read_depth_tempt(j1)/read_depth_tempt(j1-1));
            end
        elseif t_seq_vec_tempt(1) ~= 0
            r_est_l_min(1) = r_exp_l(1);
            for j1 = 2:vec_length
                r_est_l_min(j1) = r_est_l_min(j1-1)/...
                    2^(t_seq_vec_tempt(j1)-t_seq_vec_tempt(j1-1))*...
                    read_depth_tempt(j1)/read_depth_tempt(j1-1);
            end
        end
        
        % Estimate fitness of a genotype by optimization of the likelihood 
        % function, given an initial inferred value of fitness
        options = optimoptions(@fminunc, 'Algorithm', 'quasi-newton', ...
            'MaxFunEvals',1e8, 'MaxIter',1000, 'Display','off');
        [x_ini_est_opt(i)] = fminunc(@xini_opt_m2_complex, x_ini0(i), options);    
    end
end
end
