function [ x_ini_est_opt ] = xini_opt_command(t_seq_vec, BC_num_mat, ...
    x_ini0, x_mean_est, kappa_vec, read_depth)
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
%               sequencing, see Levy et al. Nature 2015 519, 81-186.
%               size = 1 * length(t_seq_vec)
%
% -- read_depth: a vector of the total read number of the population at each
%                sequencing time point,
%                size = 1 * length(t_seq_vec)
%
% -- fitness_type: 'Wrightian' or 'Malthusian'
%
% OUTPUTS
% -- x_ini_est_opt: a vector of the estimated fitness of each genotype,
%                   size = genotypes * 1
%
%--------------------------------------------------------------------------
global t_seq_vec_tempt read_depth_tempt x_mean_est_term ...
    r_exp_l kappa_vec_tempt r_est_l_min

t_seq_vec_tempt = t_seq_vec;
read_depth_tempt = read_depth;
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
        r_est_l_min(1) = r_exp_l(1);
        for j1 = 2:vec_length
            r_est_l_min(j1) = r_est_l_min(j1-1)/...
                2^(t_seq_vec_tempt(j1)-t_seq_vec_tempt(j1-1))*...
                read_depth_tempt(j1)/read_depth_tempt(j1-1);
        end
        
        t_seq_vec_left = floor(t_seq_vec);
        t_seq_vec_right = ceil(t_seq_vec);
        for j1 = 2:vec_length
            x_mean_est_part = interp1([t_seq_vec(j1-1),t_seq_vec(j1)],...
                [x_mean_est(j1-1),x_mean_est(j1)], t_seq_vec_right(j1-1):t_seq_vec_left(j1));
            x_mean_est_term(j1-1) = prod(1+x_mean_est_part(1:end-1))*...
                (1+x_mean_est_part(end))^(t_seq_vec(j1)-t_seq_vec_left(j1))*...
                (1+x_mean_est(j1-1))^(t_seq_vec_right(j1-1)-t_seq_vec(j1-1));
        end
        
        % Estimate fitness of a genotype by optimization of the likelihood
        % function, given an initial inferred value of fitness
        options = optimoptions(@fminunc, 'Algorithm', 'quasi-newton', ...
            'MaxFunEvals',1e8, 'MaxIter',1000, 'Display','off');
        [x_ini_est_opt(i)] = fminunc(@xini_opt_Wrightian, x_ini0(i), options);
        % see xini_opt_Wrightian.m
    end
end
end
