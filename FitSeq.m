function [ file_name, x_estimate_result, r_estimate_result, x_mean_est ] = ...
    FitSeq(t_seq_vec, BC_num_mat_original, effective_cell_depth, deltat, varargin)
%--------------------------------------------------------------------------
% FitSeq
% ESTIMATE FITNESS OF EACH GENOTYPE IN A COMPETETIVE POOLED GROWTH EXPERIMENT
%
%
% INPUTS
% -- t_seq_vec: a vector of all sequencing time points
%
% -- BC_num_mat_original: a matrix of read number of each genotype at each
%                         sequencing time point,
%                         size = genotypes * length(t_seq_vec) 
%
% -- effective_cell_depth: a vector of the effective cell number (number of 
%                          cells transferred at the bottleneck) of the population 
%                          at each sequencing time point,
%                          size = 1 * length(t_seq_vec)
%
% --  deltat: number of generations between two succesive cell transfers
%             This is required in addition to t_seq_vec because not every 
%             cell transfer is necessarely sequenced (e.g. t_seq_vec = [0, 3, 6,
%             9, 15])
%
% -- 'format': optional, file format of the output file, 'csv'(default) or 'mat'
%
% -- 'kappa': optional, a noise parameter that characterizes the total noise 
%             introduced by growth, cell transfer, DNA extraction, PCR, and 
%             sequencing, default value is 2.5, from Levy et al. Nature 2015 519, 
%             181-186. To measure kappa empirically, see that reference. 
%
% -- 'opt_cycle': optional, the number of cycles used when using likelihood 
%                 optimization method to estimate fitness, default value is 2
%
% -- file_name: the name of the file(s) written by the function.
%               When 'format' is set to 'mat', output will be:
%                   Fit-Seq_result_*Time*.mat' 
%               When 'format' is set to 'csv', output will be:
%                   Fit-Seq_result_EstimatedFitness_*Time*.csv'
%                   Fit-Seq_result_EstimatedReads_*Time*.csv'
%
% OUTPUTS
% -- x_estimate_result: a vector of the estimated fitness of each genotype
%
% -- r_estimate_result: a matrix of the estimated read number of each genotype 
%                       at each sequencing time point,
%                       size = genotypes * length(t_seq_vec)
%
% -- x_mean_est: a vector of the estimated mean fitness of the population
%                at each sequencing time point,
%                size = 1 * length(t_seq_vec)
% 
%--------------------------------------------------------------------------
% Parse inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numvarargs = length(varargin);
if numvarargs > 6
    error('FitSeq:TooManyInputs', 'requires at most 3 optional inputs');
end
vararg_update = {'format','csv', 'kappa',2.5, 'opt_cycle',2};
pos1 = find(strcmpi(varargin,'format'))+1;
pos2 = find(strcmpi(varargin,'kappa'))+1;
pos3 = find(strcmpi(varargin,'opt_cycle'))+1;
if numvarargs/2 == 3
    vararg_update([2,4,6]) = varargin([pos1,pos2,pos3]);
elseif numvarargs/2 == 2
    if isempty(pos1)
        vararg_update([4,6]) = varargin([pos2,pos3]);
    elseif isempty(pos2)
        vararg_update([2,6]) = varargin([pos1,pos3]);
    elseif isempty(pos3)
        vararg_update([2,4]) = varargin([pos1,pos2]);
    end
elseif numvarargs/2 == 1
    if ~isempty(pos1)
        vararg_update(2) = varargin(pos1);
    elseif ~isempty(pos2)
        vararg_update(4) = varargin(pos2);
    elseif ~isempty(pos3)
        vararg_update(6) = varargin(pos3);
    end
end
ouptput_format = vararg_update{2}; 
kappa = vararg_update{4};
opt_cycle_num = vararg_update{6};

BC_num_mat = BC_num_mat_original;
lineage = size(BC_num_mat,1);

% Replace any 0 entries between two non-zero entries in a row with 1  
% This is done becasue 0s must be removed to avoid undefined terms in the 
% likelihood optimization function.   
BN1 = BC_num_mat'>0~=0;
[~,index_1] = sort(BN1);
last_nonzero_point = (index_1(end,:).*any(BN1))';
for i = 1:lineage
    pos_modify = BC_num_mat(i,1:last_nonzero_point(i))==0;
    BC_num_mat(i,pos_modify)=1;
end

linear_fit_points = 2; % Number of time points used in initial fitness estimate
%                        based on a log-linear regression

% Replace 0s in the 2nd column of BC_num_mat_original with 1. 
% This is done to allow for a log-linear regression on all time points
BC_num_mat(BC_num_mat(:,2)==0,2:linear_fit_points)=1;

read_depth_keep = sum(BC_num_mat_original);
cell_depth_keep =  effective_cell_depth;


% Start of the main function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Make an inital guess of the fitnesses by log-linear regression using
%         the first several sequencing time points, the number of time points 
%         used is defined by linear_fit_points
BC_freq_mat_tempt = BC_num_mat(:,1:linear_fit_points)./...
    repmat(read_depth_keep(1:linear_fit_points),lineage,1);
fit_param1 = nan(lineage,2);
for i = 1:lineage
    fit_param1(i,:) = polyfit(t_seq_vec(1:linear_fit_points),...
        log(BC_freq_mat_tempt(i,:)),1);
end
% Wrighten fitness
x_ini_est_linear = exp(fit_param1(:,1))-1;
x_ini_est_linear = (1+x_ini_est_linear)/(1+nansum(BC_num_mat(:,1).*...
    x_ini_est_linear)./sum(BC_num_mat(:,1)))-1;


% Step 2: Estimate all fitnesses by optimization of the likelihood function
kappa_vec = kappa*ones(1,length(t_seq_vec));
x_ini0 = x_ini_est_linear;
x_mean_est = nansum(BC_num_mat.*repmat(x_ini_est_linear,1,length(t_seq_vec)))./sum(BC_num_mat);

[ x_ini_est_opt_tempt1 ] = xini_opt_command(t_seq_vec, BC_num_mat, ...
    read_depth_keep, cell_depth_keep, x_ini0, x_mean_est, kappa_vec, deltat);
% see xini_opt_command.m

% Step 3: Use updated fitness values to re-estimate the mean fitness and all fitnesses
%         by optimization of the likelihood function. The number of times this step 
%         will repeat is equal to opt_cycle_num - 1
opt_cycle_num_count = 1;
x_ini_est_opt_tempt = (1+x_ini_est_opt_tempt1)/(1+nansum(BC_num_mat(:,1).*...
    x_ini_est_opt_tempt1)./sum(BC_num_mat(:,1)))-1;
pos = x_ini_est_opt_tempt<-1;
x_ini_est_opt_tempt(pos) = x_ini_est_linear(pos);
    
while opt_cycle_num_count < opt_cycle_num   
    x_ini0 = x_ini_est_opt_tempt;
    x_mean_est = nansum(BC_num_mat.*repmat(x_ini0,1,length(t_seq_vec)))./sum(BC_num_mat);
    [ x_ini_est_opt_tempt ] = xini_opt_command(t_seq_vec, BC_num_mat, ...
        read_depth_keep, cell_depth_keep, x_ini0, x_mean_est, kappa_vec, deltat);
        % see xini_opt_command.m
    x_ini_est_opt_tempt = (1+x_ini_est_opt_tempt)/(1+nansum(BC_num_mat(:,1).*...
        x_ini_est_opt_tempt)./sum(BC_num_mat(:,1)))-1;  
    pos = x_ini_est_opt_tempt<-1;
    x_ini_est_opt_tempt(pos) = x_ini_est_linear(pos);
    opt_cycle_num_count = opt_cycle_num_count+1;
end
x_estimate_result = x_ini_est_opt_tempt;
% End of the main function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Report Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_estimate_result = zeros(lineage,length(t_seq_vec));
for i = 1:lineage
    x0_l = x_estimate_result(i);
    r_exp_l = BC_num_mat(i,:);
    [r_estimate_result(i,:),~] = xini_est_m2_complex(t_seq_vec,x0_l,r_exp_l,...
        read_depth_keep,cell_depth_keep,x_mean_est,kappa_vec,deltat);
        %see xini_est_m2_complex.m
end

dt = datestr(now,'yyyymmdd-HHMMSSFFF');
switch lower(ouptput_format)
    case {'mat'}
        file_name = ['Fit-Seq_result_' dt '.mat'];
        save(file_name,'file_name','x_estimate_result','r_estimate_result')
    case {'csv'}
        file_name_1 = ['Fit-Seq_result_EstimatedFitness_' dt '.csv'];
        file_name_2 = ['Fit-Seq_result_EstimatedReads_' dt '.csv'];
        csvwrite(file_name_1,x_estimate_result)
        csvwrite(file_name_2,r_estimate_result)
        file_name = ['Fit-Seq_result_' dt '.csv'];
end
end
