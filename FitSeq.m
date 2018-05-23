function [ file_name, x_estimate_result, r_estimate_result, x_mean_est ] = ...
    FitSeq(t_seq_vec, BC_num_mat_original, effective_cell_depth, deltat, varargin)
% -------------------------------------------------------------------------
% FitSeq
% ESTIMATE FITNESS OF EACH GENOTYPE IN A COMPETETIVE POOLED GROWTH EXPERIMENT
%
% INPUTS
% -- t_seq_vec: a vector of all sequencing time points
% -- BC_num_mat_original: a matrix of read number of each genotype at each
%                         sequencing time point
%                         size = genotypes * length(t_seq_vec) 
% -- effective_cell_depth: a vector of the effective cell number of population 
%                          at each sequenced time point,
%                          size = length(t_seq_vec)
% --  deltat: number of generations between two succesive cell transfers [ISN'T THIS REDUNDANT WITH t_seq_vec?]
% -- 'format': optional, file format of outputs, 'csv'(default) or 'mat'
% -- 'kappa': optional, [DEFINE WHAT THIS IS], default value is 2.5
% -- 'opt_cycle': optional, [DEFINE WHAT THIS IS], default value is 2
%
% OUTPUTS
% -- file_name: name of the file generated
% -- x_estimate_result: vector of the estimated fitness of each genotype
% -- r_estimate_result: matrix of the estimated read number of each genotype at each
%                       sequencing time point
%                       size = genotypes * length(t_seq_vec)
% -------------------------------------------------------------------------
%%
% Parse Inputs
numvarargs = length(varargin);
if numvarargs > 6
    error('FitSeq:TooManyInputs', 'requires at most 3 optional inputs');
end

optargs = {'format','csv', 'kappa',2.5, 'opt_cycle',2};
if numvarargs==6 && strcmpi(varargin{1},'format') && ...
        strcmpi(varargin{2},'kappa') && strcmpi(varargin{3},'opt_cycle')
    optargs ([2,4,6]) = varargin([2,4,6]);
elseif numvarargs==6 && strcmpi(varargin{1},'format') && ...
        strcmpi(varargin{3},'kappa') && strcmpi(varargin{2},'opt_cycle')
    optargs ([2,4,6]) = varargin([2,6,4]);
elseif numvarargs==6 && strcmpi(varargin{2},'format') && ...
        strcmpi(varargin{1},'kappa') && strcmpi(varargin{3},'opt_cycle')
    optargs ([2,4,6]) = varargin([4,2,6]);
elseif numvarargs==6 && strcmpi(varargin{2},'format') && ...
        strcmpi(varargin{3},'kappa') && strcmpi(varargin{1},'opt_cycle')
    optargs ([2,4,6]) = varargin([4,6,2]);
elseif numvarargs==6 && strcmpi(varargin{3},'format') && ...
        strcmpi(varargin{1},'kappa') && strcmpi(varargin{2},'opt_cycle')
    optargs ([2,4,6]) = varargin([6,2,4]);
elseif numvarargs==6 && strcmpi(varargin{3},'format') && ...
        strcmpi(varargin{2},'kappa') && strcmpi(varargin{1},'opt_cycle')
    optargs ([2,4,6]) = varargin([6,4,2]);
elseif numvarargs==4 && strcmpi(varargin{1},'format') && ...
        strcmpi(varargin{2},'kappa')
    optargs ([2,4]) = varargin([2,4]);
elseif numvarargs==4 && strcmpi(varargin{2},'format') && ...
        strcmpi(varargin{1},'kappa')
    optargs ([2,4]) = varargin([4,2]);
elseif numvarargs==4 && strcmpi(varargin{1},'format') && ...
        strcmpi(varargin{2},'opt_cycle')
    optargs ([2,6]) = varargin([2,4]);
elseif numvarargs==4 && strcmpi(varargin{2},'format') && ...
        strcmpi(varargin{1},'opt_cycle')
    optargs ([2,6]) = varargin([4,2]);
elseif numvarargs==4 && strcmpi(varargin{1},'kappa') && ...
        strcmpi(varargin{2},'opt_cycle')
    optargs ([4,6]) = varargin([2,4]);
elseif numvarargs==4 && strcmpi(varargin{2},'kappa') && ...
        strcmpi(varargin{1},'opt_cycle')
    optargs ([4,6]) = varargin([4,2]);
elseif numvarargs==2 && strcmpi(varargin{1},'format')
    optargs (2) = varargin(2);
elseif numvarargs==4 && strcmpi(varargin{1},'kappa')
    optargs (4) = varargin(2);
elseif numvarargs==4 && strcmpi(varargin{1},'opt_cycle')
    optargs (6) = varargin(2);
end

ouptput_format = optargs{2}; 
kappa = optargs{4};
opt_cycle_num = optargs{6};
%%
BC_num_mat = BC_num_mat_original;
lineage = size(BC_num_mat,1);

% Replace zero entries between two non-zero entries in BC_num_mat_original as 1
% [CAN YOU WRITE A JUSTIFICATION FOR WHY YOU DO THIS]
BN1 = BC_num_mat'>0~=0;
[~,index_1] = sort(BN1);
last_nonzero_point = (index_1(end,:).*any(BN1))';
for i = 1:lineage
    pos_modify = BC_num_mat(i,1:last_nonzero_point(i))==0;
    BC_num_mat(i,pos_modify)=1;
end

linear_fit_points = 2; % Number of time points used in initial fitness estimate
%                        based on a log-linear regression
% [WOULDN'T IT BE BETTER TO MOVE THIS TO A varargin THAT CAN BE DEFINED BY THE USER?]

% Replace zero entries in the 2nd column of BC_num_mat_original as 1
% This is done to allow for a log-linear regression on all time points
BC_num_mat(BC_num_mat(:,2)==0,2:linear_fit_points)=1;


read_depth_keep = sum(BC_num_mat_original);
cell_depth_keep =  effective_cell_depth;
%%
% Start of the main function
% Step 1: Make an inital guess of the fitness by log-linear regression using
%         the first several sequencing time points, the number of time points 
%         used is defined by linear_fit_points

BC_freq_mat_tempt = BC_num_mat(:,1:linear_fit_points)./...
    repmat(read_depth_keep(1:linear_fit_points),lineage,1);
fit_param1 = nan(lineage,2);
for i = 1:lineage
    fit_param1(i,:) = polyfit(t_seq_vec(1:linear_fit_points),...
        log(BC_freq_mat_tempt(i,:)),1);
end
% Malthusian fitness
% x_ini_est_linear = fit_param1(:,1);
% x_ini_est_linear = x_ini_est_linear - ...
%    nansum(BC_num_mat(:,1).*x_ini_est_linear)./sum(BC_num_mat(:,1));   

% Wrighten fitness
x_ini_est_linear = exp(fit_param1(:,1))-1;
x_ini_est_linear = (1+x_ini_est_linear)/...
    (1+nansum(BC_num_mat(:,1).*x_ini_est_linear)./sum(BC_num_mat(:,1)))-1;


%%
% Step 2: Estimate all fitnesses by optimization of the likelihood function
kappa_vec = kappa*ones(1,length(t_seq_vec));
x_ini0 = x_ini_est_linear;
x_mean_est = nansum(BC_num_mat.*repmat(x_ini_est_linear,1,...
    length(t_seq_vec)))./sum(BC_num_mat);
%%
[ x_ini_est_opt_tempt1 ] = xini_opt_command(t_seq_vec, BC_num_mat, ...
    read_depth_keep, cell_depth_keep, x_ini0, x_mean_est, kappa_vec, deltat);

%%
% Step 3: Use updated fitness values to re-estimate all fitnesses
% by optimization of the likelihood function.
% The number of times step will repeat is equal to opt_cycle_num - 1
opt_cycle_num_count = 1;
% x_ini_est_opt_tempt = x_ini_est_opt_tempt1 - ...
%         nansum(BC_num_mat(:,1).*x_ini_est_opt_tempt1)./sum(BC_num_mat(:,1));
x_ini_est_opt_tempt = (1+x_ini_est_opt_tempt1)/(1+...
        nansum(BC_num_mat(:,1).*x_ini_est_opt_tempt1)./sum(BC_num_mat(:,1)))-1;
pos = x_ini_est_opt_tempt<-1;
x_ini_est_opt_tempt(pos) = x_ini_est_linear(pos);
    
while opt_cycle_num_count < opt_cycle_num   
    x_ini0 = x_ini_est_opt_tempt;
    x_mean_est = nansum(BC_num_mat.*repmat(x_ini0,1,length(t_seq_vec)))./...
        sum(BC_num_mat);
    [ x_ini_est_opt_tempt ] = xini_opt_command(t_seq_vec, BC_num_mat, ...
        read_depth_keep, cell_depth_keep, x_ini0, x_mean_est, kappa_vec, deltat);
%     x_ini_est_opt_tempt = x_ini_est_opt_tempt - ...
%         nansum(BC_num_mat(:,1).*x_ini_est_opt_tempt)./sum(BC_num_mat(:,1));
    x_ini_est_opt_tempt = (1+x_ini_est_opt_tempt)/(1+...
        nansum(BC_num_mat(:,1).*x_ini_est_opt_tempt)./sum(BC_num_mat(:,1)))-1;
    pos = x_ini_est_opt_tempt<-1;
    x_ini_est_opt_tempt(pos) = x_ini_est_linear(pos);
    opt_cycle_num_count = opt_cycle_num_count+1;
end

x_estimate_result = x_ini_est_opt_tempt;
% End of the main function

%%
% Report Results
r_estimate_result = zeros(lineage,length(t_seq_vec));
for i = 1:lineage
    x0_l = x_estimate_result(i);
    r_exp_l = BC_num_mat(i,:);
    [r_estimate_result(i,:),~] = xini_est_m2_complex(t_seq_vec,x0_l,r_exp_l,...
        read_depth_keep,cell_depth_keep,x_mean_est,kappa_vec,deltat);
end

dt = datestr(now,'yyyymmdd-HHMMSSFFF');
switch lower(ouptput_format)
    case {'mat'}
        file_name = ['Fit-Seq_result_' dt '.mat'];
        save(file_name,'file_name','x_estimate_result','r_estimate_result')
    case {'csv'}
        file_name_1 = ['Fit-Seq_result_EstimatedFitness' dt '.csv'];
        file_name_2 = ['Fit-Seq_result_EstimatedReads' dt '.csv'];
        csvwrite(file_name_1,x_estimate_result)
        csvwrite(file_name_2,r_estimate_result)
        file_name = ['Fit-Seq_result_' dt];
end



