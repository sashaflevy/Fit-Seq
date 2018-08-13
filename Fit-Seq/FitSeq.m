function [ x_estimate_result, r_estimate_result, x_mean_estimate_result ] = ...
    FitSeq(BC_num_mat_original, t_seq_vec, cell_depth, file_name, varargin)
%--------------------------------------------------------------------------
% FitSeq
% ESTIMATE FITNESS OF EACH GENOTYPE IN A COMPETETIVE POOLED GROWTH EXPERIMENT
%
%
% INPUTS
% -- BC_num_mat_original: a matrix of read number of each genotype at each
%                         sequencing time point,
%                         size = number of genotypes * length(t_seq_vec).
%                         This is a required input.
%
% -- t_seq_vec: a vector of all sequencing time points.
%               Use [] as input if it is unknown.
%               Either unempty ([]) t_seq_vec or unempty cell_depth is a required input.
%
% -- cell_depth: a matrix of the cell number across growth cycles. The first
%                row of this matrix is the cell number after bottleneck
%                (before growth) for each time point. The second row of this
%                matrix is the cell number before bottleneck (after growth).
%                t_seq_vec is calculated from this matrix, and
%                length(t_seq_vec) = size(cell_depth,2) + 1.
%                If t_seq_vec is also supplied by the user, the t_seq_vec
%                calculated from this matrix will be used instead.
%                Use [] as input if it is unknown.
%
% -- file_name: the name of the file(s) written by the function.
%               When 'format' is set to 'mat', output will be:
%                   file_name_Fit-Seq_result_*Time*.mat'
%               When 'format' is set to 'csv', output will be:
%                   file_name_Fit-Seq_result_EstimatedFitness_*Time*.csv'
%                   file_name_Fit-Seq_result_EstimatedReads_*Time*.csv'
%               Use [] as input if it is not given.
%
% -- 'format': optional, file format of the output file, 'csv'(default) or 'mat'
%
% -- 'kappa': optional, a noise parameter that characterizes the total noise
%             introduced by growth, cell transfer, DNA extraction, PCR, and
%             sequencing, default value is 2.5, from Levy et al. Nature 2015 519,
%             181-186. To measure kappa empirically, see that reference.
%
% -- 'opt_cycle_max': optional, maximum number of cycles used when using likelihood
%                 optimization method to estimate fitness, default value is
%                 10. Increse this cycle number might increase the fitness
%                 estimation accuracy, but extend the compute time.
%
% -- 'platform': optional, the platform that run Fit-Seq, 'MATLAB (default)' or 'OCTAVE'
%
%
% OUTPUTS
% -- x_estimate_result: a vector of the estimated fitness of each genotype
%
% -- r_estimate_result: a matrix of the estimated read number of each genotype
%                       at each sequencing time point,
%                       size = number of genotypes * length(t_seq_vec)
%
% -- x_mean_estimate_result: a vector of the estimated mean fitness of the
%                            population at each sequencing time point
%                            size = 1 * length(t_seq_vec)
%--------------------------------------------------------------------------
% Parse inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
numvarargs = length(varargin);
if numvarargs > 8
    error('FitSeq:TooManyInputs', 'requires at most 4 optional inputs');
end

vararg_update = {'format','csv', 'kappa',2.5, 'opt_cycle_max',10, 'platform','MATLAB'};
pos1 = find(strcmpi(varargin,'format'))+1;
pos2 = find(strcmpi(varargin,'kappa'))+1;
pos3 = find(strcmpi(varargin,'opt_cycle_max'))+1;
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
%%
ouptput_format = vararg_update{2};
kappa = vararg_update{4};
opt_cycle_num_max = vararg_update{6};
platform_choice = vararg_update{8};

if ~isempty(cell_depth)
    t_seq_vec_interval = max(log2(cell_depth(2,:)./cell_depth(1,:)),0);
    t_seq_vec = zeros(1,size(BC_num_mat_original,2));
    for i1 = 1:length(t_seq_vec_interval)
        t_seq_vec(i1+1) = sum(t_seq_vec_interval(1:i1));
    end
end


BC_num_mat = BC_num_mat_original;
lineage = size(BC_num_mat,1);

% Replace any 0 entries between two non-zero entries in a row with 1
% This is done becasue 0s must be removed to avoid undefined terms in the
% likelihood optimization function.
BN1 = BC_num_mat'>0;
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
x_ini_est_linear = exp(fit_param1(:,1))-1;
x_ini_est_linear = (1+x_ini_est_linear)/(1+nansum(BC_num_mat(:,1).*...
    x_ini_est_linear)./sum(BC_num_mat(:,1)))-1;

%%
switch lower(platform_choice)
    case{'matlab'}
        % Step 2: Estimate all fitnesses by optimization of the likelihood function
        kappa_vec = kappa*ones(1,length(t_seq_vec));
        
        r_estimate_mat_tempt = zeros(lineage,length(t_seq_vec));
        likelihood_vec_tempt = zeros(lineage,1);
        
        x_ini0 = x_ini_est_linear;
        x_mean_est = nansum(BC_num_mat.*repmat(x_ini_est_linear,1,length(t_seq_vec)))./sum(BC_num_mat);
        
        [ x_ini_est_opt_tempt ] = xini_opt_command(t_seq_vec, BC_num_mat, x_ini0, ...
            x_mean_est, kappa_vec, read_depth_keep);
        % see xini_opt_command.m
        
        for i = 1:lineage
            x0_l = x_ini_est_opt_tempt(i);
            r_exp_l = BC_num_mat(i,:);
            [r_estimate_mat_tempt(i,:),likelihood_vec_tempt(i)] = ...
                xini_est_Wrightian(t_seq_vec,r_exp_l,x0_l,...
                x_mean_est,kappa_vec,read_depth_keep);
            % see xini_est_Wrightian.m
        end
        
        pos = (((BC_num_mat(:,1)/sum(BC_num_mat(:,1)))./...
            (BC_num_mat(:,2)/sum(BC_num_mat(:,2))))>...
            2^(t_seq_vec(2)-t_seq_vec(1))) | x_ini_est_opt_tempt<-1;
        x_ini_est_opt_tempt(pos) = x_ini_est_linear(pos);
        likelihood_all_old = sum(likelihood_vec_tempt);
        
        likelihood_all_result_linear = sum(likelihood_vec_tempt);
        x_estimate_result_linear = x_ini_est_opt_tempt;
        x_mean_est_result_linear = x_mean_est;
        r_estimate_result_linear = r_estimate_mat_tempt;
        
        
        % Step 3: Use updated fitness values to re-estimate the mean fitness and all fitnesses
        %         by optimization of the likelihood function. The number of times this step
        %         will repeat is equal to opt_cycle_num_max - 1 or when total
        %         likelihood doesn't increasing
        opt_cycle_num_count = 1;
        likelihood_all_new = likelihood_all_old;
        while (likelihood_all_new  >= likelihood_all_old) && (opt_cycle_num_count < opt_cycle_num_max)
            
            r_estimate_result = r_estimate_mat_tempt;
            if opt_cycle_num_count == 1
                x_ini0 = x_ini_est_opt_tempt;
                x_mean_est = nansum(r_estimate_mat_tempt.*...
                    repmat(x_ini0,1,length(t_seq_vec)))./sum(r_estimate_mat_tempt);
                x_mean_est(1) = 0;
            else
                x_ini0 = x_ini_est_opt_tempt;
                x_mean_est_tempt = nansum(r_estimate_mat_tempt.*...
                    repmat(x_ini0,1,length(t_seq_vec)))./sum(r_estimate_mat_tempt);
                x_mean_est = mean(x_mean_est-x_mean_est_tempt)+x_mean_est_tempt;
            end
            
            [ x_ini_est_opt_tempt ] = xini_opt_command(t_seq_vec, BC_num_mat, x_ini0, ...
                x_mean_est, kappa_vec, read_depth_keep);
            % see xini_opt_command.m
            
            for i = 1:lineage
                x0_l = x_ini_est_opt_tempt(i);
                r_exp_l = BC_num_mat(i,:);
                [r_estimate_mat_tempt(i,:),likelihood_vec_tempt(i)] = ...
                    xini_est_Wrightian(t_seq_vec,r_exp_l,x0_l,...
                    x_mean_est,kappa_vec,read_depth_keep);
                % see xini_est_Wrightian.m
            end
            
            pos = (((BC_num_mat(:,1)/sum(BC_num_mat(:,1)))./...
                (BC_num_mat(:,2)/sum(BC_num_mat(:,2))))>...
                2^(t_seq_vec(2)-t_seq_vec(1))) | x_ini_est_opt_tempt<-1;
            x_ini_est_opt_tempt(pos) = x_ini_est_linear(pos);
            
            likelihood_all_old = likelihood_all_new;
            likelihood_all_new = sum(likelihood_vec_tempt);
            
            if  (likelihood_all_new < likelihood_all_old) && opt_cycle_num_count == 1
                likelihood_all_result = likelihood_all_result_linear;
                x_estimate_result = x_estimate_result_linear;
                x_mean_estimate_result = x_mean_est_result_linear;
                r_estimate_result = r_estimate_result_linear;
                break
            end
            
            if likelihood_all_new >= likelihood_all_old
                likelihood_all_result = likelihood_all_new;
                x_estimate_result = x_ini_est_opt_tempt;
                x_mean_estimate_result = x_mean_est;
                r_estimate_result = r_estimate_mat_tempt;
                opt_cycle_num_count = opt_cycle_num_count+1;
            end
            
        end
     
        
    case{'octave'}
        % Step 2: Estimate all fitnesses by optimization of the likelihood function
        kappa_vec = kappa*ones(1,length(t_seq_vec));
        
        r_estimate_mat_tempt = zeros(lineage,length(t_seq_vec));
        likelihood_vec_tempt = zeros(lineage,1);
        
        x_ini0 = x_ini_est_linear;
        x_mean_est = nansum(BC_num_mat.*repmat(x_ini_est_linear,1,length(t_seq_vec)))./sum(BC_num_mat);
        
        [ x_ini_est_opt_tempt ] = xini_opt_command_OCTAVE(t_seq_vec, BC_num_mat, x_ini0, ...
            x_mean_est, kappa_vec, read_depth_keep);
        % see xini_opt_command_OCTAVE.m
        
        for i = 1:lineage
            x0_l = x_ini_est_opt_tempt(i);
            r_exp_l = BC_num_mat(i,:);
            [r_estimate_mat_tempt(i,:),likelihood_vec_tempt(i)] = ...
                xini_est_Wrightian(t_seq_vec,r_exp_l,x0_l,...
                x_mean_est,kappa_vec,read_depth_keep);
            % see xini_est_Wrightian.m
        end
        
        pos = (((BC_num_mat(:,1)/sum(BC_num_mat(:,1)))./...
            (BC_num_mat(:,2)/sum(BC_num_mat(:,2))))>...
            2^(t_seq_vec(2)-t_seq_vec(1))) | x_ini_est_opt_tempt<-1;
        x_ini_est_opt_tempt(pos) = x_ini_est_linear(pos);
        likelihood_all_old = sum(likelihood_vec_tempt);
        
        likelihood_all_result_linear = sum(likelihood_vec_tempt);
        x_estimate_result_linear = x_ini_est_opt_tempt;
        x_mean_est_result_linear = x_mean_est;
        r_estimate_result_linear = r_estimate_mat_tempt;
        
        % Step 3: Use updated fitness values to re-estimate the mean fitness and all fitnesses
        %         by optimization of the likelihood function. The number of times this step
        %         will repeat is equal to opt_cycle_num_max - 1 or when total
        %         likelihood doesn't increasing
        opt_cycle_num_count = 1;
        likelihood_all_new = likelihood_all_old;
        while (likelihood_all_new  >= likelihood_all_old) && (opt_cycle_num_count < opt_cycle_num_max)
            
            r_estimate_result = r_estimate_mat_tempt;
            if opt_cycle_num_count == 1
                x_ini0 = x_ini_est_opt_tempt;
                x_mean_est = nansum(r_estimate_mat_tempt.*...
                    repmat(x_ini0,1,length(t_seq_vec)))./sum(r_estimate_mat_tempt);
                x_mean_est(1) = 0;
            else
                x_ini0 = x_ini_est_opt_tempt;
                x_mean_est_tempt = nansum(r_estimate_mat_tempt.*...
                    repmat(x_ini0,1,length(t_seq_vec)))./sum(r_estimate_mat_tempt);
                x_mean_est = mean(x_mean_est-x_mean_est_tempt)+x_mean_est_tempt;
            end
            
            [ x_ini_est_opt_tempt ] = xini_opt_command_OCTAVE(t_seq_vec, BC_num_mat, x_ini0, ...
                x_mean_est, kappa_vec, read_depth_keep);
            % see xini_opt_command_OCTAVE.m
            
            for i = 1:lineage
                x0_l = x_ini_est_opt_tempt(i);
                r_exp_l = BC_num_mat(i,:);
                [r_estimate_mat_tempt(i,:),likelihood_vec_tempt(i)] = ...
                    xini_est_Wrightian(t_seq_vec,r_exp_l,x0_l,...
                    x_mean_est,kappa_vec,read_depth_keep);
                % see xini_est_Wrightian.m
            end
            
            pos = (((BC_num_mat(:,1)/sum(BC_num_mat(:,1)))./...
                (BC_num_mat(:,2)/sum(BC_num_mat(:,2))))>...
                2^(t_seq_vec(2)-t_seq_vec(1))) | x_ini_est_opt_tempt<-1;
            x_ini_est_opt_tempt(pos) = x_ini_est_linear(pos);
            
            likelihood_all_old = likelihood_all_new;
            likelihood_all_new = sum(likelihood_vec_tempt);
            
            if  (likelihood_all_new < likelihood_all_old) && opt_cycle_num_count == 1
                likelihood_all_result = likelihood_all_result_linear;
                x_estimate_result = x_estimate_result_linear;
                x_mean_estimate_result = x_mean_est_result_linear;
                r_estimate_result = r_estimate_result_linear;
                break
            end
            
            if likelihood_all_new >= likelihood_all_old
                likelihood_all_result = likelihood_all_new;
                x_estimate_result = x_ini_est_opt_tempt;
                x_mean_estimate_result = x_mean_est;
                r_estimate_result = r_estimate_mat_tempt;
                opt_cycle_num_count = opt_cycle_num_count+1;
            end
        end 
end
% End of the main function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Report Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = datestr(now,'yyyymmdd-HHMMSSFFF');
switch lower(ouptput_format)
    case {'mat'}
        file_name_full = [file_name '_Fit-Seq_result_' dt '.mat'];
        switch lower(platform_choice)
            case{'matlab'}
                save(file_name_full,'x_estimate_result','r_estimate_result',...
                    'x_mean_estimate_result','likelihood_all_result','opt_cycle_num_count')
            case{'octave'}
                save('-mat7-binary', file_name_full,'x_estimate_result','r_estimate_result',...
                    'x_mean_estimate_result','likelihood_all_result','opt_cycle_num_count')
        end
        
    case {'csv'}
        file_name_full_1 = [file_name '_Fit-Seq_result_' dt '_EstimatedFitness.csv'];
        file_name_full_2 = [file_name '_Fit-Seq_result_' dt '_EstimatedReads.csv'];
        file_name_full_3 = [file_name '_Fit-Seq_result_' dt '_EstimatedMeanFitness.csv'];
        csvwrite(file_name_full_1, x_estimate_result)
        csvwrite(file_name_full_2, r_estimate_result)
        csvwrite(file_name_full_3, x_mean_estimate_result)
end
end
