function [r_est_l,likelihood_value] = xini_est_Wrightian(t_seq_vec_tempt,...
    r_exp_l,x0_l,x_mean_est_tempt,kappa_vec_tempt,read_depth_tempt)
% -------------------------------------------------------------------------
% xini_est_Wrightian_s1
% A SUB-FUNCTION CALLED BY Fit-Seq.m TO CALCULATE THE ESTIMTED READ NUMBER 
% OF A GENOTYPE AT ALL SEQUENCING TIME POINTS AND THE LIKEHOOD VALUE OF THE
% GENOTYPE WITHOUT TOTAL CELL NUMBER AT EACH SEQUENCING TIME POINT FOR 
% WRIGHTIAN FITNESS
%
%
% INPUTS
% -- t_seq_vec_tempt: a vector of all sequencing time points
%
% -- r_exp_l: observed read number of a genotype at each sequencing time
%             point
%
% -- x0_l: fitness of a genotype
%
% -- x_mean_est_tempt: a vector of the mean fitness of the population at 
%                      each sequencing time point, 
%                      size = 1 * length(t_seq_vec_tempt)
%
% -- kappa_vec_tempt: a vector of the kappa value at each sequencing time point, 
%                     kappa is a noise parameter that characterizes the total 
%                     noise introduced by growth, cell transfer, DNA extraction, 
%                     PCR, and sequencing, see Levy et al. Nature 2015 519, 81-186.
%                     size = 1 * length(t_seq_vec_tempt)
%
% -- read_depth_tempt: a vector of the total read number of the population 
%                      at each sequencing time point, 
%                      size = 1 * length(t_seq_vec_tempt)
%
%
% OUTPUTS
% -- r_est_l: a vector of the estimated read number of a genotype at each 
%             sequencing time point, 
%             size = 1*length(t_seq_vec_tempt)
%
% -- likelihood_value: likelihood value of a genptype
%
% -------------------------------------------------------------------------
%%
vec_length = length(t_seq_vec_tempt);
r_est_l = zeros(1, vec_length); % r_est_l: estimated read number of a genotype 
                                % at each sequencing time point
r_est_l(1) = r_exp_l(1);

% Calculate the minimum read number for a genotype at each sequencing time point
% This minimum represents the number of reads that are expected for a genotype 
% that is not growing and being lost to dilution  
r_est_l_min = zeros(1, vec_length);
r_est_l_min(1) = r_exp_l(1);
for j1 = 2:vec_length
        r_est_l_min(j1) = r_est_l_min(j1-1)/...
            2^(t_seq_vec_tempt(j1)-t_seq_vec_tempt(j1-1))*...
            read_depth_tempt(j1)/read_depth_tempt(j1-1);
end

% Estimate the relative fitness (relative to mean fitness) for a genotype,
% and calculate the estimated read number of the genotype at all sequencing
% time points
t_seq_vec = t_seq_vec_tempt;
t_seq_vec_left = floor(t_seq_vec);
t_seq_vec_right = ceil(t_seq_vec);
x_mean_est = x_mean_est_tempt;
x_mean_est_term = nan(size(t_seq_vec));
for j1 = 2:vec_length
    x_mean_est_part = interp1([t_seq_vec(j1-1),t_seq_vec(j1)],...
        [x_mean_est(j1-1),x_mean_est(j1)], t_seq_vec_right(j1-1):t_seq_vec_left(j1));
    x_mean_est_term(j1-1) = prod(1+x_mean_est_part(1:end-1))*...
        (1+x_mean_est_part(end))^(t_seq_vec(j1)-t_seq_vec_left(j1))*...
        (1+x_mean_est(j1-1))^(t_seq_vec_right(j1-1)-t_seq_vec(j1-1));
end
for j1 = 2:vec_length  
    r_est_l(j1) = max(r_exp_l(j1-1)*max(1+x0_l,0)...
        ^(t_seq_vec_tempt(j1)-t_seq_vec_tempt(j1-1))/x_mean_est_term(j1-1)*...
        (read_depth_tempt(j1)/read_depth_tempt(j1-1)), r_est_l_min(j1));      
end

% Calculate the value of the likelihood function for a genotype using the 
% true and estimated read number data of the genotype
likelihood_vec_log = nan(1,vec_length);
pos1 = r_est_l>=20 & r_exp_l>0;
likelihood_vec_log(pos1) = 1/4*log(r_est_l(pos1))-...
    1/2*log(4*pi*kappa_vec_tempt(pos1)) - 3/4*log(r_exp_l(pos1))-...
    (sqrt(r_exp_l(pos1))-sqrt(r_est_l(pos1))).^2./kappa_vec_tempt(pos1);

pos2 = r_est_l>0 & r_est_l<20 & r_exp_l>0 & r_exp_l<10;
likelihood_vec_log(pos2) = r_exp_l(pos2).*log(r_est_l(pos2))-r_est_l(pos2)...
    -log(factorial(r_exp_l(pos2)));

pos3 = r_est_l>0 & r_est_l<20 & r_exp_l>=10;
likelihood_vec_log(pos3) = r_exp_l(pos3).*log(r_est_l(pos3))-r_est_l(pos3)...
    -r_exp_l(pos3).*log(r_exp_l(pos3))+r_exp_l(pos3)-...
    1/2*log(2*pi*r_exp_l(pos3));

pos4 = r_est_l>0 & r_est_l<20 & r_exp_l == 0;
likelihood_vec_log(pos4) = -r_est_l(pos4);

likelihood_value = nansum(likelihood_vec_log(2:end));
end
