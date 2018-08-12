function [output] = xini_opt_Wrightian(input)
% -------------------------------------------------------------------------
% xini_opt_Wrightian_s1
% A SUB-FUNCTION CALLED BY xini_opt_command_s1.m TO CALCULATE THE OPPOSITE 
% VALUE (-y) OF THE LIKELIHOOD FUNCTION FOR A GENOTYPE GIVEN THE FITNESS OF 
% THE GENOTYPE WITHOUT TOTAL CELL NUMBER AT EACH SEQUENCING TIME POINT FOR 
% WRIGHTIAN FITNESS
%
%
% INPUTS
% -- input: fitness of a genotype
%
%
% OUTPUTS
% -- output: opposite likelihood value (-y) of the genotype
%
% -------------------------------------------------------------------------
global t_seq_vec_tempt read_depth_tempt x_mean_est_term ...
    r_exp_l kappa_vec_tempt r_est_l_min

x0_l = input;
vec_length = length(t_seq_vec_tempt);
r_est_l = zeros(1, vec_length); % r_est_l: estimated read number of a genotype 
                                % at each sequencing time point
r_est_l(1) = r_exp_l(1);

% Estimate the relative fitness (relative to mean fitness) for a genotype,
% and calculate the estimated read number of the genotype at all sequencing
% time points
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

pos2 = r_est_l>0 & r_est_l<20 & r_exp_l>0 & r_exp_l<20;
likelihood_vec_log(pos2) = r_exp_l(pos2).*log(r_est_l(pos2))-r_est_l(pos2)...
    -log(factorial(r_exp_l(pos2)));

pos3 = r_est_l>0 & r_est_l<20 & r_exp_l>=20;
likelihood_vec_log(pos3) = r_exp_l(pos3).*log(r_est_l(pos3))-r_est_l(pos3)...
    -r_exp_l(pos3).*log(r_exp_l(pos3))+r_exp_l(pos3)-...
    1/2*log(2*pi*r_exp_l(pos3));

pos4 = r_est_l>0 & r_est_l<20 & r_exp_l == 0;
likelihood_vec_log(pos4) = -r_est_l(pos4);

output = -nansum(likelihood_vec_log(2:end));
end
