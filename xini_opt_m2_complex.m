function [output] = xini_opt_m2_complex(input)
% -------------------------------------------------------------------------
% xini_opt_m2_complex
% A SUB-FUNCTION USED IN xini_opt_command.m TO CALCULATE LIKELIHOOD VALUE
% NEEDED IN OPTIMIZATION
%
% INPUTS
% -- input: fitness of a genotype
%
% OUTPUTS
% -- output: likelihood value of a genotype given its fitness
% -------------------------------------------------------------------------

% [CAN YOU ADD A FEW MORE COMMENTS THROUGHOUT TO EXPLAIN WHAT IS HAPPENING?]
global t_seq_vec_tempt read_depth_tempt cell_depth_tempt x_mean_est_tempt ...
    r_exp_l kappa_vec_tempt r_est_l_min
%%
x0_l = input;
vec_length = length(t_seq_vec_tempt);
r_est_l = zeros(1, vec_length);
% r_est_l: estimated read number of a genotype at each sequencing time point
r_est_l(1) = r_exp_l(1);

for j1 = 2:vec_length
    x_mean_est_interp = interp1([t_seq_vec_tempt(j1-1),t_seq_vec_tempt(j1)],...
        [x_mean_est_tempt(j1-1),x_mean_est_tempt(j1)], ...
        t_seq_vec_tempt(j1-1):(t_seq_vec_tempt(j1)-1));
    x_ini_rela = 1+x_mean_est_interp;
    r_est_l(j1) = max(r_exp_l(j1-1)*max(1+x0_l,0)...
        ^(t_seq_vec_tempt(j1)-t_seq_vec_tempt(j1-1))/prod(x_ini_rela)*...
        (cell_depth_tempt(j1-1)*read_depth_tempt(j1))/...
        (read_depth_tempt(j1-1)*cell_depth_tempt(j1)), r_est_l_min(j1));      
end

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
    

