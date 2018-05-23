function [r_est_l,likelihood_value] = xini_est_m2_complex(t_seq_vec_tempt,...
    x0_l,r_exp_l,read_depth_tempt,cell_depth_tempt,x_mean_est_tempt,...
    kappa_vec_tempt,deltat_tempt)
% -------------------------------------------------------------------------
% xini_est_m2_complex
% A SUB-FUNCTION USED IN Fit-Seq.m TO ESTIMTE READ NUMBER OF A GENOTYPE AT
% ALL SEQUENCING TIME POINTS GIVEN FITNESS
%
% INPUTS
% -- t_seq_vec: a vector of all sequencing time points
% -- x0_l: fitness of a genptype
% -- r_exp_l: observed read number of a genotype at each sequencing time
%             point
% -- read_depth: vector of total read number of the population at each sequencing
%                time points, 1 * length(t_seq_vec) 
% -- cell_depth: vector of the total effective cell number of the population at 
%                each sequencing time point, 1 * length(t_seq_vec) 
% -- x_mean_est: vector of the mean fitness of the population at each sequencing 
%                time point, 1 * length(t_seq_vec) 
%
% OUTPUTS
% -- r_est_l: vector of the estimated read number of a genotype at each sequencing 
%             time point, 1*length(t_seq_vec_tempt) 
% -------------------------------------------------------------------------
%%
% [CAN YOU ADD A FEW MORE COMMENTS THROUGHOUT TO EXPLAIN WHAT IS HAPPENING?]
vec_length = length(t_seq_vec_tempt);
r_est_l = zeros(1, vec_length);
% r_est_l: estimated read number of a genotype at each sequencing time point
r_est_l(1) = r_exp_l(1);

r_est_l_min = zeros(1, vec_length);
if t_seq_vec_tempt(1) == 0
    r_est_l_min(1) = r_exp_l(1)/2^deltat_tempt;
    r_est_l_min(2) = r_est_l_min(1)/2^(t_seq_vec_tempt(2)-t_seq_vec_tempt(1)-deltat_tempt);  
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

for j1 = 2:vec_length
    x_mean_est_interp = interp1([t_seq_vec_tempt(j1-1),t_seq_vec_tempt(j1)],...
        [x_mean_est_tempt(j1-1),x_mean_est_tempt(j1)], ...
        t_seq_vec_tempt(j1-1):(t_seq_vec_tempt(j1)-1));
    x_ini_rela = 1+x_mean_est_interp;
    r_est_l(j1) = max(r_exp_l(j1-1)*max(1+x0_l,0)...
        ^(t_seq_vec_tempt(j1)-t_seq_vec_tempt(j1-1))*prod(x_ini_rela)*...
        (cell_depth_tempt(j1-1)*read_depth_tempt(j1))/...
        (read_depth_tempt(j1-1)*cell_depth_tempt(j1)),r_est_l_min(j1));
end

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
