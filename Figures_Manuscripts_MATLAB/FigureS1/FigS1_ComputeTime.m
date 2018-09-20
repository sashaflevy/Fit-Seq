%%   Simulated data for Figure S1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lineage_list = [1e2,1e3,1e4,1e5,1e6];
deltat = 10;
t_evo = 50;
read_depth_average = 100;
noise_option = [1,1,1,1,1];

file_name_list1 = cell(length(lineage_list),1);
for i1 = 1:length(lineage_list)
    lineage = lineage_list(i1);
    cell_num_ini = 1e2*ones(lineage,1);
    x_ini = random('Normal',0,0.1, [lineage,1]);
    
    [ file_name_list1{i1} ] = evolution_0mut_complex(lineage,t_evo,...
        cell_num_ini,x_ini,deltat,read_depth_average,noise_option,'format','mat');
end
dt = datestr(now,'yyyymmdd-HHMMSSFFF');
file_name = ['file_name_list1_Fig_ComputeTime_' dt '.mat'];
save(file_name,'file_name_list1')



load('file_name_list1_Fig_ComputeTime_20180805-222626794.mat')
compute_time_vec = zeros(size(file_name_list1));
for i2 = 1:size(file_name_list1,1)
    tic
    load(file_name_list1{i2})
    BC_num_mat_original = num_mat_sequencing(:,1:4);
    t_seq_vec = t_seq_vec(1:4);
    cell_depth = [];
    [ ~, ~, ~] = FitSeq(BC_num_mat_original, t_seq_vec, ...
        cell_depth, [], 'format','mat','opt_cycle_max',2);

    compute_time_vec(i2) = toc;
end

dt = datestr(now,'yyyymmdd-HHMMSSFFF');
file_name = ['compute_time_Fig_ComputeTime_' dt '.mat'];
save(file_name,'compute_time_vec')




%%   Figure S1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

line_color1 = [239,0,0]/255;
line_color2 = [0,184,229]/255;
line_color3 = [120,120,120]/255;
line_color4 = [180,180,180]/255;
line_color5 = [0,187,87]/255;

lineage_list = [1e2,1e3,1e4,1e5,1e6];
load('compute_time_Fig_ComputeTime_20180806-012535595.mat')

figure('rend','painters','pos',[1 11 260 220])
ax1 = axes('Position',[0.2 0.3 0.7 0.6]);
ax1.ActivePositionProperty = 'outerposition';
plot(lineage_list, compute_time_vec, '.-', 'MarkerFaceColor',line_color3, ...
    'LineWidth',3, 'MarkerSize',20, 'Color',line_color3);
box on
pbaspect([1 1 1])
set(gca,'xscale','log')
set(gca,'yscale','log')
xlim([1e2,1e6])
ylim([0.6,6e4])
xticks(lineage_list)
xticklabels({'10^2','10^3','10^4','10^5','10^6'})
yticks([0.6,6,60,600,6000,60000])
yticklabels({'10^{-2}','10^{-1}','10^0','10^1','10^2','10^3'})
set(gca,'FontSize',10.8,'FontName','Helvetica')

ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1)
text(0.55,0.11,{'Number of genotypes' ; 'of the population, L'},...
    'FontSize',12,'FontName','Helvetica','HorizontalAlignment','center')
text(0.125,0.6,'Compute time (minutes)','Rotation',...
    90,'FontSize',12,'FontName','Helvetica','HorizontalAlignment','center');



