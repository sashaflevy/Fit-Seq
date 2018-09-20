%% Load and organize original real data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename1 = 'ReadData.csv';
read_data_cell = textscan(fopen(filename1),'%f %f %s %s %f %f %f %f %f %f %f %f %f',...
    'Delimiter',',','headerLines', 1);

BC_num_mat_original_all = cell2mat(read_data_cell(1,6:10));

x_jamie = cell2mat(read_data_cell(1,11));

t_seq_vec = 0:3:12;
BC_num_mat_original_r1 = round(BC_num_mat_original_all(1:2500,:));
BC_num_mat_original_r2 = round(BC_num_mat_original_all(2501:5000,:));
BC_num_mat_original_r3 = round(BC_num_mat_original_all(5001:7500,:));

pos1 = sum(BC_num_mat_original_r1>0,2)<=2 | ...
    sum(BC_num_mat_original_r2>0,2)<=2 | ...
    sum(BC_num_mat_original_r3>0,2)<=2 | ...
    sum(BC_num_mat_original_r1(:,1:2)==0,2)==2 | ...
    sum(BC_num_mat_original_r2(:,1:2)==0,2)==2 | ...
    sum(BC_num_mat_original_r3(:,1:2)==0,2)==2;

pos2 = ((BC_num_mat_original_r1(:,1)/sum(BC_num_mat_original_r1(:,1)))./...
    (BC_num_mat_original_r1(:,2)/sum(BC_num_mat_original_r1(:,2))))>...
    2^(t_seq_vec(2)-t_seq_vec(1));

pos3 = ((BC_num_mat_original_r2(:,1)/sum(BC_num_mat_original_r2(:,1)))./...
    (BC_num_mat_original_r2(:,2)/sum(BC_num_mat_original_r2(:,2))))>...
    2^(t_seq_vec(2)-t_seq_vec(1));

pos4 = ((BC_num_mat_original_r3(:,1)/sum(BC_num_mat_original_r3(:,1)))./...
    (BC_num_mat_original_r3(:,2)/sum(BC_num_mat_original_r3(:,2))))>...
    2^(t_seq_vec(2)-t_seq_vec(1));

pos = ~(pos1 | pos2 | pos3 | pos4);
% pos = ~pos3;

BC_num_mat_r1 = [BC_num_mat_original_r1(pos,:); sum(BC_num_mat_original_r1(~pos,:))];
BC_num_mat_r2 = [BC_num_mat_original_r2(pos,:); sum(BC_num_mat_original_r2(~pos,:))];
BC_num_mat_r3 = [BC_num_mat_original_r3(pos,:); sum(BC_num_mat_original_r3(~pos,:))];



%%   Fitness estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cell_depth = [];
[x_estimate_result_r1, ~, ~] = ...
    FitSeq(BC_num_mat_r1, t_seq_vec, cell_depth, [],'format','mat','kappa',1);

[x_estimate_result_r2, ~, ~] = ...
    FitSeq(BC_num_mat_r2, t_seq_vec, cell_depth, [],'format','mat','kappa',1);

[x_estimate_result_r3, ~, ~] = ...
    FitSeq(BC_num_mat_r3, t_seq_vec, cell_depth, [],'format','mat','kappa',1);



%%   Figure S2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('Fit-Seq_result_20180805-221423798.mat')
x_estimate_result_r1 = x_estimate_result;

load('Fit-Seq_result_20180805-221452537.mat')
x_estimate_result_r2 = x_estimate_result;

load('Fit-Seq_result_20180805-221520190.mat')
x_estimate_result_r3 = x_estimate_result;

x_jamie_r1 = x_jamie(1:2500);
x_jamie_r2 = x_jamie(2501:5000);
x_jamie_r3 = x_jamie(5001:7500);

pos1 = sum(BC_num_mat_original_r1>0,2)<=2 | ...
    sum(BC_num_mat_original_r2>0,2)<=2 | ...
    sum(BC_num_mat_original_r3>0,2)<=2 | ...
    sum(BC_num_mat_original_r1(:,1:2)==0,2)==2 | ...
    sum(BC_num_mat_original_r2(:,1:2)==0,2)==2 | ...
    sum(BC_num_mat_original_r3(:,1:2)==0,2)==2;

pos2 = ((BC_num_mat_original_r1(:,1)/sum(BC_num_mat_original_r1(:,1)))./...
    (BC_num_mat_original_r1(:,2)/sum(BC_num_mat_original_r1(:,2))))>...
    2^(t_seq_vec(2)-t_seq_vec(1));

pos3 = ((BC_num_mat_original_r2(:,1)/sum(BC_num_mat_original_r2(:,1)))./...
    (BC_num_mat_original_r2(:,2)/sum(BC_num_mat_original_r2(:,2))))>...
    2^(t_seq_vec(2)-t_seq_vec(1));

pos4 = ((BC_num_mat_original_r3(:,1)/sum(BC_num_mat_original_r3(:,1)))./...
    (BC_num_mat_original_r3(:,2)/sum(BC_num_mat_original_r3(:,2))))>...
    2^(t_seq_vec(2)-t_seq_vec(1));

pos = ~(pos1 | pos2 | pos3 | pos4);
lineage = sum(pos)+1;

BC_num_mat = BC_num_mat_r1;
BC_num_mat_tempt = BC_num_mat;
BC_num_mat_tempt(BC_num_mat_tempt(:,1) == 0,1) = 1;
BC_num_mat = BC_num_mat_tempt;
read_depth_keep = sum(BC_num_mat);
BC_freq_mat_tempt = BC_num_mat./repmat(read_depth_keep,lineage,1);
fit_param1 = nan(lineage,2);
for i = 1:lineage
    pos = BC_freq_mat_tempt(i,:)~=0;
    if sum(pos)>=2
        fit_param1(i,:) = polyfit(t_seq_vec(pos),...
            log(BC_freq_mat_tempt(i,pos)),1);
    end
end
x_est_linear_tempt_r1 = exp(fit_param1(:,1))-1;
x_est_linear_r1 = (1+x_est_linear_tempt_r1)/...
    (1+nansum(BC_num_mat(:,1).*x_est_linear_tempt_r1)./...
    sum(BC_num_mat(:,1)))-1;
    

BC_num_mat = BC_num_mat_r2;
BC_num_mat_tempt = BC_num_mat;
BC_num_mat_tempt(BC_num_mat_tempt(:,1) == 0,1) = 1;
BC_num_mat = BC_num_mat_tempt;
read_depth_keep = sum(BC_num_mat);
BC_freq_mat_tempt = BC_num_mat./repmat(read_depth_keep,lineage,1);
fit_param1 = nan(lineage,2);
for i = 1:lineage
    pos = BC_freq_mat_tempt(i,:)~=0;
    if sum(pos)>=2
        fit_param1(i,:) = polyfit(t_seq_vec(pos),...
            log(BC_freq_mat_tempt(i,pos)),1);
    end
end
x_est_linear_tempt_r2 = exp(fit_param1(:,1))-1;
x_est_linear_r2 = (1+x_est_linear_tempt_r2)/...
    (1+nansum(BC_num_mat(:,1).*x_est_linear_tempt_r2)./...
    sum(BC_num_mat(:,1)))-1;


BC_num_mat = BC_num_mat_r3;
BC_num_mat_tempt = BC_num_mat;
BC_num_mat_tempt(BC_num_mat_tempt(:,1) == 0,1) = 1;
BC_num_mat = BC_num_mat_tempt;
read_depth_keep = sum(BC_num_mat);
BC_freq_mat_tempt = BC_num_mat./repmat(read_depth_keep,lineage,1);
fit_param1 = nan(lineage,2);
for i = 1:lineage
    pos = BC_freq_mat_tempt(i,:)~=0;
    if sum(pos)>=2
        fit_param1(i,:) = polyfit(t_seq_vec(pos),...
            log(BC_freq_mat_tempt(i,pos)),1);
    end
end
x_est_linear_tempt_r3 = exp(fit_param1(:,1))-1;
x_est_linear_r3 = (1+x_est_linear_tempt_r3)/...
    (1+nansum(BC_num_mat(:,1).*x_est_linear_tempt_r3)./...
    sum(BC_num_mat(:,1)))-1;






line_color1 = [251,131,45]/255;
line_color2 = [77,154,74]/255;
line_color3 = [120,120,120]/255;
line_color4 = [180,180,180]/255;
line_color5 = [0,187,87]/255;

figure('rend','painters','pos',[1 11 500 400]);
ax1 = axes('Position',[0.2 0.2 0.8 0.8]);
ax1.ActivePositionProperty = 'outerposition';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subaxis(2,2,1, 'MarginBottom',0.17,'MarginTop',0.08,'MarginLeft',...
        0.15,'MarginRight',0.1,'SpacingHoriz',0.1,'SpacingVert',0.18);

sz = 15*ones(size(x_estimate_result_r2));    
plot([-10,10],[-10,10],':','LineWidth',3,'Color',line_color4)
hold on
scatter(x_est_linear_r2(1:end-1), x_est_linear_r3(1:end-1), ...
    sz(1:end-1), 'MarkerEdgeColor','none', 'MarkerFaceColor',line_color1, ...
    'MarkerFaceAlpha',0.1, 'MarkerEdgeAlpha',0, 'LineWidth',.1);
hold off
box on
pbaspect([1 1 1]) 
xlim([-1,2])
ylim([-1,2])
set(gca,'FontSize',10.8,'FontName','Helvetica')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subaxis(2,2,2, 'MarginBottom',0.17,'MarginTop',0.08,'MarginLeft',...
        0.15,'MarginRight',0.1,'SpacingHoriz',0.1,'SpacingVert',0.18);

sz = 15*ones(size(x_estimate_result_r2));    
plot([-10,10],[-10,10],':','LineWidth',3,'Color',line_color4)
hold on
scatter(x_estimate_result_r2(1:end-1), x_estimate_result_r3(1:end-1), ...
    sz(1:end-1), 'MarkerEdgeColor','none', 'MarkerFaceColor',line_color1, ...
    'MarkerFaceAlpha',0.1, 'MarkerEdgeAlpha',0, 'LineWidth',.1);
hold off
box on
pbaspect([1 1 1]) 
xlim([-1,2])
ylim([-1,2])
set(gca,'FontSize',10.8,'FontName','Helvetica')


subaxis(2,2,3, 'MarginBottom',0.17,'MarginTop',0.08,'MarginLeft',...
        0.15,'MarginRight',0.1,'SpacingHoriz',0.1,'SpacingVert',0.18);

sz = 15*ones(size(x_estimate_result_r2));    
plot([-10,10],[-10,10],':','LineWidth',3,'Color',line_color4)
hold on
scatter(x_estimate_result_r2(1:end-1), x_est_linear_r2(1:end-1), ...
    sz(1:end-1), 'MarkerEdgeColor','none', 'MarkerFaceColor',line_color2, ...
    'MarkerFaceAlpha',0.1, 'MarkerEdgeAlpha',0, 'LineWidth',.1);
hold off
box on
pbaspect([1 1 1]) 
xlim([-1,2])
ylim([-1,2])
set(gca,'FontSize',10.8,'FontName','Helvetica')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subaxis(2,2,4, 'MarginBottom',0.17,'MarginTop',0.08,'MarginLeft',...
        0.15,'MarginRight',0.1,'SpacingHoriz',0.1,'SpacingVert',0.18);

sz = 15*ones(size(x_estimate_result_r2));    
plot([-10,10],[-10,10],':','LineWidth',3,'Color',line_color4)
hold on
scatter(x_estimate_result_r3(1:end-1), x_est_linear_r3(1:end-1), ...
    sz(1:end-1), 'MarkerEdgeColor','none', 'MarkerFaceColor',line_color2, ...
    'MarkerFaceAlpha',0.1, 'MarkerEdgeAlpha',0, 'LineWidth',.1);
hold off
box on
pbaspect([1 1 1]) 
xlim([-1,2])
ylim([-1,2])
set(gca,'FontSize',10.8,'FontName','Helvetica')


ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1) % sets ax1 to current axes

text(0.198-0.13,0.915+0.05,'a','FontSize',15,'FontName','Helvetica','FontWeight','bold');
text(0.622-0.13,0.915+0.05,'b','FontSize',15,'FontName','Helvetica','FontWeight','bold')
text(0.198-0.13,0.448+0.05,'c','FontSize',15,'FontName','Helvetica','FontWeight','bold')
text(0.622-0.13,0.448+0.05,'d','FontSize',15,'FontName','Helvetica','FontWeight','bold')

text(0.295+0.02,0.605-0.055,{'Estimated fitness' '(log-linear regression, R1)'},...
    'FontSize',12,'FontName','Helvetica','HorizontalAlignment','center')
text(0.72+0.02,0.605-0.055,{'Estimated fitness' '(Fit-Seq, R1)'},...
    'FontSize',12,'FontName','Helvetica','HorizontalAlignment','center')


text(0.295+0.02,0.14-0.055,{'Estimated fitness' '(Fit-Seq, R1)'},...
    'FontSize',12,'FontName','Helvetica','HorizontalAlignment','center')
text(0.72+0.02,0.14-0.055,{'Estimated fitness' '(Fit-Seq, R2)'},...
    'FontSize',12,'FontName','Helvetica','HorizontalAlignment','center')

text(0.175-0.055,0.775,{'Estimated fitness' '(log-linear regression, R2)'},...
    'Rotation',90,'FontSize',12,'FontName','Helvetica','HorizontalAlignment','center');
text(0.5985-0.055,0.775,{'Estimated fitness' '(Fit-Seq, R2)'},...
    'Rotation',90,'FontSize',12,'FontName','Helvetica','HorizontalAlignment','center');

text(0.175-0.055,0.31,{'Estimated fitness' '(log-linear regression, R1)'},...
    'Rotation',90,'FontSize',12,'FontName','Helvetica','HorizontalAlignment','center');
text(0.5985-0.055,0.31,{'Estimated fitness' '(log-linear regression, R2)'},...
    'Rotation',90,'FontSize',12,'FontName','Helvetica','HorizontalAlignment','center');

