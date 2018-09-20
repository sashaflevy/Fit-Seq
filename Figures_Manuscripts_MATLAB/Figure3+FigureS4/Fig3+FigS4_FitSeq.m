%%   Simulated data for Figure 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lineage = 5e3;
cell_num_ini = 100*ones(lineage,1);
sigma_vec = 0.01:0.01:0.4;
D_vec = 100:100:1000;
T_vec = [96:-12:12,6];
M_vec = [2,3,4,5,7];
noise_option = [1,1,1,1,1];

file_name_list1 = cell(length(sigma_vec),length(D_vec),length(T_vec),length(M_vec));

for i1 = 1:length(sigma_vec)
    x_ini = random('Normal',0,sigma_vec(i1), [lineage,1]);
    for i2 = 1:length(D_vec)
        for i3 = 1:length(T_vec)
            for i4 = 1:length(M_vec)
                if T_vec(i3)~=6 || M_vec(i4)~=5
                    t_evo = T_vec(i3);
                    M = M_vec(i4);
                    deltat = t_evo/(M-1);
                    read_depth_average = round(D_vec(i2)/M);
                    [ file_name_list1{i1,i2,i3,i4} ] = evolution_0mut_complex(...
                        lineage, t_evo, cell_num_ini, x_ini, deltat, ...
                        read_depth_average, noise_option, 'format','mat');       
                end
            end
        end
    end
end
dt = datestr(now,'yyyymmdd-HHMMSSFFF');
file_name1 = ['file_name_list1_' dt '.mat'];
save(file_name1,'file_name_list1')



load('file_name_list1_20180806-092114836.mat')
lineage = 5e3;
sigma_vec = 0.01:0.01:0.4;
D_vec = 100:100:1000;
T_vec = [96:-12:12,6];
M_vec = [2,3,4,5,7];

x_estimate_result_mat = cell(length(sigma_vec),length(D_vec),length(T_vec),length(M_vec));
x_estimate_result_mat(:) = {nan(lineage,1)};

for i1 = 1:length(sigma_vec)
    for i2 = 1:length(D_vec)
        for i3 = 1:length(T_vec)
            for i4 = 1:length(M_vec)
                if (T_vec(i3)~=6 || M_vec(i4)~=5)
                    load(file_name_list1{i1,i2,i3,i4})
                    M = M_vec(i4);
                    BC_num_mat_original = num_mat_sequencing(:,1:M);
                    t_seq_vec = t_seq_vec(1:M);
                    cell_depth = [];
                    [x_estimate_result_mat{i1,i2,i3,i4}, ~, ~] = ...
                        FitSeq(BC_num_mat_original, t_seq_vec, cell_depth, [], 'format','mat');
                end
            end
        end
    end
end

dt = datestr(now,'yyyymmdd-HHMMSSFFF');
file_name3 = ['x_estimate_result_mat_' dt '.mat'];
save(file_name3,'x_estimate_result_mat')




lineage = 5e3;
sigma_vec = 0.01:0.01:0.4;
D_vec = 100:100:1000;
T_vec = [96:-12:12,6];
M_vec = [2,3,4,5,7];

load('file_name_list1_20180806-092114836.mat')
load('x_estimate_result_mat_20180808-121852533.mat')

x_estimate_error_mat = nan(length(sigma_vec),length(D_vec),length(T_vec),length(M_vec));

sliding_window_size1 = 500;
sliding_window_size2 = 10;
tempt1 = length(1:sliding_window_size2:(lineage-sliding_window_size1+1));

x_estimate_error_sliding_window_cell = cell(length(sigma_vec),...
    length(D_vec),length(T_vec),length(M_vec));
x_estimate_error_sliding_window_cell(:) = {nan(tempt1,1)};
x_ini_sort_sliding_window_cell = x_estimate_error_sliding_window_cell;

for i1 = 1:length(sigma_vec)
    for i2 = 1:length(D_vec)
        for i3 = 1:length(T_vec)
            for i4 = 1:length(M_vec)
                if T_vec(i3)~=6 || M_vec(i4)~=5
                    load(file_name_list1{i1,i2,i3,i4})
                    x_ini_tempt = x_ini;
                    x_estimate_result = x_estimate_result_mat{i1,i2,i3,i4};
                    x_estimate_error_mat(i1,i2,i3,i4) = ...
                        nanmean(abs(x_estimate_result-x_ini_tempt));
                    
                    [x_ini_sort, x_ini_sort_pos] = sort(x_ini_tempt);
                    x_estimate_result_sort = x_estimate_result(x_ini_sort_pos);
                    for k1 = 1:tempt1
                        tempt2 = ((k1-1)*sliding_window_size2+1):...
                            ((k1-1)*sliding_window_size2+sliding_window_size1);
                        x_estimate_error_sliding_window_cell{i1,i2,i3,i4}(k1) = ...
                            nanmean(abs(x_estimate_result_sort(tempt2)...
                            - x_ini_sort(tempt2)));
                        x_ini_sort_sliding_window_cell{i1,i2,i3,i4}(k1) = ...
                            nanmean(x_ini_sort(tempt2));
                    end
                end
            end
        end
    end
end

dt = datestr(now,'yyyymmdd-HHMMSSFFF');
file_name1 = ['x_estimate_error_mat_' dt '.mat'];
save(file_name1,'x_estimate_error_mat')
file_name2 = ['x_estimate_error_sliding_window_cell_' dt '.mat'];
save(file_name2,'x_estimate_error_sliding_window_cell')
file_name3 = ['x_ini_sort_sliding_window_cell_' dt '.mat'];
save(file_name3,'x_ini_sort_sliding_window_cell')



%%   Figure 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('x_estimate_result_mat_20180808-121852533.mat')
load('x_estimate_error_mat_20180808-123618242.mat')

lineage = 5e3;
sigma_vec = 0.01:0.01:0.4;
D_vec = 100:100:1000;
T_vec = [96:-12:12,6];
M_vec = [2,3,4,5,7];

x_estimate_error_sigmaD_optdeltat = nan(length(sigma_vec),length(D_vec));
x_estimate_error_sigmaD_optM = nan(length(sigma_vec),length(D_vec));
x_estimate_error_sigmaD_optT = nan(length(sigma_vec),length(D_vec));
x_estimate_error_sigmaD_opterror = nan(length(sigma_vec),length(D_vec));
x_estimate_error_2d = nan(length(T_vec),length(M_vec));

x_estimate_error_sigmaD_optT_fixedM = nan(length(sigma_vec),length(D_vec));
x_estimate_error_sigmaD_opterror_fixedM = nan(length(sigma_vec),length(D_vec));
x_estimate_error_1d = nan(length(T_vec),1);

for i1 =  1:length(sigma_vec)
    for i2 = 1:length(D_vec)        
        x_estimate_error_2d(:) = x_estimate_error_mat(i1,i2,:,:);
        x_estimate_error_2d_min = min(min(x_estimate_error_2d));
        [r1,c1] = find(x_estimate_error_2d==x_estimate_error_2d_min);
        x_estimate_error_sigmaD_opterror(i1,i2) = x_estimate_error_2d_min;
        x_estimate_error_sigmaD_optM(i1,i2) = M_vec(c1);
        x_estimate_error_sigmaD_optT(i1,i2) = T_vec(r1);
        x_estimate_error_sigmaD_optdeltat(i1,i2) = T_vec(r1)/(M_vec(c1)-1);
        
        x_estimate_error_1d(:) = x_estimate_error_mat(i1,i2,:,3);
        x_estimate_error_1d_min = min(x_estimate_error_1d);
        r1 = find(x_estimate_error_1d==x_estimate_error_1d_min);
        x_estimate_error_sigmaD_opterror_fixedM(i1,i2) = x_estimate_error_1d_min;
        x_estimate_error_sigmaD_optT_fixedM(i1,i2) = T_vec(r1); 
    end
end




load('x_estimate_result_mat_20180808-121852533.mat')
load('x_estimate_error_sliding_window_cell_20180808-123618242.mat')
load('x_ini_sort_sliding_window_cell_20180808-123618242.mat')

line_color1 = [239,0,0]/255;
line_color2 = [0,184,229]/255;
line_color3 = [120,120,120]/255;
line_color4 = [180,180,180]/255;

figure('rend','painters','pos',[1 11 850 410])
ax1 = axes('Position',[0.1 0.4 0.16 0.8]);
ax1.ActivePositionProperty = 'outerposition';

i1 = 20; i2 = 5;
x_estimate_error_2d = nan(length(T_vec),length(M_vec));
x_estimate_error_2d(:) = x_estimate_error_mat(i1,i2,:,:);
imagesc(flip(x_estimate_error_2d',2), ...
    'AlphaData',~isnan(flip(x_estimate_error_2d',2)))
axis('xy')
hold on
p1 = 9-find(x_estimate_error_sigmaD_optT(i1,i2)==T_vec)+1;
p2 = find(x_estimate_error_sigmaD_optM(i1,i2)==M_vec);
plot(p1, p2, '.','MarkerSize',18,'Color',line_color1)
hold off
pbaspect([1 1 1])
xticks([1,2:2:9])
xticklabels({'6','12','36','60','84'})
yticks(1:5)
yticklabels({'2','3','4','5','7'})
set(gca,'FontSize',10.8,'FontName','Helvetica')
colormap(brewermap(64*8,'*Spectral'))
cb1 = colorbar('Ticks',0.03:0.04:0.15,'TickLabels',...
    {'0.03','0.07','0.11','0.15'});
set(cb1,'FontSize',10.8,'FontName','Helvetica');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax1 = axes('Position',[0.4 0.4 0.5 0.8]);
ax1.ActivePositionProperty = 'outerposition';

imagesc(sigma_vec,D_vec,x_estimate_error_sigmaD_opterror_fixedM');
sigma_step_size = sigma_vec(2)-sigma_vec(1);
for i1 = 1:length(sigma_vec)
    for i2 = 1:length(D_vec)
        text(sigma_vec(i1),D_vec(i2),...
            num2str(x_estimate_error_sigmaD_optT_fixedM(i1,i2),'%d'),...
            'FontSize',7.5,'FontName','Helvetica','HorizontalAlignment', ...
            'center')
    end
end
axis('xy')
pbaspect([3 1 1])
xticks(0.05:0.05:0.4)
xticklabels({'0.05', '0.1', '0.15', '0.2','0.25','0.3','0.35','0.4'})
yticks(D_vec(2:2:end))
yticklabels({'200','400','600','800','1000'})
set(gca,'FontSize',10.8,'FontName','Helvetica')
colormap(brewermap(64*8,'*Spectral'))
cb1 = colorbar('Ticks',0.005:0.02:0.085,'TickLabels',...
    {'0.005','0.025','0.045','0.065','0.085'});
set(cb1,'FontSize',10.8,'FontName','Helvetica');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% line_color_list = [9,129,84;190,156,46;251,131,45;1,99,146;199,46,41]/255;
line_color_list = [9,129,84;178,171,210;251,131,45;120,120,120;199,46,41]/255;
x_estimate_error_2d = nan(length(T_vec),length(M_vec));
i1_list = [5,15,30]; % sigma
i2 = 6; % D_vec
for k1 = 1:3
    i1 = i1_list(k1);
    ax1 = axes('Position',[0.0848+(k1-1)*0.2 0.15 0.15 0.3]);
    ax1.ActivePositionProperty = 'outerposition';
    hold on
    if k1 == 1
        tempt1 = 0.006;
        plot([0.005,0.005],[0,tempt1],':','LineWidth',3,'Color',line_color4)
        tempt3 = 0.018;
        text(-0.05,tempt3*0.08,'2','FontSize',10.8,'FontName','Helvetica',...
            'FontWeight','bold','Color',line_color3*0.6,'HorizontalAlignment','center');
        text(0.05,tempt3*0.08,'3','FontSize',10.8,'FontName','Helvetica',...
            'FontWeight','bold','Color',line_color3*0.6)
        
    elseif k1 == 2
        tempt1 = 0.024; tempt2 = 0.012;
        plot([-0.14,-0.14],[0,tempt1],':','LineWidth',3,'Color',line_color4)
        plot([0.12,0.12],[0,tempt2],':','LineWidth',3,'Color',line_color4)
        tempt3 = 0.045;
        text(-0.24,tempt3*0.08,'1','FontSize',10.8,'FontName','Helvetica',...
            'FontWeight','bold','Color',line_color3*0.6);
        text(0,tempt3*0.08,'2','FontSize',10.8,'FontName','Helvetica',...
            'FontWeight','bold','Color',line_color3*0.6)
        text(0.24,tempt3*0.08,'3','FontSize',10.8,'FontName','Helvetica',...
            'FontWeight','bold','Color',line_color3*0.6)
        
    elseif k1 == 3
        tempt1 = 0.045; tempt2 = 0.03;
        plot([-0.3,-0.3],[0, tempt1],':','LineWidth',3,'Color',line_color4)
        plot([0.12,0.12],[0,tempt2],':','LineWidth',3,'Color',line_color4)
        tempt3 = 0.135;
        text(-0.48,tempt3*0.08,'1','FontSize',10.8,'FontName','Helvetica',...
            'FontWeight','bold','Color',line_color3*0.6);
        text(-0.1,tempt3*0.08,'2','FontSize',10.8,'FontName','Helvetica',...
            'FontWeight','bold','Color',line_color3*0.6)
        text(0.36,tempt3*0.08,'3','FontSize',10.8,'FontName','Helvetica',...
            'FontWeight','bold','Color',line_color3*0.6)
    end
    
    x_estimate_error_2d(:) = x_estimate_error_mat(i1,i2,:,:);
    x_estimate_error_2d_min_vec = min(x_estimate_error_2d);
    
    i4 = 1;   % M = 2
    % r1_pos_2 = find(x_estimate_error_2d(:,i4)==x_estimate_error_2d_min_vec(i4));
    if k1 == 1
        r1_pos_2 = 7;
    elseif k1 == 2
        r1_pos_2 = 8;
    elseif k1 == 3
        r1_pos_2 = 9;
    end
    
    lg_str_2 = ['M = ' num2str(M_vec(i4)) ', T = ' num2str(T_vec(r1_pos_2))];
    
    h(1) = plot(x_ini_sort_sliding_window_cell{i1,i2,r1_pos_2,i4}, ...
        x_estimate_error_sliding_window_cell{i1,i2,r1_pos_2,i4},...
        '-','LineWidth',2,'Color',line_color_list(1,:));
    
    i4 = 3;   % M = 4
    r1_pos_1 = find(x_estimate_error_2d(:,i4)==x_estimate_error_2d_min_vec(i4));
    lg_str_1 = ['M = ' num2str(M_vec(i4)) ', T = ' num2str(T_vec(r1_pos_1))];
    h(2) = plot(x_ini_sort_sliding_window_cell{i1,i2,r1_pos_1,i4}, ...
        x_estimate_error_sliding_window_cell{i1,i2,r1_pos_1,i4},...
        '-','LineWidth',2,'Color',line_color_list(3,:));
    hold off
    box on
    pbaspect([1 1 1])
    tempt1 = sigma_vec(i1_list(k1))*2;
    xlim([-tempt1,tempt1])
    ylim([0,tempt3])
    xticks(-0.8*tempt1:tempt1*0.8:0.8*tempt1)
    yticks(0:(tempt3/3):tempt3)
    if k1 == 1
        xticklabels({'-0.08','0','0.08'})
        yticklabels({'0.000','0.006','0.012','0.018'})
    elseif k1 == 2
        xticklabels({'-0.24','0','0.24'})
        yticklabels({'0.000','0.015','0.030','0.045'})
    elseif k1 == 3
        xticklabels({'-0.48','0','0.48'})
        yticklabels({'0.000','0.045','0.090','0.135'})
    end
    set(gca,'FontSize',10.8,'FontName','Helvetica')
    
    text(-0.5*tempt1,tempt3-(tempt3)/9,lg_str_1,...
        'FontSize',10.8,'FontName','Helvetica','Color',line_color_list(3,:));
    text(-0.5*tempt1,tempt3-(tempt3)/4,lg_str_2,...
        'FontSize',10.8,'FontName','Helvetica','Color',line_color_list(1,:));
    
    text(0,1.12*tempt3,['\sigma = ', num2str(sigma_vec(i1_list(k1)))],...
        'FontSize',10.8,'FontName','Helvetica','HorizontalAlignment','center')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax1 = axes('Position',[0.73 0.155 0.21 0.32]);
ax1.ActivePositionProperty = 'outerposition';

i2 = 6;
i3 = find(T_vec==24);
hold on
h(5) = plot(sigma_vec,x_estimate_error_mat(:,i2,i3,5), '-', ...
    'LineWidth',3, 'Color',line_color_list(5,:));
h(4) = plot(sigma_vec,x_estimate_error_mat(:,i2,i3,4), '-', ...
    'LineWidth',3, 'Color',line_color_list(4,:));
h(3) = plot(sigma_vec,x_estimate_error_mat(:,i2,i3,3), '-', ...
    'LineWidth',3, 'Color',line_color_list(3,:));
h(2) = plot(sigma_vec,x_estimate_error_mat(:,i2,i3,2), '-', ...
    'LineWidth',3, 'Color',line_color_list(2,:));
h(1) = plot(sigma_vec,x_estimate_error_mat(:,i2,i3,1), '-', ...
    'LineWidth',3, 'Color',line_color_list(1,:));
hold off
box on
pbaspect([1.25 1 1])
xlim([sigma_vec(1)-0.5*(sigma_vec(2)-sigma_vec(1)), ...
    sigma_vec(end)+0.5*(sigma_vec(2)-sigma_vec(1))])
ylim([0,0.36])
set(gca,'FontSize',10.8,'FontName','Helvetica')
lg = legend(h,'M = 2','M = 3','M = 4','M = 5','M = 7');
set(lg, 'FontSize',10.8, 'FontName','Helvetica', ...
    'Position',[0.76   0.31    0.05   0.12])


ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1) % sets ax1 to current axes
text(0.16,0.582,{'Total number of' 'generations grown, T'},...
    'FontSize',12,'FontName','Helvetica','HorizontalAlignment','center')
text(0.05,0.798,{'Number of' 'time points, M'},'Rotation',90,'FontSize',...
    12,'FontName','Helvetica','HorizontalAlignment','center');
text(0.284,0.698,'Fit-Seq error, \epsilon','Rotation',90,'FontSize',12,...
    'FontName','Helvetica');

text(0.51,0.57,'Standard deviation of fitness, \sigma','FontSize',12,...
    'FontName','Helvetica')
text(0.342,0.62,'Total reads per genotype, D','Rotation',90,'FontSize',12,...
    'FontName','Helvetica');
text(0.936,0.798,{'Fit-Seq error, \epsilon'},...
    'Rotation',90,'FontSize',12,'FontName','Helvetica',...
    'HorizontalAlignment','center');

text(0.34,0.084,'Fitness','FontSize',12,'FontName','Helvetica')
text(0.032,0.21,'Fit-Seq error, \epsilon','Rotation',90,'FontSize',12,...
    'FontName','Helvetica');

text(0.74,0.084,'Standard deviation of fitness, \sigma',...
    'FontSize',12,'FontName','Helvetica')
text(0.692,0.31,'Fit-Seq error, \epsilon',...
    'Rotation',90,'FontSize',12,'FontName','Helvetica',...
    'HorizontalAlignment','center');


text(0.094-0.08,0.944+0.01,'a','FontSize',15,'FontName','Helvetica','FontWeight','bold')
text(0.392-0.08,0.97+0.01,'b','FontSize',15,'FontName','Helvetica','FontWeight','bold')
text(0.094-0.08,0.458+0.01,'c','FontSize',15,'FontName','Helvetica','FontWeight','bold')
text(0.738-0.08+0.01,0.488+0.01,'d','FontSize',15,'FontName','Helvetica','FontWeight','bold')




%%  Figure S4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_estimate_error_2d = nan(length(T_vec),length(M_vec));
bottom = min(min(min(min(x_estimate_error_mat(:,:,:,1:4)))));
top = max(max(max(max(x_estimate_error_mat(:,:,:,1:4))))); 
for k1 = 1 % 1:4
    figure('rend','painters','pos',[1 11 800 630])
    for i1 = (10*k1-9):(k1*10)
        for i2 = 1:length(D_vec)
            x_estimate_error_2d(:) = x_estimate_error_mat(i1,i2,:,:);
            subaxis(length(D_vec),10,(length(D_vec)-i2)*10+i1-(k1-1)*10, ...
                'MarginBottom',0.125,'MarginTop',0.05,'MarginLeft',0.15,...
                'MarginRight',0.2,'SpacingHoriz',0.005,'SpacingVert',0.005); 
            imagesc(flip(x_estimate_error_2d',2),'AlphaData',~isnan(flip(x_estimate_error_2d',2)))
            axis('xy')
            hold on
            p1 = 9-find(x_estimate_error_sigmaD_optT(i1,i2)==T_vec)+1;
            p2 = find(x_estimate_error_sigmaD_optM(i1,i2)==M_vec);
            plot(p1,p2, '.','MarkerSize',10,'Color',line_color1)
            hold off
            pbaspect([1 1 1])
            if i2 == 1
                xticks([1,2:2:9])
                xticklabels({'6', '12', '36', '60','84'})
                text(5,-2.3,['\sigma = ' num2str(sigma_vec(i1))],...
                    'FontSize',10.8,'FontName','Helvetica','HorizontalAlignment','center');
            else
                set(gca,'xtick',[])
            end
            if i1-(k1-1)*10 == 1
                yticks(1:length(M_vec))
                yticklabels({'2','3','4','5','7'});
                text(-4.5,2.5,['D = ' num2str(D_vec(i2))],'Rotation',90,...
                    'FontSize',10.8,'FontName','Helvetica','HorizontalAlignment','center');
            else
                set(gca,'ytick',[])
            end
            set(gca,'FontSize',7,'FontName','Helvetica')
            colormap(brewermap(64*8,'*Spectral'))
            caxis manual
            caxis([bottom top]);
        end
    end
    cb1 = colorbar('Ticks',0.005:0.04:0.405,'TickLabels',...
        {'0.005','0.045','0.085','0.125','0.165','0.205',...
        '0.245','0.285','0.325','0.365','0.405'});
    set(cb1,'FontSize',10.8,'FontName','Helvetica','Position',[.81,.126,.018,.823])
end
ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1) % sets ax1 to current axes
if k1 == 1
    text(0.045,0.9585+0.01,'a','FontSize',16,'FontName','Helvetica','FontWeight','bold')
elseif k1 == 2
    text(0.045,0.9585+0.01,'b','FontSize',16,'FontName','Helvetica','FontWeight','bold')
elseif k1 == 3
    text(0.045,0.9585+0.01,'c','FontSize',16,'FontName','Helvetica','FontWeight','bold')
elseif k1 == 4
    text(0.045,0.9585+0.01,'d','FontSize',16,'FontName','Helvetica','FontWeight','bold')
end

text(0.37,0.036,'Number of generations grown, T','FontSize',12,...
    'FontName','Helvetica');
text(0.075,0.445,'Number of time points, M','Rotation',90,'FontSize',12,...
    'FontName','Helvetica');
text(0.895,0.5,'Fit-Seq error, \epsilon','Rotation',90,'FontSize',10.8,...
    'FontName','Helvetica');

