%%  Simulated data for Figure 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normal distribution
lineage = 1e5;
x_ini_1 = pearsrnd(0,0.1,0,3,lineage,1);  % normal distribution

% Left skewed distribution
x_ini_tempt = pearsrnd(0,0.15,2.5,20,lineage,1);
pos_delete = x_ini_tempt>0.49 | x_ini_tempt<-0.49;
x_ini_tempt(pos_delete) = normrnd(0,0.05,sum(pos_delete),1);
x_ini_2 = -(x_ini_tempt-mean(x_ini_tempt));
sum(x_ini_2<-0.5 | x_ini_2>0.5)

% Right skewed distribution
x_ini_tempt = pearsrnd(0,0.15,2.5,20,lineage,1);
pos_delete = x_ini_tempt>0.49 | x_ini_tempt<-0.49;
x_ini_tempt(pos_delete) = normrnd(0,0.05,sum(pos_delete),1);
x_ini_3 = x_ini_tempt-mean(x_ini_tempt);
sum(x_ini_3<-0.5 | x_ini_3>0.5)


t_evo = 50;
cell_num_ini = 1e4*ones(lineage,1);
for i1 = 1:3
    if i1 == 1
        x_ini = x_ini_1;
    elseif i1 == 2
        x_ini = x_ini_2;
    elseif i1 == 3
        x_ini = x_ini_3;
    end
    
    pick_across_fitness = [0.4, 0.25, 0.12, 0, -0.2, -0.4]; % pick_across_rank
    pick_across_fitness_index = nan(size(pick_across_fitness));
    for j1 = 1:length(pick_across_fitness)
        pick_across_fitness_index(j1) = ...
            find(abs(x_ini-pick_across_fitness(j1))...
            ==min(abs(x_ini-pick_across_fitness(j1))));
        x_ini(pick_across_fitness_index(j1)) = pick_across_fitness(j1);
    end
    
    % Determinstic simulation (Simplified)
    noise_option = 0;
    [ file_name1 ] = evolution_0mut_simplified(lineage, t_evo, cell_num_ini, ...
        x_ini, noise_option, 'format','mat');
end



%%  Figure 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
line_color1 = [239,0,0]/255;
line_color2 = [0,184,229]/255;
line_color3 = [120,120,120]/255;
line_color4 = [180,180,180]/255;
line_color5 = line_color2;
line_color6 = line_color1;

figure('rend','painters','pos',[1 11 644 500]);
% Figure 1a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('data_evo_simu_0mut_simplified_20180516-111123819.mat'); x_ini_1 = x_ini;
load('data_evo_simu_0mut_simplified_20180516-111125454.mat'); x_ini_2 = x_ini;
load('data_evo_simu_0mut_simplified_20180516-111126924.mat'); x_ini_3 = x_ini;

pick_across_fitness = [0.4, 0.25, 0.12, 0, -0.2, -0.4]; % pick_across_rank
stepsize = 0.01;
x_range_dist = -0.5:stepsize:0.5;
y_simu_1 = histcounts(x_ini_1,x_range_dist)/(lineage*stepsize);
y_simu_2 = histcounts(x_ini_2,x_range_dist)/(lineage*stepsize);
y_simu_3 = histcounts(x_ini_3,x_range_dist)/(lineage*stepsize);
x_range_dist_2 = (x_range_dist(1:end-1)+x_range_dist(2:end))/2;

y_theory_1 = normpdf(x_range_dist_2,0,0.1);

x_ini_tempt = pearsrnd(0,0.15,2.5,20,1e7,1);
pos_delete = x_ini_tempt>0.5 | x_ini_tempt<-0.5;
x_ini_tempt(pos_delete) = normrnd(0,0.05,sum(pos_delete),1);
x_ini_2_theory = x_ini_tempt-mean(x_ini_tempt);
y_theory_2 = histcounts(-x_ini_2_theory,x_range_dist)/(1e7*stepsize);
y_theory_3 = histcounts(x_ini_2_theory,x_range_dist)/(1e7*stepsize);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subaxis(3,3,1, 'MarginBottom',0.05,'MarginTop',0.02,'MarginLeft',0.14,...
    'MarginRight',0.25,'SpacingHoriz',0.06,'SpacingVert',0.05);
bar(x_range_dist_2,y_simu_1,1,'EdgeColor','none','FaceAlpha',0.1,...
    'FaceColor',line_color4)
hold on
plot(x_range_dist_2,y_theory_1,'LineWidth',3, 'Color',line_color4);
bar(x_range_dist_2,y_simu_2,1,'EdgeColor','none','FaceAlpha',0.1,...
    'FaceColor',line_color3)
plot(x_range_dist_2,y_theory_2,'-','LineWidth',3,'Color',line_color3);
plot([pick_across_fitness(1)+stepsize/4,pick_across_fitness(1)+stepsize/4],...
    [0,4.5], '--','Color',line_color1, 'LineWidth',3);
plot([pick_across_fitness(2)+stepsize/4,pick_across_fitness(2)+stepsize/4],...
    [0,4.5], '-.','Color',line_color1, 'LineWidth',3);
plot([pick_across_fitness(3)+stepsize/4,pick_across_fitness(3)+stepsize/4],...
    [0,4.5], ':','Color',line_color1, 'LineWidth',3);
plot([pick_across_fitness(4)+stepsize/4,pick_across_fitness(4)+stepsize/4],...
    [0,4.5], '--','Color',line_color2, 'LineWidth',3);
plot([pick_across_fitness(5)+stepsize/4,pick_across_fitness(5)+stepsize/4],...
    [0,4.5], '-.','Color',line_color2, 'LineWidth',3);
plot([pick_across_fitness(6)+stepsize/4,pick_across_fitness(6)+stepsize/4],...
    [0,4.5], ':','Color',line_color2, 'LineWidth',3);
hold off
box on
pbaspect([1 1 1])
xlim([-0.5,0.5])
ylim([0,4.5])
xticks(-0.5:0.25:0.5)
xticklabels({'-0.5','-0.25','0','0.25','0.5'})
set(gca,'FontSize',10.8,'FontName','Helvetica')  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subaxis(3,3,2, 'MarginBottom',0.05,'MarginTop',0.02,'MarginLeft',0.14,...
    'MarginRight',0.25,'SpacingHoriz',0.06,'SpacingVert',0.05);
bar(x_range_dist_2,y_simu_1,1,'EdgeColor','none','FaceAlpha',0.1,'FaceColor',line_color3)
hold on
plot(x_range_dist_2,y_theory_1,'-','LineWidth',3,'Color',line_color3);
plot([pick_across_fitness(1)+stepsize/4,pick_across_fitness(1)+stepsize/4],...
    [0,4.5], '--','Color',line_color1, 'LineWidth',3);
plot([pick_across_fitness(2)+stepsize/4,pick_across_fitness(2)+stepsize/4],...
    [0,4.5], '-.','Color',line_color1, 'LineWidth',3);
plot([pick_across_fitness(3)+stepsize/4,pick_across_fitness(3)+stepsize/4],...
    [0,4.5], ':','Color',line_color1, 'LineWidth',3);
plot([pick_across_fitness(4)+stepsize/4,pick_across_fitness(4)+stepsize/4],...
    [0,4.5], '--','Color',line_color2, 'LineWidth',3);
plot([pick_across_fitness(5)+stepsize/4,pick_across_fitness(5)+stepsize/4],...
    [0,4.5], '-.','Color',line_color2, 'LineWidth',3);
plot([pick_across_fitness(6)+stepsize/4,pick_across_fitness(6)+stepsize/4],...
    [0,4.5], ':','Color',line_color2, 'LineWidth',3);
hold off
box on
pbaspect([1 1 1])
xlim([-0.5,0.5])
ylim([0,4.5])
xticks(-0.5:0.25:0.5)
xticklabels({'-0.5','-0.25','0','0.25','0.5'})
set(gca,'FontSize',10.8,'FontName','Helvetica')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subaxis(3,3,3, 'MarginBottom',0.05,'MarginTop',0.02,'MarginLeft',0.14,...
    'MarginRight',0.25,'SpacingHoriz',0.06,'SpacingVert',0.05);
bar(x_range_dist_2,y_simu_1,1,'EdgeColor','none','FaceAlpha',0.1,...
    'FaceColor',line_color4)
hold on
plot(x_range_dist_2,y_theory_1,'LineWidth',3, 'Color',line_color4);
bar(x_range_dist_2,y_simu_3,1,'EdgeColor','none','FaceAlpha',0.1,...
     'FaceColor',line_color3)
plot(x_range_dist_2,y_theory_3,'-','LineWidth',3,'Color',line_color3);
h = zeros(1,5);
h(1) = plot([pick_across_fitness(1)+stepsize/4,pick_across_fitness(1)+stepsize/4],...
    [0,4.5], '--','Color',line_color1, 'LineWidth',3);
h(2) = plot([pick_across_fitness(2)+stepsize/4,pick_across_fitness(2)+stepsize/4],...
    [0,4.5], '-.','Color',line_color1, 'LineWidth',3);
h(3) = plot([pick_across_fitness(3)+stepsize/4,pick_across_fitness(3)+stepsize/4],...
    [0,4.5], ':','Color',line_color1, 'LineWidth',3);
h(4) = plot([pick_across_fitness(4)+stepsize/4,pick_across_fitness(4)+stepsize/4],...
    [0,4.5], '--','Color',line_color2, 'LineWidth',3);
h(5) = plot([pick_across_fitness(5)+stepsize/4,pick_across_fitness(5)+stepsize/4],...
    [0,4.5], '-.','Color',line_color2, 'LineWidth',3);
h(6) = plot([pick_across_fitness(6)+stepsize/4,pick_across_fitness(6)+stepsize/4],...
    [0,4.5], ':','Color',line_color2, 'LineWidth',3);
hold off
box on
pbaspect([1 1 1])
xlim([-0.5,0.5])
ylim([0,4.5])
xticks(-0.5:0.25:0.5)
xticklabels({'-0.5','-0.25','0','0.25','0.5'})
set(gca,'FontSize',10.8,'FontName','Helvetica')

[hleg, hobj, hout, mout]=legend(h,'  0.4','  0.25','  0.12',...
    '  0','  -0.2','  -0.4');
hL = findobj(hobj,'type','line');
set(hL,'XData',[0.05,0.545]) % 0.6 is used to control width of line label

set(hleg,'Position',[0.78   0.737   0.12   0.21])
temp = [hleg; hleg.ItemText];
set(temp, 'FontName', 'Helvetica','FontSize',10.8)


% Figure 1b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k1 = 1:3
    if k1 == 1
        load('data_evo_simu_0mut_simplified_20180516-111125454.mat');
    elseif k1 == 2
        load('data_evo_simu_0mut_simplified_20180516-111123819.mat');
    elseif k1 == 3
        load('data_evo_simu_0mut_simplified_20180516-111126924.mat');
    end

    pick_across_fitness = [0.4, 0.25, 0.12, 0, -0.2, -0.4]; % pick_across_rank
    pick_across_fitness_index = zeros(size(pick_across_fitness));
    for i1 = 1:length(pick_across_fitness)
        if sum(x_ini==pick_across_fitness(i1)) == 1
            pick_across_fitness_index(i1) = find(x_ini==pick_across_fitness(i1));
        end
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate fold enrichment for each genotype
    fold_enrichment_mat = nan(lineage,t_evo+1);
    for i2 = 1:(t_evo+1)
        fold_enrichment_mat(:,i2) = log2(cell_num_evo(:,i2)./cell_num_evo(:,1));
    end
    fold_enrichment_mat_modify = fold_enrichment_mat;
    fold_enrichment_mat_modify(fold_enrichment_mat<=-8)=-8.5;
    fold_enrichment_mat_modify(fold_enrichment_mat>8)=8.5;
    % converge <=-8 into one row and converge >=8 into one row
    
    % Calculate distribution of fold enrichment
    x_range_1 = [-8.5,-8:0.25:8,8.5];
    fold_enrichment_dis = nan(length(x_range_1)-1,t_evo+1);
    for i2 = 1:(t_evo+1)
        fold_enrichment_dis(:,i2) ...
            = histcounts(fold_enrichment_mat_modify(:,i2), x_range_1)/lineage;
    end
    fold_enrichment_dis_log = log10(fold_enrichment_dis);
    fold_enrichment_dis_log(fold_enrichment_dis==0)=NaN;
        
    % Calculate fold enrichment trajectory of 5 picked genotypes
    bar_x_pick_across_fitness = ...
        discretize(fold_enrichment_mat_modify(pick_across_fitness_index,:),x_range_1);
    bar_y_pick_across_fitness = nan(length(pick_across_fitness_index),t_evo+1);
    for i = 1:length(pick_across_fitness_index)
        for j = 1:t_evo
            bar_y_pick_across_fitness(i,j) ...
                = fold_enrichment_dis(bar_x_pick_across_fitness(i,j),j);
        end
    end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subaxis(3,3,k1+3, 'MarginBottom',0.05,'MarginTop',0.02,'MarginLeft',0.14,...
        'MarginRight',0.25,'SpacingHoriz',0.06,'SpacingVert',0.015);
    set(gcf, 'PaperPosition', [2 1 4 2]);
    x_range_2 = -8:0.25:8;
    % ax1 = axes;
    imagesc(0:t_evo,(x_range_2(1:end-1)+x_range_2(2:end))/2, ...
        fold_enrichment_dis_log(2:end-1,:));
    axis('xy')
    xticks(0:10:40)
    xticklabels({'0','10','20','30','40'})
    yticks([-8.75+0.125,-6:3:6,8.75-0.125])
    yticklabels({'<-8','-6','-3','0','3','6','>8'})
    hold on
    imagesc(0:t_evo, -8.625, fold_enrichment_dis_log(1,:));
    imagesc(0:t_evo, 8.625, fold_enrichment_dis_log(end,:));
    caxis([-5, 0])
   
    colormap(brewermap(64*8,'*RdBu'))
    if k1 == 3   
        cb1 = colorbar;
        set(cb1,'YTick',-5:0,'YTickLabel',{'10^{-5}','10^{-4}','10^{-3}',...
            '10^{-2}','10^{-1}','10^0'},'FontSize',10.8,'FontName','Helvetica');
        set(cb1,'Position',[0.78 0.4095 0.025 0.2105])
        shading interp
    end

    fold_enrichment_mat_pick_across_fitness = ...
        fold_enrichment_mat(pick_across_fitness_index,:);
    fold_enrichment_mat_pick_across_fitness...
        (fold_enrichment_mat_pick_across_fitness>8) = 8.5+0.125;
    fold_enrichment_mat_pick_across_fitness...
        (fold_enrichment_mat_pick_across_fitness<=-8) = -8.5-0.125;
    
    x_range_4 = [-9.5+0.25,-8:0.25:8,9.5-0.25];
    % plot([-0.5,40.5], [0,0],'-', 'LineWidth',3,'Color', line_color3,'MarkerSize',15)
    plot(0:t_evo, fold_enrichment_mat_pick_across_fitness(1,:),...
        '--', 'LineWidth',3,'Color', line_color1,'MarkerSize',15)
    plot(0:t_evo, fold_enrichment_mat_pick_across_fitness(2,:),...
        '-.', 'LineWidth',3,'Color', line_color1,'MarkerSize',15)
    plot(0:t_evo, fold_enrichment_mat_pick_across_fitness(3,:),...
        ':', 'LineWidth',3,'Color', line_color1,'MarkerSize',15)
    plot(0:t_evo, fold_enrichment_mat_pick_across_fitness(4,:),...
        '--', 'LineWidth',3,'Color', line_color2,'MarkerSize',15)
    plot(0:t_evo, fold_enrichment_mat_pick_across_fitness(5,:),...
        '-.', 'LineWidth',3,'Color', line_color2,'MarkerSize',15)
    plot(0:t_evo, fold_enrichment_mat_pick_across_fitness(6,:),...
        ':', 'LineWidth',3,'Color', line_color2,'MarkerSize',15)
    hold off
    xlim([-0.5,40.5])
    ylim([-9.25+0.125,9.25-0.125])
    box on
    pbaspect([1 1 1])
    set(gca,'FontSize',10.8,'FontName','Helvetica')
    

end


% Figure 1c
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k1 = 1:3
    if k1 == 1
        load('data_evo_simu_0mut_simplified_20180516-111125454.mat');
    elseif k1 == 2
        load('data_evo_simu_0mut_simplified_20180516-111123819.mat');
    elseif k1 == 3
        load('data_evo_simu_0mut_simplified_20180516-111126924.mat');
    end
    
    subaxis(3,3,k1+6, 'MarginBottom',0.05,'MarginTop',0.02,'MarginLeft',0.14,...
        'MarginRight',0.25,'SpacingHoriz',0.06,'SpacingVert',0);
    hold on
    plot(0:40, x_mean(1:41),'Color',line_color3,'LineWidth',3)
    plot([-0.5,40.5], [pick_across_fitness(1),pick_across_fitness(1)],...
        '--', 'LineWidth',3,'Color', line_color1,'MarkerSize',15)
    plot([-0.5,40.5], [pick_across_fitness(2),pick_across_fitness(2)],...
        '-.', 'LineWidth',3,'Color', line_color1,'MarkerSize',15)
    plot([-0.5,40.5], [pick_across_fitness(3),pick_across_fitness(3)],...
        ':', 'LineWidth',3,'Color', line_color1,'MarkerSize',15)
    plot([-0.5,40.5], [pick_across_fitness(4),pick_across_fitness(4)],...
        '--', 'LineWidth',3,'Color', line_color2,'MarkerSize',15)
    hold off
    box on
    pbaspect([1 0.8 1])
    xlim([-0.5,40.5])
    ylim([0,0.55])
    xticks(0:10:40)
    xticklabels({'0','10','20','30','40'})
    yticks(0:0.2:0.6)
    yticklabels({'0','0.2','0.4','0.6'})
    set(gca,'FontSize',10.8,'FontName','Helvetica')
    
end

ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1)

text(0.0425,0.966,'a','FontSize',15,'FontName','Helvetica','FontWeight','bold');
text(0.0425,0.64,'b','FontSize',15,'FontName','Helvetica','FontWeight','bold')
text(0.0425,0.309,'c','FontSize',15,'FontName','Helvetica','FontWeight','bold')

text(0.415,0.728-0.056,'Fitness','FontSize',12,'FontName','Helvetica')
text(0.3675,0.401-0.056,'Time (generations)','FontSize',12,'FontName','Helvetica')
text(0.3675,0.112-0.056,'Time (generations)','FontSize',12,'FontName','Helvetica')

text(0.08,0.738,'Probability density','Rotation',90,'FontSize',12,'FontName','Helvetica');
text(0.08,0.4,'Fold-enrichment (log_2)','Rotation',90,'FontSize',12,'FontName','Helvetica');
text(0.08,0.138,'Mean fitness','Rotation',90,'FontSize',12,'FontName','Helvetica');      

text(0.865,0.384,'Frequency of genotypes','Rotation',90,'FontSize',12,'FontName','Helvetica');


