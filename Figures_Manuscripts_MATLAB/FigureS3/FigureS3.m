%%  Simulated data for Figure S3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('file_name_list1_20180805-200034240.mat')
read_depth_average = 2e2;
file_name_list1_FoldEnrichment = cell(size(file_name_list1));

for i1 = 1:size(file_name_list1,1)
    for i2 = 1:size(file_name_list1,2)
        load(file_name_list1{i1,i2})   
        [ file_name_list1_FoldEnrichment{i1,i2} ] = ...
            evolution_0mut_complex(lineage,t_evo,cell_num_ini,x_ini,...
            deltat,read_depth_average,noise_option,'format','mat');
    end
end
dt = datestr(now,'yyyymmdd-HHMMSSFFF');
file_name1 = ['file_name_list1_FoldEnrichment_' dt '.mat'];
save(file_name1,'file_name_list1_FoldEnrichment','x_ini_list','cell_num_ini_list')




load('file_name_list1_FoldEnrichment_20180805-204848651.mat')
lineage = 5e3;
x_estimate_result_mat_linear = cell(size(file_name_list1_FoldEnrichment));
x_estimate_result_mat_linear(:) = {nan(lineage,1)};
for i3 = 1:size(file_name_list1_FoldEnrichment,1)
    for i4 = 1:size(file_name_list1_FoldEnrichment,2)
        load(file_name_list1_FoldEnrichment{i3,i4})
        BC_num_mat_tempt = num_mat_sequencing(:,[1,4]);
        BC_num_mat_tempt(BC_num_mat_tempt(:,1) == 0,1) = 1;
        BC_num_mat = BC_num_mat_tempt;
        read_depth_keep = sum(BC_num_mat);
        t_seq_vec = t_seq_vec([1,4]);
        BC_freq_mat_tempt = BC_num_mat./repmat(read_depth_keep,lineage,1);
        fit_param1 = nan(lineage,2);
        for i = 1:lineage
            pos = BC_freq_mat_tempt(i,:)~=0;
            if sum(pos)>=2
                fit_param1(i,:) = polyfit(t_seq_vec(pos),...
                    log(BC_freq_mat_tempt(i,pos)),1);
            end
        end
        x_estimate_result_linear_tempt = exp(fit_param1(:,1))-1;
        x_estimate_result_mat_linear{i3,i4} = (1+x_estimate_result_linear_tempt)/...
            (1+nansum(BC_num_mat(:,1).*x_estimate_result_linear_tempt)./...
            sum(BC_num_mat(:,1)))-1;
    end
end
dt = datestr(now,'yyyymmdd-HHMMSSFFF');
file_name2 = ['x_estimate_result_mat_FoldEnrichment_' dt '.mat'];
save(file_name2,'x_estimate_result_mat_linear')




%%   Figure S3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('file_name_list1_FoldEnrichment_20180805-204848651.mat')
load('x_estimate_result_mat_FoldEnrichment_20180805-204914992.mat')

mean_error_mat = zeros([size(file_name_list1_FoldEnrichment),3]);
for i3 = 1:size(file_name_list1_FoldEnrichment,1)
    for i4 = 1:size(file_name_list1_FoldEnrichment,2)
        load(file_name_list1_FoldEnrichment{i3,i4})
        if i4 == 1
            mean_error_mat(i3,i4,1) = nanmean(abs(x_estimate_result_mat_linear{i3,i4}-x_ini));
        elseif i4 == 2
            pos = cell_num_ini == min(cell_num_ini);
            mean_error_mat(i3,i4,2) = nanmean(abs(x_estimate_result_mat_linear{i3,i4}(pos)-x_ini(pos)));
            mean_error_mat(i3,i4,3) = nanmean(abs(x_estimate_result_mat_linear{i3,i4}(~pos)-x_ini(~pos)));
        end
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
line_color3 = [249,145,55]/255; % Orange
line_color4 = [110,169,36]/255; % Green

line_color5 = [120,120,120]/255;
line_color6 = [180,180,180]/255;

sigma_list = [0.15,0.25,0.35, round(std(x_ini_list{4}),2),round(std(x_ini_list{5}),2)];

figure('rend','painters','pos',[1 11 950 800]);
for i1 = 1:size(file_name_list1_FoldEnrichment,1)
    load(file_name_list1_FoldEnrichment{i1,1})
    
    subaxis(1,5,i1, 'MarginBottom',0.8,'MarginTop',0.01,...
        'MarginLeft',0.3,'MarginRight',0.05,'SpacingHoriz',0.035,'SpacingVert',0.035);
    
    x_lim = max(abs(x_ini));
    if i1 <= 3
        stepsize = (sigma_list(i1)/5);
        x_range = -x_lim:stepsize:x_lim;
        y_dist = histcounts(x_ini,x_range)/(lineage*stepsize);
        y_theory = normpdf(x_range,0,sigma_list(i1));
    elseif i1 == 4
        stepsize = 0.025;
        x_range = -x_lim:stepsize:x_lim;
        y_theory_tempt1 = normpdf(x_range,0.315,0.15);
        y_theory_tempt2 = normpdf(x_range,-0.315,0.15);
        y_theory = y_theory_tempt1/2 + y_theory_tempt2/2;
        y_dist = histcounts(x_ini,x_range)/(lineage*stepsize);
        
    elseif i1 == 5
        stepsize = 0.025;
        x_range = -x_lim:stepsize:x_lim;
        y_theory_tempt1 = normpdf(x_range,0.41,0.11);
        y_theory_tempt2 = normpdf(x_range,-0.41,0.11);
        y_theory_tempt3 = normpdf(x_range,0,0.11);
        y_theory = y_theory_tempt1*(round(lineage/3)/lineage) + ...
            y_theory_tempt2*(round(lineage/3)/lineage) + ...
            y_theory_tempt3*((lineage-2*round(lineage/3))/lineage);
        y_dist = histcounts(x_ini,x_range)/(lineage*stepsize);
    end
    
    bar((x_range(1:end-1)+x_range(2:end))/2,y_dist,1,'EdgeColor','none',...
        'FaceAlpha',0.3,'FaceColor',line_color6)
    hold on
    plot(x_range, y_theory, 'Color',line_color6, 'LineWidth',3)
    hold off
    box on
    pbaspect([1 1 1])
    xlim([-1.5,1.5])
    ylim([0,3])
    xticks(-1.5:0.75:1.5)
    xticklabels({'-1.5','-0.75','0','0.75','1.5'})
    set(gca,'FontSize',10.8,'FontName','Helvetica')
    text(0.8,2.6,['\sigma = ', num2str(sigma_list(i1),'%1.2f')],...
        'FontSize',10.8,'FontName','Helvetica','HorizontalAlignment','center')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    

for i1 = 1:size(file_name_list1_FoldEnrichment,1)
    if i1 == 1
        pos1 = 0.6;
    elseif i1 == 2
        pos1 = 0.9;
    elseif i1 == 3
        pos1 = 1.2;
    elseif i1 == 4
        pos1 = 0.9;
    elseif i1 == 5
        pos1 = 0.9;
    end
    
    for i2 = 1:size(file_name_list1_FoldEnrichment,2)
        load(file_name_list1_FoldEnrichment{i1,i2})
        sz = 15*ones(size(x_ini));
        
        if i2 == 1
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            subaxis(4,5,i1, 'MarginBottom',0.1,'MarginTop',0.27,...
                'MarginLeft',0.3,'MarginRight',0.05,'SpacingHoriz',0.035,...
                'SpacingVert',0.035);
            
            x_est_tempt = x_estimate_result_mat_linear{i1,i2};
            pos_nan = isnan(x_est_tempt);
            x_est_tempt(pos_nan) = -pos1*1.1;
            
            plot([-10,10],[-10,10],':','LineWidth',3,'Color',line_color6)
            hold on
            scatter(x_ini(~pos_nan),x_est_tempt(~pos_nan),...
                sz(~pos_nan), 'MarkerEdgeColor','none', 'MarkerFaceColor',line_color3, ...
                'MarkerFaceAlpha',0.1, 'MarkerEdgeAlpha',0, 'LineWidth',.1);        
            scatter(x_ini(pos_nan),x_est_tempt(pos_nan),...
                sz(pos_nan), 'MarkerEdgeColor','none', 'MarkerFaceColor',line_color5, ...
                'MarkerFaceAlpha',0.1, 'MarkerEdgeAlpha',0, 'LineWidth',.1);        
            hold off
            box on
            pbaspect([1 1 1])
            xticks(-pos1:pos1:pos1)
            xticklabels({num2str(-pos1),'0',num2str(pos1)})
            yticks([-pos1*1.1,0,pos1])
            yticklabels({'-Inf','0',num2str(pos1)})
            xlim([-pos1*1.25,pos1*1.25])
            ylim([-pos1*1.25,pos1*1.25])
            set(gca,'FontSize',10.8,'FontName','Helvetica')
            text(-pos1,pos1,['\epsilon',' = ' num2str(mean_error_mat(i1,i2,1),'% 1.3f')],...
                    'FontSize',10.8,'FontName','Helvetica','Color',line_color3)
            text(0.65*pos1,-0.88*pos1,[ num2str(sum(pos_nan)/lineage*100,'% 2.2f'),'%'],...
                'FontSize',10.8,'FontName','Helvetica','Color',line_color5*0.5,...
                'HorizontalAlignment','center')
            
                       
        elseif i2 == 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            subaxis(4,5,i1+5, 'MarginBottom',0.1,'MarginTop',0.27,...
                'MarginLeft',0.3,'MarginRight',0.05,'SpacingHoriz',0.035,...
                'SpacingVert',0.035);
            
            x_est_tempt = x_estimate_result_mat_linear{i1,i2};
            pos_nan = isnan(x_est_tempt);
            x_est_tempt(pos_nan) = -pos1*1.1;    
            pos = cell_num_ini == min(cell_num_ini);
            plot([-10,10],[-10,10],':','LineWidth',3,'Color',line_color6)
            hold on 
  
            scatter(x_ini(pos&(~pos_nan)),x_est_tempt(pos&(~pos_nan)),...
                sz(pos&(~pos_nan)), 'MarkerEdgeColor','none', 'MarkerFaceColor',line_color4, ...
                'MarkerFaceAlpha',0.1, 'MarkerEdgeAlpha',0, 'LineWidth',.1);
            scatter(x_ini((~pos)&(~pos_nan)),x_est_tempt((~pos)&(~pos_nan)),...
                sz((~pos)&(~pos_nan)), 'MarkerEdgeColor','none', 'MarkerFaceColor',line_color3, ...
                'MarkerFaceAlpha',0.1, 'MarkerEdgeAlpha',0, 'LineWidth',.1);
            scatter(x_ini(pos_nan),x_est_tempt(pos_nan),...
                sz(pos_nan), 'MarkerEdgeColor','none', 'MarkerFaceColor',line_color5, ...
                'MarkerFaceAlpha',0.1, 'MarkerEdgeAlpha',0, 'LineWidth',.1); 
            hold off
            box on
            pbaspect([1 1 1])
            
            xticks(-pos1:pos1:pos1)
            xticklabels({num2str(-pos1),'0',num2str(pos1)})
            yticks([-pos1*1.1,0,pos1])
            yticklabels({'-Inf','0',num2str(pos1)})
            xlim([-pos1*1.25,pos1*1.25])
            ylim([-pos1*1.25,pos1*1.25])
            set(gca,'FontSize',10.8,'FontName','Helvetica')
            text(-pos1,pos1,['\epsilon',' = ' num2str(mean_error_mat(i1,i2,3),'% 1.3f')],...
                'FontSize',10.8,'FontName','Helvetica','Color',line_color3)
            text(-pos1,0.7*pos1,['\epsilon',' = ' num2str(mean_error_mat(i1,i2,2),'% 1.3f')],...
                'FontSize',10.8,'FontName','Helvetica','Color',line_color4)
            
            
            text(0.65*pos1,-0.88*pos1,[ num2str(sum(pos_nan)/lineage*100,'% 2.2f'),'%'],...
                'FontSize',10.8,'FontName','Helvetica','Color',line_color5*0.5, ...
                'HorizontalAlignment','center')    
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i2 = 1:size(file_name_list1_FoldEnrichment,2)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subaxis(4,1,i2, 'MarginBottom',0.1,'MarginTop',0.27,'MarginLeft',...
        0.07,'MarginRight',0.828,'SpacingHoriz',0.035,'SpacingVert',0.035);
            
    if i2 == 1
        bar(0.5,1,0.4,'EdgeColor','none','FaceAlpha',0.9,'FaceColor',line_color3)
    elseif i2 == 2
        bar(0,0.5,0.4,'EdgeColor','none','FaceAlpha',0.9,'FaceColor',line_color4)
        hold on
        bar(1,0.5,0.4,'EdgeColor','none','FaceAlpha',0.9,'FaceColor',line_color3)
        hold off
    end
    box on
    pbaspect([1 1 1])
    
    xlim([-0.5,1.5])
    if i2 == 1
        ylim([0,1.16])
        xticks(0.5)
        xticklabels('100')
        yticks(0:0.5:1)
        yticklabels({'0','0.5','1'})
        
    elseif i2 ==2 
        ylim([0,0.58])
        xticks([0,1])
        xticklabels({'20','180'})
        yticks(0:0.25:0.5)
        yticklabels({'0','0.25','0.5'})
    end
    set(gca,'FontSize',10.8,'FontName','Helvetica')
end


ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1)

text(0.21,0.636+0.1,'a','FontSize',15,'FontName','Helvetica','FontWeight','bold');
text(0.21,0.47+0.1,'b','FontSize',15,'FontName','Helvetica','FontWeight','bold')


text(0.299+0.033,0.838-0.045,'Fitness','FontSize',12,'FontName','Helvetica')
text(0.436+0.033,0.838-0.045,'Fitness','FontSize',12,'FontName','Helvetica')
text(0.299+0.137*2+0.033,0.838-0.045,'Fitness','FontSize',12,'FontName','Helvetica')
text(0.299+0.137*3+0.033,0.838-0.045,'Fitness','FontSize',12,'FontName','Helvetica')
text(0.299+0.137*4+0.033,0.838-0.045,'Fitness','FontSize',12,'FontName','Helvetica')

text(0.288-0.035,0.895,{'Probability' ; 'density function'},'Rotation',...
    90,'FontSize',12,'FontName','Helvetica','HorizontalAlignment','center');


text(0.299+0.02,0.44-0.045,'True fitness','FontSize',12,'FontName','Helvetica')
text(0.436+0.02,0.44-0.045,'True fitness','FontSize',12,'FontName','Helvetica')
text(0.299+0.137*2+0.02,0.44-0.045,'True fitness','FontSize',12,'FontName','Helvetica')
text(0.299+0.137*3+0.02,0.44-0.045,'True fitness','FontSize',12,'FontName','Helvetica')
text(0.299+0.137*4+0.02,0.44-0.045,'True fitness','FontSize',12,'FontName','Helvetica')


text(0.288-0.045,0.665,{'Estimated fitness'; '(fold-enrichment)'},'Rotation',...
    90,'FontSize',12,'FontName','Helvetica','HorizontalAlignment','center');
text(0.288-0.045,0.498,{'Estimated fitness' ; '(fold-enrichment)'},'Rotation',...
    90,'FontSize',12,'FontName','Helvetica','HorizontalAlignment','center');


text(0.07,0.44-0.045,'Initial cell number','FontSize',12,'FontName','Helvetica')

text(0.025,0.636+0.025,{'Frequency'},'Rotation',...
    90,'FontSize',12,'FontName','Helvetica','HorizontalAlignment','center');
text(0.025,0.47+0.025,{'Frequency'},'Rotation',...
    90,'FontSize',12,'FontName','Helvetica','HorizontalAlignment','center');


