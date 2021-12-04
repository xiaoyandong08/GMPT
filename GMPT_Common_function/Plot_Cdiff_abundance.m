function Plot_Cdiff_abundance(Num_sample,Cdiff,simulation_index,group_name,simulation_absolute_OTU_table_Cdiff,simulation_relative_OTU_table_Cdiff)

for i = 1 : length(group_name)
    simulation_absolute_abundance_Cdiff(i) = mean(simulation_absolute_OTU_table_Cdiff{i}(Cdiff,:));
    simulation_relative_abundance_Cdiff(i) = mean(simulation_relative_OTU_table_Cdiff{i}(Cdiff,:));
end

[a,b] = sort(simulation_relative_abundance_Cdiff,'ascend');
sort_group_name = group_name(b);
simulation_absolute_OTU_table_Cdiff = simulation_absolute_OTU_table_Cdiff(b);
simulation_relative_OTU_table_Cdiff = simulation_relative_OTU_table_Cdiff(b);

color = parula(length(sort_group_name)); % Generate color values

for i = 1 : size(simulation_relative_OTU_table_Cdiff,2)
    relative_abundance_Cdiff{i} = simulation_relative_OTU_table_Cdiff{i}(Cdiff,:);
end
relative_abundance_Cdiff = cell2mat(relative_abundance_Cdiff')';

for i = 1 : size(relative_abundance_Cdiff,2) - 1
    [h,p(i)] = ttest2(relative_abundance_Cdiff(:,i),relative_abundance_Cdiff(:,i+1));
end

bar_color = parula(length(sort_group_name)); % Generate color values
% bar_color = cellfun(@hex2rgb,{'#caedd8','#f7ebe7','#f4dfd4','#f4d3c3','#f4bfab','#f2a48c','#f78e72','#f76d52'},'UniformOutput',false);
% bar_color = cell2mat(bar_color');
mean_relative_abundance_Cdiff = mean(relative_abundance_Cdiff,1);
figure;hold on;box on
b = bar([1:1:length(mean_relative_abundance_Cdiff)],mean_relative_abundance_Cdiff,0.6);
b.FaceColor = 'flat';
b.CData = bar_color;

errorbar(mean(relative_abundance_Cdiff,1),std(relative_abundance_Cdiff,1,1)/sqrt(Num_sample),'.k')
set(gca,'Xlim',[0.5 length(mean_relative_abundance_Cdiff)+0.5],'XTick',[1:1:length(mean_relative_abundance_Cdiff)],'XTickLabelRotation',45,'XTickLabel',sort_group_name,'fontsize',10);
ylabel('Relative abundance of species X')
set(gca,'fontsize',14,'TickDir','out')
set(gcf,'position',[30 19 535 437])
set(gca,'yscale','log', 'YMinorTick','on')
ylim([1e-6 1e-1])
yt = get(gca, 'YTick');
ytkvct = sort(10.^(-linspace(1, 10*size(yt,2), 10*size(yt,2))));
set(gca, 'YTick', ytkvct);
set(gcf,'position',[39 467 647 239])
%p-value<0.05 (*), <0.01(**), <0.001(***); >0.05(NS)
% for i = 1 : size(relative_abundance_Cdiff,2) - 1
%     line([i i+1],[mean_relative_abundance_Cdiff(i+1)*1.1 mean_relative_abundance_Cdiff(i+1)*1.1],'Color','k')
%     if p(i)>0.05
%         significance = 'NS';
%     elseif p(i)<=0.05 & p(i)>0.01
%         significance = '*';
%     elseif p(i)<=0.01 & p(i)>0.001
%         significance = '**';
%     elseif p(i)<=0.001 
%         significance = '***';
%     end
%     text(mean([i i+1]),mean_relative_abundance_Cdiff(i+1)*1.1,significance,'HorizontalAlignment','center','fontsize',14)
% end
% title(['Sim = ' num2str(simulation_index)])
%%%%%%%

AllSteady = cell2mat(simulation_relative_OTU_table_Cdiff)';
D = squareform(My_rJSD(AllSteady));
% D = squareform(pdist(AllSteady));


figure;hold on;box on;
[Y,e] = cmdscale(D);
% Y = tsne(D);
for i = 1 : length(sort_group_name)
    scatter(Y((i-1)*Num_sample+1:i*Num_sample,1),Y((i-1)*Num_sample+1:i*Num_sample,2),...
    70,'Marker','o','MarkerEdgeColor',hex2rgb('262626'),'MarkerFaceColor',color(i,:),'LineWidth',0.25);
end
legend(sort_group_name,'location','best')
title('PcoA')
set(gca,'fontsize',14,'TickDir','out')
xlabel('PC1');ylabel('PC2')
set(gcf,'position',[566 48 535 437])
set(gcf,'position',[39 780 647 475])

end