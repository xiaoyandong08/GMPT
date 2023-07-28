function FigS1_b2g_GMPT_performance_by_ANCOM

rng(123456)
Diag_Type = 'Non_diagStable_connectivity';
Cdiff = 1;

folder_time = '20220429T222328';
Lefse_or_ANCOM = 'ANCOM'
wscore_threshold = 0.6;
simulation_index = 5

analysis_type    = 'relative';
print_file = 'false';

aa = dir(Diag_Type);
aa = {aa.name};
basf = regexp(aa,folder_time,'match');
index = find(cellfun(@(basf) ~isempty(basf),basf));
folder_name = [Diag_Type '/' aa{index}];
aa = dir(folder_name);
aa = {aa.name};
basf = regexp(aa,'Many_times','match');
index = find(cellfun(@(basf) ~isempty(basf),basf));
load([folder_name '/' aa{index(end)}])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = A{simulation_index};
r = r{simulation_index};

OTU_table = OTU_table{simulation_index};
OTU_table_Cdiff = OTU_table_Cdiff{simulation_index};
simulation_relative_OTU_table_Cdiff = relative_OTU_table_Cdiff{simulation_index};
simulation_absolute_OTU_table_Cdiff = OTU_table_Cdiff{simulation_index};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
exchange_name = [];
for i = 1 : length(exchange_rate)
    WB_name{i} = ['WB-' num2str(exchange_rate(i))];
    BW_name{i} = ['BW-' num2str(exchange_rate(i))];
    exchange_name = [exchange_name num2str(exchange_rate(i)) '_'];
end

group_name = {'W',WB_name{1:end},BW_name{1:end},'B'};
if strcmp(analysis_type,'relative')
    folder_name = ['./' folder_name '/' analysis_type '_With_Cdiff_exchange_' exchange_name '/'];
elseif strcmp(analysis_type,'absolute')
    folder_name = ['./' folder_name '/With_Cdiff_exchange_' exchange_name abundance_type '/'];
end

part_pheno_group_name = group_name;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



for i = 1 : length(Num_sample)
    for k = 1 : size(OTU_table,2)
        temp_OTU_table{k} = OTU_table{k}(:,1:Num_sample(i));
        temp_OTU_table_Cdiff{k} = OTU_table_Cdiff{k}(:,1:Num_sample(i));
    end
    Lefse_folder_name = [folder_name 'N' num2str(N) 'S' num2str(simulation_index) 'sample_' num2str(Num_sample(i))];
    [part_of_group_name{i},Fraction_present{i},List{i},RHO{i},effect{i}] = ddd(Lefse_folder_name,analysis_type,part_pheno_group_name,...
        temp_OTU_table,temp_OTU_table_Cdiff,Cdiff,Lefse_or_ANCOM,group_name,wscore_threshold);
end

if strcmp(Lefse_or_ANCOM,'ANCOM')
    [effect_TwoParts] = Analysis_FirstSecond_from_ANCOM(Lefse_folder_name,A,wscore_threshold);
elseif strcmp(Lefse_or_ANCOM,'ALDEx2')
    [effect_TwoParts] = Analysis_FirstSecond_from_ALDEx2(Lefse_folder_name,A,effect_threshold);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% case_control_pair = {'W','BW-0.1'};
case_control_pair = {part_of_group_name{1}{1} part_of_group_name{1}{end}};

case_control_effect = effect{1};
all_pair = nchoosek(group_name,2);
for i = 1 : size(all_pair,1)
    if (strcmp(all_pair(i,1),case_control_pair{1})&&strcmp(all_pair(i,2),case_control_pair{2}))||(strcmp(all_pair(i,1),case_control_pair{2})&&strcmp(all_pair(i,2),case_control_pair{1}))

        case_control_index = i;
        break;
    end
end
case_control_effect=case_control_effect{case_control_index};

[~,case_control_species_index] = sort(abs(case_control_effect(:,2)),'descend');
case_control_effect = case_control_effect(case_control_species_index,:);

index = find(case_control_species_index==Cdiff);
case_control_species_index(index) = [];
case_control_effect(index,:) = [];
interaction = A(case_control_species_index,Cdiff);
case_control_data = [case_control_effect interaction];

OTU_table_Cdiff_control = simulation_relative_OTU_table_Cdiff{strcmp(group_name,case_control_pair{1})};
OTU_table_Cdiff_case    = simulation_relative_OTU_table_Cdiff{strcmp(group_name,case_control_pair{2})};
for i = 1 : size(case_control_data,1)
    if case_control_data(i,3)~=0
        speccies_index = case_control_data(i,1);
        case_abundance = mean(OTU_table_Cdiff_case(speccies_index,:));
        control_abundance = mean(OTU_table_Cdiff_control(speccies_index,:));
        if case_abundance>control_abundance
            case_control_data(i,2) = 1 * case_control_data(i,2);
        elseif case_abundance<control_abundance
            case_control_data(i,2) = -1 * case_control_data(i,2);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,TwoParts_species_index] = sort(abs(effect_TwoParts(:,2)),'descend');
TwoParts_effect = effect_TwoParts(TwoParts_species_index,:);
index = find(TwoParts_species_index==Cdiff);
TwoParts_species_index(index) = [];
TwoParts_effect(index,:) = [];
interaction = A(TwoParts_species_index,Cdiff);
TwoParts_data = [TwoParts_effect interaction];

for i = 1 : length(group_name)
    relative_abundance_Cdiff(i) = mean(simulation_relative_OTU_table_Cdiff{i}(Cdiff,:));
end

[a,b] = sort(relative_abundance_Cdiff,'ascend');
group_name(b);
OTU_table_Cdiff_part_control = cell2mat(simulation_relative_OTU_table_Cdiff(b(1:length(b)/2)));
OTU_table_Cdiff_part_case    = cell2mat(simulation_relative_OTU_table_Cdiff(b(1+length(b)/2:end)));


for i = 1 : size(TwoParts_data,1)
    if TwoParts_data(i,3)~=0
        speccies_index = TwoParts_data(i,1);
        case_abundance = mean(OTU_table_Cdiff_part_case(speccies_index,:));
        control_abundance = mean(OTU_table_Cdiff_part_control(speccies_index,:));
        if case_abundance>control_abundance
            TwoParts_data(i,2) = 1 * TwoParts_data(i,2);
        elseif case_abundance<control_abundance
            TwoParts_data(i,2) = -1 * TwoParts_data(i,2);
        end
    end
end
% case_control_data(:,2) = case_control_data(:,2).*case_control_data(:,3);
case_control_data(:,3) = [];
% TwoParts_data(:,2) = TwoParts_data(:,2).*TwoParts_data(:,3);
TwoParts_data(:,3) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sample_size_index = length(Num_sample); % the number of sample, sample size
first_order = N;
plot_matrix = List{sample_size_index}(1:first_order,3:size(List{sample_size_index},2));
if strcmp(print_file,'true')
    filename = [Lefse_or_ANCOM '_MPT_list_for_Fig2.txt'];
    taxa_name = List{sample_size_index}(1:first_order,1);
    phenotype_name = part_of_group_name{sample_size_index};
    Write_file_for_MPT_list(filename,plot_matrix,phenotype_name,taxa_name)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

index = [1:N];%find(PVAL{sample_size_index} <= 1);
index = setdiff(index,find(List{sample_size_index}(:,1)==Cdiff));
data = [setdiff(List{sample_size_index}(:,1),Cdiff,'stable') RHO{sample_size_index}(index)' A(List{sample_size_index}(index,1),Cdiff)];
GMPT_data = data;
if strcmp(print_file,'true')
    filename = [Lefse_or_ANCOM '_Spearman_correlation_for_Fig2.txt'];
    taxa_name = GMPT_data(:,1);
    Write_file_for_Spearman_correlation(filename,GMPT_data,taxa_name)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
set(gcf,'position',[30 142 1049 451])
set(gcf,'Position',[40 90 996 1042])
set(gcf,'position',[31 534 799 840])
set(gcf,'position',[31 612 810 733])
% subplot('position',[0.1 0.75 0.76 0.2])
% subaxis(3,1,1,'SpacingVertical', 0.05, 'PaddingTop', 0.02)
subplot(3,1,1)
hold on;box on;
% rectangle('Position',[0 0.5 N-1+0.5 1.1-0.5],'EdgeColor', 'n', 'FaceColor', [hex2rgb('F74461') 0.3])
% rectangle('Position',[0 -1.1 N-1+0.5 1.1-0.5],'EdgeColor', 'n', 'FaceColor', [hex2rgb('0072BD') 0.3])
% plot([0 N-1+0.5],[0.5 0.5],'k--')
% plot([0 N-1+0.5],[-0.5 -0.5],'k--')
for i = 1 : size(GMPT_data,1)
    if GMPT_data(i,3) < 0
        stem(i ,GMPT_data(i,2),'MarkerSize',12,'Color',hex2rgb('#0072BD'),'Marker','o','MarkerFaceColor',hex2rgb('#0072BD'),'LineStyle','-');
    elseif GMPT_data(i,3) > 0
        stem(i,GMPT_data(i,2),'MarkerSize',12,'Color',hex2rgb('#F74461'),'Marker','o','MarkerFaceColor',hex2rgb('#F74461'),'LineStyle','-');
    elseif GMPT_data(i,3) == 0
        stem(i,GMPT_data(i,2),'MarkerSize',12,'Color',hex2rgb('#CCCCCC'),'Marker','o','MarkerFaceColor',hex2rgb('#CCCCCC'),'LineStyle','-');
    end
end
ylabel('Spearman correlation');
xlabel('Species (in descending order of their frequency in all pair-wise DAA)')
xlim([0 N-1+0.5])
ylim([-1.1 1.1])
set(gca,'fontsize',14,'TickDir','out')
set(gca,'XTick',[1:first_order-1],'XTickLabel',GMPT_data(:,1),'XTickLabelRotation',0);
title([Lefse_or_ANCOM '-GMPT, Sim index = ' num2str(simulation_index)])

% subplot('position',[0.1 0.1 0.76 0.33])
% subaxis(3,1,2,'SpacingVertical', 0.05, 'PaddingTop', 0.02)
subplot(3,1,2)
hold on;box on;
for i = 1 : size(case_control_data,1)
    if case_control_data(i,3) < 0
        stem(i ,case_control_data(i,2),'MarkerSize',12,'Color',hex2rgb('#0072BD'),'Marker','o','MarkerFaceColor',hex2rgb('#0072BD'),'LineStyle','-');
    elseif case_control_data(i,3) > 0
        stem(i,case_control_data(i,2),'MarkerSize',12,'Color',hex2rgb('#F74461'),'Marker','o','MarkerFaceColor',hex2rgb('#F74461'),'LineStyle','-');
    elseif case_control_data(i,3) == 0
        stem(i,case_control_data(i,2),'MarkerSize',12,'Color',hex2rgb('#CCCCCC'),'Marker','o','MarkerFaceColor',hex2rgb('#CCCCCC'),'LineStyle','-');
    end
end
ylim([-N-5 N+5])
plot([0 N-1+0.5],[N*wscore_threshold N*wscore_threshold],'k--')
plot([0 N-1+0.5],[-N*wscore_threshold -N*wscore_threshold],'k--')
ylabel('W score');
xlabel('Species (in descending order of absolute values of W score)')
xlim([0 N-1+0.5])
set(gca,'fontsize',14,'TickDir','out','ytick',[-N:10:N])
set(gca,'XTick',[1:first_order-1],'XTickLabel',case_control_data(:,1),'XTickLabelRotation',0);
title([Lefse_or_ANCOM '-MWAS_1, Sim index = ' num2str(simulation_index) ', ' case_control_pair{1} ', ' case_control_pair{2}])

% subaxis(3,1,3,'SpacingVertical', 0.05, 'PaddingTop', 0.02)
subplot(3,1,3)
hold on;box on;
for i = 1 : size(TwoParts_data,1)
    if TwoParts_data(i,3) < 0
        stem(i ,TwoParts_data(i,2),'MarkerSize',12,'Color',hex2rgb('#0072BD'),'Marker','o','MarkerFaceColor',hex2rgb('#0072BD'),'LineStyle','-');
    elseif TwoParts_data(i,3) > 0
        stem(i,TwoParts_data(i,2),'MarkerSize',12,'Color',hex2rgb('#F74461'),'Marker','o','MarkerFaceColor',hex2rgb('#F74461'),'LineStyle','-');
    elseif TwoParts_data(i,3) == 0
        stem(i,TwoParts_data(i,2),'MarkerSize',12,'Color',hex2rgb('#CCCCCC'),'Marker','o','MarkerFaceColor',hex2rgb('#CCCCCC'),'LineStyle','-');
    end
end
ylim([-N-5 N+5])
plot([0 N-1+0.5],[N*wscore_threshold N*wscore_threshold],'k--')
plot([0 N-1+0.5],[-N*wscore_threshold -N*wscore_threshold],'k--')
ylabel('W score');
xlabel('Species (in descending order of absolute values of W score)')
xlim([0 N-1+0.5])
set(gca,'fontsize',14,'TickDir','out','ytick',[-N:10:N])
set(gca,'XTick',[1:first_order-1],'XTickLabel',TwoParts_data(:,1),'XTickLabelRotation',0);
title([Lefse_or_ANCOM '-MWAS_2, Sim index = ' num2str(simulation_index)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case_control_data(abs(case_control_data(:,2))<N*wscore_threshold,2) = 0;
TwoParts_data(abs(TwoParts_data(:,2))<N*wscore_threshold,2) = 0;

dot_color = [hex2rgb('#0072BD');hex2rgb('#F74461');hex2rgb('#1AA85A');hex2rgb('#262626')];

num_non_zero = sum(GMPT_data(:,3)~=0);
TwoKinds_Perfect_accuracy = [1:num_non_zero repmat(num_non_zero,1,N-1-num_non_zero)]/num_non_zero;
for kk = 1 : N-1
    
    temp_data1 = GMPT_data(1:kk,:);
    temp_data1 = temp_data1(temp_data1(:,3)~=0,:);
    TwoKinds_accuracy_GMPT(1,kk) = sum(sign(temp_data1(:,2))-sign(temp_data1(:,3))==0)/num_non_zero;
    
    temp_data2 = case_control_data(1:kk,:);
    temp_data2 = temp_data2(temp_data2(:,3)~=0,:);
    TwoKinds_accuracy_case_control(1,kk) = sum(sign(temp_data2(:,2))-sign(temp_data2(:,3))==0)/num_non_zero;
    
    temp_data3 = TwoParts_data(1:kk,:);
    temp_data3 = temp_data3(temp_data3(:,3)~=0,:);
    TwoKinds_accuracy_TwoParts(1,kk) = sum(sign(temp_data3(:,2))-sign(temp_data3(:,3))==0)/num_non_zero;
end


TwoKinds_accuracy=[TwoKinds_accuracy_GMPT;TwoKinds_accuracy_case_control;TwoKinds_accuracy_TwoParts;TwoKinds_Perfect_accuracy];

save('N30_ANCOM_GMPT.txt','GMPT_data','-ascii');
save('N30_ANCOM_CaseControl.txt','case_control_data','-ascii');
save('N30_ANCOM_TwoParts.txt','TwoParts_data','-ascii');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Binary = readtable('N30_Binary_classification_metrics_by_ANCOM.txt');

methods = {'GMPT','Case_control','Two_part'};

for m = 1 : length(methods)
        
    index = find(strcmp(Binary.Type,methods{m}));
    if length(index)~=N-1
        aa = 1;
    end
    f1(m,:) = Binary.f1(index);
    gm(m,:) = Binary.gm(index);
    
end

dot_color = [hex2rgb('#FF00FF');hex2rgb('#FF8C00');hex2rgb('#26A9E1');hex2rgb('#A8A8AC')];
figure;
set(gcf,'position',[31 534 799 840])
set(gcf,'position',[31 612 810 733])
subplot(3,1,1)
hold on;box on;
plot([1:N-1],f1(1,:),'linewidth',4,'color',dot_color(1,:))
plot([1:N-1],f1(2,:),'linewidth',4,'color',dot_color(2,:))
plot([1:N-1],f1(3,:),'linewidth',4,'color',dot_color(3,:))
set(gca,'fontsize',14,'TickDir','out','XMinorTick','on')
%legend('GMPT','Tradition','Two parts','Location','best')
ylabel('F1-score')
xlabel('Top-K species')
title([Lefse_or_ANCOM '-GMPT, Sim index = ' num2str(simulation_index)])
set(gca,'Ytick',[0:0.2:1])

subplot(3,1,2)
hold on;box on;
plot([1:N-1],gm(1,:),'linewidth',4,'color',dot_color(1,:))
plot([1:N-1],gm(2,:),'linewidth',4,'color',dot_color(2,:))
plot([1:N-1],gm(3,:),'linewidth',4,'color',dot_color(3,:))
set(gca,'fontsize',14,'TickDir','out','XMinorTick','on')
%legend('GMPT','Tradition','Two parts','Location','best')
ylabel('Geometry mean')
xlabel('Top-K species')
title([Lefse_or_ANCOM '-MWAS_1, Sim index = ' num2str(simulation_index) ', ' case_control_pair{1} ', ' case_control_pair{2}])
set(gca,'Ytick',[0:0.2:1])

subplot(3,1,3)
hold on;box on;
plot([1:N-1],TwoKinds_accuracy_GMPT,'linewidth',4,'color',dot_color(1,:))
plot([1:N-1],TwoKinds_accuracy_case_control,'linewidth',4,'color',dot_color(2,:))
plot([1:N-1],TwoKinds_accuracy_TwoParts,'linewidth',4,'color',dot_color(3,:))
plot([1:N-1],TwoKinds_Perfect_accuracy,'linewidth',4,'color',dot_color(4,:))
set(gca,'fontsize',14,'TickDir','out','XMinorTick','on')
legend('GMPT','MWAS_1','MWAS_2','Perfect inference','Location','best')
ylabel('Accuracy')
xlabel('Top-K species')
title([Lefse_or_ANCOM '-MWAS_2, Sim index = ' num2str(simulation_index)])
set(gca,'Ytick',[0:0.2:1])
end


% function [part_of_group_name,Fraction_present,List,RHO,PVAL,effect] = ddd(Lefse_folder_name,analysis_type,group_name,OTU_table,OTU_table_Cdiff,A,r,Cdiff,Num_sample,Lefse_or_ANCOM,effect_threshold)
function [part_of_group_name,Fraction_present,List,RHO,effect]      = ddd(Lefse_folder_name,analysis_type,part_pheno_group_name,OTU_table,OTU_table_Cdiff,Cdiff,Lefse_or_ANCOM,group_name,wscore_threshold)

N = size(OTU_table{1},1);

for j = 1 : size(OTU_table,2)
    relative_OTU_table{j} = bsxfun(@times,OTU_table{j},1./sum(OTU_table{j},1));
    relative_OTU_table_Cdiff{j} = bsxfun(@times,OTU_table_Cdiff{j},1./sum(OTU_table_Cdiff{j},1));
    
    absolute_abundance_Cdiff(j) = mean(OTU_table_Cdiff{j}(Cdiff,:));
    relative_abundance_Cdiff(j) = mean(relative_OTU_table_Cdiff{j}(Cdiff,:));
end
new_severity = [1:1:size(OTU_table,2)];

if strcmp(analysis_type,'relative')
    [aa,severity]=sort(relative_abundance_Cdiff,'ascend');
    b = severity(new_severity);
    part_of_OTU_table_Cdiff = relative_OTU_table_Cdiff(b);

elseif strcmp(analysis_type,'absolute')
    [aa,severity]=sort(absolute_abundance_Cdiff,'ascend');
    b = severity(new_severity);
    part_of_OTU_table_Cdiff = OTU_table_Cdiff(b);

end
part_of_group_name = group_name(b);

for j = 1 : size(part_of_OTU_table_Cdiff,2)
    part_of_abundance_Cdiff(j) = mean(part_of_OTU_table_Cdiff{j}(Cdiff,:));
    part_of_abundance_error_Cdiff(j) = std(part_of_OTU_table_Cdiff{j}(Cdiff,:),1);
end

%% show C diff abundance among selected phenotype
% figure;hold on
% bar([1:1:length(part_of_group_name)],part_of_abundance_Cdiff,0.6)
% errorbar(part_of_abundance_Cdiff,part_of_abundance_error_Cdiff,'.r')
% set(gca,'Xlim',[0.5 length(part_of_group_name)+0.5],'XTick',[1:1:length(part_of_group_name)],'XTickLabel',part_of_group_name,'XticklabelRotation',45,'fontsize',10);
% ylim([0 max(ylim)])
% ylabel('abundance of C diff')
% title([analysis_type 'abundance'])




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Lefse_folder_name = [folder_name 'N' num2str(N) 'C' num2str(C) 'S' num2str(simulation_index) 'sample_' num2str(Num_sample)];


if strcmp(Lefse_or_ANCOM,'ANCOM')
    [Diff_taxa,ALL_present,Fraction_present,effect] = Analysis_results_from_ANCOM(Lefse_folder_name,N,part_pheno_group_name,group_name,wscore_threshold);
elseif strcmp(Lefse_or_ANCOM,'ALDEx2')
    [Diff_taxa,ALL_present,Fraction_present,effect] = Analysis_results_from_ALDEx2(Lefse_folder_name,N,part_pheno_group_name,group_name,effect_threshold);
end


% [a,b] = sort(Fraction_present,'descend');
% List = [b a];
% for i = 1 : size(List,1)
%     b = List(i,1);
%     for j = 1 : size(part_of_OTU_table_Cdiff,2)
%         temp = mean(part_of_OTU_table_Cdiff{j},2);
%         List(i,j+2) = temp(b);
%     end
% end
% for i = 1 : size(List,1)
%     [RHO(i)] = corr(List(i,3:size(List,2))',[1:length(group_name)]','Type','Spearman');
%     PVAL(i) = 0;
%     %[RHO(i),PVAL(i)] = corr(List(i,3:size(List,2))',List(1,3:size(List,2))','Type','Spearman');
% end

[a,b] = sort(Fraction_present,'descend');
List = [b a];

all_table = cell2mat(part_of_OTU_table_Cdiff);
for i = 1 : size(List,1)
    b = List(i,1);
    List(i,3:size(all_table,2)+2) = all_table(b,:);
end

index = find(List(:,1)==Cdiff);
for i = 1 : size(List,1)
    a = List(i,3:size(List,2))';
    b = List(index,3:size(List,2))';
    [RHO(i)] = corr(a,b,'Type','Spearman');
    PVAL(i) = 0;
    %[RHO(i),PVAL(i)] = corr(List(i,3:size(List,2))',List(1,3:size(List,2))','Type','Spearman');
end

end



