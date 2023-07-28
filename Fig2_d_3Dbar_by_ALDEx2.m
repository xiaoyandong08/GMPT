function Fig2_d_3Dbar_by_ALDEx2

rng(123456)
Diag_Type = 'Non_diagStable_connectivity';
Cdiff = 1;

folder_time = '20220429T222328';
Lefse_or_ANCOM = 'ALDEx2'
effect_threshold = 0.5;
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

OTU_table = OTU_table{simulation_index};
OTU_table_Cdiff = OTU_table_Cdiff{simulation_index};

for i = 1 : length(Num_sample)
    for k = 1 : size(OTU_table,2)
        temp_OTU_table{k} = OTU_table{k}(:,1:Num_sample(i));
        temp_OTU_table_Cdiff{k} = OTU_table_Cdiff{k}(:,1:Num_sample(i));
    end
    Lefse_folder_name = [folder_name 'N' num2str(N) 'S' num2str(simulation_index) 'sample_' num2str(Num_sample(i))];
    [part_of_group_name{i},Fraction_present{i},List{i},RHO{i},effect{i}] = ddd(Lefse_folder_name,analysis_type,part_pheno_group_name,...
        temp_OTU_table,temp_OTU_table_Cdiff,Cdiff,Lefse_or_ANCOM,group_name,effect_threshold);
end

if strcmp(Lefse_or_ANCOM,'ANCOM')
    [effect_TwoParts] = Analysis_FirstSecond_from_ANCOM(Lefse_folder_name,A);
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

[~,case_control_species_index] = sort(abs(case_control_effect),'descend');
case_control_effect = case_control_effect(case_control_species_index);

index = find(case_control_species_index==Cdiff);
case_control_species_index(index) = [];
case_control_effect(index) = [];
interaction = A(case_control_species_index,Cdiff);
case_control_data = [case_control_species_index' case_control_effect' interaction];

[~,TwoParts_species_index] = sort(abs(effect_TwoParts),'descend');
TwoParts_effect = effect_TwoParts(TwoParts_species_index);
index = find(TwoParts_species_index==Cdiff);
TwoParts_species_index(index) = [];
TwoParts_effect(index) = [];
interaction = A(TwoParts_species_index,Cdiff);
TwoParts_data = [TwoParts_species_index' TwoParts_effect' interaction];

% case_control_data(abs(case_control_data(:,2))<effect_threshold,2) = 0;
% TwoParts_data(abs(TwoParts_data(:,2))<effect_threshold,2) = 0;

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

figure;
h = bar3(plot_matrix',0.7);
for i = 1 : first_order
    if find(List{sample_size_index}(i,1)==Cdiff)
        h(i).FaceColor = hex2rgb('#1AA85A');
    else
        if A(List{sample_size_index}(i,1),1) > 0
            h(i).FaceColor = hex2rgb('#F74461');
        elseif A(List{sample_size_index}(i,1),1) < 0
            h(i).FaceColor = hex2rgb('#0072BD');
        elseif A(List{sample_size_index}(i,1),1) == 0
            h(i).FaceColor = hex2rgb('#CCCCCC');
        end
    end
end
set(gca,'YTick',[1:1:length(part_of_group_name{sample_size_index})],'YTickLabel',part_of_group_name{sample_size_index},'YTickLabelRotation',45,...
    'XTick',[1:first_order],'XTickLabel',{List{sample_size_index}(1:first_order,1)},'XTickLabelRotation',0,'fontsize',14);
set(gca,'fontsize',14);
ylim([0 1+length(part_of_group_name{sample_size_index})])
xlim([0.1 first_order+0.9])
set(gcf,'position',[648 67 999 911])
set(gcf,'position',[648 67 1474 911])
view(-12,29)
ylabel('Phenotypes')
xlabel('Species (in descending order of their frequency in all pair-wise DAA)')
zlabel('Relative abundance')
title([Lefse_or_ANCOM ', Sim index = ' num2str(simulation_index)])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end


% function [part_of_group_name,Fraction_present,List,RHO,PVAL,effect] = ddd(Lefse_folder_name,analysis_type,group_name,OTU_table,OTU_table_Cdiff,A,r,Cdiff,Num_sample,Lefse_or_ANCOM,effect_threshold)
function [part_of_group_name,Fraction_present,List,RHO,effect]      = ddd(Lefse_folder_name,analysis_type,part_pheno_group_name,OTU_table,OTU_table_Cdiff,Cdiff,Lefse_or_ANCOM,group_name,effect_threshold)

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
    [Diff_taxa,ALL_present,Fraction_present,effect] = Analysis_results_from_ANCOM(Lefse_folder_name,N,part_pheno_group_name,group_name);
elseif strcmp(Lefse_or_ANCOM,'ALDEx2')
    [Diff_taxa,ALL_present,Fraction_present,effect] = Analysis_results_from_ALDEx2(Lefse_folder_name,N,part_pheno_group_name,group_name,effect_threshold);
end


[a,b] = sort(Fraction_present,'descend');
List = [b a];
for i = 1 : size(List,1)
    b = List(i,1);
    for j = 1 : size(part_of_OTU_table_Cdiff,2)
        temp = mean(part_of_OTU_table_Cdiff{j},2);
        List(i,j+2) = temp(b);
    end
end
for i = 1 : size(List,1)
    [RHO(i)] = corr(List(i,3:size(List,2))',[1:length(group_name)]','Type','Spearman');
    PVAL(i) = 0;
    %[RHO(i),PVAL(i)] = corr(List(i,3:size(List,2))',List(1,3:size(List,2))','Type','Spearman');
end

% [a,b] = sort(Fraction_present,'descend');
% List = [b a];
% 
% all_table = cell2mat(part_of_OTU_table_Cdiff);
% for i = 1 : size(List,1)
%     b = List(i,1);
%     List(i,3:size(all_table,2)+2) = all_table(b,:);
% end
% 
% index = find(List(:,1)==Cdiff);
% for i = 1 : size(List,1)
%     a = List(i,3:size(List,2))';
%     b = List(index,3:size(List,2))';
%     [RHO(i)] = corr(a,b,'Type','Spearman');
%     PVAL(i) = 0;
%     %[RHO(i),PVAL(i)] = corr(List(i,3:size(List,2))',List(1,3:size(List,2))','Type','Spearman');
% end

end



