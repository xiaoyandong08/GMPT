function Plot_GMPT_pipeline
addpath('../GMPT_Common_function')

rng(123456)
Diag_Type = 'Non-diagStable';
Cdiff = 1;

folder_time = '20210907T175417'; %20200813T112758
Lefse_or_ANCOM = 'ALDEx2'
simulation_index = 2

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

Plot_Cdiff_abundance(Num_sample,Cdiff,simulation_index,group_name,simulation_absolute_OTU_table_Cdiff,simulation_relative_OTU_table_Cdiff)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read data
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OTU_table = OTU_table{simulation_index};
OTU_table_Cdiff = OTU_table_Cdiff{simulation_index};

for i = 1 : length(Num_sample)
    for k = 1 : size(OTU_table,2)
        temp_OTU_table{k} = OTU_table{k}(:,1:Num_sample(i));
        temp_OTU_table_Cdiff{k} = OTU_table_Cdiff{k}(:,1:Num_sample(i));
    end
    aa = dir(folder_name);
    basf = regexp({aa.name},'sample','match');
    if find(cellfun(@(basf) ~isempty(basf),basf))
        Lefse_folder_name = [folder_name 'N' num2str(N) 'C' num2str(C) 'S' num2str(simulation_index) 'sample_' num2str(Num_sample(i))];
    else
        Lefse_folder_name = [folder_name 'N' num2str(N) 'C' num2str(C) 'S' num2str(simulation_index)];
    end
    [part_of_group_name{i},Fraction_present{i},List{i},RHO{i},PVAL{i}] = ddd(Lefse_folder_name,analysis_type,group_name,...
        temp_OTU_table,temp_OTU_table_Cdiff,A,r,Cdiff,Num_sample(i),Lefse_or_ANCOM);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sample_size_index = length(Num_sample); % the number of sample, sample size
first_order = N;
plot_matrix = List{sample_size_index}(1:first_order,3:size(List{sample_size_index},2));
index = find(List{sample_size_index}(:,1)==Cdiff);

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
view(-12,29)
ylabel('Phenotypes')
xlabel('Species (in descending order of their frequency in all pair-wise DAA)')
zlabel('Relative abundance')
title([Lefse_or_ANCOM ', Sim index = ' num2str(simulation_index)])

label = List{sample_size_index}(1:first_order,1);
label = label(label~=Cdiff);
figure;
subplot('position',[0.1 0.45 0.75 0.48])
h = imagesc(plot_matrix(2:end,:)');
set(gca,'YTick',[1:1:length(part_of_group_name{sample_size_index})],'YTickLabel',part_of_group_name{sample_size_index},'YTickLabelRotation',45,...
    'XTick',[1:first_order-1],'XTickLabel',label,'XTickLabelRotation',0,'fontsize',14,'TickDir','out');
set(gca,'fontsize',14);
ylim([0.5 0.5+length(part_of_group_name{sample_size_index})])
xlim([0.5 first_order-1+0.5])
ylabel('Phenotypes')
xlabel('Species (in descending order of their frequency in all pair-wise DAA)')
title([Lefse_or_ANCOM ', Sim index = ' num2str(simulation_index)])
set(gcf,'position',[30 142 1049 451])

figure;
subplot('position',[0.15 0.45 0.80 0.48])
h = imagesc(plot_matrix(List{sample_size_index}(1:first_order,1)==Cdiff,:)');
set(gca,'YTick',[1:1:length(part_of_group_name{sample_size_index})],'YTickLabel',part_of_group_name{sample_size_index},'YTickLabelRotation',45,...
    'TickDir','out');
set(gca,'fontsize',14);
ylim([0.5 0.5+length(part_of_group_name{sample_size_index})])
xlim([0.5 1+0.5])
ylabel('Phenotypes')
title([Lefse_or_ANCOM ', Sim index = ' num2str(simulation_index)])
set(gcf,'position',[30 142 35 451])
caxis([min(plot_matrix(:)) max(plot_matrix(:))])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
index = find(PVAL{sample_size_index} <= 1);
index = setdiff(index,find(List{sample_size_index}(:,1)==Cdiff));
data = [setdiff(List{sample_size_index}(:,1),Cdiff,'stable') RHO{sample_size_index}(index)' A(List{sample_size_index}(index,1),Cdiff)];

if strcmp(print_file,'true')
    filename = [Lefse_or_ANCOM '_Spearman_correlation_for_Fig2.txt'];
    taxa_name = data(:,1);
    Write_file_for_Spearman_correlation(filename,data,taxa_name)
end

figure
set(gcf,'position',[30 327 535 420])
subplot('position',[0.1 0.55 0.76 0.33])
hold on;box on;
set(gca,'YColor',hex2rgb('#262626'))
for i = 1 : size(data,1)
    if data(i,3) < 0
        stem(i ,data(i,2),'MarkerSize',12,'Color',hex2rgb('#0072BD'),'Marker','o','MarkerFaceColor',hex2rgb('#0072BD'),'LineStyle','-');
    elseif data(i,3) > 0
        stem(i,data(i,2),'MarkerSize',12,'Color',hex2rgb('#F74461'),'Marker','o','MarkerFaceColor',hex2rgb('#F74461'),'LineStyle','-');
    elseif data(i,3) == 0
        stem(i,data(i,2),'MarkerSize',12,'Color',hex2rgb('#CCCCCC'),'Marker','o','MarkerFaceColor',hex2rgb('#CCCCCC'),'LineStyle','-');
    end
end
ylabel('Spearman correlation');
xlabel('Species')
xlim([0 N-1+0.5])
ylim([-1.1 1.1])
set(gca,'fontsize',14,'TickDir','out')
set(gca,'XTick',[1:first_order-1],'XTickLabel',data(:,1),'XTickLabelRotation',0);

threshold = [1:1:N-1];
true_promoter = find(data(:,3)>0);
true_inhibitor = find(data(:,3)<0);
true_neutral = find(data(:,3)==0);


for i = 1 : length(threshold)

    true_label = zeros(3,N-1);
    true_label(1,find(sign(data(:,3))==-1)) = 1;
    true_label(2,find(sign(data(:,3))==0)) = 1;
    true_label(3,find(sign(data(:,3))==1)) = 1;
    
    infer = data(1:threshold(i),:);
    infer_promoter = find(infer(:,2)>0);
    infer_inhibitor = find(infer(:,2)<0);
    infer_neutral = setdiff([1:N-1],unique([infer_promoter' infer_inhibitor']));
  
    infer_label = zeros(3,N-1);
    infer_label(1,find(sign(infer(:,2))==-1)) = 1;
    infer_label(3,find(sign(infer(:,2))==1)) = 1;
    infer_label(2,setdiff([1:N-1],unique([infer_promoter' infer_inhibitor']))) = 1;
    
    if find(sum(infer_label,2) == 0)
        true_label = true_label(find(sum(infer_label,2) ~= 0),:);
        infer_label = infer_label(find(sum(infer_label,2) ~= 0),:);
    end
    [out] = evaluation(true_label,infer_label);
    
    f1_macro(i) = out.fscoreMacro/100;
    f1_micro(i) = out.fscoreMicro/100;

end
title([Lefse_or_ANCOM ', Sim index = ' num2str(simulation_index)])

subplot('position',[0.1 0.1 0.76 0.33])
hold on;box on;
plot(threshold,f1_macro,'-s','MarkerSize',12,'color',hex2rgb('#262626'),'MarkerFaceColor',hex2rgb('#262626'),'LineWidth',1)
plot(threshold,f1_micro,'-o','MarkerSize',12,'color',hex2rgb('#1AA85A'),'MarkerFaceColor',hex2rgb('#1AA85A'),'LineWidth',1)
legend('f1-marco','f1-mirco')
set(gca,'fontsize',14,'TickDir','out')
xlabel('Threshold');ylabel('F1 score')
xlim([0 N-1+0.5])
ylim([0.6 1.03])
% title(analysis_type)
set(gcf,'position',[30 142 1049 451])

if strcmp(print_file,'true')
    filename = [Lefse_or_ANCOM '_f1_macro_f1_micro_for_Fig2.txt'];
    data = [f1_macro' f1_micro'];
    Write_file_for_f1_macro_f1_micro(filename,data,threshold+1)
end


end


function [part_of_group_name,Fraction_present,List,RHO,PVAL] = ddd(Lefse_folder_name,analysis_type,group_name,OTU_table,OTU_table_Cdiff,A,r,Cdiff,Num_sample,Lefse_or_ANCOM)
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

if strcmp(Lefse_or_ANCOM,'Lefse')
    [Diff_taxa,ALL_present,Fraction_present] = Analysis_results_from_LefSe(Lefse_folder_name,A,part_of_group_name);
elseif strcmp(Lefse_or_ANCOM,'ANCOM')
    [Diff_taxa,ALL_present,Fraction_present] = Analysis_results_from_ANCOM(Lefse_folder_name,A,part_of_group_name);
elseif strcmp(Lefse_or_ANCOM,'ALDEx2')
    [Diff_taxa,ALL_present,Fraction_present] = Analysis_results_from_ALDEx2(Lefse_folder_name,A,part_of_group_name);
end

% [Diff_severity,diff_num_taxa,fraction_sign]=Output_table_of_differential_pair(part_of_group_name,part_of_OTU_table,Diff_taxa,new_severity,A,Cdiff);
% fraction_sign(isnan(fraction_sign))=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

% figure;hold on;box on
% x = [];y = [];
% x = [x;Fraction_present(setdiff(1:N,Cdiff))];
% y = [y;(A(setdiff(1:N,Cdiff),Cdiff))];
% scatter(x(find(y<0)),y(find(y<0)),'s','MarkerEdgeColor',hex2rgb('#0072BD'),'MarkerFaceColor',hex2rgb('#0072BD'))
% scatter(x(find(y>0)),y(find(y>0)),'s','MarkerEdgeColor',hex2rgb('#F74461'),'MarkerFaceColor',hex2rgb('#F74461'))
% scatter(x(find(y==0)),y(find(y==0)),'s','MarkerEdgeColor',hex2rgb('#CCCCCC'),'MarkerFaceColor',hex2rgb('#CCCCCC'))
% set(gca,'fontsize',14);
% xlabel('# of a species present in pair')
% ylabel('interaction strength of species i to C. diff')
% xlim([0 nchoosek(length(part_of_group_name),2)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


% first_order = 100;
% plot_matrix = List(1:first_order,3:size(List,2));
% figure;
% h = bar3(plot_matrix');
% for i = 1 : first_order
%     if A(List(i,1),1) > 0
%         h(i).FaceColor = hex2rgb('#F74461');
%     elseif A(List(i,1),1) < 0
%         h(i).FaceColor = hex2rgb('#0072BD');
%     elseif A(List(i,1),1) == 0
%         h(i).FaceColor = hex2rgb('#CCCCCC');
%     end
% end
% h(1).FaceColor = hex2rgb('#262626');
% set(gca,'YTick',[1:1:length(group_name)],'YTickLabel',group_name(severity),'YTickLabelRotation',45,...
%     'XTick',[1:first_order],'XTickLabel',{List(1:first_order,1)},'XTickLabelRotation',0,'fontsize',14);
% set(gca,'fontsize',14);
% ylim([0 1+length(group_name)])
% xlim([0.1 first_order+0.9])
% set(gcf,'position',[521 67 1194 911])
% view(-38,38)
% % xlabel('Phenotypes (from less severe to more severe)')
% % ylabel('Species index (frequency in pair comparasions from high to low)')



% figure;hold on;box on
% yyaxis left
% set(gca,'YColor',hex2rgb('#262626'))
% for i = 1 : N
%     if A(List(i,1),Cdiff) < 0
%         stem(i,RHO(i),'Color',hex2rgb('#0072BD'),'Marker','o','MarkerFaceColor',hex2rgb('#0072BD'),'LineStyle','-');
%     elseif A(List(i,1),Cdiff) > 0
%         stem(i,RHO(i),'Color',hex2rgb('#F74461'),'Marker','o','MarkerFaceColor',hex2rgb('#F74461'),'LineStyle','-');
%     elseif A(List(i,1),Cdiff) == 0
%         stem(i,RHO(i),'Color',hex2rgb('#CCCCCC'),'Marker','o','MarkerFaceColor',hex2rgb('#CCCCCC'),'LineStyle','-');
%     end
% end
% ylabel('Spearman correlation');
% xlabel('Species index (sorted by descending present count)')
% xlim([0 N+1])
% set(gcf,'position',[10 608 400 316])
% set(gca,'fontsize',14);
% window_size = 5;
% moving_var=  zeros(1,length(RHO) - window_size);
% for i = 1 : length(RHO) - window_size
%     moving_var(i) = var(abs(RHO(i:i+window_size-1)));
% end
% yyaxis right
% set(gca,'YColor',hex2rgb('#1AA85A'))
% plot([1:length(moving_var)],moving_var,'-','color',hex2rgb('#1AA85A'),'LineWidth',3)
% % plot([1:length(PVAL)],PVAL,'o','color',hex2rgb('#1AA85A'))
% ylabel('Moving variance')





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



end

function [Diff_taxa,ALL_present,Fraction_present] = Analysis_results_from_ANCOM(folder_path,A,part_of_group_name)

dirData = dir([folder_path '/ANCOM_res' '/*_res.txt']);
All_file_name = {dirData.name};
if length(All_file_name)<66
    aa = 1;
end
pair_name = nchoosek(part_of_group_name,2);
for i = 1 : length(All_file_name)
    temp = regexp(All_file_name{i}, '_', 'split');
    index = str2num(temp{2});
    Diff_taxa{index} = Read_ANCOM_results([folder_path '/ANCOM_res/' All_file_name{i}]);
    
end

ALL_present = [];
Fraction_present = zeros(size(A,1),1);
for i = 1 : size(A,1)
    count = 0;
    for j = 1 :  length(Diff_taxa)
        if length(intersect(i,Diff_taxa{1,j}(:)))>0
            count = count + 1;
        end
    end
    Fraction_present(i) = count;
    if count == length(Diff_taxa)
        ALL_present = [ALL_present i];
    end
end
end


function Diff_taxa = Read_ANCOM_results(filename)
Diff_taxa = zeros(0,1);
fid = fopen(filename);
nextline = fgetl(fid);
nextline = fgetl(fid);

count = 1;
while isstr(nextline) & length(nextline)>0
    temp = regexp(nextline, '\t', 'split');
    if strcmp(temp{7},'TRUE')
        temp1 = regexp(temp{2}, '"', 'split');
        temp1 = regexp(temp1{2}, 'OTU', 'split');
        Diff_taxa(count,1) = str2num(temp1{2});
        count = count + 1;
    end
    nextline = fgetl(fid);
end
fclose(fid);
end

function [Diff_taxa,ALL_present,Fraction_present] = Analysis_results_from_ALDEx2(folder_path,A,part_of_group_name)

dirData = dir([folder_path '/ALDEx2_res' '/*_res.txt']);
All_file_name = {dirData.name};
pair_name = nchoosek(part_of_group_name,2);
for i = 1 : length(All_file_name)
    temp = regexp(All_file_name{i}, '_', 'split');
    index = str2num(temp{2});
    Diff_taxa{index} = Read_ALDEx2_results([folder_path '/ALDEx2_res/' All_file_name{i}]);
    
end

ALL_present = [];
Fraction_present = zeros(size(A,1),1);
for i = 1 : size(A,1)
    count = 0;
    for j = 1 :  length(Diff_taxa)
        if length(intersect(i,Diff_taxa{1,j}(:)))>0
            count = count + 1;
        end
    end
    Fraction_present(i) = count;
    if count == length(Diff_taxa)
        ALL_present = [ALL_present i];
    end
end
end


function Diff_taxa = Read_ALDEx2_results(filename)
Diff_taxa = zeros(0,1);
fid = fopen(filename);
nextline = fgetl(fid);
nextline = fgetl(fid);


count = 1;
while isstr(nextline) & length(nextline)>0
    temp = regexp(nextline, '\t', 'split');
    if str2num(temp{15})>0.5||str2num(temp{15})<-0.5
        temp1 = regexp(temp{1}, '"', 'split');
        Diff_taxa(count,1) = str2num(temp1{2});
        count = count + 1;
    end
    nextline = fgetl(fid);
end
fclose(fid);

end


function [Diff_taxa,ALL_present,Fraction_present] = Analysis_results_from_LefSe(folder_path,A,part_of_group_name)

dirData = dir([folder_path '/*.res']);
All_file_name = {dirData.name};
if length(All_file_name)<66
    aa = 1;
end
pair_name = nchoosek(part_of_group_name,2);
for i = 1 : length(All_file_name)
    temp = regexp(All_file_name{i}, '_', 'split');
    temp1 = regexp(temp{end}, '.res', 'split');
    temp{end} =  temp1{1};
    temp_name = {temp{4} temp{6}};
    for j = 1 : size(pair_name,1)
        if sum(strcmp(temp_name,pair_name(j,:)))==2 || sum(strcmp(temp_name([2 1]),pair_name(j,:)))==2
            Diff_taxa{j} = Read_LefSe_results([folder_path '/' All_file_name{i}]);
            break
        end
    end
    
end

ALL_present = [];
Fraction_present = zeros(size(A,1),1);
for i = 1 : size(A,1)
    count = 0;
    for j = 1 :  length(Diff_taxa)
        if length(intersect(i,Diff_taxa{1,j}(:)))>0
            count = count + 1;
        end
    end
    Fraction_present(i) = count;
    if count == length(Diff_taxa)
        ALL_present = [ALL_present i];
    end
end
end


function Diff_taxa = Read_LefSe_results(filename)
Diff_taxa = zeros(0,1);
fid = fopen(filename);
nextline = fgetl(fid);

count = 1;
while isstr(nextline) & length(nextline)>0
    temp = regexp(nextline, '\t', 'split');
    if ~strcmp(temp{3},'')
        temp1 = regexp(temp{1}, '_', 'split');
        Diff_taxa(count,1) = str2num(temp1{2});
        count = count + 1;
    end
    nextline = fgetl(fid);
end
fclose(fid);

end