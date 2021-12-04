function Many_times_different_sample_size(N,C,delta,diag,VarianceType,exchange_rate,Num_sample,simulation)
% N = [100];
% C = [0.4];
% delta = [0.2];
% diag = -1;
% VarianceType = 2; % 1: Travis; 2: Variance
% exchange_rate = [0.05 0.1 0.2 0.3 0.5];
% Num_sample = [5 10 15 20];
% simulation = 5;

FunctionType = 1;
abundance_type = 'absolute';
if VarianceType == 1
    Cdiff_disease_abundance = 0.2;
else
    Cdiff_disease_abundance = 0.2;
end
Cdiff_health_abundance = 1e-3;
Strongest_inhibitor_present = 0;
select_white_black_mixed = 'mixed';

folder_subname = datestr(datetime('now'),30);
if VarianceType == 1
    father_folder = ['DiagStable/' 'N' num2str(N) '_C' num2str(C) '_Pheno' num2str(length(exchange_rate)*2+2) '_' folder_subname];
else
    father_folder = ['Non-diagStable/' 'N' num2str(N) '_C' num2str(C) '_Pheno' num2str(length(exchange_rate)*2+2) '_' folder_subname];
end
for i = 1 : length(exchange_rate)
    WB_name{i} = ['WB-' num2str(exchange_rate(i))]; 
    BW_name{i} = ['BW-' num2str(exchange_rate(i))];
end
group_name = {'W',WB_name{1:end},BW_name{1:end},'B'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1 : simulation
    i
    [A{i},r{i},OTU_table{i},OTU_table_Cdiff{i},relative_OTU_table{i},relative_OTU_table_Cdiff{i}]=...
        Microbe_Phenotype_Triangulation_v3(N,C,delta,diag,VarianceType,FunctionType,Num_sample,i,exchange_rate,abundance_type,father_folder,...
        Cdiff_disease_abundance,Cdiff_health_abundance,Strongest_inhibitor_present,select_white_black_mixed);
end
save(['./' father_folder '/Many_times_' abundance_type '_and_exchange_rate_N' num2str(N) '_C' num2str(C) '_' folder_subname '.mat'])

end


function [A,r,OTU_table,OTU_table_Cdiff,relative_OTU_table,relative_OTU_table_Cdiff]=Microbe_Phenotype_Triangulation_v3...
    (N,C,delta,diag,VarianceType,FunctionType,Num_sample,simulation,exchange_rate,abundance_type,father_folder,...
    Cdiff_disease_abundance,Cdiff_health_abundance,Strongest_inhibitor_present,select_white_black_mixed)

h1 = 0.5;
h2 = 0.1;
time = [0:0.01:50];
Cdiff = 1;

exchange_name = [];
for i = 1 : length(exchange_rate)
    WB_name{i} = ['WB-' num2str(exchange_rate(i))]; 
    BW_name{i} = ['BW-' num2str(exchange_rate(i))];
    exchange_name = [exchange_name num2str(exchange_rate(i)) '_'];
end
group_name = {'W',WB_name{1:end},BW_name{1:end},'B'};
perfect_group_name = {'W',WB_name{1:end},BW_name{end:-1:1},'B'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while 1
    while 1
        %%% generate ecological network
        [A,r] = Generate_Network_A_of_BandW(N,C,delta,diag,VarianceType,time,FunctionType,h1,h2,Cdiff,abundance_type,Cdiff_disease_abundance,Cdiff_health_abundance,select_white_black_mixed);
%         load 'special-A.mat';
        %%% generate B, W
        [X_W,X_B,X_W_Cdiff,X_B_Cdiff,regenerate_A]=Generate_Samples_W_and_B(A,r,time,FunctionType,h1,h2,max(Num_sample),Cdiff,abundance_type,Cdiff_disease_abundance,Cdiff_health_abundance,Strongest_inhibitor_present,select_white_black_mixed);
        if regenerate_A == 0
            break;
        end
    end
    %%% generate WB, BW
    for j = 1 : length(exchange_rate)
        [X_WB{j},X_BW{j},X_WB_Cdiff{j},X_BW_Cdiff{j},X_W_before_cohouse{j},X_B_before_cohouse{j}] = Generate_Samples_BW_and_WB_with_exchange_rate...
            (A,r,time,FunctionType,h1,h2,Cdiff,max(Num_sample),exchange_rate(j),abundance_type,Cdiff_disease_abundance,Cdiff_health_abundance,Strongest_inhibitor_present,select_white_black_mixed);
    end
    %%% get OTU table
    OTU_table = {X_W X_WB{1:length(exchange_rate)} X_BW{1:length(exchange_rate)} X_B};
    OTU_table_Cdiff = {X_W_Cdiff X_WB_Cdiff{1:length(exchange_rate)} X_BW_Cdiff{1:length(exchange_rate)} X_B_Cdiff};
    for i = 1 : size(OTU_table,2)
        relative_OTU_table{i} = bsxfun(@times,OTU_table{i},1./sum(OTU_table{i},1));
        relative_OTU_table_Cdiff{i} = bsxfun(@times,OTU_table_Cdiff{i},1./sum(OTU_table_Cdiff{i},1));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i = 1 : 2+2*length(exchange_rate)
        relative_abundance_Cdiff(:,i) = relative_OTU_table_Cdiff{i}(Cdiff,:)';
    end
    [a,b] = sort(mean(relative_abundance_Cdiff,1),'ascend');
    sort_group_name = group_name(b);
    relative_abundance_Cdiff = relative_abundance_Cdiff(:,b);
    
    for i = 1 : size(relative_abundance_Cdiff,2) - 1
        [h,p(i)] = ttest2(relative_abundance_Cdiff(:,i),relative_abundance_Cdiff(:,i+1));
    end
    
    if sum(strcmp(sort_group_name,perfect_group_name))==length(sort_group_name)&sum(p<0.05)>=size(relative_abundance_Cdiff,2) - 1 - 2
        color = parula(length(sort_group_name)); % Generate color values
        mean_relative_abundance_Cdiff = mean(relative_abundance_Cdiff,1);
        figure;hold on;box on
        b = bar([1:1:length(mean_relative_abundance_Cdiff)],mean_relative_abundance_Cdiff,0.6);
        b.FaceColor = 'flat';
        b.CData = color;
        
        errorbar(mean(relative_abundance_Cdiff,1),std(relative_abundance_Cdiff,1,1)/sqrt(Num_sample),'.r')
        set(gca,'Xlim',[0.5 length(mean_relative_abundance_Cdiff)+0.5],'XTick',[1:1:length(mean_relative_abundance_Cdiff)],'XTickLabelRotation',45,'XTickLabel',sort_group_name,'fontsize',10);
        ylim([0 max(ylim)])
        ylabel('Relative abundance of C. difficile')
        set(gca,'fontsize',14)
        set(gcf,'position',[30 19 535 437])
        set(gca,'yscale','log')
        %p-value<0.05 (*), <0.01(**), <0.001(***); >0.05(NS)
        for i = 1 : size(relative_abundance_Cdiff,2) - 1
            line([i i+1],[mean_relative_abundance_Cdiff(i+1)*1.1 mean_relative_abundance_Cdiff(i+1)*1.1],'Color','k')
            if p(i)>0.05
                significance = 'NS';
            elseif p(i)<=0.05 & p(i)>0.01
                significance = '*';
            elseif p(i)<=0.01 & p(i)>0.001
                significance = '**';
            elseif p(i)<=0.001
                significance = '***';
            end
            text(mean([i i+1]),mean_relative_abundance_Cdiff(i+1)*1.1,significance,'HorizontalAlignment','center','fontsize',14)
        end
        break;
    end
end

% figure;hold on
% [~,severity]=sort(abundance_Cdiff,'ascend');
% bar([1:1:size(abundance_Cdiff,2)],abundance_Cdiff(severity),0.6)
% errorbar(abundance_Cdiff(severity),abundance_error_Cdiff(severity),'.r')
% set(gca,'Xlim',[0.5 size(abundance_Cdiff,2)+0.5],'XTick',[1:1:size(abundance_Cdiff,2)],'XTickLabel',group_name(severity),'fontsize',14);
% ylim([0 max(ylim)])
% ylabel('abundance of C diff')
% set(gcf,'position',[680 608 558 370])
% set(gcf,'position',[680 608 558 370])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% generate text file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% for LefSe
% folder = ['./' father_folder '/No_Cdiff_exchange_'   exchange_name abundance_type '/N' num2str(N) 'C' num2str(C) 'S' num2str(simulation)];
% Generate_file_for_LefSe(folder,OTU_table,N,group_name,Num_sample);

folder = ['./' father_folder '/absolute' '_With_Cdiff_exchange_' exchange_name '/N' num2str(N) 'C' num2str(C) 'S' num2str(simulation)];
Generate_file_for_LefSe(folder,OTU_table_Cdiff,N,group_name,Num_sample);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% generate text file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% for LefSe

% folder = ['./' father_folder '/Convert_absolute_to_relative/No_Cdiff_exchange_'   exchange_name abundance_type '/N' num2str(N) 'C' num2str(C) 'S' num2str(simulation)];
% Generate_file_for_LefSe(folder,relative_OTU_table,N,group_name,Num_sample);

folder = ['./' father_folder '/relative' '_With_Cdiff_exchange_' exchange_name '/N' num2str(N) 'C' num2str(C) 'S' num2str(simulation)];
Generate_file_for_LefSe(folder,relative_OTU_table_Cdiff,N,group_name,Num_sample);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function Generate_file_for_LefSe(folder,OTU_table,N,group_name,Num_sample)

all_pair = nchoosek(1:1:size(OTU_table,2),2);
taxa_name = num2cell([1:1:N]);
for j = 1 : length(Num_sample)
    temp_folder = [folder 'sample_' num2str(Num_sample(j))];
    if exist(temp_folder)==0
        mkdir(temp_folder);
    end
    for k = 1 : size(OTU_table,2)
        temp_OTU_table{k} = OTU_table{k}(:,1:Num_sample(j));
    end
    for i = 1 : size(all_pair,1)
        Temp_table = cell2mat(temp_OTU_table(:,all_pair(i,:)));
        subject_ID = [repmat(group_name(all_pair(i,1)),1,Num_sample(j)) repmat(group_name(all_pair(i,2)),1,Num_sample(j))];
        filename = [temp_folder '/For_Lefse_' num2str(i) '_' group_name{all_pair(i,1)} '_and_' group_name{all_pair(i,2)} '.txt'];
        Write_file_for_LEFSe(filename,Temp_table,[repmat(group_name(all_pair(i,1)),1,Num_sample(j)) repmat(group_name(all_pair(i,2)),1,Num_sample(j))],subject_ID,taxa_name);
    end
end
end


