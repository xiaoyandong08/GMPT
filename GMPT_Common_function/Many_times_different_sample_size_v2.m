function Many_times_different_sample_size_v2(dynamical_parameter,N,C,delta,diag,VarianceType,exchange_rate,Num_sample,simulation,noise)
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
    if strcmp(dynamical_parameter,'connectivity')
        father_folder = ['Non_diagStable_' dynamical_parameter '/' 'N' num2str(N) '_C' num2str(C) '_Pheno' num2str(length(exchange_rate)*2+2) '_' folder_subname];
    elseif strcmp(dynamical_parameter,'delta')
        father_folder = ['Non_diagStable_' dynamical_parameter '/' 'N' num2str(N) '_delta' num2str(delta) '_Pheno' num2str(length(exchange_rate)*2+2) '_' folder_subname];
    elseif strcmp(dynamical_parameter,'noise')
        father_folder = ['Non_diagStable_' dynamical_parameter '/' 'N' num2str(N) '_noise' num2str(noise) '_Pheno' num2str(length(exchange_rate)*2+2) '_' folder_subname];
    end
end
for i = 1 : length(exchange_rate)
    WB_name{i} = ['WB-' num2str(exchange_rate(i))]; 
    BW_name{i} = ['BW-' num2str(exchange_rate(i))];
end
group_name = {'W',WB_name{1:end},BW_name{1:end},'B'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1 : simulation
    i
    [A{i},r{i},OTU_table{i},OTU_table_Cdiff{i},relative_OTU_table{i},relative_OTU_table_Cdiff{i},abundance_Cdiff{i},abundance_error_Cdiff{i}]=...
        Microbe_Phenotype_Triangulation_v4(N,C,delta,diag,VarianceType,FunctionType,Num_sample,i,exchange_rate,abundance_type,father_folder,...
        Cdiff_disease_abundance,Cdiff_health_abundance,Strongest_inhibitor_present,select_white_black_mixed,noise);
end
save(['./' father_folder '/Many_times_' abundance_type '_and_exchange_rate_N' num2str(N) '_' folder_subname '.mat'])

end


function [A,r,OTU_table,OTU_table_Cdiff,relative_OTU_table,relative_OTU_table_Cdiff,abundance_Cdiff,abundance_error_Cdiff]=Microbe_Phenotype_Triangulation_v4...
    (N,C,delta,diag,VarianceType,FunctionType,Num_sample,simulation,exchange_rate,abundance_type,father_folder,...
    Cdiff_disease_abundance,Cdiff_health_abundance,Strongest_inhibitor_present,select_white_black_mixed,noise)

h1 = 0.5;
h2 = 0.1;
time = [0:0.01:50];
Cdiff = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while 1
    %%% generate ecological network
    [A,r] = Generate_Network_A_of_BandW(N,C,delta,diag,VarianceType,time,FunctionType,h1,h2,Cdiff,abundance_type,Cdiff_disease_abundance,Cdiff_health_abundance,select_white_black_mixed);
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
    OTU_table{i} = OTU_table{i} + OTU_table{i}.*(2*noise*rand(N,max(Num_sample))-noise);
    OTU_table_Cdiff{i} = OTU_table_Cdiff{i} + OTU_table_Cdiff{i}.*(2*noise*rand(N,max(Num_sample))-noise);
end

for i = 1 : size(OTU_table,2)
    relative_OTU_table{i} = bsxfun(@times,OTU_table{i},1./sum(OTU_table{i},1));
    relative_OTU_table_Cdiff{i} = bsxfun(@times,OTU_table_Cdiff{i},1./sum(OTU_table_Cdiff{i},1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1 : 2+2*length(exchange_rate)
    abundance_Cdiff(i) = mean(OTU_table_Cdiff{i}(Cdiff,:));
    abundance_error_Cdiff(i) = std(OTU_table_Cdiff{i}(Cdiff,:),1);
end
exchange_name = [];
for i = 1 : length(exchange_rate)
    WB_name{i} = ['WB-' num2str(exchange_rate(i))]; 
    BW_name{i} = ['BW-' num2str(exchange_rate(i))];
    exchange_name = [exchange_name num2str(exchange_rate(i)) '_'];
end
group_name = {'W',WB_name{1:end},BW_name{1:end},'B'};


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

folder = ['./' father_folder '/absolute' '_With_Cdiff_exchange_' exchange_name '/N' num2str(N) 'S' num2str(simulation)];
Generate_file_for_LefSe(folder,OTU_table_Cdiff,N,group_name,Num_sample);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% generate text file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% for LefSe

% folder = ['./' father_folder '/Convert_absolute_to_relative/No_Cdiff_exchange_'   exchange_name abundance_type '/N' num2str(N) 'C' num2str(C) 'S' num2str(simulation)];
% Generate_file_for_LefSe(folder,relative_OTU_table,N,group_name,Num_sample);

folder = ['./' father_folder '/relative' '_With_Cdiff_exchange_' exchange_name '/N' num2str(N) 'S' num2str(simulation)];
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


