function Plot_cohousing_process
addpath('../GMPT_Common_function')

N = 30;
diff_C = 0.05; 
delta = 0.2;
diag = -1;
VarianceType = 2; % 1: Travis; 2: Variance
exchange_rate = [0.1 0.3 0.5]; % phenotype = 2 * length(exchange_rate) + 2
Num_sample = 20;
simulation = 10;

% for i = 1 : length(diff_C)
%     Many_times_different_sample_size(N,diff_C(i),delta,diag,VarianceType,exchange_rate,Num_sample,simulation)
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Diag_Type = 'Non-diagStable';
folder_time = '20210907T175417';
analysis_type    = 'relative';
Cdiff = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

simulation_index = 2

A = A{simulation_index};
r = r{simulation_index};

simulation_relative_OTU_table_Cdiff = relative_OTU_table_Cdiff{simulation_index};
simulation_absolute_OTU_table_Cdiff = OTU_table_Cdiff{simulation_index};

% Plot_Cdiff_abundance(Num_sample,Cdiff,simulation_index,group_name,simulation_absolute_OTU_table_Cdiff,simulation_relative_OTU_table_Cdiff)

togml(A,r,'A30_20210907.gml')
%%%%%%%%%%%%%%%% W and B
time = [0:0.01:50];
h1 = 0.5;
h2 = 0.1;
FunctionType = 1;
abundance_type = 'absolute';

%% generate B, W
         
sample_index = 3;      
figure('position',[238 365 610 202]);%[206 491 897 308]
count=1;
for i = [1 8];
    ini = zeros(N,1);
    ini(find(simulation_relative_OTU_table_Cdiff{i}(:,sample_index)~=0)) = 0.2;
    species_index = find(ini~=0);
    [X,X_ini]=glv_Euler_type(ini,A,r,time,FunctionType,h1,h2,abundance_type);
    X_ini(Cdiff) = 0.1;
    [X_Cdiff,~]=glv_Euler_type(X_ini,A,r,time,FunctionType,h1,h2,abundance_type);
    all_data = [X';X_Cdiff'];
    
    subplot(1,2,count);hold on;
    %all_data(:,Cdiff) = [];
    plot([time time(end)+time],all_data)
    plot(time(end)+time,X_Cdiff(Cdiff,:),'r','LineWidth',3)
    scatter(time(end)+time(1),X_Cdiff(Cdiff,1),60,'ro','filled')
    xlabel('time');
    ylabel('Abundance of Species X')
    set(gca,'fontsize',14)
    title(['Phenotype-' group_name{i}],'fontsize',14)
    count=count+1;
    
    y_lim = get(gca,'Ylim');
    v = [[0 time(end) time(end) 0]' [y_lim(1) y_lim(1) y_lim(2) y_lim(2)]'];
    f = [1 2 3 4];
    patch('Faces',f,'Vertices',v,'FaceColor',hex2rgb('00AEEF'),'EdgeColor','none','FaceAlpha',.25)
    v = [[time(end) 2*time(end) 2*time(end) time(end)]' [y_lim(1) y_lim(1) y_lim(2) y_lim(2)]'];
    f = [1 2 3 4];
    patch('Faces',f,'Vertices',v,'FaceColor',hex2rgb('F74461'),'EdgeColor','none','FaceAlpha',.25)
end

local_W_index = find(simulation_relative_OTU_table_Cdiff{1}(:,sample_index)~=0);
local_B_index = find(simulation_relative_OTU_table_Cdiff{8}(:,sample_index)~=0);

length(local_W_index)
length(local_B_index)
%%%%%%%%%%%%%%%%WB-0.1 WB 0.5 WB-1

mycolor = jet(N);
exchange_rate = [0.5];
for j = 1 : length(exchange_rate)
    
    initial = zeros(N,1);
    initial(local_W_index) = 0.2;
    [XX_W_before_cohouse,X_W_before_cohouse]=glv_Euler_type(initial,A,r,time,FunctionType,h1,h2,abundance_type);
    
    WB_intro = unique(randsample(local_B_index,ceil(exchange_rate(j)*length(local_B_index))));
    WB_intro = [1 2 4 15 19 20 22 23]
    X_W_before_cohouse(WB_intro) = X_W_before_cohouse(WB_intro)+0.1;
    [XX_WB_Cohouse,X_WB_Cohouse]=glv_Euler_type(X_W_before_cohouse,A,r,time,FunctionType,h1,h2,abundance_type);
    
    X_WB_Cohouse(Cdiff) = 0.2;
    [XX_WB_Cdiff,~]=glv_Euler_type(X_WB_Cohouse,A,r,time,FunctionType,h1,h2,abundance_type);
    all_data = [XX_W_before_cohouse';XX_WB_Cohouse';XX_WB_Cdiff'];

    %subplot(1,3,j);hold on;
    figure;hold on;
    h=plot([time time(end)+time time(end)*2+time],all_data(:,local_W_index));
    set(h,{'color'},num2cell(mycolor(local_W_index,:),2))
    if length(WB_intro)>0
        h=plot([time(end)+time time(end)*2+time],all_data(1+size(all_data,1)/3:size(all_data,1),WB_intro));
        set(h,{'color'},num2cell(mycolor(WB_intro,:),2))

    end
    plot(time(end)*2+time,XX_WB_Cdiff(Cdiff,:),'r','LineWidth',3)
    scatter(time(end)*2+time(1),XX_WB_Cdiff(Cdiff,1),60,'ro','filled')
    xlabel('time');
    ylabel('Abundance of Species X')
    set(gca,'fontsize',14)
    title(['WB-' num2str(exchange_rate(j))],'fontsize',14)
    ylim([0 1.2])
    
    y_lim = get(gca,'Ylim');
    v = [[0 time(end) time(end) 0]' [y_lim(1) y_lim(1) y_lim(2) y_lim(2)]'];
    f = [1 2 3 4];
    patch('Faces',f,'Vertices',v,'FaceColor',hex2rgb('00AEEF'),'EdgeColor','none','FaceAlpha',.25)
    v = [[time(end) 2*time(end) 2*time(end) time(end)]' [y_lim(1) y_lim(1) y_lim(2) y_lim(2)]'];
    f = [1 2 3 4];
    patch('Faces',f,'Vertices',v,'FaceColor',hex2rgb('2BB673'),'EdgeColor','none','FaceAlpha',.25)
    v = [[2*time(end) 3*time(end) 3*time(end) 2*time(end)]' [y_lim(1) y_lim(1) y_lim(2) y_lim(2)]'];
    f = [1 2 3 4];
    patch('Faces',f,'Vertices',v,'FaceColor',hex2rgb('F74461'),'EdgeColor','none','FaceAlpha',.25)
    set(gcf,'position',[616 801 461 186])
end






end