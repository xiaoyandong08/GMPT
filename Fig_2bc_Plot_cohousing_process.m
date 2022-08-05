function Fig_2bc_Plot_cohousing_process

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = A{simulation_index};
r = r{simulation_index};

simulation_relative_OTU_table_Cdiff = relative_OTU_table_Cdiff{simulation_index};
simulation_absolute_OTU_table_Cdiff = OTU_table_Cdiff{simulation_index};

% Plot_Cdiff_abundance(Num_sample,Cdiff,simulation_index,group_name,simulation_absolute_OTU_table_Cdiff,simulation_relative_OTU_table_Cdiff)

% togml(A,r,'A30_20220429T091254_simulation3.gml')
%%%%%%%%%%%%%%%% W and B
time = [0:0.01:50];
h1 = 0.5;
h2 = 0.1;
FunctionType = 1;
abundance_type = 'absolute';

%% generate B, W
% Cdiff_health_abundance = 1e-5;
% Num_sample = 5;
% [X_W,X_B,X_W_Cdiff,X_B_Cdiff,regenerate_A]=Generate_Samples_W_and_B(A,r,time,FunctionType,h1,h2,Num_sample,Cdiff,abundance_type,Cdiff_disease_abundance,Cdiff_health_abundance,Strongest_inhibitor_present,select_white_black_mixed);
%         
sample_index = 13; 

local_W_index = find(simulation_relative_OTU_table_Cdiff{1}(:,sample_index)~=0);
local_B_index = find(simulation_relative_OTU_table_Cdiff{8}(:,sample_index)~=0);

A(local_W_index,1)
A(local_B_index,1)
length(local_W_index);
length(local_B_index);

mycolor = jet(N);
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
    h = plot([time time(end)+time],all_data(:,species_index));
    set(h,{'color'},num2cell(mycolor(species_index,:),2))
    plot(time(end)+time,X_Cdiff(Cdiff,:),'r','LineWidth',3)
    scatter(time(end)+time(1),X_Cdiff(Cdiff,1),60,'ro','filled')
    xlabel('time');
    ylabel('Abundance of Species X')
    set(gca,'fontsize',14,'TickDir','out')
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


%%%%%%%%%%%%%%%%WB-0.1 WB 0.5 WB-1

% figure('position',[74 559 1191 303])
mycolor = jet(N);
exchange_rate = [0.5];
for j = 1 : length(exchange_rate)
    
    initial = zeros(N,1);
    initial(local_W_index) = 0.2;
    [XX_W_before_cohouse,X_W_before_cohouse]=glv_Euler_type(initial,A,r,time,FunctionType,h1,h2,abundance_type);
    
    WB_intro = unique(randsample(local_B_index,ceil(exchange_rate(j)*length(local_B_index))));
    WB_intro';
    WB_intro = [6    12    14    15    17    19    22];
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
    set(gca,'fontsize',14,'TickDir','out')
    title(['WB-' num2str(exchange_rate(j))],'fontsize',14)
    %ylim([0 1])
    
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
    set(gcf,'position',[616 801 461 202])
end


end