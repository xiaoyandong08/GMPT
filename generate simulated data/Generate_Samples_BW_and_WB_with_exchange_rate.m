function [X_WB,X_BW,X_WB_Cdiff,X_BW_Cdiff,X_W_before_cohouse,X_B_before_cohouse] = Generate_Samples_BW_and_WB_with_exchange_rate(A,r,time,FunctionType,h1,h2,Cdiff,Num_sample,change,abundance_type,Cdiff_disease_abundance,Cdiff_health_abundance,Strongest_inhibitor_present,select_white_black_mixed)
N = size(A,1);
for i = 1 : Num_sample
    while 1
        
        if strcmp(select_white_black_mixed,'mixed')
            %         local_W_index = unique([Cdiff;randsample([1:1:N],randi(N))']);
            %         local_B_index = unique([Cdiff;randsample([1:1:N],randi(N))']);
            promoters  = find(A(:,Cdiff)>0);
            inhibitors = find(A(:,Cdiff)<0);
            neutrals    = find(A(:,Cdiff)==0);
            if isempty(neutrals)
                neutrals = Cdiff;
            end
            local_W_index = unique([Cdiff;randsample(inhibitors,randi(length(inhibitors))); randsample(promoters,randi(ceil(0.1*length(promoters)))); randsample(neutrals,randi(length(neutrals)))]);
            local_B_index = unique([Cdiff;randsample(promoters,randi(length(promoters))); randsample(inhibitors,randi(ceil(0.1*length(inhibitors)))); randsample(neutrals,randi(length(neutrals)))]);
        else
            local_W_index = unique([Cdiff;randsample(find(A(:,Cdiff)<=0),randi(length(find(A(:,Cdiff)<=0))))]);
            local_B_index = unique([Cdiff;randsample(find(A(:,Cdiff)>=0),randi(length(find(A(:,Cdiff)>=0))))]);
        end
        
        if Strongest_inhibitor_present == 1
            Strongest_inhibitor = find(A(:,Cdiff) ==  min(A(setdiff(1:N,Cdiff),Cdiff)));
            local_W_index = unique([local_W_index;Strongest_inhibitor]);
        end
        
        initial = zeros(N,1);
        initial(local_W_index) = 0.2;
        [dx_W,dxx_W]=glv_Euler_type(initial,A,r,time,FunctionType,h1,h2,abundance_type);
        
        initial = zeros(N,1);
        initial(local_B_index) = 0.2;
        [dx_B,dxx_B]=glv_Euler_type(initial,A,r,time,FunctionType,h1,h2,abundance_type);
        
        if dxx_W(Cdiff)<Cdiff_health_abundance && dxx_B(Cdiff)>Cdiff_disease_abundance
            
            initial = zeros(N,1);
            initial(setdiff(local_W_index,Cdiff)) = 0.2;
            [~,X_W_before_cohouse(:,i)]=glv_Euler_type(initial,A,r,time,FunctionType,h1,h2,abundance_type);
            
            initial = zeros(N,1);
            initial(setdiff(local_B_index,Cdiff)) = 0.2;
            [~,X_B_before_cohouse(:,i)]=glv_Euler_type(initial,A,r,time,FunctionType,h1,h2,abundance_type);
            
            WB_intro = unique(randsample(local_B_index,ceil(change*length(local_B_index))));
            ini_WB = zeros(N,1);
            ini_WB(setdiff([local_W_index;WB_intro],Cdiff)) = 0.2;
            [~,X_WB(:,i)]=glv_Euler_type(ini_WB,A,r,time,FunctionType,h1,h2,abundance_type);
            ini_WB(Cdiff) = 0.2;
            [~,X_WB_Cdiff(:,i)]=glv_Euler_type(ini_WB,A,r,time,FunctionType,h1,h2,abundance_type);
            
            BW_intro = unique(randsample(local_W_index,ceil(change*length(local_W_index))));
            ini_BW = zeros(N,1);
            ini_BW(setdiff([local_B_index;BW_intro],Cdiff)) = 0.2;
            [~,X_BW(:,i)]=glv_Euler_type(ini_BW,A,r,time,FunctionType,h1,h2,abundance_type);
            ini_BW(Cdiff) = 0.2;
            [~,X_BW_Cdiff(:,i)]=glv_Euler_type(ini_BW,A,r,time,FunctionType,h1,h2,abundance_type);
            
            if sum(isnan(X_WB(:,i)))==0 &&sum(isnan(X_BW(:,i)))==0&&sum(isnan(X_WB_Cdiff(:,i)))==0 &&sum(isnan(X_BW_Cdiff(:,i)))==0
                break;
            end
        end
    end
    
end
end