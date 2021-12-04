function [X_W,X_B,X_W_Cdiff,X_B_Cdiff,regenerate_A]=Generate_Samples_W_and_B(A,r,time,FunctionType,h1,h2,Num_sample,Cdiff,abundance_type,Cdiff_disease_abundance,Cdiff_health_abundance,Strongest_inhibitor_present,select_white_black_mixed)
N = size(A,1);
count = 0;
regenerate_A = 0;
for i = 1 : Num_sample
    while 1
        if count == 1e3
            X_W= [];X_B=[];X_W_Cdiff=[];X_B_Cdiff=[];
            regenerate_A = 1;
            continue; 
        end
       
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
            
            X_W_Cdiff(:,i) = dxx_W;
            X_B_Cdiff(:,i) = dxx_B;
            
            initial = zeros(N,1);
            initial(setdiff(local_W_index,Cdiff)) = 0.2;
            [dx_W,dxx_W]=glv_Euler_type(initial,A,r,time,FunctionType,h1,h2,abundance_type);
            
            initial = zeros(N,1);
            initial(setdiff(local_B_index,Cdiff)) = 0.2;
            [dx_B,dxx_B]=glv_Euler_type(initial,A,r,time,FunctionType,h1,h2,abundance_type);
            
            X_W(:,i) = dxx_W;
            X_B(:,i) = dxx_B;
            if sum(isnan(X_W(:,i)))==0 &&sum(isnan(X_B(:,i)))==0 &&sum(isnan(X_W_Cdiff(:,i)))==0 &&sum(isnan(X_B_Cdiff(:,i)))==0
                count = 0;
                break;
            else 
                count = count + 1;
            end
        end
    end
end
end