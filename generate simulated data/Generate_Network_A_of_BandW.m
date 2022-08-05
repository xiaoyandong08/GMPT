function [A,r] = Generate_Network_A_of_BandW(N,C,delta,diag,VarianceType,time,FunctionType,h1,h2,Cdiff,abundance_type,Cdiff_disease_abundance,Cdiff_health_abundance,select_white_black_mixed,shuffle_SF_Amatrix)
gama = 2.5;
while 1
    while 1
%         [A,r] = generate_A(N,C,delta,diag,VarianceType);
        [A,r]=generate_SF_A(N,C,gama,delta,diag,VarianceType);
        if shuffle_SF_Amatrix == 1
            index = randperm(N);
            A = A(index,index);
            r = r(index);
        end
        promoters  = find(A(:,Cdiff)>0);
        inhibitors = find(A(:,Cdiff)<0);
        neutrals    = find(A(:,Cdiff)==0);
        if length(promoters) >0 & length(inhibitors)>0 & length(neutrals)>0 & sum(sum(A~=0,2)'+sum(A~=0,1)==2)==0
            break;
        end
    end
    if strcmp(select_white_black_mixed,'mixed')
%         local_W_index = unique([Cdiff;randsample([1:1:N],randi(N))']);
%         local_B_index = unique([Cdiff;randsample([1:1:N],randi(N))']);
        
        if isempty(neutrals)
            neutrals = Cdiff;
        end
        local_W_index = unique([Cdiff;randsample(inhibitors,randi(length(inhibitors))); randsample(promoters,randi(ceil(0.1*length(promoters)))); randsample(neutrals,randi(length(neutrals)))]);     
        local_B_index = unique([Cdiff;randsample(promoters,randi(length(promoters))); randsample(inhibitors,randi(ceil(0.1*length(inhibitors)))); randsample(neutrals,randi(length(neutrals)))]);
    else
        local_W_index = unique([Cdiff;randsample(find(A(:,Cdiff)<=0),randi(length(find(A(:,Cdiff)<=0))))]);
        local_B_index = unique([Cdiff;randsample(find(A(:,Cdiff)>=0),randi(length(find(A(:,Cdiff)>=0))))]);
    end
    initial = zeros(N,1);
    initial(local_W_index) = 0.2;
    [dx_W,dxx_W]=glv_Euler_type(initial,A,r,time,FunctionType,h1,h2,abundance_type);
    partial1 = sum(A(find(dxx_W~=0),Cdiff)>0)/sum(A(find(dxx_W~=0),Cdiff)<0);
    
    initial = zeros(N,1);
    initial(local_B_index) = 0.2;
    [dx_B,dxx_B]=glv_Euler_type(initial,A,r,time,FunctionType,h1,h2,abundance_type);
    partial2 = sum(A(find(dxx_B~=0),Cdiff)>0)/sum(A(find(dxx_B~=0),Cdiff)<0);

    if dxx_W(Cdiff)<Cdiff_health_abundance && dxx_B(Cdiff)>Cdiff_disease_abundance
        aa = 1;
    end
    if dxx_W(Cdiff)<Cdiff_health_abundance && dxx_B(Cdiff)>Cdiff_disease_abundance && partial1 >=0.5 && partial1 <= 2 && partial2 >=0.5 && partial2 <= 2
        aa = 1;
    end
    if dxx_W(Cdiff)<Cdiff_health_abundance && dxx_B(Cdiff)>Cdiff_disease_abundance && partial1 >=0.5 && partial1 <= 2 && partial2 >=0.5 && partial2 <= 2 && sum(A(find(dxx_W~=0),Cdiff)>0)>=2 && sum(A(find(dxx_B~=0),Cdiff)>0)>=2
    %if dxx_W(Cdiff)<Cdiff_health_abundance && dxx_B(Cdiff)>Cdiff_disease_abundance
        break;
    end
end
    
end