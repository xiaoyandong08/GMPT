function [Nodes,r]=generate_SF_A(N,m,gama,delta,diag,VarianceType)
% m is sparsity
epsilon = 0;
if VarianceType == 1
    zigma = 1/sqrt(N*(2+delta));
elseif VarianceType == 2
    zigma = delta;
end

Nodes=zeros(N,N);
pd=zeros(N,1);
pp=zeros(N,1);
alfa=1/(gama-1);
edge_num=N*(N-1)*m;

%1.get nodes' weight
for i=1:N
    pd(i)=1/(i^alfa);
end
%2.normalized node weights
sum_p=sum(pd);
pp(1)=pd(1)/sum_p;
for i=2:N
    pp(i)=pp(i-1)+pd(i)/sum_p;
end

%3.add edge
for i=1:edge_num
    %i
    ADD_ONE_EDGE=0;
    while ADD_ONE_EDGE==0
        node_1_Len=find(pp>rand(1));
        node_1=node_1_Len(1);
        node_2_Len=find(pp>rand(1));
        node_2=node_2_Len(1);
        while node_2==node_1%avoid self-loop
            node_2_Len=find(pp>rand(1));
            node_2=node_2_Len(1);
        end
        if Nodes(node_1,node_2)==0
            flag = 1;
            while flag == 1
                Nodes(node_1,node_2) = zigma * randn(1);
                if abs(Nodes(node_1,node_2))>epsilon
                    flag = 0;
                end
            end
            ADD_ONE_EDGE=1;
        else
            ADD_ONE_EDGE=0;
        end%Nodes(node_1,node_2)
    end%while ADD_ONE_EDGE==0
end% for i=1:edge_num
Nodes(1:N+1:N^2)=diag;
r = rand(1,N);
end