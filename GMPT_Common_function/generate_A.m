function [A,r] = generate_A(N,C,delta,diag,VarianceType)
r = rand(1,N);
%r = 0.5*ones(1,N);

while 1
    A = zeros(N,N);
    if VarianceType == 1
        zigma = 1/sqrt(N*(2+delta));
    elseif VarianceType == 2
        zigma = delta;
    end
    for i = 1 : N
        for j = 2 : N
            if rand(1) < C
                A(i,j) =  ((zigma * randn(1)));
            end
        end
        A(i,i) = diag;
    end
    index1 = randsample(2:N,5);
    A(index1,1) = abs(zigma * randn(1,5));
    A(randsample(setdiff(2:N,index1),5),1) = -abs(zigma * randn(1,5));
    if sum(sum(A~=0,1)+sum(A~=0,2)'==2)==0
        break;
    end
end
end