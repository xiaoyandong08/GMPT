function Final_M0 = Contribution_Matrix(A,r,x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Final_M0(i,j) is the contribution of species i on species j
%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time = [0:0.1:100];

N = size(A,1);

[dx1,dx,time]=glv_ODE45(A,r,x,1);
% Extinction_index = find(dx1<1e-3);

for i = 1 : size(dx,1)
    [f,goodness] = fit(time(end-600:end)',log(dx(i,end-600:end))','poly1');
    slope(i) = f.p1;
end

count = 0;
[~,sort_index] = sort(abs(slope),'descend');
while 1
    
    part_A = A(setdiff([1:N],sort_index(1:count)),setdiff([1:N],sort_index(1:count)));
    part_r = r(setdiff([1:N],sort_index(1:count)));
    part_x = -inv(part_A')*part_r;
    x_start(count+1,1:N) = 0;
    x_start(count+1,setdiff([1:N],sort_index(1:count))) = part_x';
    count = count + 1;
    if isempty(find(part_x<0))
        break;
    end
    
end
Extinction_index = sort_index(1:count-1);


Remain_index = setdiff([1:N],Extinction_index);
new_A = A(Remain_index,Remain_index);
new_r = r(Remain_index);
for i = 1 : size(new_A,1)
    for j = 1 : size(new_A,2)
        temp = new_A;
        temp(i,:) = [];
        temp(:,j) = [];
        M0(j,i) = -(-1)^(i+j) * det(temp)*new_r(j) /det(new_A);
    end
end
Final_M0 = zeros(N,N);
Final_M0(Remain_index,Remain_index) = M0;
% if length(Extinction_index)>0
%     Final_M0(setdiff([1:size(A,1)],Extinction_index),Extinction_index) = -Inf;  
% end
% Final_M0(Extinction_index,:) = NaN;
% Final_M0(:,Extinction_index) = NaN;

end

function [X_set,XX,time]=glv_ODE45(A,r,x,Type)
N = size(A,1);

f = @(t,x) Lotka_Volterra_ode(x,A,r,Type); % define the system
options = odeset('AbsTol',1e-14,'RelTol',1e-12,'NonNegative',ones(1,N)); %options for the system
sol = ode45(f,[0 100],x,options);

X_set = sol.y(:,end); 
XX = sol.y;
time = sol.x;
end

function dx = Lotka_Volterra_ode(x,A,r,Type)
dx = zeros(size(A,1),1);

switch Type
    case 1  % Functional type 1
        dx = x.*(A'*x + r);
    case 2 % Functional type 2
        dx = x.*(A'*(x./(1+x)) + r);
    case 3 % Functional type 2
        dx = x.*(A'*(x.^2) + r);
end
end