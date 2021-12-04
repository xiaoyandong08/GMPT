function togml(A,r,filename)
N=size(A,1);
AA = zeros(N,N);
[a]=find(A~=0);
AA(a) = 1;AA([1:N+1:N*N])=0;

% [D] = all_shortest_paths(sparse(AA));
% dis = D(:,1);
fib=fopen(filename,'wt');

fprintf(fib,'%s\n','graph');
fprintf(fib,'%s\n','[');
fprintf(fib,'%s\n','directed 1');

for i=1:N
    fprintf(fib,'%s\n','node');
    fprintf(fib,'%s\n','[');
    fprintf(fib,'%s','id ');
    fprintf(fib,'%d\n',i);
    fprintf(fib,'%s','size ');
    fprintf(fib,'%f\n',r(i));
%     fprintf(fib,'%s','dis ');
%     if isfinite(dis(i))
%         fprintf(fib,'%d\n',dis(i));
%     else
%         fprintf(fib,'%d\n',0);
%     end
    fprintf(fib,'%s\n',']');
end
for i=1:N
    for j=1:N
        if(A(i,j)>0) && (i~=j)
            fprintf(fib,'%s\n','edge');
            fprintf(fib,'%s\n','[');
            fprintf(fib,'%s','source ');
            fprintf(fib,'%d\n',i);
            fprintf(fib,'%s','target ');
            fprintf(fib,'%d\n',j);
            fprintf(fib,'%s','Width ');
            fprintf(fib,'%f\n',abs(A(i,j)));
            fprintf(fib,'%s','label ');
            fprintf(fib,'%s\n','"1"');
            fprintf(fib,'%s\n',']');
        end
        if(A(i,j)<0) && (i~=j)
            fprintf(fib,'%s\n','edge');
            fprintf(fib,'%s\n','[');
            fprintf(fib,'%s','source ');
            fprintf(fib,'%d\n',i);
            fprintf(fib,'%s','target ');
            fprintf(fib,'%d\n',j);
            fprintf(fib,'%s','Width ');
            fprintf(fib,'%f\n',abs(A(i,j)));
            fprintf(fib,'%s','label ');
            fprintf(fib,'%s\n','"-1"');
            fprintf(fib,'%s\n',']');
        end
%         if (i==j)
%             fprintf(fib,'%s\n','edge');
%             fprintf(fib,'%s\n','[');
%             fprintf(fib,'%s','source ');
%             fprintf(fib,'%d\n',i);
%             fprintf(fib,'%s','target ');
%             fprintf(fib,'%d\n',j);
%             fprintf(fib,'%s','Width ');
%             fprintf(fib,'%f\n',100*abs(A(i,j)));
%             fprintf(fib,'%s','label ');
%             fprintf(fib,'%s\n','"3"');
%             fprintf(fib,'%s\n',']');
%         end
    end
end
fprintf(fib,'%s',']');
fclose(fib);