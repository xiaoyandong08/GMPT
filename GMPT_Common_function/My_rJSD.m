function Y = My_rJSD(X)

% KLD1 = @(x,y) myprod(x,log(x./y));
% KLD = @(x,y) sum(bsxfun(KLD1,repmat(x,size(y,1),1),y),2);
KLD = @(x,y) nansum(x.*log(x./y),2);
% KLD = @(x,y) myprod(x,log(x./y));
rJSD=@(x,y)(sqrt(0.5 * KLD(repmat(x,size(y,1),1), (repmat(x,size(y,1),1)+y)/2) + 0.5 * KLD(y, (repmat(x,size(y,1),1)+y)/2)));

[n,p] = size(X);
Y = zeros(1,n*(n-1)./2);
k = 1;
for i = 1:n-1
    XI = X(i,:);
    XJ = X((i+1):n,:);
    XI(XI==0) = nan;
    XJ(XJ(:)==0) = nan;
    Y(k:(k+n-i-1)) = feval(rJSD,XI,XJ)';

    %Y(k:(k+n-i-1)) = feval(rJSD,X(i,:),X((i+1):n,:))';
    k = k + (n-i);
end

end

function z = myprod(x,y)
for i = 1 : size(x,1)

    temp1 = find(x(i,:)~=0);
    temp2 = find(y(i,:)~=0);
    index = intersect(temp1,temp2);
    XI = x(index);
    XJ = y(index);
    z(i) = x(i,find(x(i,:)~=0)&find(y(i,:)>0))*y(i,find(x(i,:)>0)&find(y(i,:)>0))';

end


% z = x.*y;
% z(isnan(z)) = 0;
% z(~isfinite(z)) = 0;
end