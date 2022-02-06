function [yy,forward,fwdDiff,derivative]=hermiteInterpolation(x,y,xx)
% Ensure column input
x=x(:);
y=y(:);
xx=xx(:);

% Sort data
[x,IX]=sort(x);
y=y(IX);

% Preallocate
n=size(x,1);

% Define h
h=x(2:n)-x(1:n-1);

% a=y
a=y;

% m
m=(a(2:n)-a(1:n-1))./h;

% Preallocate b,c
b=nan(n,1);
c=nan(n,1);
d=nan(n,1);

% b
b(1)=1/(x(3)-x(1))*((x(3)+x(2)-2*x(1))*(y(2)-y(1))/(x(2)-x(1))...
        -(x(2)-x(1))*(y(3)-y(2))/(x(3)-x(2)));
% for i=2:n-1
%     b(i)=1/(x(i+1)-x(i-1))*((x(i+1)-x(i))*(y(i)-y(i-1))/(x(i)-x(i-1))...
%         +(x(i)-x(i-1))*(y(i+1)-y(i))/(x(i+1)-x(i)));
% end
b(2:end-1)=1./(x(3:end)-x(1:end-2)).*((x(3:end)-x(2:end-1)).*(y(2:end-1)-y(1:end-2))./(x(2:end-1)-x(1:end-2))...
    +(x(2:end-1)-x(1:end-2)).*(y(3:end)-y(2:end-1))./(x(3:end)-x(2:end-1)));

b(n)=-1/(x(n)-x(n-2))*((x(n)-x(n-1))*(y(n-1)-y(n-2))/(x(n-1)-x(n-2))...
    -(2*x(n)-x(n-1)-x(n-2))*(y(n)-y(n-1))/(x(n)-x(n-1)));

% Generate c and d
% for i=1:n-1
%     c(i)=(3*m(i)-b(i+1)-2*b(i))/h(i);
%     d(i)=(b(i+1)+b(i)-2*m(i))/h(i)^2;
% end
c(1:n-1)=(3*m(1:n-1)-b(2:n)-2*b(1:n-1))./h(1:n-1);
d(1:n-1)=(b(2:n)+b(1:n-1)-2*m(1:n-1))./h(1:n-1).^2;

% Generate output
n2=size(xx,1);
yy=nan(n2,1);

% If two outputs then compute forward
fwdTrue=false;
fwdDiffTrue=false;
derTrue=false;
if nargout>1
    fwdTrue=true;
    forward=nan(n2,1);
end
if nargout>2
    fwdDiffTrue=true;
    fwdDiff=nan(n2,1);
end
if nargout>3
    derTrue=true;
    derivative=nan(n2,1);
end

for j=1:n2
    % isnan? Then skip
    if isnan(xx(j))
        continue
    end
    
    % Find index. Will not work for extrapolation!
    i=sum(x<=xx(j));
    if xx(j)>=x(n)
        i=i-1;
    end
    
    % Flat extrapolation
    if i==0
        i=1;
        xx(j)=x(1);
    elseif xx(j)>x(n)
        xx(j)=x(n);
    end
    
    % Calculate value
    yy(j)=a(i)+b(i)*(xx(j)-x(i))+c(i)*(xx(j)-x(i))^2+d(i)*(xx(j)-x(i))^3;
    
    % Calculate Forward
    if fwdTrue
        forward(j)=a(i)+b(i)*(2*xx(j)-x(i))+c(i)*(xx(j)-x(i))*(3*xx(j)-x(i))...
                   +d(i)*(xx(j)-x(i))^2*(4*xx(j)-x(i));
    end
    
    % Calculate forward diff
    if fwdDiffTrue
        fwdDiff(j)=2*b(i)+2*c(i)*(3*xx(j)-2*x(i))+2*d(i)*(xx(j)-x(i))*(6*xx(j)-3*x(i));
    end
    
    % calculate derivative
    if derTrue
        derivative(j)=b(i)+2*c(i)*(xx(j)-x(i))+3*d(i)*(xx(j)-x(i))^2;
    end
end
end