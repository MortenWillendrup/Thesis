
b=0.01;
a=0.5;
sig=0.01;
r0=0.02;
B=@(a,t,T)1/a*(1-exp(-a*(T-t)));
A=@(a,b,sig,t,T)((B(a,t,T)-(T-t))*(a*b-0.5*sig^2)/(a^2)-(sig^2*(B(a,b,T)).^2)./(4*a));
P=@(a,b,r,t,T)exp(A(a,b,sig,t,T)-B(a,t,T)*r);
dBdt=@(a,t,T)(B(a,t+0.0001,T)-B(a,t,T))/0.0001;
dAdt=@(a,b,sig,t,T)(A(a,b,sig,t+0.0001,T)-A(a,b,sig,t,T))/0.0001;
dPdt=@(a,b,r,t,T,sig)(dAdt(a,b,sig,t,T)-dBdt(a,t,T)*r)*P(a,b,r,t,T);
dPdr=@(a,b,r,t,T,sig)-B(a,t,T)*P(a,b,r,t,T);
dP2dr2=@(a,b,r,t,T,sig)(B(a,t,T))^2*P(a,b,r,t,T);


T=10;
N=1000;
dt=T/N;
r=nan(N,1);
bond=r;
B=r;
V=r;
r(1)=r0;
V(1)=100;
dW=randn(N,1)*sqrt(dt);
t=(0:dt:T)';
mat=1:T;
cashFlows=(t==round(t)&t>0)*4;
cashFlows(end)=cashFlows(end)+100;
for i=1:N
    if i>1
        r(i)=r(i-1)+(b-a*r(i-1))*dt+sig*dW(i);
        if i==2
            B(i-1)=bond(i-1);
        end
        alpha=(dPdt(a,b,r(i),t(i),T,sig)+(b-a*r(i))*dPdr(a,b,r(i),t(i),T,sig)+0.5*sig^2*dP2dr2(a,b,r(i),t(i),T,sig))/P(a,b,r(i),t(i),T);
        beta=(sig*dPdr(a,b,r(i),t(i),T,sig))/P(a,b,r(i),t(i),T);
        B(i)=B(i-1)+alpha*B(i-1)*dt+beta*B(i-1)*dW(i)-cashFlows(i);
    end
    disk=P(a,b,r(i),t(i),(t(i)+dt:dt:T)');
    bond(i)=disk'*cashFlows(t>t(i));
    
    if i>1
        V(i)=bond(i)/bond(i-1)*V(i-1)+cashFlows(i);
    end
end
plot(bond)
hold on
plot(B)