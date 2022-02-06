function price=bondOptionHullWhite(Bs,Bt,K,t,T,S,kappa,sigma,type)
typeFlag='call';
if nargin>8
    if any(strcmpi({'call','c'},type))
        typeFlag='call';
    elseif any(strcmpi({'put','p'},type))
        typeFlag='put';
    else
        error('Type must be call or put.')
    end
end

% Calculate vol term
V=1/kappa*(1-exp(-kappa*(S-T)))*sqrt(sigma^2/(2*kappa)*(1-exp(-2*kappa*(T-t))));

% Calculate d1 and d2
d1=1/V*log(Bs/(K*Bt))+1/2*V;
d2=d1-V;

% Calc price
if strcmpi(typeFlag,'call')
    price=Bs*normcdf(d1)-Bt*K*normcdf(d2);
else
    price=Bt*K*normcdf(-d2)-Bs*normcdf(-d1);
end

end