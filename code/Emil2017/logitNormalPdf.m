function val=logitNormalPdf(x,mu,sigma)
x=max(x,0);
    val=1/(sigma*sqrt(2*pi)).*exp(-(logit(x)-mu).^2./(2*sigma^2)).*1./(x.*(1-x));
    %val=1/2*(1+erf((logit(x)-mu)./(sqrt(2*sigma^2))));
    
end

function val=logit(x)
    val=log(x./(1-x));
end