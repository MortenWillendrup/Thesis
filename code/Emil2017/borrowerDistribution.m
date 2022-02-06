function [pdf,cdf]=borrowerDistribution(x,Xstar,lambda,A,B,n)
% If Xstar is empty then we are just in the beta distribution
if isempty(Xstar)
    pdf=betapdf(x,A,B);
    cdf=betacdf(x,A,B);
else
    % preallocation
    pdf=nan(size(x));
    cdf=nan(size(x));
    
    % Calculate prepayment probabilities
    prepayProbs=nan(size(Xstar));
    for i=1:size(Xstar,1)
        prepayProbs(i)=prepaymentProbability(Xstar(i),Xstar(1:i-1),lambda,A,B,n);
    end

    % Prepayment probability
    P=1-exp(-1/n*lambda);
        
    % Loop over x
    for j=1:max(size(x))

        % CDF
        cdf(j)=1/(1-sum(prepayProbs))*1/P*prepaymentProbability(x(j),Xstar,lambda,A,B,n);

        % PDF
        N=size(Xstar,1)-sum(Xstar<=x(j));
        pdf(j)=1/(1-sum(prepayProbs))*betapdf(x(j),A,B)*(1-P)^N;
    end
end
end