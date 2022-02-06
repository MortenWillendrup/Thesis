function [val,prepayProbs]=modelPrepayments(Alpha,Beta,lambda,Xstar,n)
% Calculate prepayment probabilities
prepayProbs=nan(size(Xstar));
for i=1:size(Xstar,1)
    prepayProbs(i)=prepaymentProbability(Xstar(i),Xstar(1:i-1),lambda,Alpha,Beta,n);
end

% % % If market prepayments are provided then use these
% % if nargin>5
% %     marketRemain=cumprod(1-[0;marketPrepayments(1:end-1)]);
% %     marketProbs=marketPrepayments.*marketRemain;
% %     
% %     % Find realised lambda
% %     N=size(marketPrepayments,1);
% %     
% %     % bisection
% %     tempProbs=2;
% %     ub=lambda+1;
% %     lb=max(lambda-1,0);
% %     while abs(sum(marketProbs)-sum(tempProbs))>0.0001
% %         tempProbs=nan(N,1);
% %         mb=0.5*(ub+lb);
% %         for i=1:N
% %             tempProbs(i)=prepaymentProbability(Xstar(i),Xstar(1:i-1),mb,Alpha,Beta,n);
% %         end
% %         if sum(marketProbs)>sum(tempProbs)
% %             lb=mb;
% %         else
% %             ub=mb;
% %         end
% %     end
% %     
% %     prepayProbs=nan(size(Xstar));
% %     LAMBDA=nan(size(Xstar));
% %     LAMBDA(1:N)=mb;
% %     LAMBDA(N+1:end)=lambda;
% %     for i=1:size(Xstar,1)
% %         prepayProbs(i)=prepaymentProbability(Xstar(i),Xstar(1:i-1),LAMBDA(1:i),Alpha,Beta,n);
% %     end
% %     
% %     %
% %     %  scaling=sum(prepayProbs(1:size(marketProbs,1)))/sum(marketProbs);
% %     %  prepayProbs(1:size(marketProbs,1))=marketProbs;
% %     %  prepayProbs(size(marketProbs,1)+1:end)=prepayProbs(size(marketProbs,1)+1:end)*scaling;
% % end

% Calculate prepayment rates
remaining=1-cumsum([0;prepayProbs(1:end-1)]);
val=prepayProbs./remaining;
end