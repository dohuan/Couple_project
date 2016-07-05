function y = betalik(vP, mX, vy)
k = length(vP);
eta = mX*vP(1:k-1); 
mu = exp(eta) ./ (1+exp(eta)); 
phi = vP(k);
y = -sum( gammaln(phi) - gammaln(mu*phi)- gammaln((1-mu)*phi) + (mu*phi-1) .* log(vy) + ( (1-mu)*phi-1 ) .* log(1-vy) );

% g = zeros(k,1);
% ynew = log( vy ./ (1-vy) );
% munew = psi(mu*phi) - psi((1-mu)*phi);
% T = diag( exp(eta) ./ (1+exp(eta)) .^2 );
% g(1:k-1) = -phi*mX'*T*(ynew-munew); 
% g(k) = sum( psi(phi) - mu .* psi(mu*phi) - (1-mu) .* psi((1-mu)*phi) + mu .* log(vy) + (1-mu) .* log(1-vy) );