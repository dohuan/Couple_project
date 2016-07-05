%  PROGRAM: betareg.m 
% 
%  USE: Estimates the parameters of a beta regression model
% 
%  AUTHOR: Francisco Cribari & Silvia Ferrari. 
%
%  TRANSLATED TO MATLAB: Willem-Jan de Goeij
% 
%  DATE: August 11, 2009. 
%  
%  VERSION: 1.0. 

%  BASED ON: BETA REGRESSION FOR MODELLING RATES AND PROPORTIONS, Francisco Cribari & Silvia Ferrari. 


function betareg(vy, mX);
format short g;
n = length(vy);
p = size(mX,2);

if(max(vy) >= 1 || min(vy) <= 0) 
    error(sprintf('\n\nERROR: DATA OUT OF RANGE (0,1)!\n\n')); 
end

if(p >= n) 
     error(sprintf('\n\nERROR: NUMBER OF COVARIATES CANNOT EXCEED NUMBER OF OBSERVATIONS!\n\n'));
end

ynew = log( vy ./ (1-vy) );

if(p > 1) 
     betaols = mX \ ynew; 
elseif(p==1) 
     betaols = mean(ynew);
end

olsfittednew = mX*betaols; 
olserrorvar = sum((ynew-olsfittednew).^2)/(n-p); 
olsfitted = exp(olsfittednew) ./ (1 + exp(olsfittednew)); 

ybar = mean(vy); 
yvar = var(vy);   

% starting values
vps = [betaols;(mean( 1 ./ (olserrorvar*(olsfitted .* (1-olsfitted)))) - 1)];

disp(sprintf('\nBETA REGRESSION ESTIMATION'));
disp(sprintf('--------------------------')); 

disp(sprintf('\nMEAN AND VARIANCE OF Y: %10.5f', [ybar yvar]));  
disp(sprintf('\nINITIAL VALUES FOR THE ML ESTIMATION: %16.5f', vps));

%options = optimset('GradObj','on');
%vP = fminunc(@(vP) betalik(vP, mX, vy), vps, options);
vP = fminsearch(@(vP) betalik(vP, mX, vy), vps);
 
etahat = mX*vP(1:p); 
muhat = exp(etahat) ./ (1+exp(etahat)); 
phihat = vP(p+1); 
psi1 = psi(1, muhat*phihat); 
psi2 = psi(1,(1-muhat)*phihat); 
T = diag( exp(etahat) ./ (1+exp(etahat)) .^2 );
W = diag(phihat*(psi1+psi2)) * (T .^2); 
vc = phihat*(psi1.*muhat-psi2.*(1-muhat)); 
D = diag(psi1.*(muhat.^2)+psi2.*(1-muhat).^2-psi(1,phihat));
tempinv = inv(mX'*W*mX); 
g = trace(D)-(1/phihat)*vc'*T'*mX*tempinv*mX'*T*vc; 
K1 = tempinv*(g*eye(p)+(1/phihat)*mX'*T*vc*vc'*T'*mX*tempinv);
K2 = -tempinv*mX'*T*vc; 
fisherinv = (1/(phihat*g)) * ( [K1,K2; -vc'*T'*mX*tempinv,phihat]);
disp(sprintf('\nPARAMETER ESTIMATES AND ASYMPTOTIC STANDARD ERRORS: '));
stderrors = sqrt(diag(fisherinv)); 
zstats = vP ./ stderrors; 
disp('estimates   std. errors    z stats    p-values'); 
mOutput = [vP stderrors zstats 2*(1-normcdf(abs(zstats)))];
for (i=1:size(mOutput,1))
        disp(sprintf('%9.5f %13.5f %10.5f %11.5f',mOutput(i,:)));
end

disp(sprintf('\nASYMPTOTIC COVARIANCE MATRIX OF ML ESTIMATES:'));
disp(fisherinv);      
pseudoR2 = corr(etahat,ynew)^2;
disp(sprintf('\nPSEUDO R2: %8.4f', pseudoR2));           

plot(1:n, vy, 'b'); hold on;
plot(1:n, muhat, 'r'); hold off;
title('Simulated values (blue) vs fitted values (red)');