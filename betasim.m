function [vy, mX] = betasim(n, vBeta, phi)
rand('state', 0);
randn('state', 0);
mX = [ones(n,1) randn(n,length(vBeta)-1)];
mu = exp(mX*vBeta)./(1+exp(mX*vBeta));
vy = betarnd(mu*phi, (1-mu)*phi);