n = 1000;
vBeta = [1;0.5];
phi = 100;
[vy, mX] = betasim(n, vBeta, phi);
[R2, vP] = betareg(vy, mX);