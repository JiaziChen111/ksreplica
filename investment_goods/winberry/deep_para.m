bbeta = 0.95;
ttau = 0.1; % search cost
aalpha = 0.3; % y = a*k^aalpha*l^v
v = 0.6;
aalpha0 = 2; % search elasticity
ddelta = 0.02;
pphi = 0.0000001; % price of disinvestment relative to investment
rrhox = 0.98;
ssigmax_low = 0.01;
ssigmax_high= 0.011;
ssigmax_grid = [ssigmax_low,ssigmax_high];
Pssigmax = [0.9 0.1; 0.1 0.9];
ppsi = 0; % quadratic cost of investment adjustment
