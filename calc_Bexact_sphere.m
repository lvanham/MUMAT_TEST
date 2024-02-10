%% Load points
xyz = importdata('points.dat');
%% Define constants
R = 5; % Major radius [m]
Rcm = 100*R;
B0 = [0 0 1]; % Applied B-field [T]
mu_r = 500; % Relative magnetic permeability 
chi_m = mu_r-1; % Magnetic susceptibility
%% Function for calculating B_exact
fnc_Bexact = @(r) (chi_m)/(1+chi_m/3)*(R/norm(r))^3*( dot(B0,r)*r/norm(r)^2 - B0/3 ) + B0;
%% Allocate array for B_exact then calculate
B_exact = zeros(size(xyz));
for i = 1:length(xyz)
    B_exact(i,:) = fnc_Bexact(xyz(i,:)/100); %divide by 100 since xyz in cm
end
%% Export
save B_exact.dat B_exact -ascii