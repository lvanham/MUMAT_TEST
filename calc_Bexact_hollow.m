clear;
%% Load points
xyz = importdata('Hollow sphere\points_hollow.dat');
%% 
mu_0 = 4*pi*1e-7;
zhat = [0 0 1];

%% Define system constants

R = 5; % Major radius [m]
a = R; % Radius inner sphere
b = 2*R; % Radius outer sphere
mu_prime = 500; % Relative magnetic permeability 
%mu_prime = mu_r/mu_0; % Helper mu
B0 = 1.0; % External B field strength
H0 = B0/mu_0; % External H field strength

%% Calculate delta
num = 9*mu_prime;
denom = (2*mu_prime+1)*(mu_prime+2)-2*(a^3/b^3)*(mu_prime-1)^2;
delta = -num/denom*H0;

%% Function for calculating B_exact
fnc_Bexact = @(r) -mu_0*delta*zhat;

%% Allocate array for B_exact then calculate
B_exact = zeros(size(xyz));
for i = 1:length(xyz)
    B_exact(i,:) = fnc_Bexact(xyz(i,:)/100); %divide by 100 since xyz in cm
end
%% Export
save B_exact_hollow.dat B_exact -ascii