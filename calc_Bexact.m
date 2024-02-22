clc; clear;

%% Setup
spheretype = 'hollow'; 
%spheretype = 'solid';
p = '\\wsl.localhost\Debian\home\lvh\dist_test_hollow\';
xyz = importdata([ p 'points.dat']);
%% 
mu_0 = 4*pi*1e-7;
zhat = [0 0 1];

%% System constants

R = 0.5; % Major radius [m]
a = R; % Radius inner sphere (hollow)
b = 2*R; % Radius outer sphere (hollow)
mu_prime = 500; % Relative magnetic permeability 
%mu_prime = mu_r/mu_0; % Helper mu
B0 = [0 0 1.0]; % External B field strength

switch spheretype
    case 'hollow' % Hollow sphere
        %% Calculate delta
        num = 9*mu_prime;
        denom = (2*mu_prime+1)*(mu_prime+2)-2*(a^3/b^3)*(mu_prime-1)^2;
        H0 = norm(B0)/mu_0; % External H field strength
        delta = -num/denom*H0;

        fnc_Bexact = @(r) -mu_0*delta*zhat;

    case 'solid' % Solid sphere
        chi_m = mu_prime-1; % Magnetic susceptibility

        fnc_Bexact = @(r) (chi_m)/(1+chi_m/3)*(R/norm(r))^3*( dot(B0,r)*r/norm(r)^2 - B0/3 ) + B0;

end
%% Allocate array for B_exact then calculate
B_exact = zeros(size(xyz));
for i = 1:length(xyz)
    B_exact(i,:) = fnc_Bexact(xyz(i,:)/1000); %divide by 1000 since xyz in mm
end

%% Export
save([p 'B_exact.dat'], 'B_exact', '-ascii')