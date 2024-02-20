clc; clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%        Adjust ranges and points here         %%
% Range and points for radial coordinate
lims_r = [0.0 500.0]; 
num_r = 501;
% Range and points for polar angle
lims_theta = [0 pi];
num_theta = 37;
% Range and points for azimuthal angle
lims_phi = [0 2*pi];
num_phi = 1;
%% Do not touch anything else, you silly goose. %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Spacing
d_r = (lims_r(2)-lims_r(1))/(num_r-1);
d_theta = (lims_theta(2)-lims_theta(1))/(num_theta-1);
d_phi = (lims_phi(2)-lims_phi(1))/(num_phi-1);

% Generate spherical arrays from inputs
arr_r       = [lims_r(1):d_r:lims_r(2)];
arr_theta   = [lims_theta(1):d_theta:lims_theta(2)];
arr_phi     = [lims_phi(1):d_phi:lims_phi(2)];

% Intialise cartesian coordinate arrays
num_pts = num_r*num_theta*num_phi;
new_x = zeros(num_pts,1);
new_y = zeros(num_pts,1);
new_z = zeros(num_pts,1);

% Loop to create xyz
jj = 1;
for i_r = 1:num_r
    r = arr_r(i_r);
    for i_theta = 1:num_theta
        theta = arr_theta(i_theta);
        for i_phi = 1:num_phi
            phi = arr_phi(i_phi);
        
            new_x(jj) = r*sin(theta)*cos(phi);
            new_y(jj) = r*sin(theta)*sin(phi);
            new_z(jj) = r*cos(theta);
            
            jj = jj + 1;
        end
    end
end

new_xyz = [new_x, new_y, new_z];


%% Compare to known poins
comp_xyz = importdata('\\wsl.localhost\Debian\home\lvh\dist_test_hollow\points.dat');
comp_x = comp_xyz(:,1);
comp_y = comp_xyz(:,2); 
comp_z = comp_xyz(:,3);

close all; figure;
plot(new_x-comp_x)
%% Output
formatSpec = '%15.7f,%15.7f,%15.7f\n';
fileID = fopen('points.dat','w');
fprintf(fileID,formatSpec,new_xyz');
fclose(fileID);
%fprintf(formatSpec,new_xyz(15012,:))