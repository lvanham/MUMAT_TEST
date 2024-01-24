%% Clean up any leftovers from earlier run
close all; clear;
%% Load simulation data
xyz = importdata('points.dat');
B = importdata('B.dat');
num_theta = 51;

%% Calculate B_exact
% Constants
R = 5; % Major radius [m]
Rcm = 100*R;
B0 = [0 0 1]; % Applied B-field [T]
mu_r = 500; % Relative magnetic permeability 
chi_m = mu_r-1; % Magnetic susceptibility

% Function for calculating B_exact
fnc_Bexact = @(r) (chi_m)/(1+chi_m/3)*(R/norm(r))^3*( dot(B0,r)*r/norm(r)^2 - B0/3 ) + B0;

% Allocate array for B_exact then calculate
B_exact = zeros(size(xyz));
for i = 1:length(xyz)
    B_exact(i,:) = fnc_Bexact(xyz(i,:)/100); % /100 as xyz is in cm
end
%B_exact = importdata('Bexact.dat'); %c.f. this if necessary
%% scatter locations [2D]
x = xyz(:,1); z = xyz(:,3);
if exist('ptsfig')
    delete(ptsfig)
end
ptsfig = figure('Name','sampled points','NumberTitle','off');
scatter(x,z,2); hold on;
h = rectangle( 'Position',Rcm*[-1 -1 2 2],'Curvature',[1 1],'FaceColor',[0.5 0.5 0.5],'LineStyle','--');
xlim([min(x) max(x)]); ylim([min(z) max(z)])
xlabel('x [cm]'); ylabel('y [cm]')
axis square

%%
rho = sqrt(sum(xyz.^2,2));
modB = sqrt(sum(B.^2,2));
modBact = sqrt(sum(B_exact.^2,2));

pts = length(xyz);
B(1) = 0;

%% Comparing |B| to |Bexact| 
if exist('compfig')
    delete(compfig)
end

compfig = figure('Name','|B| vs |Bexact|','NumberTitle','off');
h1 = gobjects(1,num_theta); h2 = gobjects(1,num_theta);
for i = 1:num_theta
    h1(i) = plot(rho(i:num_theta:pts)-Rcm,modB(i:num_theta:pts),'.');
    hold on;
    h2(i) = plot(rho(i:num_theta:pts)-Rcm,modBact(i:num_theta:pts), 'Color', h1(i).Color);
    xlabel('dist from sphere [cm]'); ylabel('mod B [T]')
end
cols = zeros(num_theta,3);
for i = 1:num_theta
    cols(i,:) = h1(i).Color;
end
xlim(Rcm*[-0.02 1])
%% Error in |B|
errB = (modB-modBact)./modBact;
if exist('errfig') %#ok<*EXIST> 
    delete(errfig)
end

errfig = figure('Name','Error in |B|','NumberTitle','off');
h3 = gobjects(1,num_theta);
for i = 1:num_theta
    h3(i) = plot(rho(i:num_theta:pts)-Rcm,errB(i:num_theta:pts),'LineStyle','-','Marker','.','Color',cols(i,:));
    hold on
    xlabel('dist from sphere [cm]'); ylabel('error in mod B')

end

ylim([-0.05 0.05])
xlim(Rcm*[-0.02 1])

%% Error with polar angle
theta = [0:pi/50:pi];

if exist('thetafig') %#ok<*EXIST> 
    delete(thetafig)
end
r_indices = [4 16 31 151 301 1501 3001];
thetafig = figure('Name','Error with polar angle','NumberTitle','off'); 
for i = 1:length(r_indices)
    r_ind = r_indices(i);
    plot(theta,errB(1+(num_theta*(r_ind-1)):num_theta*r_ind),'-o','Color',cols(i,:),'MarkerFaceColor',cols(i,:),'DisplayName',['dr=' char(string(norm(xyz(r_ind*num_theta,:))-500)) ' cm' ]);
    hold on
end
xlim([0 pi]); ylim([-0.1 0.1])
legend('box','off')
%% Component-wise error in B
Bexact_x = B_exact(:,1); Bexact_z = B_exact(:,3);
B_x = B(:,1); B_z = B(:,3);
errBx = (B_x-Bexact_x)./Bexact_x;
errBz = (B_z-Bexact_z)./Bexact_z;
if exist('errfigx') %#ok<*EXIST> 
    delete(errfigx)
end

errfigx = figure('Name','Error in Bx','NumberTitle','off');
h4 = gobjects(1,num_theta);
for i = 1:num_theta
    h4(i) = plot(rho(i:num_theta:pts)-Rcm,errBx(i:num_theta:pts),'LineStyle','-','Marker','.','Color',cols(i,:));
    hold on
    xlabel('dist from sphere [cm]'); ylabel('error in Bx')
end

ylim([-0.05 0.05])
xlim(Rcm*[-0.02 1])

if exist('errfigz') %#ok<*EXIST> 
    delete(errfigz)
end

errfigz = figure('Name','error in Bz','NumberTitle','off');
h5 = gobjects(1,num_theta);
for i = 1:num_theta
    h5(i) = plot(rho(i:num_theta:pts)-Rcm,errBz(i:num_theta:pts),'LineStyle','-','Marker','.','Color',cols(i,:));
    hold on
    xlabel('dist from sphere [cm]'); ylabel('error in Bz')
end

ylim([-0.05 0.05])
xlim(Rcm*[-0.02 1])