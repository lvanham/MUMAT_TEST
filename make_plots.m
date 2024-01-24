%% Clean up any leftovers from earlier run
close all; clear;
%% Load simulation data
xyz = importdata('points.dat');
B = importdata('B.dat');
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
xlim([0 max(x)]); ylim([0 max(z)])
xlabel('x [cm]'); ylabel('y [cm]')
axis square

%%
rho = sqrt(sum(xyz.^2,2));
modB = sqrt(sum(B.^2,2));
modBact = sqrt(sum(B_exact.^2,2));

intv = 7;
pts = length(xyz);
B(1) = 0;

%% Comparing |B| to |Bexact| 
if exist('compfig')
    delete(compfig)
end
l = ["0$" "\pi/12$" "\pi/6$" "\pi/4$" "\pi/3$" "5\pi/6$" "\pi/2$"];
compfig = figure('Name','|B| vs |Bexact|','NumberTitle','off');
h1 = gobjects(1,intv); h2 = gobjects(1,intv);
for i = 1:intv
    h1(i) = plot(rho(i:intv:pts)-Rcm,modB(i:intv:pts),'.');
    hold on;
    h2(i) = plot(rho(i:intv:pts)-Rcm,modBact(i:intv:pts),'DisplayName',['$\theta = ' char(l(i))], 'Color', h1(i).Color);
    xlabel('dist from sphere [cm]'); ylabel('mod B [T]')
end
legend(h2,'interpreter', 'latex','box','off')
xlim(Rcm*[-0.02 1])
%% Error in |B|
errB = (modB-modBact)./modBact;
if exist('errfig') %#ok<*EXIST> 
    delete(errfig)
end

errfig = figure('Name','Error in |B|','NumberTitle','off');
h3 = gobjects(1,intv);
for i = 1:intv
    h3(i) = plot(rho(i:intv:pts)-Rcm,errB(i:intv:pts),'LineStyle','-','Marker','.','DisplayName',['$\theta = ' char(l(i))],'Color',h1(i).Color);
    hold on
    xlabel('dist from sphere [cm]'); ylabel('error in mod B')

end
legend(h3,'interpreter', 'latex','box','off')
ylim([-0.05 0.05])
xlim(Rcm*[-0.02 1])

%% Component-wise error in B
Bexact_x = B_exact(:,1); Bexact_z = B_exact(:,3);
B_x = B(:,1); B_z = B(:,3);
errBx = (B_x-Bexact_x)./Bexact_x;
errBz = (B_z-Bexact_z)./Bexact_z;
if exist('errfigx') %#ok<*EXIST> 
    delete(errfigx)
end

errfigx = figure('Name','Error in Bx','NumberTitle','off');
h4 = gobjects(1,intv);
for i = 1:intv
    h4(i) = plot(rho(i:intv:pts)-Rcm,errBx(i:intv:pts),'LineStyle','-','Marker','.','DisplayName',['$\theta = ' char(l(i))],'Color',h1(i).Color);
    hold on
    xlabel('dist from sphere [cm]'); ylabel('error in Bx')
end
legend(h4,'interpreter', 'latex','box','off')

ylim([-0.05 0.05])
xlim(Rcm*[-0.02 1])

if exist('errfigz') %#ok<*EXIST> 
    delete(errfigz)
end

errfigz = figure('Name','error in Bz','NumberTitle','off');
h5 = gobjects(1,intv);
for i = 1:intv
    h5(i) = plot(rho(i:intv:pts)-Rcm,errBz(i:intv:pts),'LineStyle','-','Marker','.','DisplayName',['$\theta = ' char(l(i))],'Color',h1(i).Color);
    hold on
    xlabel('dist from sphere [cm]'); ylabel('error in Bz')
end
legend(h5,'interpreter', 'latex','box','off')

ylim([-0.05 0.05])
xlim(Rcm*[-0.02 1])