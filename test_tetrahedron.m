clear;
%%

tet = read_mumat('Tetrahedron\tetrahedron_mu.dat');
coords = tet.coords';

tet_vert_x = coords(:,1); tet_cen_x = sum(tet_vert_x,1)/4;
tet_vert_y = coords(:,2); tet_cen_y = sum(tet_vert_y,1)/4;
tet_vert_z = coords(:,3); tet_cen_z = sum(tet_vert_z,1)/4;
tet_cen = [tet_cen_x tet_cen_y tet_cen_z];

%% Import grid and Bexact
xyz = importdata('Tetrahedron\points_tetrahedron.dat');
B = importdata('Tetrahedron\B_tetrahedron.dat');
thetas = [0:pi/5:pi]; ntheta =  length(thetas);
thetas_str = ["0"; "\frac{\pi}{5}"; "\frac{2\pi}{5}"; "\frac{3\pi}{5}"; "\frac{4\pi}{5}";"\pi"];
phis = [0:2*pi/5:2*pi]; nphi = length(phis);
phis_str =["0"; "\frac{2\pi}{5}"; "\frac{4\pi}{5}"; "\frac{6\pi}{5}"; "\frac{8\pi}{5}";"2\pi"];
%linestyles = ["-"; "--"; "-."; ":";"-";"--"];
markers = [ "square"; "diamond"; "*"; "o";"+"; "^"];
fs = 50;

%%
x = xyz(:,1); y = xyz(:,2); z = xyz(:,3);

x_sub = zeros(1,length(x)/36);
y_sub = zeros(1,length(x)/36);
z_sub = zeros(1,length(x)/36);
B_dat = zeros(1,length(x)/36);

%% Color gradient 
col_num = 6;
col_arr = zeros(col_num,3);
col_arr(1,:) = [1 0 0];
col_arr(col_num,:) = [0 0.2 1];

for i = 2:col_num-1
    col_arr(i,:) = col_arr(1,:) + (col_arr(col_num,:)-col_arr(1,:))/(col_num-1)*(i-1);
end


close all
scatterfig = figure('Position',1e3*[  -0.038200000000000   1.376200000000000   0.560000000000000   0.420000000000000]);
%scatter3(tet_vert_x, tet_vert_y, tet_vert_z,'k',...
%    'LineWidth',2,'MarkerFaceColor','k');
scatter3(tet_cen_x, tet_cen_y, tet_cen_z,'r','MarkerFaceColor','r')
grid off
hold on

plot3(tet_vert_x, tet_vert_y, tet_vert_z,'k-',...
    'LineWidth',2)
plot3([tet_vert_x(1) tet_vert_x(3)], ...
    [tet_vert_y(1) tet_vert_y(3)], ...
    [tet_vert_z(1) tet_vert_z(3)],'k-',...
    'LineWidth',2)
plot3([tet_vert_x(1) tet_vert_x(4)], ...
    [tet_vert_y(1) tet_vert_y(4)], ...
    [tet_vert_z(1) tet_vert_z(4)],'k-',...
    'LineWidth',2)
plot3([tet_vert_x(2) tet_vert_x(4)], ...
    [tet_vert_y(2) tet_vert_y(4)], ...
    [tet_vert_z(2) tet_vert_z(4)],'k-',...
    'LineWidth',2)

xlim(100*[-2 2]+tet_cen_x);
ylim(100*[-2 2]+tet_cen_y);
zlim(100*[-2 2]+tet_cen_z);
xl = xlim;
xlim([xl(1) 100]);

Bfig = figure('outerposition',[1 1 1920 1440]);
h = gobjects(1,ntheta+1);
l = 1;
for theta_ind = 1:ntheta
    for phi_ind = 4:4

        j = 1;
        for i = phi_ind+(theta_ind-1)*6:36:length(x)
            x_sub(j) = x(i);
            y_sub(j) = y(i);
            z_sub(j) = z(i);
            B_dat(j) = B(i);
            j = j + 1;
        end
        B_dat = abs(B_dat)/max(abs(B_dat));
        r_sub = vecnorm([x_sub-tet_cen_x; y_sub-tet_cen_y; z_sub-tet_cen_z]);
        figure(Bfig)
        h(l) = plot(0.01*r_sub, B_dat, ...
            'LineWidth',12,...%'Color',col_arr(theta_ind,:), ...
            'DisplayName',['$\theta= ' char(thetas_str(theta_ind)) '$']);% ', \phi= ' char(phis_str(phi_ind)) '$']);
        h(l).Marker = markers(phi_ind);
        h(l).MarkerIndices = 1:10:length(x_sub);
        h(l).MarkerSize = 18;
        h(l).MarkerFaceColor = h(l).Color;
        hold on
        figure(scatterfig)
        %scatter3(x_sub(1:10:length(x_sub)), y_sub(1:10:length(y_sub)), z_sub(1:10:length(z_sub)), ...
        %    'MarkerEdgeColor',h(l).Color,...
        %    'MarkerFaceColor',h(l).Color);
        l = l + 1;

    end
end

figure(Bfig)
x2=[1.2:0.001:5];
y2=1./x2.^3;
h(l) = plot(x2,y2/max(y2),'k-.','LineWidth',8,'DisplayName','$y=1/x^3$');
%legend(h,'box','off','FontSize',20)
xlabel('Distance from tetrahedron center [m]','interpreter','none','FontSize',fs)
ylabel('Normalised magnetic field strength','interpreter','none','FontSize',fs)
title('Dipole behavior of a magnetised tetrahedron','interpreter','none','FontSize',fs)
xlim([1.2 3]);

%
% figure(scatterfig)
% caz =  -1.506756898778590e+02;
% cel =    0.187850749372330;
% 
% view(caz,cel);
% axis off
% 
% saveas(gcf,'Figures\tetrahedron.png');
% close(scatterfig)
% 
% figure(Bfig)
% img = imread('Figures\tetrahedron.png');
% image([1.9 2.99],[0.3 0.99],img);
% ylim([0 1]);
%
% x1 = 2.31; x2 = 2.6; x3 = 2.76;
% text(x1,0.94,'$\theta=0$','FontSize',fs,'interpreter', 'latex');
% text(x2,0.94,'$\theta=\pi/5$','FontSize',fs,'interpreter', 'latex');
% text(x3,0.73,'$\theta=2\pi/5$','FontSize',fs,'interpreter', 'latex');
% text(x3,0.53,'$\theta=3\pi/5$','FontSize',fs,'interpreter', 'latex');
% text(x2,0.32,'$\theta=4\pi/5$','FontSize',fs,'interpreter', 'latex');
% text(x1,0.32,'$\theta=\pi$','FontSize',fs,'interpreter', 'latex');
%
legend(h(1:end-1),'box','off','FontSize',fs,'Location','northeast','NumColumns',1)
        Bfig.CurrentAxes.FontSize=fs;
ylim([0 1.05])