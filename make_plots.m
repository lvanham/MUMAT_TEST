%% clean up
close all; clear;

%% load data
xyz = importdata('points.dat');
B = importdata('B.dat');
B_act = importdata('Bexact.dat');

%% scatter locations [2D]
x = xyz(:,1); z = xyz(:,3);
if exist('ptsfig')
    delete(ptsfig)
end
ptsfig = figure('Name','sampled points','NumberTitle','off');
scatter(x,z,2); hold on;
h = rectangle( 'Position',[-500 -500 1000 1000],'Curvature',[1 1],'FaceColor',[0.5 0.5 0.5],'LineStyle','--');
xlim([0 max(x)]); ylim([0 max(z)])
xlabel('x [cm]'); ylabel('y [cm]')
axis square
%%
rho = sqrt(sum(xyz.^2,2));
modB = sqrt(sum(B.^2,2));
modBact = sqrt(sum(B_act.^2,2));

intv = 7;
pts = length(xyz);
B(1) = 0;

%%
if exist('compfig')
    delete(compfig)
end
l = ["0$" "\pi/12$" "\pi/6$" "\pi/4$" "\pi/3$" "5\pi/6$" "\pi/2$"];
compfig = figure('Name','|B| vs |Bexact|','NumberTitle','off');
h = gobjects(1,intv); g = gobjects(1,intv);
for i = 1:intv
    g(i) = plot(rho(i:intv:pts),modB(i:intv:pts),'.');
    hold on;
    h(i) = plot(rho(i:intv:pts),modBact(i:intv:pts),'DisplayName',['$\theta = ' char(l(i))], 'Color', g(i).Color);
end
legend(h,'interpreter', 'latex','box','off')
%%
errB = (modB-modBact)./modBact;
if exist('errfig') %#ok<*EXIST> 
    delete(errfig)
end

errfig = figure('Name','Error in |B|','NumberTitle','off');
h2 = gobjects(1,intv);
for i = 1:intv
    h2(i) = plot(rho(i:intv:pts),errB(i:intv:pts),'LineStyle','-','Marker','.','DisplayName',['$\theta = ' char(l(i))],'Color',h(i).Color);
    hold on
end
legend(h2,'interpreter', 'latex','box','off')
ylim([-0.05 0.05])
xlim([490 2000])

