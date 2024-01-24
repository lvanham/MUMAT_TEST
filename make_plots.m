%% clean up
close all; clear;

%% load data
xyz = importdata('points.dat');
B = importdata('B.dat');
B_act = importdata('Bexact.dat');

%% scatter locations [2D]
x = xyz(:,1); z = xyz(:,3);
if exist('ptsfig')
    close(ptsfig)
end
ptsfig = figure('Name','sampled points','NumberTitle','off','Position',[10 10 800 800]);
scatter(x,z,2); hold on;
h = rectangle( 'Position',[-500 -500 1000 1000],'Curvature',[1 1],'FaceColor',[0.5 0.5 0.5]);
xlim([0 max(x)]); ylim([0 max(z)])
axis square
%%
rho = sqrt(sum(xyz.^2,2));
modB = sqrt(sum(B.^2,2));
modBact = sqrt(sum(B_act.^2,2));

colors={'r','g','b','c','m','y','k'};
intv = 7;
pts = length(xyz);
B(1) = 0;

for i = 1:intv
    plot(rho(i:intv:pts),modB(i:intv:pts),'o','Color',colors{i});
    hold on;
    plot(rho(i:intv:pts),modBact(i:intv:pts),'-','Color',colors{i});
end



xlim([0 rho(end)])