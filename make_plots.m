%% clean up
close all; clear;

%% load data
x = importdata('points.dat');
B = importdata('B.dat');
B_act = importdata('Bexact.dat');

rho = sqrt(sum(x.^2,2));
modB = sqrt(sum(B.^2,2));
modBact = sqrt(sum(B_act.^2,2));

colors={'r','g','b','c','m','y','k'};
intv = 7;
pts = length(x);
B(1) = 0;

for i = 1:intv
    plot(rho(i:intv:pts),modB(i:intv:pts),'o','Color',colors{i});
    hold on;
    plot(rho(i:intv:pts),modBact(i:intv:pts),'-','Color',colors{i});
end



xlim([0 rho(end)])