intv = 1;
Rcm = 500;

xyz = importdata(char("\\wsl.localhost\Debian\home\lvh\dist_test_hollow\points.dat"));

x = xyz(1:intv:end,1); z = xyz(1:intv:end,3);

figure;
rectangle( 'Position',Rcm*[-2 -2 4 4],'Curvature',[1 1],'FaceColor',[0.7 0.7 0.7],'LineStyle','--');
rectangle( 'Position',Rcm*[-1 -1 2 2],'Curvature',[1 1],'FaceColor',[1 1 1],'LineStyle','--');
hold on; %quiver(x,z,Bx_exact,Bz_exact,'DisplayName','exact')
xlim(Rcm*[0 2])
ylim(Rcm*[-2 2])
xlabel('x [cm]');
ylabel('z [cm]');
%legend('box','off','interpreter','none')

scatter(x,z)