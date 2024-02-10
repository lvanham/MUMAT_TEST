function [] = quiver_B_hollow(varargin)
%QUIVER_B_HOLLOW Summary of this function goes here
%   Detailed explanation goes here

intv = 40;
sf = 2000;
Rcm = 500;
%Bexact = importdata("B_exact.dat");
xyz = importdata("C:\Users\lucas\Desktop\MUMAT_TEST\Hollow sphere\points_hollow.dat");

x = xyz(1:intv:end,1); z = xyz(1:intv:end,3);
%Bx_exact = Bexact(1:intv:end,1); Bz_exact = Bexact(1:intv:end,3);

figure;
rectangle( 'Position',Rcm*[-2 -2 4 4],'Curvature',[1 1],'FaceColor',[0.7 0.7 0.7],'LineStyle','--');
rectangle( 'Position',Rcm*[-1 -1 2 2],'Curvature',[1 1],'FaceColor',[1 1 1],'LineStyle','--');

hold on; %quiver(x,z,Bx_exact,Bz_exact,'DisplayName','exact')
xlim(Rcm*[0 1.5])
ylim(Rcm*[-1.5 1.5])
legend('box','off','interpreter','none')
for i = 1:nargin
    try
        B = importdata(varargin{i});
    catch
        disp(['Could not load file ' varargin{i}])
       
        return
    end
    Bx = B(1:intv:end,1); Bz = B(1:intv:end,3);
    quiver(x,z,Bx*sf,Bz*sf,'AutoScale','off','LineWidth',1,'DisplayName',varargin{i})
end

end
