function [] = plot_polerr_hollow(varargin)
%PLOT_POLAR_HOLLOW Plots error in B from mumat vs polar angle
%   All data must use the grid of points defined in points.dat
%   Call using plot_polerr_hollow(ir, path to .dat files )
%   ir: radial index, 1 <= ir <= # of r values
%   path to .dat files: character array, e.g. 'file.dat'

if nargin <= 1
    disp('Not enough arguments')
    return
end

try
    ir = varargin{1};
    if ~isa(ir, 'double')
        disp('Radial index not of type double')
        return
    end
catch
    disp('Could not read radial index')
    return
end

ir = max(ir,1);

thetas = [0:pi/50:pi]; 
ntheta = length(thetas);

B_exact = importdata("Hollow sphere\B_exact_hollow.dat"); 
B_exact_x = B_exact(:,1); B_exact_z = B_exact(:,3);
mod_B_exact = sqrt(sum(B_exact.^2,2));

len = length(B_exact); len = len/ntheta;

if ir > len
    disp(['Radial index out of bounds, clamping to ' num2str(len)])
    ir = len;
end

xyz = importdata("Hollow sphere\points_hollow.dat");
dist = norm(xyz(1 + ntheta*(ir-1),:));

filec = nargin-1;
hx = gobjects(1,filec);
hz = gobjects(1,filec);
hn = gobjects(1,filec);
for i = 1:filec
    file = varargin{i+1};
    try
        B_data = importdata(file);
        
    catch
        disp(['Failed to read file ' file])
        return
    end
    B_data_x = B_data(:,1); B_data_z = B_data(:,3);
    mod_B_data = sqrt(sum(B_data.^2,2));
    error_Bx = (B_data_x-B_exact_x);%./abs(B_exact_x+0.0001);
    error_Bz = (B_data_z-B_exact_z)./abs(B_exact_z);
    error_Bn = (mod_B_data-mod_B_exact)./mod_B_exact;

%     if i == 1
%         xfig = figure('Name','Error in x component','NumberTitle','off');
%     else
%         figure(xfig)
%     end
%     hx(i) = plot(thetas, error_Bx(1 + (ntheta*(ir-1) ):ntheta*ir), ...
%         '-o','DisplayName',file);
%     hold on

    if i == 1
        zfig = figure('Name','Error in z component','NumberTitle','off');
    else
        figure(zfig)
    end
    L = strfind(file,'\');
    n = file(L(end)+1:end);
    hz(i) = plot(thetas, error_Bz(1 + (ntheta*(ir-1) ):ntheta*ir), ...
        '-o','DisplayName',n,'LineWidth',2,'MarkerSize',8);
    hz(i).MarkerFaceColor = hz(i).Color;
    hold on    
% 
%     if i == 1
%         nfig = figure('Name','Error in norm','NumberTitle','off');
%     else
%         figure(nfig)
%     end
%     
%     hn(i) = plot(thetas, error_Bn(1 + (ntheta*(ir-1) ):ntheta*ir), ...
%         '-o','DisplayName',file);
%     hold on 

end

%     figure(nfig)
%     yline(0,'k--');
%     legend(hn,'box','off','interpreter','none')
%     xlabel('polar angle [rad]')
%     ylabel('$\epsilon\ |B|$','interpreter','latex')
%     xlim([0 pi])
%     ylim([-0.06 0.06])    
%     set(gca, 'xtick',[0 pi/4 pi/2 3*pi/4 pi], ...
%     'XTickLabel',{'$0$', '$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'}, ...
%     'TickLabelInterpreter','latex')
%     title(['dx = ' char(string(dist*10)) ' mm'])

    figure(zfig)
    yline(0,'k--','LineWidth',2);
    legend(hz,'box','off','interpreter','none')
    xlabel('Polar angle [rad]','interpreter','none','FontSize',20)
    ylabel('Error in magnetic field strength','interpreter','none','FontSize',20)
    xlim([0 pi])
    %ylim([-0.06 0.06])    
    set(gca, 'xtick',[0 pi/4 pi/2 3*pi/4 pi], ...
    'XTickLabel',{'$0$', '$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'}, ...
    'TickLabelInterpreter','latex')
    set(gca,'FontSize',20)
    title(['Relative error in magnetic field strength for hollow sphere at ' char(string(5000-dist*10)) ' mm'], ...
        'interpreter','none');
%     
%     figure(xfig)
%     yline(0,'k--');
%     legend(hx,'box','off','interpreter','none')
%     xlabel('polar angle [rad]')
%     ylabel('$\Delta\ B_x$','interpreter','latex')
%     xlim([0 pi])
%     ylim([-2e-4 1e-4]) 
%     set(gca, 'xtick',[0 pi/4 pi/2 3*pi/4 pi], ...
%     'XTickLabel',{'$0$', '$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'}, ...
%     'TickLabelInterpreter','latex')
%     title(['dx = ' char(string(dist*10)) ' mm'])


end

