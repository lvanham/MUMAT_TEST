function [] = MUMATplot(ir, plottype, ptsfile, Bfile, cells, lexclmid)
%MUMATplot Make different plots for verification of MUMAT
%   Call using MUMATplot(ir, plottype, ptsfile, Bfile, cells )
%   ir: radial index, 1 <= ir <= # of r values
%   ptsfile: path to coordinate file
%   Bfile: path to file with exact B field
%   cells: cell object with paths to simulation B fields
%   lexclmid: bool for disabling (1) datapoint at theta=pi/2
%   plottypes: ([theta] = function of theta)
%       'comp_x' 
%           plot x-component of B-field [theta]
%       'comp_z' 
%           plot z-component of B-field [theta]
%       'mod_B'
%           plot |B| [theta]
%       'diff_x' 
%           diff between x-comps of B_sim and B_exact [theta]
%       'diff_z' 
%           diff between z-comps of B_sim and B_exact [theta]
%       'diff_mod_B' 
%           diff between |B_sim| and |B_exact| [theta]
%       'error_x'
%           plot relative error between x-comps of B_sim and B_exact[theta]
%       'error_z'
%           plot relative error between z-comps of B_sim and B_exact[theta]
%       'error_mod_B'
%           plot relative error between B_sim and B_exact [theta]
%       'error_mod_B_avg'
%           plot rel. error between B_sim and B_exact averaged over theta

%% Cases
cases = {'comp_x'; 'comp_z'; 'mod_B'; 'diff_x'; 'diff_z'; 'diff_mod_B';...
    'error_x'; 'error_z';  'error_mod_B'; 'error_mod_B_avg'}; 
if ~contains(plottype,cases)
    disp([plottype ' not recognized as a valid plot type.'])
    return
end

%% Setup: Reading
if ~isa(ir, 'double')
    disp('Radial index not of type double')
    return
end
ir = max(ir,1);

R = 500; % Radius in mm
theta = [0:pi/36:pi]; % array of angles
ntheta = length(theta);

%% Read points file
try
    xyz = importdata(ptsfile);
catch
   disp(['Failed to read ptsfile ' ptsfile])
   return
end

try
    B_exact = importdata(Bfile); 
    B_exact_x = B_exact(:,1); 
    B_exact_z = B_exact(:,3); 
    mod_B_exact = sqrt(sum(B_exact.^2,2));

    len = length(B_exact)/ntheta;
    if ir > len
        disp(['Radial index out of bounds, clamping to ' num2str(len)])
        ir = len;
    end
    dist = norm(xyz(1 + ntheta*(ir-1),:))-R;
    slice = [1 + (ntheta*(ir-1) ):ntheta*ir];
catch
   disp(['Failed to read Bfile ' Bfile])
   return
end 

B_struct = struct;
filec = length(cells);

for i = 1:filec
    file = cells{i};
    try
        B_temp = importdata(file);
        B_temp_x = B_temp(:,1); 
        B_temp_z = B_temp(:,3);

        L = strfind(file,'\');
        n = file(L(end)+1:end);
       
        B_struct.(['file' char(string(i)) '_name']) = n;
        B_struct.(['file' char(string(i)) '_path']) = file;
        B_struct.(['file' char(string(i)) '_B' ]) = B_temp  ;
        B_struct.(['file' char(string(i)) '_Bx']) = B_temp_x;
        B_struct.(['file' char(string(i)) '_Bz']) = B_temp_z;
        B_struct.(['file' char(string(i)) '_modB']) = sqrt(sum(B_temp.^2,2));
     
    catch
        disp(['Failed to read file ' file])
        return
    end
end

%% Plot options
lw = 4;

%% Case selection
switch plottype
   
    case 'comp_x'

        hx = gobjects(1,filec+1);
        figure;
        hx(1) = plot(theta,B_exact_x(slice),'k-', ...
            'DisplayName','Exact', ...
            'LineWidth',lw);
        hold on
        for i = 1:filec
            stp = ['file' char(string(i))];
            data = B_struct.([stp '_Bx'])(slice);
            hx(i+1) = plot(theta, data,'o', ...
                'DisplayName',B_struct.([stp '_name']), ...
                'LineWidth',lw, ...
                'MarkerFaceColor','auto');
            hx(i+1).MarkerFaceColor = hx(i+1).Color;
        end
        legend(hx,'box','off','interpreter','none')
        xlabel('polar angle [rad]')
        ylabel('$B_x$ [T]','interpreter','latex')
        xlim([0 pi])
        ylim([-1.5 1.5]*max(abs(B_exact_x)));
        set(gca, 'xtick',[0 pi/4 pi/2 3*pi/4 pi], ...
        'XTickLabel',{'$0$', '$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'}, ...
        'TickLabelInterpreter','latex')
        title(['distance = ' char(string(dist)) ' mm'])

    case 'comp_z'

        hz = gobjects(1,filec+1);
        figure;
        hz(1) = plot(theta,B_exact_z(slice),'k-', ...
            'DisplayName','Exact', ...
            'LineWidth',lw);
        hold on
        for i = 1:filec
            stp = ['file' char(string(i))];
            data = B_struct.([stp '_Bz'])(slice);
            hz(i+1) = plot(theta, data,'o', ...
                'DisplayName',B_struct.([stp '_name']), ...
                'LineWidth',lw);
            hz(i+1).MarkerFaceColor = hz(i+1).Color;
        end
        legend(hz,'box','off','interpreter','none')
        xlabel('polar angle [rad]')
        ylabel('$B_z$ [T]','interpreter','latex')
        xlim([0 pi])
        ylim([-1.5 1.5]*max(abs(B_exact_z)));
        set(gca, 'xtick',[0 pi/4 pi/2 3*pi/4 pi], ...
        'XTickLabel',{'$0$', '$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'}, ...
        'TickLabelInterpreter','latex')
        title(['distance  = ' char(string(dist)) ' mm'])
    
    case 'mod_B'

        hz = gobjects(1,filec+1);
        figure;
        hz(1) = plot(theta,mod_B_exact(slice),'k-', ...
            'DisplayName','Exact', ...
            'LineWidth',lw);
        hold on
        for i = 1:filec
            stp = ['file' char(string(i))];
            data = B_struct.([stp '_modB'])(slice);
            hz(i+1) = plot(theta, data,'-o', ...
                'DisplayName',B_struct.([stp '_name']), ...
                'LineWidth',lw);
            hz(i+1).MarkerFaceColor = hz(i+1).Color;
        end
        legend(hz,'box','off','interpreter','none')
        xlabel('polar angle [rad]')
        ylabel('$|B|$ [T]','interpreter','latex')
        xlim([0 pi])
        ylim([-1.5 1.5]*max(abs(B_exact_z)));
        set(gca, 'xtick',[0 pi/4 pi/2 3*pi/4 pi], ...
        'XTickLabel',{'$0$', '$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'}, ...
        'TickLabelInterpreter','latex')
        title(['distance  = ' char(string(dist)) ' mm'])

    case 'diff_x'

        hx = gobjects(1,filec+1);
        figure; hx(1) = yline(0,'k--','LineWidth',lw);
        hold on
        for i = 1:filec
            stp = ['file' char(string(i))];
            data = B_struct.([stp '_Bx'])(slice)-B_exact_x(slice);
            hx(i+1) = plot(theta, data,'o', ...
                'DisplayName',B_struct.([stp '_name']), ...
                'LineWidth',lw, ...
                'MarkerFaceColor','auto');
            hx(i+1).MarkerFaceColor = hx(i+1).Color;
        end
        legend(hx(2:end),'box','off','interpreter','none')
        xlabel('polar angle [rad]')
        ylabel('$\Delta B_x$ [T]','interpreter','latex')
        xlim([0 pi])
        set(gca, 'xtick',[0 pi/4 pi/2 3*pi/4 pi], ...
        'XTickLabel',{'$0$', '$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'}, ...
        'TickLabelInterpreter','latex')
        title(['distance = ' char(string(dist)) ' mm'])
      
    case 'diff_z'

        hz = gobjects(1,filec+1);
        figure; hz(1) = yline(0,'k--','LineWidth',lw);
        hold on
        for i = 1:filec
            stp = ['file' char(string(i))];
            data = B_struct.([stp '_Bz'])(slice)-B_exact_z(slice);
            hz(i+1) = plot(theta, data,'-o', ...
                'DisplayName',B_struct.([stp '_name']), ...
                'LineWidth',lw);
            hz(i+1).MarkerFaceColor = hz(i+1).Color;
        end
        legend(hz(2:end),'box','off','interpreter','none')
        xlabel('polar angle [rad]')
        ylabel('$\Delta B_z$ [T]','interpreter','latex')
        xlim([0 pi])
        set(gca, 'xtick',[0 pi/4 pi/2 3*pi/4 pi], ...
        'XTickLabel',{'$0$', '$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'}, ...
        'TickLabelInterpreter','latex')
        title(['distance  = ' char(string(dist)) ' mm'])
    
    case 'diff_mod_B'

        hz = gobjects(1,filec+1);
        figure; hz(1) = yline(0,'k--','LineWidth',lw);
        hold on
        for i = 1:filec
            stp = ['file' char(string(i))];
            data = B_struct.([stp '_modB'])(slice)-mod_B_exact(slice);
            hz(i+1) = plot(theta, data,'-o', ...
                'DisplayName',B_struct.([stp '_name']), ...
                'LineWidth',lw);
            hz(i+1).MarkerFaceColor = hz(i+1).Color;
        end
        legend(hz(2:end),'box','off','interpreter','none')
        xlabel('polar angle [rad]')
        ylabel('$\Delta |B|$ [T]','interpreter','latex')
        xlim([0 pi])
        set(gca, 'xtick',[0 pi/4 pi/2 3*pi/4 pi], ...
        'XTickLabel',{'$0$', '$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'}, ...
        'TickLabelInterpreter','latex')
        title(['distance  = ' char(string(dist)) ' mm'])
    
    case 'error_x'

        hx = gobjects(1,filec+1);
        figure; hx(1) = yline(0,'k--','LineWidth',lw);
        hold on
        for i = 1:filec
            stp = ['file' char(string(i))];
            data = B_struct.([stp '_Bx'])(slice)-B_exact_x(slice);
            data = data./abs(B_exact_x(slice));
            hx(i+1) = plot(theta, data,'-o', ...
                'DisplayName',B_struct.([stp '_name']), ...
                'LineWidth',lw, ...
                'MarkerFaceColor','auto');
            hx(i+1).MarkerFaceColor = hx(i+1).Color;
        end
        legend(hx(2:end),'box','off','interpreter','none')
        xlabel('polar angle [rad]')
        ylabel('error in $B_x$','interpreter','latex')
        xlim([0 pi])
        set(gca, 'xtick',[0 pi/4 pi/2 3*pi/4 pi], ...
        'XTickLabel',{'$0$', '$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'}, ...
        'TickLabelInterpreter','latex')
        title(['distance = ' char(string(dist)) ' mm'])
      
    case 'error_z'

        hz = gobjects(1,filec+1);
        figure; hz(1) = yline(0,'k--','LineWidth',lw);
        hold on
        for i = 1:filec
            stp = ['file' char(string(i))];
            data = B_struct.([stp '_Bz'])(slice)-B_exact_z(slice);
            data = data./abs(B_exact_z(slice));
            hz(i+1) = plot(theta, data,'-o', ...
                'DisplayName',B_struct.([stp '_name']), ...
                'LineWidth',lw);
            hz(i+1).MarkerFaceColor = hz(i+1).Color;
        end
        legend(hz(2:end),'box','off','interpreter','none')
        xlabel('polar angle [rad]')
        ylabel('error in $B_z$','interpreter','latex')
        xlim([0 pi])
        set(gca, 'xtick',[0 pi/4 pi/2 3*pi/4 pi], ...
        'XTickLabel',{'$0$', '$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'}, ...
        'TickLabelInterpreter','latex')
        title(['distance  = ' char(string(dist)) ' mm'])

    case 'error_mod_B'

        hx = gobjects(1,filec+1);
        figure; hx(1) = yline(0,'k--','LineWidth',lw);
        hold on

        if lexclmid
            theta(ceil(end/2)) = [];
            ntheta = length(theta);
            rectangle( 'Position',[pi/2-pi/(2*ntheta) -10 pi/ntheta  20],'FaceColor',[0.7 0.7 0.7],'LineStyle','--');
        end

        yl1 = 100000; yl2 = -100000;
        for i = 1:filec
            stp = ['file' char(string(i))];
            data = B_struct.([stp '_modB'])(slice)-mod_B_exact(slice);
            data = data./(mod_B_exact(slice));
            if lexclmid
                data(ceil(end/2)) = [];
                yl1 = min(yl1,min(data)*1.5);
                yl2 = max(yl2,max(data)*1.5);
                ylim([yl1 yl2])
            end
            hx(i+1) = plot(theta, data,'-o', ...
                'DisplayName',B_struct.([stp '_name']), ...
                'LineWidth',lw, ...
                'MarkerFaceColor','auto');
            hx(i+1).MarkerFaceColor = hx(i+1).Color;
        end
        legend(hx(2:end),'box','off','interpreter','none')
        xlabel('polar angle [rad]')
        ylabel('relative error in $|B|$','interpreter','latex')
        xlim([0 pi])
        set(gca, 'xtick',[0 pi/4 pi/2 3*pi/4 pi], ...
        'XTickLabel',{'$0$', '$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'}, ...
        'TickLabelInterpreter','latex')
        title(['distance = ' char(string(dist)) ' mm'])
       
    case 'error_mod_B_avg'

        if lexclmid
            theta(ceil(end/2)) = [];
            ntheta = length(theta);
        end

        hx = gobjects(1,filec);
        figure; 
        hold on
        for i = 1:filec
            stp = ['file' char(string(i))];
            data = B_struct.([stp '_modB'])(slice)-mod_B_exact(slice);
            data = abs(data)./abs(mod_B_exact(slice));
            if lexclmid
                data(ceil(end/2)) = [];
            end
            data = sum(data)/ntheta;

            hx(i) = bar(i,data, ...
                'DisplayName',B_struct.([stp '_name']));
        end
        legend(hx,'box','off','interpreter','none')
        %xlabel('polar angle [rad]')
        ylabel('average relative error in $|B|$','interpreter','latex')
        set(gca, 'XTick',[])
        %'TickLabelInterpreter','latex')
        title(['distance = ' char(string(dist)) ' mm'])
    otherwise
        disp('Plottype invalid.')
       
end  
