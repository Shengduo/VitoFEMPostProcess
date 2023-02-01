function plotFrontSurf_function(videoprefix, times)
    % Read results from hdf5 files.
    load('BRColorScale.mat');
    frontFileName = strcat('../frontsurfFiles/', videoprefix, '-frontsurf.h5');
    fontsize = 20;   

    %% Get the video of fault-parallel and fault-vertical velocity
    % Read time
    time = h5read(frontFileName, '/time');
    time = reshape(time, [1, size(time, 3)]);
    time = time - 10e-6;
    
    % Read node geometry
    nodalXYZ = h5read(frontFileName, '/geometry/vertices');
    nodalXYZ = nodalXYZ(1:2, :);
    
    % Read nodal slip
    % Input the first wire position
    WirePos1 = [-0.025657; -0.014222; 0];

    VSstart = [0.006354, 0.003522, 0]';
    VSend = [0.063204, 0.035034, 0]';
    % VSregion = 1e3 * ([norm(VSstart - WirePos1, 2), norm(VSend - WirePos1, 2)] - norm(VSstart - WirePos1, 2));
    
    % Read nodal slip, slip rate and traction
    velocity = h5read(frontFileName, '/vertex_fields/velocity');
    velocity = velocity(1:2, :, :);
    Connection = h5read(frontFileName, '/topology/cells');
    
    % Rotate into fault coordinate system
    alpha = 29 / 180 * pi;
    Q = [cos(alpha), sin(alpha); -sin(alpha), cos(alpha)];
    
    % Rotate the velocity field
    for i = 1:1:size(velocity, 3)
        velocity(:, :, i) = Q * velocity(:,:,i);
    end
    
    % Rotate nodal XYZ
    QnodalXYZ = Q * nodalXYZ;
    
    % Set plot range
    xEnlarge = 20;
    yEnlarge = 15;
    xWinRange = [0, 47];
    yWinRange = [-15, 15];
    xRange = [xWinRange(1) - xEnlarge, xWinRange(2) + xEnlarge];
    yRange = [yWinRange(1) - yEnlarge, yWinRange(2) + yEnlarge];
    cRange = [-1, 1];
    
    % Initiate figures
    figNo = 1;
    %% Save a image of fault-parallel velocity field
    for i = 1 : 1 : length(times)
        [~, t_idx] = min(abs(time - times(i) * 1e-6));
        fig = figure(figNo);
        fig.Position = [1000, 597, 2240/2.5, 1680/2.5];
        
        % Initialize names
        plotname = strcat(pwd, '/../PaperPlots/', videoprefix, '_frontsurfV_prl_t', num2str(times(i)), '.eps');
        
        % Save the figure
        print('-vector', figure(figNo) ,plotname, '-depsc', '-r500');
        
        % Shoot the video
        p = patch(1e3 * QnodalXYZ(1, :), 1e3 * QnodalXYZ(2, :), ...
            velocity(1, :, t_idx), 'EdgeColor', 'none', 'FaceColor', 'interp');
        p.Faces = Connection' + 1;
        colormap((black_rainbow_plus_long));
        
        yline(0, 'LineWidth', 2.0, 'color', 'w');
        
        c = colorbar;
        caxis(cRange);
        ylabel(c,'Fault-parallel Velocity [m/s]', 'FontName', 'Avenir', 'FontSize', fontsize, 'interpreter', 'latex');
        axis equal;
        grid on;
        xlim(xRange);
        ylim(yRange);
        
        % colorStr = '#0072BD';
        colorStr = '#4DBEEE';
        % Draw the boundary of observation window
        hold on; 
        plot([xWinRange(1), xWinRange(2)], [yWinRange(1), yWinRange(1)], '--', 'linewidth', 2.0, 'color', colorStr);
        plot([xWinRange(2), xWinRange(2)], [yWinRange(1), yWinRange(2)], '--', 'linewidth', 2.0, 'color', colorStr);
        plot([xWinRange(2), xWinRange(1)], [yWinRange(2), yWinRange(2)], '--', 'linewidth', 2.0, 'color', colorStr);
        plot([xWinRange(1), xWinRange(1)], [yWinRange(2), yWinRange(1)], '--', 'linewidth', 2.0, 'color', colorStr);
        text(12, 18, "Field of view", 'color', 'w', 'Fontsize', fontsize);
        text(54, 25, strcat(num2str(time(t_idx) * 1e6, '%.0f'), ' $\mathrm{\mu s}$'), 'color', 'w', 'Fontsize', fontsize);
        
        % Label the figures
        xlabel('$x_1$ [mm]');
        ylabel('$x_2$ [mm]');
        
        % title(strcat('Time = { }', num2str(time(t_idx) * 1e6, '%.2f'), ' [$\mathrm{\mu s}$]'));
        set(gca, 'FontSize', fontsize);
        set(gcf,'color','w');
        
        % Save the figure
        print('-vector', figure(figNo) ,plotname, '-depsc', '-r500');
        figNo = figNo + 1;
    end
    
    % Close all figures
    close all;
    
end