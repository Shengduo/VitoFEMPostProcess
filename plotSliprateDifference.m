function plotSliprateDifference(X, Y, sliprates, videoprefix, fmt)
    Diff = squeeze(sliprates(2, :, :) - sliprates(1, :, :));
    
    fontsize = 25;  
    WirePos1 = [-0.025657; -0.014222; 0];
    FaultStart = -[0.100000, 0.055430, 0]';
    FaultEnd = [0.100000, 0.055430, 0]';
    VSstart = [0.006354, 0.003522, 0]';
    VSend = [0.063204, 0.035034, 0]';

    %%
    load('BRColorScale.mat');
    fig = figure(1);
    % Trange = [0, 150];
    Trange = [0, 110];
    Xrange = [-121.599801068916, 107.070134864957];
    VSregion = [0, 64.9994511048604];
    fig.Position(3:4) = 1.5 * fig.Position(3:4);

    % Initialize names
    plotname = strcat(pwd, '/../PaperPlots_interp_1209/', videoprefix, '_sliprateDiff', '.eps');

    % Plot sliprate on X-T
    h = pcolor(X, Y, Diff);
    shading interp;

    colormap(black_rainbow_plus_long);
    hold on;
    xline(VSregion(1), 'r' ,'linewidth', 2.0);
    xline(VSregion(2), 'r' ,'linewidth', 2.0);
    text(VSregion(1)+ 8, 38, 'Region 2', 'color', 'r', 'Fontsize', fontsize);
    
    % Add x_3
    % text(-110, 100, strcat("$x_3$ = ", string(1000 * Distance_To_Surface - 5), "$\ \mathrm{mm}$"), 'color', 'w', 'Fontsize', fontsize);
    % Add the wave speeds
    
    % Several wave speed to show on the aplot
    cp = 2662.4;
    cs = 1279;
    nu = 0.35;
    cr = (0.874 + 0.196 * nu - 0.043 * nu^2 - 0.055 * nu^3) * cs;
    cX = [60, 80] - 15;
    cX = cX - norm(VSstart - WirePos1, 2) * 1e3;
    crY = [10, (cX(2) - cX(1)) * 1e3 / cr + 10];
    csY = [10, (cX(2) - cX(1)) * 1e3 / cs + 10];
    cpY = [10, (cX(2) - cX(1)) * 1e3 / cp + 10];
    
    % Add the wave speeds
    plot(cX, crY, 'w', 'linewidth', 2.0);
    text(cX(2) + 1, crY(2)+5, strcat('$c_r$ = 1.20 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 5, 'interpreter', 'latex');
    plot(cX, csY, 'w', 'linewidth', 2.0);
    text(cX(2) + 1, crY(2) - 1, strcat('$c_s$ = 1.28 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 5, 'interpreter', 'latex');
    plot(cX, cpY, 'w', 'linewidth', 2.0);
    text(cX(2) + 1, cpY(2) - 2, strcat('$c_p$ = 2.66 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 5, 'interpreter', 'latex');
    hold off;
    
    set(h, 'EdgeColor', 'None');
    c = colorbar;
    c2 = 2;
    clim([-c2, c2]);
    ylabel(c,'Slip rate difference [m/s]','FontName','Avenir','FontSize',fontsize, 'interpreter', 'latex');
    xlim(Xrange);
    ylim(Trange);
    title('X-T diagram of slip rate difference')
    xlabel('$x_1\ \mathrm{[mm]}$', 'interpreter', 'latex');
    ylabel('Time [$\mathrm{\mu s}$]', 'interpreter', 'latex');
    % title(['X-T diagram of Slip rate at Depth = ', ' ', num2str(5 - 1e3 * Distance_To_Surface, '%.0f'), ' ', '[mm]'], 'interpreter', 'latex');
    set(gca, 'FontSize', fontsize);
    % daspect([yToxRatio, 1, 1]);
    if c2 == 15
        plotname = strcat(pwd, '/../PaperPlots_interp_1209/', videoprefix, '_X-TofSlipRate_whole_Vito_surface_', num2str(Distance_To_Surface), 'full.eps');
    end
    % Save the figure
    if fmt == 'png'
        plotname = strcat(pwd, '/../PaperPlots_interp_1209/', videoprefix, '_X-TofSlipRate_whole_Vito_surface_', num2str(Distance_To_Surface), '.png');
        print(figure(1), plotname, '-dpng', '-r500');
    else
        print('-vector', figure(1) ,plotname, '-depsc', '-r500');
    end
end

