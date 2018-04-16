function processStatus(readDirectory, writeDirectory)

    % extract status for all runs
    status = extractStatus(readDirectory);
    
    % plot data
    set(0, 'DefaultTextInterpreter', 'latex');
    set(0, 'DefaultLegendInterpreter', 'latex');
    set(0, 'defaultAxesTickLabelInterpreter','latex');
    set(0, 'defaultAxesFontSize',12);
    
    fig = figure
    bar(status(:,1)',status(:,2:3),'stacked')
    xlabel('Duration T')
    ylabel('Number of runs')
    legend('Succes','Failed')
    
    % save figure 
    figname = [writeDirectory 'RunStatus'];
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    print(fig,figname,'-dpdf','-bestfit') 
end

