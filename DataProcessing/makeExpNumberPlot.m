function fig = makeExpNumberPlot(expNData)

    set(0, 'DefaultTextInterpreter', 'latex');
    set(0, 'DefaultLegendInterpreter', 'latex');
    set(0, 'defaultAxesTickLabelInterpreter','latex');
    set(0, 'defaultAxesFontSize',12);

    ti      = expNData(:,1);
    ni      = 1:size(expNData,2)-1;
    expN    = expNData(:,2:end);

    fig     = figure;
    imagesc('XData',ti,'YData',ni,'CData',expN');
    xlim([ti(1) , ti(end)])
    ylim([ni(1)-0.5 , ni(end)+0.5])
    xlabel('Time')
    ylabel('Lattice Site')
    yticks(ni)
    colorbar
end