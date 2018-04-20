function fig = makeControlPlot(controlData)

    %  ---- data format -----
    % time | control init | Fidelity init | control final | Fidelity final |

    set(0, 'DefaultTextInterpreter', 'latex');
    set(0, 'DefaultLegendInterpreter', 'latex');
    set(0, 'defaultAxesTickLabelInterpreter','latex');
    set(0, 'defaultAxesFontSize',12);

    time    = controlData(:,1);
    ui      = controlData(:,2);
    Fi      = controlData(:,3);
    uf      = controlData(:,4);
    Ff      = controlData(:,5);

    co = [    0    0.4470    0.7410
         0.8500    0.3250    0.0980
         0.9290    0.6940    0.1250
         0.4940    0.1840    0.5560
         0.4660    0.6740    0.1880
         0.3010    0.7450    0.9330
         0.6350    0.0780    0.1840];
    
    fig = figure;
    box on
    hold on

    p1(1) = plot(time,uf,'k-');
    p1(2) = plot(time,ui,'k--');
    xlabel('Time')
    ylabel('Control U')

    yyaxis right
    hold on
    p2(1) = plot(time,(1-Ff),'-','Color',co(1,:));
    p2(2) = plot(time,(1-Fi),'--','Color',co(1,:));

    xlim([0 time(end)])
    ylabel('Infidelity')
    set(gca,'YScale','log')
    legend(p1,'Final','Initial','Location','west')

end