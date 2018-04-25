function fig = makeRampPlot(expNData,controlData)

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
    ni      = 1:size(expNData,2)-1;
    expN    = expNData(:,2:end);
    
    co = [    0    0.4470    0.7410
         0.8500    0.3250    0.0980
         0.9290    0.6940    0.1250
         0.4940    0.1840    0.5560
         0.4660    0.6740    0.1880
         0.3010    0.7450    0.9330
         0.6350    0.0780    0.1840];
    
    fig = figure;
    
    subplot(2,1,1)
    hold on
    box on
    
    imagesc('XData',time,'YData',ni,'CData',expN');
    xlim([time(1) , time(end)])
    ylim([ni(1)-0.5 , ni(end)+0.5])
    ylabel('Lattice Site')
    yticks(ni)
    
    x = parula(211);
    %     index = true(1, size(x, 1));
    %     index([15:2:23 , 25:27 , 29:32 , 34:37  , 39:41 , 43:2:51]) = false;
    %     custommap = x(index, :);

    x1 = x(1:100,:);
    x2 = x(112:211,:);
    x1int = interp1(1:100,x1,logspace(0,2,150));
    x2int = interp1(1:100,flipud(x2),logspace(0,2,150));
    custommap = [ x1int ; x(106,:) ;  flipud(x2int)];

    colormap(custommap)

    c = colorbar;
    caxis([0.8 1.2])

  
    
    subplot(2,1,2)
    box on
    hold on

    p1(1) = plot(time,uf,'k-','LineWidth',2);
    p1(2) = plot(time,ui,'k--','LineWidth',2);
    xlabel('Time')
    ylabel('Control $U$')

    xlim([time(1) , time(end)])
    ylim([0 , 55])
    
    annotation(gcf,'textbox',...
    [0.022304084595696 0.873777279000283 0.0305186246418337 0.0547703180212016],...
    'String','(\textbf{a})',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FitBoxToText','off');
    
    annotation(gcf,'textbox',...
    [0.022536928251765 0.468197879858652 0.0305186246418337 0.0547703180212017],...
    'String','(\textbf{b})',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FitBoxToText','off');

    samexaxis('xmt','on','ytac','join','yld',1,'YTickAntiClash')
    
    
    yyaxis right
    hold on
    p2(1) = plot(time,(1-Ff),'-','Color',co(1,:),'LineWidth',2);
    p2(2) = plot(time,(1-Fi),'--','Color',co(1,:),'LineWidth',2);
    ylabel('Infidelity')
    set(gca,'YScale','log')
    
    yl = ylim;
    
    ytickspace = logspace(-10,-1,10);
    ytickvals  = [];
    for i = 1:length(ytickspace)
        if yl(1) <= ytickspace(i) && yl(2) >= ytickspace(i)
            ytickvals = [ytickvals , ytickspace(i)];
        end
    end
    
    if length(ytickvals) >= 1
        yticks( ytickvals )
    end
    
    legend(p1,'Final','Initial','Location','west')
end