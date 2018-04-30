function main()

    clear all; close all; clc;

    readDirectory = '../../../mnt/LinSigSeedN5fix3/';
    writeDirectory = '../../DataProcessing/Plots/LinSigSeedN5fix3/';

    Durations = {'1.0','2.0','2.5','3.0','4.0'};
    
    for i = 1:length(Durations)
        expNData = dlmread( [readDirectory 'ExpectationN_extendedT' Durations{i} '.txt'] );
        rampData = dlmread( [readDirectory 'BHrampInitialFinal_extendedT' Durations{i} '.txt'] );
        
        fig = makeExtendedRampPlot(expNData,rampData);
        figname = [ writeDirectory 'ExtendedRampPlot_T' Durations{i} ];
        fig.PaperPositionMode = 'auto';
        fig_pos = fig.PaperPosition;
        fig.PaperSize = [fig_pos(3) fig_pos(4)];
        print(fig,figname,'-dpdf','-bestfit')
    end
   
end

function fig = makeExtendedRampPlot(expNData,controlData)

    T           = controlData(end-100,1);
    trunc       = ceil(6/5 *(size(expNData,1) - 100));
    expNData    = expNData(1:trunc,:);
    controlData = controlData(1:trunc,:);


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
    

    x = parula(301);
    offset = 5;
    mid = ceil(length(x)/2);
    x1 = x(1:mid-offset,:);
    x2 = x(mid+offset:end,:);

    
    custommap = [x1 ; x(mid,:) ; x2 ];
    colormap(custommap)

    cticklab    = 0.8:0.1:1.2;
   
    
    scaledData  = symmetricLogScale(expN , 1);
    scaledcaxis = symmetricLogScale([0.8 , 1.2] , 1);
    scaledctick = symmetricLogScale(cticklab , 1);
    
    imagesc('XData',time,'YData',ni,'CData',scaledData');
    xlim([time(1) , time(end)])
    ylim([ni(1)-0.5 , ni(end)+0.5])
    ylabel('Lattice Site')
    yticks(ni)
    
    caxis(scaledcaxis)
    h_c = colorbar('YTick',scaledctick,'YTickLabel',cticklab);
    h_c.Ruler.MinorTick = 'on';
    
    
    plot([T,T], [ni(1)-0.5 , ni(end)+0.5],'--','Color',co(2,:),'LineWidth',1.5);
  
    
    subplot(2,1,2)
    box on
    hold on

    plot([T,T], [0, 55],'--','Color',co(2,:),'LineWidth',1.5);
    p1(1) = plot(time,uf,'k-','LineWidth',2);
    p1(2) = plot(time,ui,'k--','LineWidth',2);
    
    xlabel('Time')
    ylabel('Control $U$')

    xlim([time(1) , time(end)])
    ylim([0 , 55])
    
    yticks(0:10:50)
    
    
    annotation(gcf,'textbox',...
    [0.022304084595696 0.873777279000283 0.0305186246418337 0.0547703180212016],...
    'String','(\textbf{a})',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FitBoxToText','off');
    
    annotation(gcf,'textbox',...
    [0.022304084595696 0.468197879858652 0.0305186246418337 0.0547703180212017],...
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
    plt = gca;
    plt.YAxis(2).Color = co(1,:);
    
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

function scaledData = symmetricLogScale(data,center)
    
    below       = data.*(data < center);
    above       = data.*(data >= center);
    
    below       = below + ones(size(below)).*(below == 0);
    above       = above + ones(size(above)).*(above == 0);
    
    below       = abs(below-center)+center;
    
    abovelog    = log10(above);
    belowlog    = log10(below);
    scaledData  = abovelog-belowlog;
end