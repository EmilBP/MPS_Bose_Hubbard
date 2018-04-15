function fig = makeFidelityPlot(fidelityData)

    %  ---- data format -----
    % T | median | 25 perc | 75 perc | best | 

    set(0, 'DefaultTextInterpreter', 'latex');
    set(0, 'DefaultLegendInterpreter', 'latex');
    set(0, 'defaultAxesTickLabelInterpreter','latex');
    set(0, 'defaultAxesFontSize',12);

    T       = fidelityData(:,1);
    Fmed    = fidelityData(:,2);
    F25     = fidelityData(:,3);
    F75     = fidelityData(:,4);
    Fmax    = fidelityData(:,5);
    
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
    set(gca, 'YScale', 'log')
    
    color = co(2,:);


    xx = [T' , fliplr(T')];
    yy = [F25' , fliplr(F75')];
    fill(xx,yy,color,'FaceAlpha',0.4,'EdgeColor','none');
    semilogy(T,F25,'Linewidth',2,'Color',color);
    semilogy(T,F75,'Linewidth',2,'Color',color);
    semilogy(T,Fmed,'LineStyle',':','Linewidth',2,'Color',color);
    semilogy(T,Fmax,'LineStyle','--','Linewidth',2,'Color',color);
    
    limsy=get(gca,'YLim');
    set(gca,'Ylim',[limsy(1) 1]);
    
    xlabel('Duration T')
    ylabel('Fidelity F')
    
    %% plot inset
    
    axes('Position',[.57 .2 .3 .3])
    box on
    plot(T,1-Fmax,'LineStyle','-','Linewidth',1.5,'Color',color);
    set(gca, 'YScale', 'log')
    ylim([1e-4 , 1e-2])
    grid on

end