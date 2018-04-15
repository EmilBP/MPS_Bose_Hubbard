function fig = makeCompFidelityPlot(data,legendEntries)

    %  ---- data format -----
    % T | median (1) | 25 perc (1) | 75 perc (1) | best (1) | median (2) |
    % ...

    set(0, 'DefaultTextInterpreter', 'latex');
    set(0, 'DefaultLegendInterpreter', 'latex');
    set(0, 'defaultAxesTickLabelInterpreter','latex');
    set(0, 'defaultAxesFontSize',12);
    
    co = [    0    0.4470    0.7410
         0.8500    0.3250    0.0980
         0.9290    0.6940    0.1250
         0.4940    0.1840    0.5560
         0.4660    0.6740    0.1880
         0.3010    0.7450    0.9330
         0.6350    0.0780    0.1840];
    
    fig = figure;  
    
    T = data(:,1);
    
    % plot fidelity data for each method in subplot
    sub(1) = subplot(2,1,1);
    hold on
    box on
    for i = 1:(size(data,2)-1)/4
        algdata = data(:,(2+(i-1)*4):(1+i*4));
        color   = co(i,:);
        
        xx = [T' , fliplr(T')];
        yy = [algdata(:,2)' , fliplr(algdata(:,3)')];
        fill(xx,yy,color,'FaceAlpha',0.4,'EdgeColor','none');
        p(i) = plot(T,algdata(:,2),'Linewidth',2,'Color',color);
        plot(T,algdata(:,3),'Linewidth',2,'Color',color);
        plot(T,algdata(:,1),'LineStyle',':','Linewidth',2,'Color',color);
    end
       
    limsy=get(gca,'YLim');
    set(gca,'Ylim',[limsy(1) 1]);
    %set(gca, 'YScale', 'log')
    
%     xlabel('Duration T')
    ylabel('Fidelity F')
    
    legend(p,legendEntries,'Location','SouthEast')

    
    % plot best infidelity achieved for each method in subplot
    sub(2) = subplot(2,1,2);
    hold on
    box on
    for i = 1:(size(data,2)-1)/4
        algdata = data(:,(2+(i-1)*4):(1+i*4));
        color   = co(i,:);
        
        plot(T,1-algdata(:,4),'Linewidth',2,'Color',color);
    end
    
    xlabel('Duration T')
    ylabel('Infidelity 1-F')
    limsy=get(gca,'YLim');
    set(gca,'Ylim',[1e-4 0.8*1e-1]);
    yticks([1e-4 1e-3 1e-2])
    set(gca, 'YScale', 'log')
    ax = gca;
    ax.YGrid = 'on';
    ax.YMinorGrid = 'off';
    
    
    % adjust subplot size and stack them
    subpos = get(sub(1), 'Position');
    set(sub(1), 'position', [subpos(1), subpos(2)-subpos(4)*0.5, subpos(3), subpos(4)*1.5] );
    subpos = get(sub(2), 'Position');
    set(sub(2), 'position', [subpos(1), subpos(2), subpos(3), subpos(4)*0.5] );
    samexaxis('xmt','on','ytac','join','yld',1,'YTickAntiClash')
    
end