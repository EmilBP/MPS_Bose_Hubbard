function fig = makeIterationPlot(data,legendEntries)

    %  ---- data format -----
    % iter | median (1) | 25 perc (1) | 75 perc (1) | median (2) | ...

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
    
    iter = data(:,1);
  
    fig = figure;
    box on
    hold on
    
    for i = 1:(size(data,2)-1)/3
        algdata = data(:,(2+(i-1)*3):(1+i*3));
        color   = co(i,:);
        
        xx = [iter' , fliplr(iter')];
        yy = [algdata(:,2)' , fliplr(algdata(:,3)')];
        fill(xx,yy,color,'FaceAlpha',0.4,'EdgeColor','none');
        p(i) = plot(iter,algdata(:,2),'Linewidth',2,'Color',color);
        plot(iter,algdata(:,3),'Linewidth',2,'Color',color);
        plot(iter,algdata(:,1),'LineStyle',':','Linewidth',2,'Color',color);
    end

    xlabel('Iteration')
    ylabel('Cost J')
    set(gca,'YScale','log')
    
    xlim([0 , 1250])
    ylim([0.5*min(algdata(:,2)) , 1])
    
    legend(p,legendEntries)

end