function processCacheData(readDirectory, writeDirectory, duration)

    % extract cache data at duration T
    [iter , cost, ~] = extractChacheAtDuration(readDirectory,duration);
    
    % prepare data for iteration vs cost plot
    iterData(:,1) = iter;
    iterData(:,2) = median(cost,2);
    iterData(:,3) = prctile(cost,25,2);
    iterData(:,4) = prctile(cost,75,2);
    
    % plot iterData
    legtext = {'GROUP'};
    fig = makeIterationPlot(iterData,legtext);
    
    % save figure 
    figname = [writeDirectory 'CostProgressT' num2str(duration)];
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    print(fig,figname,'-dpdf','-bestfit') 
end

