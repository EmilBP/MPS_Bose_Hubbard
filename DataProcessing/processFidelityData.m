function PlotData = processFidelityData(readDirectory, writeDirectory)

    % extract data
    [Tlist, Flist, ~] = extractFinalFidelity(readDirectory);
    TFdata = [Tlist' , Flist'];
    
    % sort data with respect to T
    sortedData = sortrows(TFdata,1);
    
    % calculate median and percentiles for each unique T
    % data format T:median:25percentile:75percentile
    [uv,idy,~] = unique(sortedData(:,1));
    nu = numel(uv);
    PlotData = zeros(nu,5);
    
    idy = [idy ; size(sortedData,1)+1];
    for i = 1:nu
       PlotData(i,1) = uv(i);
       % group data for given T in x
       x = sortedData(idy(i):idy(i+1)-1,2);
       PlotData(i,2) = median(x);
       PlotData(i,3) = prctile(x,25);
       PlotData(i,4) = prctile(x,75);
       PlotData(i,5) = max(x);
    end   
    
    % plot data
    %fig = makeFidelityPlot(PlotData);
    legtext = {'GROUP'};
    fig = makeCompFidelityPlot(PlotData,legtext);
    
    % save figure 
    figname = [writeDirectory 'FidelityDuration'];
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    print(fig,figname,'-dpdf','-bestfit') 
    savefig(fig,figname)
end

