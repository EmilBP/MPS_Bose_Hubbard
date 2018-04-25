function processBestSolutions(readDirectory, writeDirectory)
 
    % extract data
    [Tlist, Flist, fnlist] = extractFinalFidelity(readDirectory);
    TFdata = [Tlist' , Flist'];
    
    % sort data with respect to T
    sortedData = sortrows(TFdata,1);
    
    
    % find filename corresponding to best F for each T
    [uv,idy,~] = unique(sortedData(:,1));
    nu = numel(uv);
    
    idy = [idy ; size(sortedData,1)+1];
    
    for i = 1:nu
        % group F for given T in x
        x = sortedData(idy(i):idy(i+1)-1,2);
        [maxF,~]     = max(x);
       
        % find index of filename of highest fidelity
        [~,index]    = ismember([uv(i) , maxF],TFdata,'rows');
        
        % plot <N> for best solutions
        searchstr    = [fnlist{index} 'ExpectationN.txt'];
        Files = dir( searchstr ); 
        for k = 1:length(Files)
            filename    = [readDirectory Files(k).name];
            expNData    = dlmread(filename);
%             hN(i)       = makeExpNumberPlot(expNData);
        end
        
        % plot control and infidelity for best solutions
        searchstr    = [fnlist{index} 'BHrampInitialFinal.txt'];
        Files = dir( searchstr ); 
        for k = 1:length(Files)
            filename    =  [readDirectory Files(k).name];
            controlData = dlmread(filename);
%             hC(i)       = makeControlPlot(controlData);
        end
        
        hNC(i) = makeRampPlot(expNData,controlData);
        disp(['Best solution for T = ' num2str(uv(i)) ' has ID ' erase(fnlist{index},readDirectory)])
        pause(0.3)
    end    

    % save figures 
    for i = 1:nu
%         fig = hN(i);
%         figname = [ writeDirectory 'ExpectedN_T' num2str(uv(i)) ];
%         fig.PaperPositionMode = 'auto';
%         fig_pos = fig.PaperPosition;
%         fig.PaperSize = [fig_pos(3) fig_pos(4)];
%         print(fig,figname,'-dpdf','-bestfit')
%         savefig(fig,figname)
        
%         fig = hC(i);
%         figname = [ writeDirectory 'OptimalRamp_T' num2str(uv(i)) ];
%         fig.PaperPositionMode = 'auto';
%         fig_pos = fig.PaperPosition;
%         fig.PaperSize = [fig_pos(3) fig_pos(4)];
%         print(fig,figname,'-dpdf','-bestfit')
%         savefig(fig,figname)
        
        fig = hNC(i);
        figname = [ writeDirectory 'RampPlot_T' num2str(uv(i)) ];
        fig.PaperPositionMode = 'auto';
        fig_pos = fig.PaperPosition;
        fig.PaperSize = [fig_pos(3) fig_pos(4)];
        print(fig,figname,'-dpdf','-bestfit')
%         savefig(fig,figname)
    end
    
end