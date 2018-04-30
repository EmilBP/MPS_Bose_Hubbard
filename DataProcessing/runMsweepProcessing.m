function main()
    clear all; close all; clc;
    
    Mlist = 5:40;
    prefix = 640480:640480+length(Mlist);
    
    readDirectory = '../../../mnt/LinSigSeedN5MSweepT2.5/';
    writeDirectory = '../../DataProcessing/Plots/LinSigSeedN5MSweepT2.5/';
    
    processData(readDirectory,writeDirectory,Mlist,prefix);
end

function processData(readDirectory, writeDirectory,BasisSize,prefix)

    PlotData = [BasisSize', zeros(length(BasisSize),4)];

    for i = 1:length(BasisSize)
        % find data for M = BasisSize(i)
        searchstr = [readDirectory num2str(prefix(i)) '*BHrampInitialFinal.txt'];
        Files = dir( searchstr );
        FData = zeros(1,length(Files));
        for k = 1:length(Files)
            filename     = [readDirectory Files(k).name];
            fidData      = dlmread(filename);
            FData(k)     = fidData(end,end);
        end
        
        PlotData(i,2) = median(FData);
        PlotData(i,3) = prctile(FData,25);
        PlotData(i,4) = prctile(FData,75);
        PlotData(i,5) = max(FData);
        
        [BasisSize(i) , length(Files)]
    end    
    
    
    % plot data
    legtext = {'GROUP'};
    fig2 = makeCompBasisFidelityPlot(PlotData,legtext);
    fig1 = makeBasisFidelityPlot(PlotData,legtext);
    
    % save figure 
    figname = [writeDirectory 'FidelityBasisSize'];
    fig1.PaperPositionMode = 'auto';
    fig_pos = fig1.PaperPosition;
    fig1.PaperSize = [fig_pos(3) fig_pos(4)];
    print(fig1,figname,'-dpdf','-bestfit') 
    savefig(fig1,figname)
    
     % save figure 
    figname = [writeDirectory 'BestFidelityBasisSize'];
    fig2.PaperPositionMode = 'auto';
    fig_pos = fig2.PaperPosition;
    fig2.PaperSize = [fig_pos(3) fig_pos(4)];
    print(fig2,figname,'-dpdf','-bestfit') 
    savefig(fig2,figname)
end

function [BasisSize , prefix] = extractBasisSize(directory)
 
    searchstr = [directory '*ProgressCache.txt'];
    Files = dir( searchstr ); 
    for k = 1:length(Files)
        filename    = [directory Files(k).name];
        CacheData   = dlmread(filename);
    
        BasisSize(k)= length(CacheData(1,:))-3;
        prefix{k}   = strtok(filename,'_');
    end
end

function fig = makeCompBasisFidelityPlot(data,legendEntries)

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
    
    M = data(:,1);
    
    % plot fidelity data for each method in subplot
    sub(1) = subplot(2,1,1);
    hold on
    box on
    for i = 1:(size(data,2)-1)/4
        algdata = data(:,(2+(i-1)*4):(1+i*4));
        color   = co(i,:);
        
        xx = [M' , fliplr(M')];
        yy = [algdata(:,2)' , fliplr(algdata(:,3)')];
        fill(xx,yy,color,'FaceAlpha',0.4,'EdgeColor','none');
        p(i) = plot(M,algdata(:,2),'Linewidth',2,'Color',color);
        plot(M,algdata(:,3),'Linewidth',2,'Color',color);
        plot(M,algdata(:,1),'LineStyle',':','Linewidth',2,'Color',color);
    end
       
    limsy=get(gca,'YLim');
    set(gca,'Ylim',[limsy(1) 1]);
    %set(gca, 'YScale', 'log')
    
%     xlabel('Duration T')
    ylabel('Fidelity $F$')
    
    legend(p,legendEntries,'Location','SouthEast')

    
    % plot best infidelity achieved for each method in subplot
    sub(2) = subplot(2,1,2);
    hold on
    box on
    for i = 1:(size(data,2)-1)/4
        algdata = data(:,(2+(i-1)*4):(1+i*4));
        color   = co(i,:);
        
        plot(M,1-algdata(:,4),'Linewidth',2,'Color',color);
    end
    
    xlabel('Basis Size')
    ylabel('Infidelity $1-F$')
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
    
    
    annotation(gcf,'textbox',...
    [0.122304084595696 0.873777279000283 0.0305186246418337 0.0547703180212016],...
    'String','(\textbf{a})',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FitBoxToText','off');
    
    samexaxis('xmt','on','ytac','join','yld',1,'YTickAntiClash')
    
    annotation(gcf,'textbox',...
    [0.122536928251765 0.262536430729644 0.0305186246418337 0.0547703180212017],...
    'String','(\textbf{b})',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FitBoxToText','off');
    
end

function fig = makeBasisFidelityPlot(data,legendEntries)

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
    
    M = data(:,1);
    
    hold on
    box on
    for i = 1:(size(data,2)-1)/4
        algdata = data(:,(2+(i-1)*4):(1+i*4));
        color   = co(i,:);
        
        xx = [M' , fliplr(M')];
        yy = [algdata(:,2)' , fliplr(algdata(:,3)')];
        fill(xx,yy,color,'FaceAlpha',0.4,'EdgeColor','none');
        p(i) = plot(M,algdata(:,2),'Linewidth',2,'Color',color);
        plot(M,algdata(:,3),'Linewidth',2,'Color',color);
        plot(M,algdata(:,1),'LineStyle',':','Linewidth',2,'Color',color);
    end
       
    limsy=get(gca,'YLim');
    set(gca,'Ylim',[limsy(1) 1]);

    xlim([M(1) M(end)])
    ylabel('Fidelity F')
    xlabel('Basis Size')
    
    legend(p,legendEntries,'Location','SouthEast')

       
    
end