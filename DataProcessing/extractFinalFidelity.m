function [Tlist , Flist, fnlist] = extractFinalFidelity(directory)
 
    searchstr = [directory '*BHrampInitialFinal.txt'];
    Files = dir( searchstr ); 
    for k = 1:length(Files)
        filename    = [directory Files(k).name];
        rampData    = dlmread(filename);
    
        Tlist(k)    = rampData(end,1);
        Flist(k)    = rampData(end,end);
        fnlist{k}   = strtok(filename,'B');
    end
end