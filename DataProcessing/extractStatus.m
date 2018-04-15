function status = extractStatus(directory)
 
    % status format: | T | Nsucc | Nfail |

    searchstr = [directory '*Status.txt'];
    Files = dir( searchstr ); 
    for k = 1:length(Files)
        filename        = [directory Files(k).name];
        statData(k,:)   = dlmread(filename);
    end
    
    sortedData = sortrows(statData,1);
    
    [uv,idy,~] = unique(sortedData(:,1));
    nu = numel(uv);
    status = zeros(nu,3);
    
    idy = [idy ; size(sortedData,1)+1];
    for i = 1:nu
       status(i,1) = uv(i);
       % group data for given T in x
       x = sortedData(idy(i):idy(i+1)-1,2);
       status(i,2) = sum(x);
       status(i,3) = length(x)-sum(x);
    end   
    
end