function [iter , cost] = extractChacheAtDuration(directory,duration)
 
    searchstr = [directory '*BHrampInitialFinal.txt'];
    Files = dir( searchstr );
    i = 1;
    for k = 1:length(Files)
        filename    = [directory Files(k).name];
        rampData    = dlmread(filename);
    
        % scan for ID of given duration
        if abs(rampData(end,1)-duration) < 1e-4
           fileID{i} = strtok(filename,'_');
           i = i+1;
        end  
    end
    
    l = 1;
    for j = 1:i-1
        searchstr = [fileID{j} '*ProgressCache.txt'];
        Files = dir( searchstr ); 
        for k = 1:length(Files)
            filename    = [directory Files(k).name];
            cacheData   = dlmread(filename);
    
            cost_ar{l}       = cacheData(:,2);  
            iter_lng(l)      = size(cacheData,1);
            l = l+1;
        end
    end
    
    
    % extend data such it has uniform length
    max_lng = max(iter_lng);
    iter = (1:max_lng)';
    
    cost    = zeros(max_lng,length(iter_lng));
    
    for k = 1:length(iter_lng)      
       % extend cost
       co   = cost_ar{k};
       co( iter_lng(k):max_lng ) = co(iter_lng(k));
       cost(:,k) = co;
       
%        % extend control
%        u = control_ar{k};
%        ext = max_lng-iter_lng(k)+1;
%        u(iter_lng(k):max_lng,:) = repmat(u(iter_lng(k),:),ext,1);
%        control(:,:,k) = u;
    end
    
end