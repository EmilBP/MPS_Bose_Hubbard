function [iter , cost, control] = extractChacheAtDuration(directory,duration)
 
    searchstr = [directory '*ProgressCache.txt'];
    Files = dir( searchstr ); 
    i = 1;
    for k = 1:length(Files)
        filename    = [directory Files(k).name];
        cacheData    = dlmread(filename);
    
        % save data from given duration in cells
        if abs(cacheData(1,3)-duration) < 1e-4
           cost_ar{i}       = cacheData(:,2);
           control_ar{i}    = cacheData(:,4:end);
           
           iter_lng(i)      = size(cacheData,1);
           i = i+1;
        end
        
    end
    
    % extend data such it has uniform length
    max_lng = max(iter_lng);
    iter = (1:max_lng)';
    
    for k = 1:i-1      
       % extend cost
       co   = cost_ar{k};
       co( iter_lng(k):max_lng ) = co(iter_lng(k));
       cost(:,k) = co;
       
       % extend control
       u = control_ar{k};
       ext = max_lng-iter_lng(k)+1;
       u(iter_lng(k):max_lng,:) = repmat(u(iter_lng(k),:),ext,1);
       control(:,:,k) = u;
    end
    
end