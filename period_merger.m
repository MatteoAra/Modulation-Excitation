function merged_int = period_merger(allscans,numframes)

merged_int=[];
for index=1:numframes
    merged_int(:,index)=mean(allscans(:,[index:numframes:size(allscans,2)]),2); 
end

end