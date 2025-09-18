function Es=sortcell(En_set)
Es=[];
for cc=1:length(En_set)
Es=[Es;En_set{cc}];
end
Es=sort(Es);


end