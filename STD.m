function Targets = STD(StaRgns, StaMap, SalMap)

Targets = [];

MasterMap = im2double(StaMap).*im2double(SalMap);
if ~isempty(StaRgns)
	% compute the average saliency of each stable region
	tmpmean = zeros(StaRgns.RegionNums,1);
	assignin('base','StaRngs',StaRgns);
	for k = 1:StaRgns.RegionNums
		pixellist = StaRgns.Regions(k).Props.PixelList;
		idx = sub2ind(size(StaMap), pixellist(:,2), pixellist(:,1));
		tmpmean(k) = mean(MasterMap(idx));
	end
	% the average saliency of all stable regions
	mean_tmpmean = mean(tmpmean);
	
	% targets
	targetNo = 1;
	for k = 1:StaRgns.RegionNums
		if tmpmean(k) >= mean_tmpmean
			Targets.Regions(targetNo) = StaRgns.Regions(k);
			targetNo = targetNo+1;
		end
	end
	Targets.RegionNums = targetNo-1;
end

end