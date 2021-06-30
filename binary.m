function bw = binary(Rgns,height,width)

bw = false(height,width);

if ~isempty(Rgns)
	for k = 1:size(Rgns.Regions,2)
		pixellist = Rgns.Regions(k).Props.PixelList;
		for p = 1:size(pixellist,1)
			bw(pixellist(p,2), pixellist(p,1)) = true;
		end
	end
end

end