function StabilityRegions = rgnstabmap(gray, param)

StabilityRegions = [];
[height, width] = size(gray);

% sequential segmentation
ThreshNo = 1;
ClusterNo = 1;	% cluster ID
for thresh = param.delta/2 : param.delta : 256-param.delta/2
	binary = im2bw(gray,thresh/255);
	binary = imcomplement(binary);
	regionno = 1;
	
	AllRgns = regionprops(binary,'FilledArea','BoundingBox','Centroid','PixelList');
	existRgns = false;
	
	% limit target size, and remove non-closed targets
	for k = 1:length(AllRgns)
		if AllRgns(k).FilledArea >= param.minarea && AllRgns(k).FilledArea <= height*width*param.maxarea
			tmpPixelList = AllRgns(k).PixelList;
			isclosed = true;
			if ~isempty(find(tmpPixelList(:,1)==1)) || ~isempty(find(tmpPixelList(:,1)==width)) || ...
					~isempty(find(tmpPixelList(:,2)==1)) || ~isempty(find(tmpPixelList(:,2)==height))
				isclosed = false;
			end
			% candidate regions
			if isclosed == true
				Candidate(ThreshNo).Threshold = thresh;
				Candidate(ThreshNo).Regions(regionno).Props = AllRgns(k);
				Candidate(ThreshNo).Regions(regionno).ClusterNo = 0; % cluster No.
				Candidate(ThreshNo).Regions(regionno).HitCount  = 1; % hit count
				regionno = regionno + 1;
				existRgns = true;
			end
		end
	end
	if existRgns == true
		% clustering in current segmented image
		[Candidate(ThreshNo).Regions,ClusterNo] = clustering(Candidate(ThreshNo).Regions,ClusterNo);
		if ThreshNo >= 2
			for k = 1:size(Candidate(ThreshNo).Regions,2)
				for m = 1:ThreshNo-1
					for n = 1:size(Candidate(m).Regions,2)
						if member(Candidate(ThreshNo).Regions(k),Candidate(m).Regions(n)) == true
							Candidate(ThreshNo).Regions(k).HitCount = Candidate(ThreshNo).Regions(k).HitCount+1;
							Candidate(m).Regions(n).HitCount = Candidate(m).Regions(n).HitCount+1;
							if Candidate(ThreshNo).Regions(k).Props.FilledArea > Candidate(m).Regions(n).Props.FilledArea
								Candidate(m).Regions(n).ClusterNo = Candidate(ThreshNo).Regions(k).ClusterNo;
							else
								Candidate(ThreshNo).Regions(k).ClusterNo = Candidate(m).Regions(n).ClusterNo;
							end
						end
					end
				end
			end
		end
		ThreshNo = ThreshNo + 1;
	end
end


if exist('Candidate','var') == 1
	ClusterNo = ClusterNo-1;
	
	% compute the average hit count of each candidate region
	avgHitCount = zeros(1,ClusterNo);
	appearCount = zeros(1,ClusterNo);
	for m = 1:size(Candidate,2)
		for n = 1:size(Candidate(m).Regions,2)
			tmpClusterNo = Candidate(m).Regions(n).ClusterNo;
			avgHitCount(tmpClusterNo) = avgHitCount(tmpClusterNo) + Candidate(m).Regions(n).HitCount;
			appearCount(tmpClusterNo) = appearCount(tmpClusterNo)+1;
		end
	end
	for k = 1:ClusterNo
		if appearCount(k)~=0
			avgHitCount(k) = avgHitCount(k)/appearCount(k);
		end
	end
	
	% only retain the candidate region which satisfies hit count >=2
	HitRgnsMask = zeros(1,ClusterNo);
	for k = 1:ClusterNo
		if avgHitCount(k) >= 2
			HitRgnsMask(k) = 1;
		end
	end
	existHitRgn = false;
	if max(HitRgnsMask) ~= 0
		HitRgnFiltered = filterrgn(Candidate,HitRgnsMask);
		existHitRgn = true;
	end
	
	if existHitRgn == true
		[FillRate,AspectRatio] = combine(HitRgnFiltered,ClusterNo);
	end
	
	if existHitRgn == true
		% Algorithm 2
		SimMeasure = resemb(FillRate,AspectRatio);
		PoStaRgns = [];				% potential candidate regions
		PoStaRgns.RegionNums = 0;	% the number of all potential candidate regions
		for k = 1:size(SimMeasure,2)
			if SimMeasure(1,k)>=1 && SimMeasure(4,k)>=1
				for p = 1:size(HitRgnFiltered(SimMeasure(8,k)).Regions,2)
					region = HitRgnFiltered(SimMeasure(8,k)).Regions(p);
					PoStaRgns = otsu(gray, param, PoStaRgns, region, k);
				end
			end
		end
		
		% post-processing
		if PoStaRgns.RegionNums > 0
			% remove the region which is enclosed in another bigger region
			tmpRgnsNeeded = ones(1,size(PoStaRgns.Regions,2));
			for m = 1:size(tmpRgnsNeeded,2)
				tmpNeedContinue = false;
				for n = 1:size(tmpRgnsNeeded,2)
					if m == n
						tmpNeedContinue = true;
						continue;
					end
					[isEnclose,idx] = enclose(PoStaRgns.Regions(m),PoStaRgns.Regions(n));
					if isEnclose == true
						switch idx
							case 1
								tmpRgnsNeeded(n) = 0;
							case 2
								tmpRgnsNeeded(m) = 0;
						end
					end
				end
				if tmpNeedContinue == true
					continue;
				end
			end
			if PoStaRgns.RegionNums > 1
				for m = 1:size(tmpRgnsNeeded,2)
					if tmpRgnsNeeded(m) == 1
						tmpBoundingBoxM = PoStaRgns.Regions(m).Props.BoundingBox;
						for n = 1:size(tmpRgnsNeeded,2)
							tmpBoundingBoxN = PoStaRgns.Regions(n).Props.BoundingBox;
							if tmpRgnsNeeded(n)==1 && m~=n && sum(tmpBoundingBoxM==tmpBoundingBoxN)==4
								tmpRgnsNeeded(n) = 0;
							end
						end
					end
				end
			end
			
			% detected stable regions
			tmpRgnNo = 1;
			for k = 1:size(tmpRgnsNeeded,2)
				if tmpRgnsNeeded(k) == 1
					StabilityRegions.Regions(tmpRgnNo).Props.Boundary = PoStaRgns.Regions(k).Props.Boundary;
					StabilityRegions.Regions(tmpRgnNo).Props.PixelList = PoStaRgns.Regions(k).Props.PixelList;
					StabilityRegions.Regions(tmpRgnNo).Props.BoundingBox = PoStaRgns.Regions(k).Props.BoundingBox;
					StabilityRegions.Regions(tmpRgnNo).ClusterNo = PoStaRgns.Regions(k).ClusterNo;
					tmpRgnNo = tmpRgnNo+1;
				end
				StabilityRegions.RegionNums = tmpRgnNo-1;	% the number of all stable regions
			end
		end
	end
end
end

function [ClusterRgns, ClusterNo] = clustering(Regions, ClusterNo)
%calls the function member to partition REGIONS into several groups
% by taking account of their spatial relationships, and labels the
% corresponding CLUSTER IDs.

ClusterRgns = Regions;
startInd = 1;
ClusterRgns(1).ClusterNo = ClusterNo;

if size(ClusterRgns,2) == 1		% if only one candidate region
	ClusterNo = ClusterNo+1;
else							% if exist multiple candidate regions
	for m = 2:size(ClusterRgns,2)
		if member(ClusterRgns(m-1),ClusterRgns(m)) == false
			endInd = m-1;
			for n = startInd:endInd
				ClusterRgns(n).ClusterNo = ClusterNo;
			end
			startInd = m;
			ClusterNo = ClusterNo+1;
		end
	end
	if m == size(ClusterRgns,2) % the last cluster
		for n = startInd:m
			ClusterRgns(n).ClusterNo = ClusterNo;
		end
		ClusterNo = ClusterNo+1;
	end
end

end

function [FillRate,AspectRatio] = combine(Regions,ClusterNo)
%combines the candidate REGIONS having the same CLUSTERNO (i.e., in the
% same cluster) in each segmented image, and returns the FILL RATE, and the
% ASPECT RATIO of this cluster.

FillRate = zeros(size(Regions,2),ClusterNo);
AspectRatio = zeros(size(Regions,2),ClusterNo);

tmpArea   = zeros(size(Regions,2),ClusterNo);
tmpLeft   = zeros(size(Regions,2),ClusterNo);
tmpTop    = zeros(size(Regions,2),ClusterNo);
tmpRight  = zeros(size(Regions,2),ClusterNo);
tmpBottom = zeros(size(Regions,2),ClusterNo);
tmpWidth  = zeros(size(Regions,2),ClusterNo);
tmpHeight = zeros(size(Regions,2),ClusterNo);
for m = 1:size(Regions,2)
	for n = 1:size(Regions(m).Regions,2)
		tmpClusterNo = Regions(m).Regions(n).ClusterNo;
		tmpArea(m,tmpClusterNo) = tmpArea(m,tmpClusterNo)+Regions(m).Regions(n).Props.FilledArea;
		if tmpLeft(m,tmpClusterNo)==0
			tmpLeft(m,tmpClusterNo) = Regions(m).Regions(n).Props.BoundingBox(1);
		elseif tmpLeft(m,tmpClusterNo) > Regions(m).Regions(n).Props.BoundingBox(1)
			tmpLeft(m,tmpClusterNo) = Regions(m).Regions(n).Props.BoundingBox(1);
		end
		if tmpTop(m,tmpClusterNo)==0
			tmpTop(m,tmpClusterNo) = Regions(m).Regions(n).Props.BoundingBox(2);
		elseif tmpTop(m,tmpClusterNo) > Regions(m).Regions(n).Props.BoundingBox(2)
			tmpTop(m,tmpClusterNo) = Regions(m).Regions(n).Props.BoundingBox(2);
		end
		if tmpRight(m,tmpClusterNo)==0
			tmpRight(m,tmpClusterNo) = Regions(m).Regions(n).Props.BoundingBox(1)+Regions(m).Regions(n).Props.BoundingBox(3);
		elseif tmpRight(m,tmpClusterNo) < Regions(m).Regions(n).Props.BoundingBox(1)+Regions(m).Regions(n).Props.BoundingBox(3)
			tmpRight(m,tmpClusterNo) = Regions(m).Regions(n).Props.BoundingBox(1)+Regions(m).Regions(n).Props.BoundingBox(3);
		end
		if tmpBottom(m,tmpClusterNo)==0
			tmpBottom(m,tmpClusterNo) = Regions(m).Regions(n).Props.BoundingBox(2)+Regions(m).Regions(n).Props.BoundingBox(4);
		elseif tmpBottom(m,tmpClusterNo) < Regions(m).Regions(n).Props.BoundingBox(2)+Regions(m).Regions(n).Props.BoundingBox(4)
			tmpBottom(m,tmpClusterNo) = Regions(m).Regions(n).Props.BoundingBox(2)+Regions(m).Regions(n).Props.BoundingBox(4);
		end
	end

	for n = 1:ClusterNo
		tmpWidth(m,n) = tmpRight(m,n)-tmpLeft(m,n);
		tmpHeight(m,n) = tmpBottom(m,n)-tmpTop(m,n);
		if tmpWidth(m,n)*tmpHeight(m,n)~=0
			FillRate(m,n) = tmpArea(m,n)/(tmpWidth(m,n)*tmpHeight(m,n));
		end
		if tmpHeight(m,n)~=0
			AspectRatio(m,n) = tmpWidth(m,n)/tmpHeight(m,n);
		end
	end
end

end

function HitRgnsFiltered = filterrgn(Candidate,HitRgnsMask)

tmpThreshNo = 1;
for m = 1:size(Candidate,2)
	tmpRgnNo = 1;
	needMod = true;
	for n = 1:size(Candidate(m).Regions,2)
		tmpClusterNo = Candidate(m).Regions(n).ClusterNo;
		for k = 1:size(HitRgnsMask,2)
			if k==tmpClusterNo && HitRgnsMask(k)==1
				if needMod == true
					HitRgnsFiltered(tmpThreshNo).Threshold = Candidate(m).Threshold;
					tmpThreshNo = tmpThreshNo+1;
					needMod = false;
				end
				HitRgnsFiltered(tmpThreshNo-1).Regions(tmpRgnNo) = Candidate(m).Regions(n);
				tmpRgnNo = tmpRgnNo+1;
			end
		end
	end
end
end

function [isEnclose,idx] = enclose(region1,region2)
% checks whether REGION1 is enclosed in REGION2 or not by exploiting
% their bounding boxes, and returns the corresponding ID (1 or 2) of the bigger one.

isEnclose = false;

area1 = region1.Props.BoundingBox(3)*region1.Props.BoundingBox(4);
area2 = region2.Props.BoundingBox(3)*region2.Props.BoundingBox(4);
if area1 < area2
	regionSmall = region1;
	regionBig = region2;
	idx = 2;
else
	regionSmall = region2;
	regionBig = region1;
	idx = 1;
end

rectSmall   = regionSmall.Props.BoundingBox;
leftSmall   = rectSmall(1);
topSmall    = rectSmall(2);
widthSmall  = rectSmall(3);
heightSmall = rectSmall(4);

rectBig   = regionBig.Props.BoundingBox;
leftBig   = rectBig(1);
topBig    = rectBig(2);
widthBig  = rectBig(3);
heightBig = rectBig(4);

if leftSmall>=leftBig && topSmall>=topBig ...
		&& leftSmall+widthSmall<=leftBig+widthBig ...
		&& topSmall+heightSmall<=topBig+heightBig ...
		&& ~(leftSmall==leftBig && topSmall==topBig && leftSmall+widthSmall==leftBig+widthBig && topSmall+heightSmall==topBig+heightBig)
	isEnclose = true;
end

end

function PoStaRgns = otsu(gray, param, PoStaRgns, region, ClusterNo)
% computes the Otsu's thresholding value of the subimage in a
% gray-scale image GRAY. The subimage is extracted from GRAY by exploiting
% the bounding box of the maximally stable region REGION. If REGION
% satisfies the constraint of area variation, it is added to POSTARNGS with its properties.

[height,width] = size(gray);

isClusterMatch = false;
if region.ClusterNo == ClusterNo
	rect    = region.Props.BoundingBox;
	rleft   = uint16(rect(1));
	rtop    = uint16(rect(2));
	rwidth  = uint16(rect(3));
	rheight = uint16(rect(4));
	rright  = uint16(min(rleft+rwidth,width));
	rbottom = uint16(min(rtop+rheight,height));
	cropimg = gray(rtop:rbottom, rleft:rright);		% subimage
	isClusterMatch = true;
end

if isClusterMatch == true
	% segment using Otsu's method
	OtsuThresh = graythresh(cropimg);
	[hasRgnLow,bwLow] = convertrgn(cropimg, OtsuThresh-param.delta/255/2);
	[hasRgnUpp,bwUpp] = convertrgn(cropimg, OtsuThresh+param.delta/255/2);
	
	% measure area variation \Delta_r
	if hasRgnLow==true && hasRgnUpp==true
		[areaLow,~,~] = longrgn(bwLow);
		[areaUpp,~,~] = longrgn(bwUpp);
		if areaLow>=param.size && abs(areaLow-areaUpp)<=param.Delta_r*areaUpp ...
				|| areaLow<param.size && abs(areaLow-areaUpp)<=param.Delta_r*param.size
			PoStaRgns.RegionNums = PoStaRgns.RegionNums+1;
			tmpNums = PoStaRgns.RegionNums;
			[~,bw] = convertrgn(cropimg,OtsuThresh);
			[~,boundary,pixellist] = longrgn(bw);
			PoStaRgns.Regions(tmpNums).ClusterNo = ClusterNo;
			
			% Boundary
			for p = 1:size(boundary,1)
				tmpx = boundary(p,2)+rleft-1;
				if tmpx<1
					tmpx = 1;
				elseif tmpx>width
					tmpx = width;
				end
				tmpy = boundary(p,1)+rtop-1;
				if tmpy<1
					tmpy = 1;
				elseif tmpy>height
					tmpy = height;
				end
				boundary(p,2) = tmpx;
				boundary(p,1) = tmpy;
			end
			PoStaRgns.Regions(tmpNums).Props.Boundary = boundary;
			
			% Pixel List
			for p = 1:size(pixellist,1)
				tmpx = pixellist(p,1)+rleft-1;
				if tmpx<1
					tmpx = 1;
				elseif tmpx>width
					tmpx = width;
				end
				tmpy = pixellist(p,2)+rtop-1;
				if tmpy<1
					tmpy = 1;
				elseif tmpy>height
					tmpy = height;
				end
				pixellist(p,1) = tmpx;
				pixellist(p,2) = tmpy;
			end
			PoStaRgns.Regions(tmpNums).Props.PixelList = pixellist;
			
			% Bounding Box
			isFirst = true;
			for p = 1:size(boundary,1)
				if isFirst == true
					minx = boundary(p,2);
					maxx = boundary(p,2);
					miny = boundary(p,1);
					maxy = boundary(p,1);
					isFirst = false;
				else
					if minx > boundary(p,2); minx = boundary(p,2); end
					if maxx < boundary(p,2); maxx = boundary(p,2); end
					if miny > boundary(p,1); miny = boundary(p,1); end
					if maxy < boundary(p,1); maxy = boundary(p,1); end
				end
			end
			if minx>=2; minx=minx-1; end
			if maxx<=width-1; maxx=maxx+1; end
			if miny>=2; miny=miny-1; end
			if maxy<=height-1; maxy=maxy+1; end
			rect = [minx,miny,maxx-minx,maxy-miny];
			PoStaRgns.Regions(tmpNums).Props.BoundingBox = rect;
		end
	end
end

end

function [hasRgn,bw] = convertrgn(gray,threshold)
%converts a gray-scale image GRAY to a binary image BW by exploiting THRESHOLD, 
% and checks whether BW contains foreground regions.

bw = im2bw(gray,threshold);
bw = imcomplement(bw);
bw = imfill(bw,'holes');
tmpBoundaries = bwboundaries(bw);

if size(tmpBoundaries,1) >= 1
	hasRgn = true;
else
	hasRgn = false;
end

end

function [area,boundary,pixellist] = longrgn(bw)
%returns the properties (AREA, BOUNDARY, and PIXELLIST) of the
% region which has the longest boundary in the binary image BW.

tmpBoundaries = bwboundaries(bw);
isFirst = true;
for t = 1:size(tmpBoundaries,1)
	if isFirst == true
		boundMax = size(tmpBoundaries{t},1);
		indexMax = t;
		isFirst = false;
	elseif boundMax < size(tmpBoundaries{t},1)
		boundMax = size(tmpBoundaries{t},1);
		indexMax = t;
	end
end
boundary = tmpBoundaries{indexMax};		% Boundary

% regional properties
bw = false(size(bw,1),size(bw,2));
for p = 1:size(boundary,1)
	tmpX = boundary(p,2);
	tmpY = boundary(p,1);
	bw(tmpY,tmpX) = true;
end
Props = regionprops(bw,'FilledArea');
area = Props.FilledArea;				% Area
bw = imfill(bw,'holes');
Props = regionprops(bw,'PixelList');
pixellist = Props.PixelList;			% Pixel List

end

function measure = resemb(FillRate,AspectRatio)

measure = zeros(10,size(FillRate,2));

for m = 1:size(FillRate,1)-1
	n = m+1;
	for k = 1:size(FillRate,2)
		
		% similarity of Fill Rate
		tmpFillRateMax = max(FillRate(m,k),FillRate(n,k));
		tmpFillRateMin = min(FillRate(m,k),FillRate(n,k));
		if tmpFillRateMin ~= 0
			tmpFillFactorSta = tmpFillRateMax / tmpFillRateMin;
		else
			tmpFillFactorSta = 0;
		end
		if tmpFillFactorSta ~= 0
			if measure(1,k) == 0
				measure(1,k) = tmpFillFactorSta;
				measure(2,k) = m;
				measure(3,k) = n;
			elseif measure(1,k) >= tmpFillFactorSta
				measure(1,k) = tmpFillFactorSta;
				measure(2,k) = m;
				measure(3,k) = n;
			end
		end
		
		% similarity of Aspect Ratio
		tmpAspectRatioMax = max(AspectRatio(m,k),AspectRatio(n,k));
		tmpAspectRatioMin = min(AspectRatio(m,k),AspectRatio(n,k));
		if tmpAspectRatioMin ~= 0
			tmpAspectRatioSta = tmpAspectRatioMax / tmpAspectRatioMin;
		else
			tmpAspectRatioSta = 0;
		end
		if tmpAspectRatioSta ~= 0
			if measure(4,k) == 0
				measure(4,k) = tmpAspectRatioSta;
				measure(5,k) = m;
				measure(6,k) = n;
			elseif measure(4,k) >= tmpAspectRatioSta
				measure(4,k) = tmpAspectRatioSta;
				measure(5,k) = m;
				measure(6,k) = n;
			end
		end
		if measure(1,k)~=0 && measure(4,k)~=0
			measure(7,k) = min([measure(2,k), measure(3,k), measure(5,k), measure(6,k)]);
			measure(8,k) = max([measure(2,k), measure(3,k), measure(5,k), measure(6,k)]);
		end
	end
end

end

function isSameCluster = member(region1,region2)
%member checks whether the regions REGION1 and REGION2 are belonged to the
% same cluster based on the center distance of them.

isSameCluster = false;

center1 = region1.Props.Centroid;
width1  = region1.Props.BoundingBox(3);
height1 = region1.Props.BoundingBox(4);

center2 = region2.Props.Centroid;
width2  = region2.Props.BoundingBox(3);
height2 = region2.Props.BoundingBox(4);

widthMin = min(width1,width2);
heightMin = min(height1,height2);

dist = (center2(1)-center1(1))^2 + (center2(2)-center1(2))^2;
if dist <= widthMin^2/4 + heightMin^2/4
	isSameCluster = true;
end

end
