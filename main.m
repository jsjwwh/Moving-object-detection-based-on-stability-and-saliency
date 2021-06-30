clc; clear; close all;

%% Datasets
data = 'Data3';				% Data1, Data2, or Data3

%% Parameters
param.size     = 100;		% small target size
param.delta    = 16;		% sample step \delta
param.Delta_r  = 0.2;		% area variation \Delta_r
param.sigma_s  = 16;		% filter parameter \sigma_s
% target area (fixed value)
param.minarea  = 4;			% pixels
param.maxarea  = 0.2;		% ratio

% make folders
if exist([data,'\StaMaps'], 'dir') ~= 7		% stability maps (.png)
	system(['md ',data,'\StaMaps']);
end
if exist([data,'\SalMaps'], 'dir') ~= 7		% saliency maps (.png)
	system(['md ',data,'\SalMaps']);
end
if exist([data,'\STD'], 'dir') ~= 7			% detected small targets (.mat and .png)
	system(['md ',data,'\STD']);
end

imgs = dir([data,'\Image\*.png']);
for num = 1:length(imgs)
	fprintf('  %3d/%3d : ', num, length(imgs));
	
	tic;
	
	rgb = imread([data,'\Image\',int2str(num),'.png']);
	if ndims(rgb) == 2		% for gray-scale image
		rgb = repmat(rgb, [1 1 3]);
	end
	gray = rgb2gray(rgb);
	
	% (1) extract stability regions, and generate stability map
	StabilityRegions = rgnstabmap(gray, param);
	StabilityMap = binary(StabilityRegions,size(gray,1),size(gray,2));
	imwrite(StabilityMap, [data,'\StaMaps\',int2str(num),'_Sta.png']);
	
	% (2) detect saliency, and generate saliency map
	SaliencyMap = rgnsalmap(rgb, param);
	imwrite(SaliencyMap, [data,'\SalMaps\',int2str(num),'_Sal.png']);
	
	% (3) obtain the final detection result by combining stability map and saliency map
	Targets = STD(StabilityRegions, StabilityMap, SaliencyMap);
	save([data,'\STD\',int2str(num),'_STD.mat'],'Targets');
	TargetsMap = binary(Targets,size(gray,1),size(gray,2));
	imwrite(TargetsMap,[data,'\STD\',int2str(num),'_STD.png']);
	
	toc;
	
end