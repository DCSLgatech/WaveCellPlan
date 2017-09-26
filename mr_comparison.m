%{
Copyright (c) 2017 Raghvendra V. Cowlagi. All rights reserved.

Copyright notice: 
=================
No part of this work may be reproduced without the written permission of
the copyright holder, except for non-profit and educational purposes under
the provisions of Title 17, USC Section 107 of the United States Copyright
Act of 1976. Reproduction of this work for commercial use is a violation of
copyright.


Disclaimer:
===========
This software program is intended for educational and research purposes.
The authors and the institutions with which the authors are affiliated are not
liable for damages resulting the application of this program, or any
section thereof, which may be attributed to any errors that may exist in
this program.


Author information:
===================
1. Raghvendra V. Cowlagi
Assistant Professor,
Aerospace Engineering Program,
Worcester Polytechnic Institute, Worcester, MA.
Email: rvcowlagi@wpi.edu
Website: http://www.wpi.edu/~rvcowlagi

2. Panagiotis Tsiotras
College of Engineering Dean's Professor,
School of Aerospace Engineering,
Georgia Institute of Technology, Atlanta, GA.
Email: tsiotras@gatech.edu
Website: http://soliton.ae.gatech.edu/labs/ptsiotra/

Description: Comparisons between different multiresolution methods
%}

clear all; close all; clc;

%% Resolution, cost parameters
simData = [];
allWindws = [1 1 2 2 2 2 2 2 3 3 4 4; 1 2 2 3 4 5 6 7 7 7 8 8; ...
	3 3 5 5 6 7 8 8 9 9 10 10];
nObsList = [zeros(1, 6) 100 150 200 250 350 450];
for jmin = -10
	
%% Environment 3: Cluttered
% fprintf('\nConfiguring environment map ...\t');
% tic
% nColors = 256;
% map		= gray(nColors);
% envSize = 2^(-jmin);
% 
% X		= nColors*ones(envSize);
% maxObsSize = 0.05*envSize;
% nObs	= nObsList(-jmin);
% for m1 = 1:nObs
% 	leftBottomX = 1 + round((0.95*envSize - 1)*rand);
% 	leftBottomY = 1 + round((0.95*envSize - 1)*rand);
% 	obsDimX		= round(maxObsSize*rand);
% 	obsDimY		= round(maxObsSize*rand);
% 	X(leftBottomX:(leftBottomX + obsDimX), ...
% 		leftBottomY:(leftBottomY + obsDimY)) = 0;
% end


% Xtrue	= X;
% toc

% figure('Position', [-1400 0 800 600]); endAxes = axes; hold on; grid on; axis equal;

% plot(5, 5, 'ko', 'MarkerFaceColor', 'k'); hold on;

end

	for qThres = 5 %[5 20 30 50]
% 	for winFcnNo = 1:3
% 		fprintf('jmin = %i, window = %i\n', jmin, winFcnNo)
% jmin= -7;				% Coarsest possible cell has dimension 2^(-jmin)
jmax= 0;				% Finest res cell has dimension 2^(-jmax)
nV	= 2^(jmax-jmin);	% Total number of pixels in each row

%% Start and Goal
rerun = 0;%input('Rerun with same cell decomposition?		');
if rerun
	load Data/WaveletResults/mrComparisonSimV01Data.mat
	rerun = 1;
else
	pInit	= (2^(-jmin))*rand(2,1);
	pGoal	= (2^(-jmin))*rand(2,1);
	rerun2	= 0;	
end
% save Data/WaveletResults/mrComparisonSimV01DataV10a.mat

%% Environment 1: Terrain
% load('map_q1632.mat');
% % fprintf('\nConfiguring environment map ...\t');
% % tic
% X		= double(Img.cdata);
% X		= imresize(X,[nV nV],'bilinear');
% nColors = 256;
% map		= gray(nColors);
% 
% Xrgb	= 0.2990*X(:,:,1)+0.5870*X(:,:,2)+0.1140*X(:,:,3);
% Xdisp	= wcodemat(Xrgb, nColors);
% Xtrue	= Xdisp;
% % toc

%% Environment 2: Large obstacles
% fprintf('\nConfiguring environment map ...\t');
% tic
nColors = 256;
envSize = 2^(-jmin);
map		= gray(nColors);

X = nColors*ones(2^(-jmin));

for m1 = round(0.1*envSize):round(0.3*envSize)
	for m2 = 1:round(0.5*envSize)
		if (m2 >= -0.5*m1 + 0.1*envSize) && (m2 >= 0.5*m1) ...
				&& (m2 <= -1.2*m1 + 0.5*envSize)
			X(m1, m2) = 0;
		end
	end
end

for m1 = round(0.5*envSize):round(0.7*envSize)
	for m2 = 1:envSize
		if (m2 >= -0.5*m1 + 0.5*envSize) && (m2 <= 0.5*m1 + 0.5*envSize)
			X(m1, m2) = 0;
		end
	end
end

X(round(0.8*envSize):envSize, round(0.2*envSize):round(0.3*envSize)) = 0;
X(round(0.73*envSize):round(0.92*envSize), round(0.7*envSize):round(0.8*envSize)) = 0;

for m1 = round(0.2*envSize):round(0.4*envSize)
	for m2 = round(-sqrt((0.1*envSize)^2 - (m1 - 0.3*envSize)^2) ...
			+ 0.5*envSize):round(sqrt((0.1*envSize)^2 - (m1 - 0.3*envSize)^2) + 0.5*envSize)
		X(m1, m2) = 0;
	end
end

X(round(0.1*envSize):round(0.4*envSize), round(0.75*envSize):round(0.9*envSize)) = 0;

Xtrue = X;
% toc

figure; endAxes = axes; hold on; grid on; axis equal;
image(Xtrue, 'XData', 0.5, 'YData', 0.5); colormap(map);

%% Environment 3: Cluttered
% % fprintf('\nConfiguring environment map ...\t');
% % tic
% nColors = 256;
% map		= gray(nColors);
% envSize = 2^(-jmin);
% 
% X		= nColors*ones(envSize);
% maxObsSize = 0.05*envSize;
% nObs	= round((-jmin)*5*rand) + round((-jmin)*15*rand);
% for m1 = 1:nObs
% 	leftBottomX = 1 + round((0.95*envSize - 1)*rand);
% 	leftBottomY = 1 + round((0.95*envSize - 1)*rand);
% 	obsDimX		= round(maxObsSize*rand);
% 	obsDimY		= round(maxObsSize*rand);
% 	X(leftBottomX:(leftBottomX + obsDimX), ...
% 		leftBottomY:(leftBottomY + obsDimY)) = 0;
% end
% 
% 
% Xtrue	= X;
% % toc
% 
% figure('Position', [-1400 0 800 600]); endAxes = axes; hold on; grid on; axis equal;
% plot(pInit(1), pInit(2), 'ko', 'MarkerFaceColor', 'k'); hold on;
% 
% image(Xtrue, 'XData', 0.5, 'YData', 0.5); colormap(map);
% return

%% Initial wavelet decomposition and MR approximation
% % fprintf('Wavelet decomposition ...\t\t');
% tic
% N =	-jmin;																	% Coarsest level of decomposition
% [Cforig, Sz]	= wavedec2(Xtrue, N, 'db1');								% Compute wavelet decomposition of original map
% wavtmp1 = toc;
% 
% % ------ Specify window function
% windw = allWindws(winFcnNo,:);
% 
% 
% % fprintf('First decomposition ...\t\t\t');
% tic
% [Cf, NzrData]	= mrDecV05(Cforig, Sz, jmax, pInit, windw);					% MR approximation
% wavtmp2 = toc;

%% Cells at finest resolution
% nVFine	= Sz(end,1)^2;
% VFine	= zeros(nVFine, 4);
% m		= 0;
% for cX = 1:Sz(end,1)
% 	for cY = 1:Sz(end,1)
% 		m		= m + 1;
% 		VFine(m,:)	= [cX-1 cY-1 1 Xtrue(cY, cX)];
% 	end
% end

%% Initial cell and adjacency computation
% fprintf('Cell and adjacency, v17 ...\t\t');
% tic
% % [G1, VCell]	= adjMat17(NzrData.ANzr, Cf, Sz);								% V from Cf, only once; G is a matrix
% [G, VCell]	= adjMat18(NzrData.ANzr, Cf, Sz);								% V from Cf, only once; G is a matrix
% wavtmp3 = toc;

%% Quadtree decomposition
fprintf('Quadtree decomposition ...\t\t');
tic
S = qdecomp(Xtrue, qThres);
tmp1 = toc;
fprintf('Threshold = %f, \t Number of vertices = %i, \t EnvSize = %i.\n\n', qThres, numel(S.locs), 2^(-2*jmin))

% fprintf('Quadtree adjacency ...\t\t\t');
tic
% C = qlink(S);
tmp2 = toc;

simData = cat(1, simData, [2^(-jmin) qThres numel(S.locs) tmp1 tmp2]);
% simData = cat(1, simData, [2^(-jmin) winFcnNo size(VCell, 1) wavtmp1 wavtmp2 wavtmp3]);
% 	end
	end
% end