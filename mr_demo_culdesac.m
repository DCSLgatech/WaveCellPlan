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
1. Raghvendra V. Cowlagi, Ph.D,
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

Description: "Demonstration" of completeness
%}
clear variables; close all; clc;
format short g

%% Resolution, cost parameters
jmin= -7;				% Coarsest possible cell has dimension 2^(-jmin)
jmax= 0;				% Finest res cell has dimension 2^(-jmax)
nV	= 2^(jmax-jmin);	% Total number of pixels in each row

%% Start and Goal
pInit	= [16.496; 2.3193];
pGoal	= [127; 64];

%% Terrain map
load('map_q1632.mat');
fprintf('\nConfiguring environment map ...\t');
tic
nColors = 256;
map		= gray(nColors);

X = nColors*ones(128);
for nY = 1:128
	tmp = 7 + round(abs(nY - 64)*10/64);
% 	tmp = 10;
	X(nY, 100:(100+tmp)) = 0;
end
X(125:128, :)		= nColors;
X(62:64, :)			= nColors;
X(65:66, 101:120)	= 0;
X(62:63, 100:120)	= 0;
X(62:66, 121:122)	= 0;
Xtrue = X;
toc

%% Initial wavelet decomposition and MR approximation
fprintf('Wavelet decomposition ...\t\t');
tic
N =	-jmin;																	% Coarsest level of decomposition
[Cforig, Sz]	= wavedec2(Xtrue, N, 'db1');								% Compute wavelet decomposition of original map
toc

% ------ Specify window function
windw = [1 1 2 2 2 2 2 2];													% Coarse to fine

fprintf('First decomposition ...\t\t\t');
tic
[Cf, NzrData]	= mrDecV05(Cforig, Sz, jmax, pInit, windw);					% MR approximation
toc

%% Cells and graph at finest resolution
fprintf('Finest resolution graph ...\t\t');
tic
nVFine	= Sz(end,1)^2;
VFine	= zeros(nVFine, 4);
m		= 0;
for cX = 1:Sz(end,1)
	for cY = 1:Sz(end,1)
		m		= m + 1;
		VFine(m,:)	= [cX-1 cY-1 1 Xtrue(cY, cX)];
	end
end

JFine	= Inf(nVFine,1);
visitRec= zeros(nVFine, 1);
backPtr	= zeros(nVFine, 1);

nEdges	= 0;
nExpEd	= nVFine*4;
edgeList= zeros(nExpEd, 3);
for m = 1:nVFine
	if (m + 1 <= nVFine) && (mod(m, Sz(end,1)) ~= 0)
		nEdges				= nEdges + 1;
		edgeList(nEdges, :) = [m (m + 1) 1];
		nEdges				= nEdges + 1;
		edgeList(nEdges, :) = [(m + 1) m 1];
	end

	if (m + Sz(end,1)) <= nVFine
		nEdges				= nEdges + 1;
		edgeList(nEdges, :) = [m (m + Sz(end,1)) 1];
		nEdges				= nEdges + 1;
		edgeList(nEdges, :) = [(m + Sz(end,1)) m 1];
	end
end
GFine = sparse(edgeList(1:nEdges,1), edgeList(1:nEdges,2), edgeList(1:nEdges,3));

pInitCellFine	= floor(pInit);
pGoalCellFine	= floor(pGoal);
[~,nodeSFine]	= ismember(pInitCellFine', VFine(:,1:2), 'rows');
[~,nodeGFine]	= ismember(pGoalCellFine', VFine(:,1:2), 'rows');
GFine			= adj2cost07(GFine, VFine, nodeGFine);

toc

%% Initial cell and adjacency computation
fprintf('Cell and adjacency, v17 ...\t\t');
tic
[G, VCell]	= adjMat17(NzrData.ANzr, Cf, Sz);								% V from Cf, only once; G is a matrix
toc

%------ Define heuristic, transition cost
heur= zeros(size(VCell,1),1);

for j = (-N):jmax
	pInitCell	= floor(pInit*(2^(j)))*(2^(-j));
	[tmp, nodeS]= ismember([pInitCell' 2^(-j)], VCell(:,1:3), 'rows');
	if tmp, break; end
end
for j = (-N):jmax
	pGoalCell	= floor(pGoal*(2^(j)))*(2^(-j));
	[tmp, nodeG]= ismember([pGoalCell' 2^(-j)], VCell(:,1:3), 'rows');
	if tmp, break; end
end
nVCell	= size(VCell,1);
G		= adj2cost07(G, VCell, nodeG);

figure('Position', [300 18 1200 950]); trueAxes = axes; axis equal; hold on;
plot(trueAxes, pInit(1), pInit(2), 'wo', 'MarkerFaceColor', 'w'); hold on;
image(Xtrue, 'XData', 0.5, 'YData', 0.5); colormap(map);
set(trueAxes, 'XLim', [0 2^(-jmin)], 'YLim', [0 2^(-jmin)], 'XTick', ...
	[], 'YTick', [], 'XColor', 'w', 'YColor', 'w');
plot(pInit(1), pInit(2), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
plot(pGoal(1), pGoal(2), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);

%----- A* to plan first path
fprintf('Initial path ...\t\t\t\t');
tic
nodeData	= aStarV03(G, nodeS, nodeG, heur);
% [nodeData, ~] = aStar(nodeS, nodeG, heur, [], G, [], [], 1);
mrSPath			= traceGreedyV05(nodeS, nodeG, nodeData);
visitRec(nodeSFine) = 1;
toc

figure('Position', [300 18 1200 950]); envAxes = axes; axis equal; hold on;
% figure; envAxes = axes; axis equal; hold on;
drawCdV05(VCell, envAxes, 'k', 1, 0, [], 0, 2^(-jmin), map)
drawCdV05(VCell, envAxes, 'b', 3, 0, mrSPath)
drawCdV05(VCell, envAxes, 'r', 3, 0, nodeG)
drawCdV05(VCell, envAxes, 'r', 2, 0, nodeS)
set(envAxes, 'XLim', [0 2^(-jmin)], 'YLim', [0 2^(-jmin)], 'XTick', ...
	[], 'YTick', [], 'XColor', 'w', 'YColor', 'w');

%% Further "incremental" planning
fprintf('Further path planning ...\t\t\n');
VAct	= VCell;	GAct	= G;	NzrAct	= NzrData;						% Multi-resolution graph stuff
pPtr	= 2;		JAct	= nodeData(nodeG).d;
nIter	= 0;		backStep= 0;											% No. iterations, time taken
XUp		= Xtrue;	CfUp	= Cforig;
VActDraw= VCell;


% envAvi = avifile('mrCompleteMovie01.avi');
envAvi = VideoWriter('mrCompleteMovie01.avi');
open(envAvi)
% trueAvi = avifile('mrCompleteMovie02.avi');

fineGraphPath	= [];
while VAct(nodeG, 3) > 1
% while nIter < 20
	nIter = nIter + 1;
	fprintf('\tIteration %i ...\t\t', nIter);
	
% % 	% -------- Movie stuff -------
	envFrame			= getframe(envAxes);
% 	trueFrame			= getframe(trueAxes);
% 	envMovie(nIter)		= envFrame;
% 	fineMovie(nIter)	= fineFrame;
% 	envAvi				= addframe(envAvi,envFrame);
    
    writeVideo(envAvi, envFrame);
% 	trueAvi				= addframe(trueAvi,trueFrame);
	
	% -------- Step forward, note dlta = (dx, dy) --------
	if ~backStep
		pAct		= (VAct(mrSPath(pPtr-1), 1:2))';
		pNext		= (VAct(mrSPath(pPtr), 1:2))';
		pNextSize	= VAct(mrSPath(pPtr), 3);								% Size of next cell (can force to be one)	
	else
		fprintf('backstepping... \t');
	end
	dlta		= pNext - pAct;
	
	% -------- For recording cost-to-go estimate --------
	[~,nodeActFine]	= ismember(pAct', VFine(:,1:2), 'rows');
	[~,nodeNextFine]= ismember(pNext', VFine(:,1:2), 'rows');
	backPtr(nodeNextFine) = nodeActFine;
	
	% -------- Record path and control --------
	fineGraphPath	= cat(2, fineGraphPath, nodeActFine);
	
	% -------- Update set of nzr detail coeffs --------
	pCells = [-N 0 0];
	for j = (-N+1):(jmax-1)													% From coarse to fine
		n = j + N + 1;
		posn	= floor((2^j)*pAct);
		pCells	= cat(1, pCells, [j posn']);
	end
	pCells = cat(1, pCells, [0 pAct']);
	NzrAct.pCells = pCells;	

	[NzrNext, inOnN, inNnO]	= getNewNzrV02(NzrAct, Sz, jmax, dlta, windw);
	% inOnN : set of nzr coefficients inOldnotNew
	% inNnO : set of nzr coefficients inNewnotOld

	if ~backStep && (numel(inOnN) == 0) && (numel(inNnO) == 0)				% In case decomposition doesn't change
		fprintf('continuing previous trajectory.\n');
		fprintf('Current location: %i\n', nodeNextFine)
		pPtr	= pPtr + 1;	
		continue;															% Basically, don't do anything, just step fwd
	end
	fprintf('\nCurrent location: %i\n', nodeNextFine)
	pPtr = 2;

	% -------- Update graph and set of cells --------
	if backStep
		[CfNext,NzrNext]= mrDecV05(CfUp, Sz, jmax, pNext, windw);			% MR approximation
		[GNext, VNext]	= adjMat17(NzrNext.ANzr, CfNext, Sz);				% V from Cf, only once; G is a matrix
		
		[CfDraw,~]	= mrDecV05(Cforig, Sz, jmax, pNext, windw);			% MR approximation
		[~, VDraw]	= adjMat17(NzrNext.ANzr, CfDraw, Sz);				% V from Cf, only once; G is a matrix
	else
		[GNext, VNext]	= getNewGraphV02(NzrAct.ANzr, GAct, ...
			VAct, inOnN, inNnO, CfUp, Sz);
		[~, VDraw]	= getNewGraphV02(NzrAct.ANzr, GAct, ...
			VActDraw, inOnN, inNnO, Cforig, Sz);
	end
	nVNext	= size(VNext, 1);
	
	% -------- Locate current cell and goal in new decomposition --------
	[~, nodeNext]	= ismember([pNext' pNextSize], VNext(:, 1:3), 'rows');
	for j = (-N):jmax
		pGoalCell	= floor(pGoal*(2^(j)))*(2^(-j));
		[tmp, nodeG]= ismember([pGoalCell' 2^(-j)], VNext(:,1:3), 'rows');
		if tmp, break; end
	end
	
	[~, nodeNextDraw]	= ismember([pNext' pNextSize], VDraw(:, 1:3), 'rows');
	for j = (-N):jmax
		pGoalCellDraw	= floor(pGoal*(2^(j)))*(2^(-j));
		[tmp, nodeGDraw]= ismember([pGoalCellDraw' 2^(-j)], VDraw(:,1:3), 'rows');
		if tmp, break; end
	end
	
	
 	% -------- Redefine heuristic for changed nodes --------
	heur = zeros(nVNext,1);
	
	% -------- Draw things --------
	cla(envAxes); hold on; grid on; axis equal;
	set(envAxes, 'XLim', [0 2^(-jmin)], 'YLim', [0 2^(-jmin)], 'XTick', ...
		[], 'YTick', [], 'XColor', 'w', 'YColor', 'w');	
% 	set(envAxes, 'XLim', [0 2^(-jmin)], 'YLim', [0 2^(-jmin)], 'XTick', ...
% 		[], 'YTick', []);	
	drawCdV05(VDraw, envAxes, 'k', 1, 0, [], 0, 2^(-jmin), map)
	drawCdV05(VDraw, envAxes, 'r', 3, 0, nodeGDraw)
	drawCdV05(VDraw, envAxes, 'r', 3, 0, nodeNextDraw)
	drawnow;
	
	
	% -------- Draw things --------
% 	cla(trueAxes); hold on; grid on; axis equal;	
% 	plot(trueAxes, pInit(1), pInit(2), 'wo', 'MarkerFaceColor', 'w'); hold on;
% 	image(Xtrue, 'XData', 0.5, 'YData', 0.5); colormap(map);
% 	set(trueAxes, 'XLim', [0 2^(-jmin)], 'YLim', [0 2^(-jmin)], 'XTick', ...
% 		[], 'YTick', [], 'XColor', 'w', 'YColor', 'w');
% 	plot(trueAxes,pInit(1), pInit(2), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
% 	plot(trueAxes,pGoal(1), pGoal(2), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
% 	plot(trueAxes,pAct(1), pAct(2), 'bs', 'MarkerFaceColor', 'b', 'MarkerSize', 10);
% 	drawnow;
	
	
% 	cla(trueAxes); hold on; grid on; axis equal;	
% 	set(trueAxes, 'XLim', [0 2^(-jmin)], 'YLim', [0 2^(-jmin)], 'XTick', ...
% 		[], 'YTick', [], 'XColor', 'w', 'YColor', 'w');
% 	plot(trueAxes, pInit(1), pInit(2), 'ko', 'MarkerFaceColor', 'k');
% % 	image(Xtrue, 'XData', 0.5, 'YData', 0.5); colormap(map);
% 	drawCdV05(VDraw, trueAxes, 'r', 2, 0, nodeNextDraw);
% 	plot(trueAxes, pInit(1), pInit(2), 'ko', 'MarkerFaceColor', 'k');
% 	plot(trueAxes, pGoal(1), pGoal(2), 'ko', 'MarkerFaceColor', 'k');
% 	drawnow;
	
	
	% -------- Run search again --------
	GNext			= adj2cost07(GNext, VNext, nodeG);
	pBackPtrNextFine= VFine(backPtr(nodeNextFine), 1:2)';
	[~, backPtrNext]= ismember([pBackPtrNextFine' 1], VNext(:, 1:3), 'rows');
	GNext(nodeNext, backPtrNext) = 0;	GNext(backPtrNext, nodeNext) = 0;
	
	nodeData= aStarV03(GNext, nodeNext, nodeG, heur);
	mrSPath = traceGreedyV05(nodeNext, nodeG, nodeData);
% 	JFine(nodeNextFine)	= nodeData(nodeG).d;
% 	JAct				= cat(1, JAct, nodeData(nodeG).d);
	mrSPathCost = nodeData(nodeG).d;
	if mrSPathCost >= 1e8 % % % obstacle cost
		pAct	= (VNext(nodeNext, 1:2))';
		pNext	= (VNext(backPtrNext, 1:2))';
% 		CfUp	= wvlCfUpdateV01(CfUp, Sz, pAct-[1 1]', -nColors, jmax);
		tmp1 = pAct + [1; 1];
		XUp(tmp1(2), tmp1(1)) = 0;
		[CfUp, Sz]	= wavedec2(XUp, N, 'db1');						% Compute wavelet decomposition of updated map
		backStep= 1;
	else
		backStep= 0;
	end
	
	visitRec(nodeNextFine)= visitRec(nodeNextFine) + 1;
	if visitRec(nodeNextFine) > 1
		tmp2 = VFine(nodeNextFine, 1:2)' + [1;1];
		tmp3 = min(XUp(tmp2(2), tmp2(1)), max(nColors - visitRec(nodeNextFine)*50, 10));
		
		XUp(tmp2(2), tmp2(1)) = tmp3;
		[CfUp, Sz]	= wavedec2(XUp, N, 'db1');
		VNext(nodeNext, 4) = tmp3;
	end
	
	% -------- Trace path and draw
	if ~backStep
% 		drawCdV05(VFine, envAxes, 'k', 2, 0, fineGraphPath);
		drawCdV05(VNext, envAxes, 'b', 3, 0, mrSPath(2:(end-1)))
% 		drawCdV05(VNext, trueAxes, 'b', 2, 0, mrSPath(2:end))
		drawnow;
	end
	
	% -------- Update everything --------
	GAct	= GNext;
	VAct	= VNext;
	NzrAct	= NzrNext;
	VActDraw= VDraw;
end
close(envAvi);
% trueAvi = close(trueAvi);

%%
% % figure('Position', [300 18 1200 950]); endAxes = axes;
% figure; endAxes = axes;
% hold on; grid on; axis equal;
% % drawCdV05(VFine(:, [1:3 8]), endAxes, 'k', 0.01, 0, [], 0, 2^(-jmin), map)
% plot(pInit(1), pInit(2), 'ko', 'MarkerFaceColor', 'k'); hold on;
% 
% image(Xtrue, 'XData', 0.5, 'YData', 0.5); colormap(map);
% drawCdV05(VFine, envAxes, 'k', 2, 0, fineGraphPath)

