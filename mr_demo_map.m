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

Description: Multi-resolution, multi-history. Target sets. MPC.
First cut of complete simulation. Direct H-cost search instead of explicit
lifted graph construction; aircraft model, bank only control
%}
clear variables; close all; clc;
% format short g
format long g
global envAxes fineAxes g vCruise
g = 9.8;																	% gravitational acceleration, m/s^2

%% Aircraft model parameters
altitude	= 1e3;															% m
temperature = 15.04 - 0.00649*altitude + 273.1;								% Kelvin
pressure	= 101.29*((temperature/288.08)^5.256);							% kPa
stdrhoD		= pressure/(0.2869*temperature);								% kg/m3

%----- Aircraft and flight parameters
% Example 2.11 from B. N. Pamadi
modelParams.mass	= 5e4/g;												% mass, kg
modelParams.CD0		= 0.02;													% skin friction drag coefficient, dim-less
modelParams.K		= 0.04;													% induced drag coefficient, dim-less
modelParams.S		= 30;													% wing area, m2
modelParams.rhoD	= stdrhoD;												% air density, kg/m3
modelParams.CLmax	= 1.2;													% max lift coefficient, dim-less

mass	= modelParams.mass;
CD0		= modelParams.CD0;
K		= modelParams.K;
S		= modelParams.S;
rhoD	= modelParams.rhoD;
CLmax	= modelParams.CLmax;

vStall	= sqrt(2*mass*g/(rhoD*S*CLmax));									% Stall speed, m/s
vCruise = 85;

%% Resolution, cost parameters
jmin= -7;				% Coarsest possible cell has dimension 2^(-jmin)
jmax= 0;				% Finest res cell has dimension 2^(-jmax)
nV	= 2^(jmax-jmin);	% Total number of pixels in each row

%% Start and Goal
% rerun = input('Rerun with same cell decomposition?		');
% if rerun
% 	load Data/mrMotionPlanTestV05Data.mat
% 	rerun = 1;
% else
% 	pInit	= (2^(-jmin))*rand(2,1);
% 	pGoal	= (2^(-jmin))*rand(2,1);
% 	rerun2	= 0;	
% end
pInit	= (2^(-jmin))*[0; 0.999999];
pGoal	= [80.1; 10]; %(2^(-jmin))*[0.99; 0];
% save Data/mrMotionPlanTestV05Data.mat

%% Terrain map
load('map_q1632.mat');
fprintf('\nConfiguring environment map ...\t');
tic
X		= double(Img.cdata);
X		= imresize(X,[nV nV],'bilinear');
nColors = 256;
map		= gray(nColors);

Xrgb	= 0.2990*X(:,:,1)+0.5870*X(:,:,2)+0.1140*X(:,:,3);
Xdisp	= wcodemat(Xrgb, nColors);
Xtrue	= Xdisp;
toc

%% Initial wavelet decomposition and MR approximation
fprintf('Wavelet decomposition ...\t\t');
tic
N =	-jmin;																	% Coarsest level of decomposition
[Cforig, Sz]	= wavedec2(Xtrue, N, 'db1');								% Compute wavelet decomposition of original map
toc

% ------ Specify window function
windw = [1 1 2 2 2 2 2 3];													% Coarse to fine
% windw = [1 1 2 2 2 2 3 4];													% Coarse to fine

fprintf('First decomposition ...\t\t\t');
tic
[Cf, NzrData]	= mrDecV05(Cforig, Sz, jmax, pInit, windw);					% MR approximation
toc

%% Cells at finest resolution
nVFine	= Sz(end,1)^2;
VFine	= zeros(nVFine, 4);
m		= 0;
for cX = 1:Sz(end,1)
	for cY = 1:Sz(end,1)
		m		= m + 1;
		VFine(m,:)	= [cX-1 cY-1 1 Xtrue(cY, cX)];
	end
end

%% Initial cell and adjacency computation
fprintf('Cell and adjacency, v17 ...\t\t');
tic
[G, VCell]	= adjMat17(NzrData.ANzr, Cf, Sz);								% V from Cf, only once; G is a matrix
toc

% Define heuristic, transition cost
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
nVCell		= size(VCell,1);
G	= adj2cost06(G, VCell, nodeG);

% Y = waverec2(Cf, Sz, 'db1');
figure('Position', [300 18 1200 950]); envAxes = axes; axis equal; hold on;
% plot(pInit(1), pInit(2), 'ko', 'MarkerFaceColor', 'k');

% image(Y, 'XData', 0.5, 'YData', 0.5); colormap(map);
drawCdV05(VCell, envAxes, 'k', 1, 0, [], 0, 2^(-jmin), map)
drawCdV05(VCell, envAxes, 'r', 3, 0, nodeG)
drawCdV05(VCell, envAxes, 'r', 2, 0, nodeS)

%% Multi-history lifted graph parameters
H		= 3;
nHMem	= 15;
dMax	= 2;
varDim.nStates = 4;
varDim.nInputs = 2;
 
vLim	= [vStall; 5*vStall]*dMax/2;
VFine	= cat(2, VFine, [ones(nVFine,2) kron(ones(nVFine,1), vLim') VFine(:,4)]);
VFine(:, 4) = [];

%% Boundary nodes and terminal penalty (for subsequent lifted search)
fprintf('Boundary nodes ...\t\t\t\t');
tic
bdData = fineResBoundaryV01(G, VCell, dMax);
toc

fprintf('Boundary terminal penalty ...\t');
tic
mrBdNodes = [];
for m = 1:numel(bdData)
	mrBdNodes	= cat(2, mrBdNodes, bdData(m).node);
end
mrNodeData = aStarV04(G', nodeG, 1:size(G,1), heur);							% A* to find cost-to-go from each boundary node to goal

for m = 1:numel(bdData)
	bdData(m).cost	= mrNodeData(bdData(m).node).d;
	tmp1			= traceGreedyV04(nodeG, bdData(m).node, mrNodeData);
	bdData(m).optP	= fliplr(tmp1);
end

toc
drawCdV05(VCell, envAxes, 'r', 2, 0, mrBdNodes)

%% Initial Condition Selection
cellScale = 1000;

VCell	= cat(2, VCell, [ones(nVCell,2) kron(ones(nVCell,1), vLim') VCell(:,4)]);
VCell(:, 4) = [];

fprintf('\n');
% rerun = input('Rerun with same initial condition?	');
% if rerun
% 	load Data/mrMotionPlanTestV06Data.mat z0
% else
% % 	z0	= [1; 127.5; -5*pi/180; (vLim(1) + (vLim(2) - vLim(1))*rand)];
	z0	= [1*cellScale; 127.5*cellScale; -5*pi/180; vCruise];
% end
% save Data/mrMotionPlanTestV08Data.mat

plot(envAxes, z0(1)/cellScale, z0(2)/cellScale, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 4);
quiver(envAxes, z0(1)/cellScale, z0(2)/cellScale, 0.2*cos(z0(3)), 0.2*sin(z0(3)), ...
	'Color', 'b', 'MaxHeadSize', 2, 'LineWidth', 2);

figure('Position', [100 18 800 600]); fineAxes = axes; hold on; grid on; axis equal
tmp1	= 1:size(VCell,1);
drawCdV05(VCell, fineAxes, 'k', 2, 0, tmp1(VCell(:,3) < 2))


drawCdV05(VCell, fineAxes, 'r', 3, 0, mrBdNodes)
plot(fineAxes, z0(1)/cellScale, z0(2)/cellScale, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 4);
quiver(fineAxes, z0(1)/cellScale, z0(2)/cellScale, 0.35*cos(z0(3)), 0.35*sin(z0(3)), ...
	'Color', 'b', 'MaxHeadSize', 4, 'LineWidth', 2);
drawnow;

%% Just for displaying costs-to-go of boundary nodes
% fprintf('\n');
% bdSortedCosts = [];
% for m = 1:numel(bdData)
% 	bdSortedCosts = cat(1, bdSortedCosts, [bdData(m).node bdData(m).cost]);
% end
% tttmp = sortrows(bdSortedCosts, 2);
% disp(tttmp)


%% Search!
%------ Manhattan heuristic
mnhtHeur = zeros(size(G,1), 1);
for m = 1:size(G,1)
	mnhtHeur(m) = 1*mrNodeData(m).d;
end

% disp(traceGreedyV04(nodeG, nodeS, mrNodeData))

fprintf('Initial path ...\t\t\t\t\n');

tic
[optCost, optPath, ~, nodeExp, nodeData] = hsmrWStatesV06(G, VCell, ...
	mnhtHeur, nodeS, z0, mrBdNodes, bdData, H, nHMem, dMax, varDim, ...
	@tilePlanV07, cellScale);
% [optCost, optPath, ~, nodeExp, nodeData] = hsmrWStatesV05(G, VCell, ...
% 	mnhtHeur, nodeS, z0, mrBdNodes, bdData, H, nHMem, dMax, varDim, @tilePlanV04d);
toc

if optCost < Inf
	zPlot	= z0; sysTraj = z0; uMPC = optPath.ip;
	for m = 1:size(uMPC, 2)	
		[tTmp, zTmp] = ode45(@(t,z) mpcPlant09(t, z, uMPC(1:2, m), ...
			modelParams), [0 uMPC(3,m)], zPlot);
		zPlot = zTmp(end,:)';
		sysTraj = cat(2, sysTraj, zPlot);
	end
	
	mrSPath = optPath.nodes;
	
	plot(fineAxes, sysTraj(1, :)/cellScale, sysTraj(2, :)/cellScale, 'b');
	drawCdV05(VCell, envAxes, 'b', 2, 0, mrSPath)
	drawCdV05(VCell, envAxes, 'b', 2, 0, (bdData(mrSPath(end) == mrBdNodes).optP));
end

%% Further "incremental" planning
fprintf('Further path planning ...\t\t\n');
VAct	= VCell;	GAct	= G;	NzrAct	= NzrData;						% Multi-resolution graph stuff
pPtr	= 2;
nIter	= 0;																% No. iterations, time taken

fineGraphPath	= [];		controlInp = [];								% Path in graph and the control input
stateTraj		= [z0; 0];													% State trajectory; the 0 is for time

% envAvi = avifile('mrMovie03.avi', 'fps', 5, 'compression', 'None');
% fineAvi= avifile('mrMovie04.avi', 'fps', 5, 'compression', 'None');
%%
% while VAct(nodeG, 3) > 1
while nIter < 45
	fprintf('Current speed: %f units/s\n', stateTraj(4,end));	
	nIter = nIter + 1;
	fprintf('\tIteration %i ...\t\t', nIter);
	
	% -------- Movie stuff -------
	envFrame			= getframe(envAxes);
	fineFrame			= getframe(fineAxes);
	envMovie(nIter)		= envFrame;
	fineMovie(nIter)	= fineFrame;
% 	envAvi				= addframe(envAvi,envFrame);
% 	fineAvi				= addframe(fineAvi,fineFrame);
	
	% -------- Step forward, note dlta = (dx, dy) --------
	pAct		= (VAct(mrSPath(pPtr-1), 1:2))';
	pNext		= (VAct(mrSPath(pPtr), 1:2))';
	pNextSize	= VAct(mrSPath(pPtr), 3);									% Size of next cell (can force to be one)	
	dlta		= pNext - pAct;
	
	% -------- For recording cost-to-go estimate --------
	[~,nodePrvFine] = ismember(pAct', VFine(:,1:2), 'rows');
	[~,nodeActFine] = ismember(pNext', VFine(:,1:2), 'rows');
	
	% -------- Record path and control --------
	fineGraphPath	= cat(2, fineGraphPath, nodePrvFine);
	thisNode		= optPath.nodes(pPtr + H);
	thisIndx		= optPath.indx(pPtr - 1);
	ipAct			= nodeData(thisNode).ipSeq(thisIndx).ip;
	controlInp		= cat(2, controlInp, ipAct);
	zk				= stateTraj(1:4, end);		tk = stateTraj(5, end);		% <-- HARD CODED FOR 4-state MODEL
	
	for m = 1:size(ipAct, 2)
		[~, zTmp]	= ode45(@(t,z) mpcPlant09(t, z, ipAct(1:2, m), ...
			modelParams), [0 ipAct(3,m)], zk);
		zk			= zTmp(end, :)';
% 		zk			= pdDiscreteModel(ipAct(1:2, m), zk, ipAct(3, m));
		tk			= tk + ipAct(3, m);
		stateTraj	= cat(2, stateTraj, [zk; tk]);
	end
	
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
	
	if (numel(inOnN) == 0) && (numel(inNnO) == 0)							% In case decomposition doesn't change
		fprintf('continuing previous trajectory.\n');
		pPtr	= pPtr + 1;	
		continue;															% Basically, don't do anything, just step fwd
	end
	pPtr = 2;
	
	% -------- Update graph and set of cells --------
	[GNext, VNext]	= getNewGraphV02(NzrAct.ANzr, GAct, ...
		VAct(:, [1:3 8]), inOnN, inNnO, Cforig, Sz);
	nVNext	= size(VNext, 1);
	VNext	= cat(2, VNext, [ones(nVNext,2) kron(ones(nVNext,1), vLim') VNext(:,4)]);
	VNext(:, 4) = [];
	% MORE GENERALLY: the previous line will have to be done by associating
	% the cells in VNext with those in VFine and getting the "color" and
	% vlim information from there.
	
	% -------- Locate current cell and goal in new decomposition --------
	[~, nodeAct]	= ismember([pNext' pNextSize], VNext(:, 1:3), 'rows');
% 	disp(nodeAct)
	for j = (-N):jmax
		pGoalCell	= floor(pGoal*(2^(j)))*(2^(-j));
		[tmp, nodeG]= ismember([pGoalCell' 2^(-j)], VNext(:,1:3), 'rows');
		if tmp, break; end
	end
		
 	% -------- Redefine heuristic for changed nodes --------
	heur = zeros(nVNext,1);
	
	% -------- Boundary points and penalty
	GNext	= adj2cost06(GNext, VNext(:,[1:3 8]), nodeG);
	bdData	= fineResBoundaryV01(GNext, VNext, dMax);
	mrBdNodes = [];	
	for m = 1:numel(bdData)
		mrBdNodes	= cat(2, mrBdNodes, bdData(m).node);
	end
	mrNodeData = aStarV04(GNext', nodeG, 1:nVNext, heur);					% A* to find cost-to-go from each boundary node to goal

	for m = 1:numel(bdData)
		bdData(m).cost	= mrNodeData(bdData(m).node).d;
		tmp1			= traceGreedyV04(nodeG, bdData(m).node, mrNodeData);
		bdData(m).optP	= fliplr(tmp1);
	end
	
	% -------- Draw things --------
	cla(envAxes); hold on; grid on; axis equal;
	drawCdV05(VNext(:,[1:3 8]), envAxes, 'k', 1, 0, [], 0, 2^(-jmin), map)
	drawCdV05(VNext, envAxes, 'r', 3, 0, nodeG)
	drawCdV05(VNext, envAxes, 'r', 2, 0, nodeAct)
	plot(envAxes, stateTraj(1,1)/cellScale, stateTraj(2,1)/cellScale, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 5);
	plot(envAxes, stateTraj(1,:)/cellScale, stateTraj(2,:)/cellScale, 'b', 'LineWidth', 2.5);
	drawnow;
	
	cla(fineAxes); hold on; grid on;
	tmp1	= 1:nVNext;
	drawCdV05(VNext, fineAxes, 'k', 2, 0, tmp1(VNext(:,3) < 2))
	for m = 1:numel(bdData)
		drawCdV05(VNext, fineAxes, 'r', 3, 0, bdData(m).node)
	end
	zPlot	= stateTraj(1:4, end); 
	plot(fineAxes, zPlot(1)/cellScale, zPlot(2)/cellScale, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6); 
	quiver(fineAxes, zPlot(1)/cellScale, zPlot(2)/cellScale, 0.65*cos(zPlot(3)), 0.65*sin(zPlot(3)), ...
		'Color', 'r', 'MaxHeadSize', 8, 'LineWidth', 2); axis equal;
	drawnow;
	
% 	fprintf('\n');
% 	bdSortedCosts = [];
% 	for m = 1:numel(bdData)
% 		bdSortedCosts = cat(1, bdSortedCosts, [bdData(m).node bdData(m).cost]);
% 	end
% 	disp(sortrows(bdSortedCosts, 2))
	
	% -------- Redefine heuristic for changed nodes --------
	mnhtHeur = zeros(nVNext,1);
	for m = 1:nVNext
		mnhtHeur(m) = 1*mrNodeData(m).d;
	end
	
	% -------- Run search again --------
% 	[optCost, optPath, ~, nodeExp, nodeData] = hsmrWStatesV05(GNext, ...
% 		VNext, mnhtHeur, nodeAct, stateTraj(1:4, end), mrBdNodes, ...
% 		bdData, H, nHMem, dMax, varDim, @tilePlanV04d);
	
	[optCost, optPath, ~, nodeExp, nodeData] = hsmrWStatesV06(GNext, ...
		VNext, mnhtHeur, nodeAct, stateTraj(1:4, end), mrBdNodes, ...
		bdData, H, nHMem, dMax, varDim, @tilePlanV07, cellScale);
	
	% -------- Trace path and draw
	if optCost < Inf
		fprintf('found a trajectory.\n');
		mrSPath = optPath.nodes;
		trajPlot = zPlot; uMPC = optPath.ip;
		for m = 1:size(uMPC, 2)
			[tTmp, zTmp] = ode45(@(t,z) mpcPlant09(t, z, uMPC(1:2, m), ...
				modelParams), [0 uMPC(3,m)], zPlot);
			zPlot	= zTmp(end,:)';
% 			tPlot	= cat(2, tPlot, tPlot(end) + tTmp(end));
			trajPlot= cat(2, trajPlot, zPlot);
		end
		plot(fineAxes, trajPlot(1, :)/cellScale, trajPlot(2, :)/cellScale, 'b', 'LineWidth', 2);
		drawCdV05(VNext, fineAxes, 'b', 3, 0, mrSPath);
		drawnow;

		drawCdV05(VNext, envAxes, 'b', 2, 0, mrSPath)
		drawCdV05(VNext, envAxes, 'b', 3, 0, (bdData(mrSPath(end) == mrBdNodes).optP));
		drawnow;		
	else
		fprintf('no trajectory, stopped.\n');
		break;
	end
	
	% -------- Update everything --------
	GAct	= GNext;
	VAct	= VNext;
	NzrAct	= NzrNext;
	
% 	procd = input('Procced?		');
% 	if procd == 5, break; end	
end
% envAvi	= close(envAvi);
% fineAvi = close(fineAvi);

save Data/mrMotionPlanTestV08ResultsAsym04

%%
figure('Position', [300 18 1200 950]); endAxes = axes;
hold on; grid on; axis equal;
% drawCdV05(VFine(:, [1:3 8]), endAxes, 'k', 0.01, 0, [], 0, 2^(-jmin), map)
plot(pInit(1)/cellScale, pInit(2)/cellScale, 'ko', 'MarkerFaceColor', 'k');

image(Xtrue, 'XData', 0.5, 'YData', 0.5); colormap(map);
plot(endAxes, stateTraj(1,:)/cellScale, stateTraj(2,:)/cellScale, 'b');

%%
% close all; clear all; clc
% load Data/mrMotionPlanTestV08Results02Asym
% 
% figure('Position', [300 18 1200 950]); endAxes = axes;
