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

Description: Target set analysis for a given tile
%}
function [targetSetData, tileData] = targetSetAnalysisV03(tileData, VCell)
%--------------------------------------------------------------------------
% Notes
% 1. tileData struct:
%		tileData.nodes	= nodes that make up the tile, ordered first to last
%		tileData.rad	= min radius of turn within each cell in the tile
%		tileData.cVert	= coordinates of A,B,C,D vertices of second cell of
%			the tile (through which traj is to be planned)
%		tileData.pType	= whether second cell has traversal across parallel
%			(1) or adjacent edges (2)
% 2. targetSetData struct:
% 		targetSetData.alfaSmp	= reachable sets for cells 2 through H+1
% 		targetSetData.btaSmp	= target sets for cells 2 through H+1; for
%			cell C = H+1, target set is ([-pi/2, pi/2] x entire exit segment)
% 		targetSetData.nSol		= number of points on each entry segment
%			for which reachable orientation set is computed
% 		targetSetData.nSmp		= number of grid points on each entry and
%			exit segments
% 		targetSetData.wL		= lower boundary of entry segments
% 		targetSetData.wU		= upper boundary of entry segments
% 		targetSetData.xSmp		= grid points on exit segments
% 		targetSetData.wSmp		= grid points on entry segments
% 3. Code matches new notation
% 4. Functions called:
%		targetOrtSetP1V01 - parallel edges, large radius of turn
%		targetOrtSetP2V01 - adjacent edges, large radius of turn

% --------- Constants
faceRef = ...																% 0 = identity; 1 = r-90; 2 = r90; 3 = FH; 4 = FV
	  [	 1	-1	1	0	0;
		-2	 2	1	2	0;
		 2	-2	1	1	0;
		-1	 1	1	4	0;
		 1	-2	2	0	0;
		-2	-1	2	2	0;
		 2	 1	2	1	0;
		-1	-2	2	4	0;
		 1	 2	2	3	0;
		 2	-1	2	4	2;
		-2	 1	2	4	1;
		-1	 2	2	4	3	];

vPermut	= [	1	2	3	4;													% Row 1 = Base: no transformation
			4	1	2	3;													% Row 2 = Rotate -90
			2	3	4	1;													% Row 3	= Rotate 90
			4	3	2	1;													% Row 4 = Flip horizontal
			2	1	4	3];													% Row 5 = Flip vertical
I		= eye(4);
invTf	= [0; 2; 1; 3; 4];													% Inverse transformations in the same order as above rows
nSmp	= 100;																% Number of sample points
nSol	= 20;																% Number of points where aU and aL are solved for

% --------- Channel Data
V	= VCell(tileData.nodes, :);
C	= numel(tileData.nodes) - 2;											% Number of rectangles in channel
chan(:, 1:2)	= (diag(V(2:(end-1), 3)))*ones(C, 2);						% Rectangle dimensions (dx, dy)
chan(:, 3:4)	= V(2:(end-1), 1:2) + chan(:, 1:2)/2;						% Rectangle center coordinates

fChan(:, 1:2)	= (diag(V(:, 3)))*ones(C+2, 2);								% Full channel, needed only to define
fChan(:, 3:4)	= V(:, 1:2) + fChan(:, 1:2)/2;								% entry and exit segments

yEx	= zeros(C,1);		zEx	= zeros(C,1);
wU	= zeros(C,1);		wL	= zeros(C,1);

wSol	= zeros(C, nSol);
wSmp	= zeros(C, nSmp);
xSmp	= zeros(C, nSmp);
btaSmp	= zeros(2*C, nSmp);													% aSmp stores relative angles, not transformed angles
btaSmp((2*C-1):(2*C), :) = [pi/2*ones(1, nSmp); -pi/2*ones(1, nSmp)];

alfaSol	= zeros(2*C, nSol);
alfaSmp(1:2:(2*C-1), :) = -Inf(C, nSmp);
alfaSmp(2:2:2*C, :)		= Inf(C, nSmp);

% --------- Find entry segments
ul		= zeros(C, 1);		exSeg	= zeros(C+1, 1);
cEdge	= zeros(C, 4);

pType	= zeros(C, 1);														% Problem 1 or 2
rTran	= zeros(C, 2);														% Transformation to bring to standard

for n = 2:(C+2)
	delX = fChan(n, 3) - fChan(n-1, 3);	
	tol = 1e-6;
	
	if abs(abs(delX) - sum(fChan(n-1:n, 1))/2) < tol
		if fChan(n, 3) > fChan(n-1, 3)										% X transition (left or right)
			exSeg(n-1)	= -1;												% Right
		else
			exSeg(n-1)	= 1;												% Left
		end
		ul(n-1) = 0;
		xLim	= fChan(n-1, 3) - exSeg(n-1)*fChan(n-1, 1)/2;
		yLim	= [min(fChan(n, 4) + fChan(n, 2)/2, fChan(n-1, 4) + fChan(n-1, 2)/2) ...
			max(fChan(n, 4) - fChan(n, 2)/2, fChan(n-1, 4) - fChan(n-1, 2)/2)];
		cEdge(n-1, :) = [xLim yLim(1) xLim yLim(2)];
	else		
		if fChan(n, 4) > fChan(n-1, 4)										% Y transition (up or down)
			exSeg(n-1)	= 2;												% Up
		else
			exSeg(n-1)	= -2;												% Down
		end
		ul(n-1) = 1;
		yLim	= fChan(n-1, 4) + (exSeg(n-1)/2)*fChan(n-1, 2)/2;
		xLim	= [min(fChan(n, 3) + fChan(n, 1)/2, fChan(n-1, 3) + fChan(n-1, 1)/2) ...
			max(fChan(n, 3) - fChan(n, 1)/2, fChan(n-1, 3) - fChan(n-1, 1)/2)];
		cEdge(n-1, :) = [xLim(1) yLim xLim(2) yLim];
	end	
end
fromFace	= -exSeg(1:C, 1);
toFace		= exSeg(2:(C+1), 1);

% --------- Get type of transition and (y, z)
cVert = zeros(4, 2*C);														% Vertices of each rectangle, coordinates in each column
for n = 1:C
	rVert(1,:) = [(chan(n, 3) - chan(n, 1)/2) (chan(n, 4) + chan(n, 2)/2)];	% Vertices of the rectangle
	rVert(2,:) = [(chan(n, 3) + chan(n, 1)/2) (chan(n, 4) + chan(n, 2)/2)];	% in CW order starting
	rVert(3,:) = [(chan(n, 3) + chan(n, 1)/2) (chan(n, 4) - chan(n, 2)/2)];	% with top-left
	rVert(4,:) = [(chan(n, 3) - chan(n, 1)/2) (chan(n, 4) - chan(n, 2)/2)];
	
	[~, loc]	= ismember([fromFace(n) toFace(n)], faceRef(:, 1:2), 'rows');
	rTran(n, :) = faceRef(loc, 4:5);										% Existing transformation on current rectangle
	tf1	= invTf(faceRef(loc, 5) + 1);	tf2 = invTf(faceRef(loc, 4) + 1);	% Get inverse transforms
	p1	= vPermut(tf1 + 1, :);			p2	= vPermut(tf2 + 1, :);
	
	cVert(:, (2*n - 1):(2*n)) = I(p2,:)*I(p1,:)*rVert;
	pType(n) = (ul(n) ~= ul(n+1)) + 1;
	
	switch pType(n)
		case 1																% Measure distance from C for y, z
			y1	= norm(cVert(3, (2*n - 1):(2*n)) - cEdge(n+1, 1:2));
			y2	= norm(cVert(3, (2*n - 1):(2*n)) - cEdge(n+1, 3:4));
		case 2																% Measure distance from D for y, z
			y1	= norm(cVert(4, (2*n - 1):(2*n)) - cEdge(n+1, 1:2));
			y2	= norm(cVert(4, (2*n - 1):(2*n)) - cEdge(n+1, 3:4));
	end
	yEx(n)	= (min(y1, y2));	zEx(n)	= (max(y1, y2));
end
for n = 1:C
	w1	= norm(cVert(4, (2*n - 1):(2*n)) - cEdge(n, 1:2));					% Measure distance from D for wL, wU
	w2	= norm(cVert(4, (2*n - 1):(2*n)) - cEdge(n, 3:4));
	wL(n)	= (min(w1, w2));
	wU(n)	= (max(w1, w2));
end
for n = 1:C
	if pType(n) == 2		
		wSolMin	= max(wL(n), 1e-3);
	else
		wSolMin	= wL(n);
	end
	wSol(n,:)	= linspace(wSolMin, wU(n), nSol);
	wSmp(n, :)	= linspace(wL(n), wU(n), nSmp);
	xSmp(n, :)	= linspace(yEx(n), zEx(n), nSmp);
end

tileData.cVert	= cVert(:, 1:2);
tileData.pType	= pType;
tileData.rTran	= rTran;
tileData.toFace = toFace(1);

% fprintf('\tProblem types: \t');		disp(pType');
% fprintf('\tTransformations: \n');	disp(rTran);
% fprintf('\tEntry segments: \n');	disp([wL'; wU']);
% fprintf('\tExit segments: \n');		disp([yEx'; zEx']);
% fprintf('\tFaces: \n');		disp([fromFace toFace]);

% --------- Cone Analysis
if pType(C) == 1
	coneFhandle = @targetOrtSetP1V03;
else
	coneFhandle = @targetOrtSetP2V03;
end
for q = 1:nSol
	w = wSol(C, q);	
	[alfa, ~,~,~]	= coneFhandle(w, chan(C,1), tileData.rad(C+1), xSmp(C,:), btaSmp((2*C-1):(2*C), :));
	alfaSol((2*C-1):(2*C), q)	= alfa';		
end

% --------- Interpolate alfaSol over the whole grid wSmp
alfaSmp(((2*C-1):(2*C)), :) = interpWInf(wSol(C, :), ...
	alfaSol(((2*C-1):(2*C)), :), wSmp(C, :), 5);

for p = (C-1):-1:1
	% ------- Transform alfa of previous to get bta for current
	if (ismember(3, rTran(p + 1, :)) && ~ismember(4, rTran(p + 1, :))) ...
			|| (~ismember(3, rTran(p + 1, :)) && ismember(4, rTran(p + 1, :)))	% If flipped once
		btaSmp((2*p-1):(2*p), :) = fliplr(-flipud(alfaSmp((2*p+1):(2*p+2), :)));
	else
		btaSmp((2*p-1):(2*p), :) = alfaSmp((2*p+1):(2*p+2), :);
	end
	if (ismember(3, rTran(p, :)) && ~ismember(4, rTran(p, :))) ...
			|| (~ismember(3, rTran(p, :)) && ismember(4, rTran(p, :)))		% If flipped once

		btaSmp((2*p-1):(2*p), :) = fliplr(-flipud(btaSmp((2*p-1):(2*p), :)));
	end

	% ------- Solve for target sets
	if pType(p) == 1
		coneFhandle = @targetOrtSetP1V03;
	else
		coneFhandle = @targetOrtSetP2V03;
	end

	for q = 1:nSol		
		w = wSol(p, q);	
		[alfa, ~,~,~]	= coneFhandle(w, chan(p,1), tileData.rad(p+1), ...
			xSmp(p,:), btaSmp((2*p-1):(2*p), :));
		alfaSol((2*p-1):(2*p), q)	= alfa';
	end

	% ------- Interpolate over all x
% 	alfaSol(((2*p-1):(2*p)), :)
% 	fprintf('***************************\n')
	alfaSmp(((2*p-1):(2*p)), :) = interpWInf(wSol(p, :), ...
		alfaSol(((2*p-1):(2*p)), :), wSmp(p, :), 5);
end

targetSetData.alfaSmp	= alfaSmp;
targetSetData.btaSmp	= btaSmp;
targetSetData.nSol		= nSol;
targetSetData.nSmp		= nSmp;
targetSetData.wL		= wL;
targetSetData.wU		= wU;
targetSetData.xSmp		= xSmp;
targetSetData.wSmp		= wSmp;
targetSetData.wSol		= wSol;
targetSetData.alfaSol	= alfaSol;
% -------------------------------------------------------------------------
%**************************************************************************
function yInter = interpWInf(x, y, xInter, polyOrder)

[~, notInfIndx] = remInf(y(1,:));
polyOrder		= min(polyOrder, (numel(notInfIndx)-1));
if numel(notInfIndx)
	indx1		= find(xInter >= x(notInfIndx(1)), 1, 'first');
	indx2		= find(xInter <= x(notInfIndx(end)), 1, 'last');
	yInter		= [-Inf(1, numel(xInter)); Inf(1, numel(xInter))];	
	yPolyApp1	= polyfit(x(notInfIndx), y(1, notInfIndx), polyOrder);
	yPolyApp2	= polyfit(x(notInfIndx), y(2, notInfIndx), polyOrder);
	yInter(:, indx1:indx2) = [polyval(yPolyApp1, xInter(indx1:indx2)); ...
		polyval(yPolyApp2, xInter(indx1:indx2))];	
else	
	yInter		= kron([-Inf; Inf], ones(1, numel(xInter)));
end
% -------------------------------------------------------------------------