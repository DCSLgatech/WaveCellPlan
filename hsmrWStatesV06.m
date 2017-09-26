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

Description: Efficient graph search for history-based cost
%}
function [optCost, optPath, nData, nodeExp, nodeData] = ...
	hsmrWStatesV06(G, VCell, heur, nodeS, initState, mrBdNodes, bdData, ...
	H, nHMem, dMax, varDim, tpHandle, cellScale)
global envAxes fineAxes g vCruise

% Descriptions of main identifiers
	% G		: transition cost matrix of original graph
	% A		: original adjacency structure (neighhbors, cost(?))
	% VCell	: locations and dimensions of cells
	% heur	: heuristic
	% nodeS	: start node in original graph
	% initState: initial state
	% nodeG	: goal node in original graph	
	% H		: length of history cells	
	% nHMem	: number of histories to be remembered per node, Inf =>
	%		remember everything
	% dMax	: Largest cell size to consider for history, Inf => all
	% varDim: variable dimensions
	%	varDim.nStates, varDim.nInputs
	% tpHandle: function handle to tile planner, i/o format below
	%	 [newCost, newSt, newIp] = tilePlannerName(nodeSeq, stAct, VCell)
% ------- Variables Initialization -------
N		= size(G, 1);														% Number of nodes in original graph
nOpen	= 0;
openL	= [];																% Sorted OPEN lists

nHSeries= [4 12 36 100 284 780 2172 5916 16268 44100]';
nHWorst	= min(nHMem, nHSeries(H));

noUse	= struct('ip', []);
ipStr	= repmat(noUse, 1, nHWorst);										% Input sequence, one per history
noUse	= struct('hist', [], 'hLen', [], 'd', [], 'h', [], 'mk', [], ...
	'nHX', [], 'st', [], 'ipSeq', ipStr, 'ot', [], 'ohst', [], 'ohL', []);
nodeData= repmat(noUse, 1, N);

for Vnode = 1:(N+1)
	nodeData(Vnode).hist= zeros(nHWorst, H+1);								% List of histories known so far
	nodeData(Vnode).hLen= zeros(nHWorst, 1);								% Length of each history, needed for "multi-history"
	nodeData(Vnode).d	= zeros(nHWorst, 1);								% Label, one per history
	nodeData(Vnode).h	= zeros(nHWorst, 1);								% 'Tail', one per history
	nodeData(Vnode).mk	= zeros(nHWorst, 1);								% Marker, 0 = NEW; 1 = OPEN; 2 = CLOSED
	nodeData(Vnode).nHX	= 0;												% Number of histories known so far (X for explored)	
	nodeData(Vnode).st	= zeros(varDim.nStates, nHWorst);					% State, one per history
	nodeData(Vnode).ot	= zeros(nHWorst, 1);								% Original 'Tail', one per history
	nodeData(Vnode).ohst= zeros(nHWorst, H+1);								% List of original histories known so far
	nodeData(Vnode).ohL = zeros(nHWorst, 1);								% Original length of each history, needed for "multi-history"
	% For some reason, initializing with zeros instead of sparse results in
	% much faster execution.
end

% ------- Algorithm Initialization -------
nodeSHist = sizeConstrHistoryV02(G, VCell, dMax, H+1, nodeS, []);
for L = H:-1:1	
	allnSHist	= sizeConstrHistoryV02(G, VCell, dMax, L, nodeS, nodeSHist);	% Get "single" histories of start node of length L+1
	nodeSHist	= cat(1, nodeSHist, [allnSHist zeros(size(allnSHist, 1), H+1-L)]);
end

VCell(:, 1:3) = VCell(:, 1:3)*cellScale;

nodeSnhbrs = find(G(nodeS, :));
% disp(initState')
relvSNhbr = 0;
for m = nodeSnhbrs
% 	abs(initState(1:2, :) - (VCell(m, 1:2) + VCell(m, 3)/2)')
	if max(abs(initState(1:2, :) - (VCell(m, 1:2) + VCell(m, 3)/2)')) ...
			<= (0.5 + 1e-5)*cellScale
		relvSNhbr = m;
		break;
	end
end
% disp(relvSNhbr)
% disp(nodeSHist)

VCell(:, 1:3) = VCell(:, 1:3)./cellScale;

if relvSNhbr
	nodeSHist = nodeSHist(nodeSHist(:, 2) == relvSNhbr, :);
end

% disp(nodeSHist)
% 
% optCost = [];
% optPath = [];
% nData = [];
% nodeExp = [];
% return

nodeExp = [];
for count = 1:size(nodeSHist, 1)
	tHistLen	= sum(nodeSHist(count, :) > 0);								% In next line, tHistLen replaces (H+2) from standard history search
	tNode		= nodeSHist(count, tHistLen);								% Last node in H+2 history, will go to OPEN
	[tHistCost, tHistSt, tHistIp] = tpHandle(nodeSHist(count, ...
		1:(tHistLen)), initState, VCell, cellScale);							% Cost of this history
	
	tnHX		= nodeData(tNode).nHX;
	if tnHX < nHMem
		tnHX				= tnHX + 1;
		nodeData(tNode).nHX	= tnHX;
		tIndx				= tnHX;
	else																	% What to do when enough histories are known
		[worstCost, tIndx]	= max(nodeData(tNode).d(1:tnHX));
		if tHistCost >= worstCost, continue; end							% If this history is worse than all known ones, ignore it
		
		% If not, remove the worst known history from OPEN and replace it
		% with this one
		[~, loc]		= ismember([tNode tIndx], openL(1:nOpen, 1:2), 'rows');
		openL(loc, :)	= [];
		nOpen			= nOpen - 1;
	end
	
% 	[tNode tIndx  tHistLen (nodeSHist(count, 2:tHistLen))]
	nodeData= updateNodeData(nodeData, tNode, tIndx, ...
		(nodeSHist(count, 2:tHistLen)), nodeS, (tHistLen - 1), ...
		tHistCost, tHistSt, tHistIp, 1, 1);
	
	nodeExp = cat(1, nodeExp, tNode);
	openL	= cat(1, openL, [tNode tIndx (tHistCost + heur(tNode))]);		% Add to OPEN list
	nOpen	= nOpen + 1;
end
openL(1:nOpen, :) = sortrows(openL(1:nOpen, :), 3);

% ------- Algorithm Iterative Steps -------
fprintf('Searching v5.0... \n');
nIter	= 0; nStExp = nOpen; nHExp = nOpen; tConn = 0;
goalClosed = 0;

while (nOpen > 0) && (~goalClosed) && (openL(1, 3) < Inf)
	nIter	= nIter + 1;
	
	fprintf('******************** Iteration %i\n', nIter);
% 	tmp5 = min(2, nOpen);
% 	disp(openL((openL(1:tmp5,3) < Inf), :))
% 	disp(openL((openL(:,3) < Inf), :))
		
	nodeAct	= openL(1, 1);		histAct	= openL(1, 2);
	
	if (nodeAct == (N+1)), goalClosed = 1; continue; end					% Goal is a dummy node
	
% 	if dispOn, fprintf('*** nodeAct = %i\n', nodeAct); end	
	fprintf('*** nodeAct = %i\n', nodeAct);
	disp(nodeData(nodeAct).hist(histAct, :))

% 	nhbrs		= A(nodeAct).nhbrs;											% Note all neighbors of nodeAct
	nhbrs		= find(G(nodeAct, :));										% Note all neighbors of nodeAct
	newOpen		= [];

	newHistLen	= nodeData(nodeAct).hLen(histAct, 1);
	newTail		= nodeData(nodeAct).hist(histAct, 1);						% Possible tail for nodeNew
	actCost		= nodeData(nodeAct).d(histAct);								% Cost to come to (nodeAct, histAct)
	actState	= nodeData(nodeAct).st(:, histAct);
	disp(actState')
	
	if (VCell(nodeAct,3) >= dMax) && (newHistLen == 2)						% Adjacent to the goal
% 		fprintf('Almost done... \n')
		[isKnownBd, locBd] = ismember(nodeAct, mrBdNodes);
		
		if isKnownBd
			goalCost = (actCost + bdData(locBd).cost);			
		else
			fprintf('***** WARNING: Oversize cell of node (%i) with no boundary cost associated. *****\n', nodeAct);
			goalCost = Inf;
		end
		
		if nodeData(N+1).mk == 0											% Goal not already in OPEN
			newOpen	= cat(1, newOpen, [(N+1) 0 goalCost]);					% Add goal to OPEN, with known terminal cost
			nodeData(N+1).mk			= 1;
			nodeData(N+1).hist(1, 1:2)	= [nodeAct histAct];
			nodeData(N+1).d				= goalCost;
		elseif goalCost < nodeData(N+1).d									% If c2come to goal can be improved
			[~, locGoal]		= ismember((N+1), openL(1:nOpen, 1));
			nOpen				= nOpen - 1;
			openL(locGoal, :)	= [];
			
			newOpen	= cat(1, newOpen, [(N+1) 0 goalCost]);
			nodeData(N+1).hist(1, 1:2)	= [nodeAct histAct];
			nodeData(N+1).d				= goalCost;
		end
		
		openL(1, :)	= [];
		nOpen		= nOpen - 1;
		
		tmpOpen = binSort(openL(1:nOpen, :), newOpen, 3);					% Add newly opened pairs to sorted open list
		if numel(tmpOpen) == 0
			nOpen	= 0;
			openL	= [];
		else
			nOpen	= size(tmpOpen, 1);
			openL(1:nOpen, :)	= tmpOpen;
		end		
		continue;
	end
	
	newHistAct	= nodeData(nodeAct).hist(histAct, 2:newHistLen);			% History corresponding to (nodeAct, histAct), newHistLen replaces (H+1) here
	
	if (VCell(nodeAct,3) >= dMax)											% Technically not adjacent to the goal, do not explore further neighbors
		newHIndx= histAct;
		nodeNew	= nodeAct;
		
		[tHistCost, newSt, newIp] = tpHandle([newTail newHistAct], ...
			actState, VCell, cellScale);
		newCost	= actCost + tHistCost;
		
		nodeData= updateNodeData(nodeData, nodeAct, histAct, ...
			newHistAct, newTail, (newHistLen-1), newCost, newSt, ...
			[nodeData(nodeNew).ipSeq(newHIndx).ip newIp], 1, 0);
				
		newOpen	= cat(1, newOpen, ...
			[nodeNew newHIndx (newCost + heur(nodeNew))]);					% Set of newly opened (node, hist) pairs		
		nhbrs	= [];		
	end

	% ---- Explore each neighbor
	for nodeNew = nhbrs
% 		if dispOn, fprintf('\t--- nodeNew = %i\n', nodeNew); end
		fprintf('\t--- nodeNew = %i\n', nodeNew);
		
		if nodeNew == nodeData(nodeAct).h(histAct), continue; end
		if nodeNew == newTail, continue; end
		if any(nodeNew == newHistAct), continue; end
		if any(G(nodeNew, newHistAct(1:(end-1)))), continue; end
		% First (H-1) elements of each row of newHistAct are "indices" of
		% histories of nodeNew through nodeAct unique up to all but one
		% element (tail of one node) of histories of nodeAct. This cost for
		% nodeNew through nodeAct for each of these "indices" is to be
		% minimized over that tail element.	
		
		nHXNew		= nodeData(nodeNew).nHX;
		[isKnown, loc]= ismember(newHistAct,  ...
			nodeData(nodeNew).hist(1:nHXNew, 1:(newHistLen-1)), 'rows');
		
% 		disp([isKnown loc])
		if isKnown
			newHIndx = loc;
			
			if (nodeData(nodeNew).mk(newHIndx) == 2), continue; end
			% If a history is already CLOSED, it cannot be improved. Comes
			% from the proposition that says every pair enters open list
			% only once. 
		else
			if nHXNew < nHMem
				newHIndx				= nHXNew + 1;
				nodeData(nodeNew).nHX	= newHIndx;
			else
				nodeNewOpenHist			= find(nodeData(nodeNew).mk(1:nHXNew) == 1);
				if numel(nodeNewOpenHist) == 0
					continue;
				else
					[worstNewCost, newHIndx]= max(nodeData(nodeNew).d(nodeNewOpenHist));
					newHIndx				= nodeNewOpenHist(newHIndx);
					if actCost >= worstNewCost, continue; end				% new history of nodeNew cannot improve cost over its known ones
				end
			end
		end
		
% 		tic
		[tHistCost, newSt, newIp] = tpHandle([newTail newHistAct ...
			nodeNew], actState, VCell, cellScale);
% 		ttTmp = toc;
		
		newCost	= actCost + tHistCost;
		
		if ~isKnown
			nodeData= updateNodeData(nodeData, nodeNew, newHIndx, ...
				[newHistAct nodeNew], newTail, newHistLen, newCost, ...
				newSt, newIp, 1, 1);
			
			newOpen	= cat(1, newOpen, ...
				[nodeNew newHIndx (newCost + heur(nodeNew))]);				% Set of newly opened (node, hist) pairs
			nStExp	= nStExp + 1;
			nHExp	= nHExp + 1;
% 			tConn = tConn + ttTmp;
			
		elseif (nodeData(nodeNew).d(newHIndx) > newCost)
			[~, t2]	= ismember([nodeNew newHIndx], openL(1:nOpen, 1:2), 'rows');
			nOpen	= nOpen - 1; 
			openL(t2, :)= [];
			
			nodeData= updateNodeData(nodeData, nodeNew, newHIndx, ...
				[newHistAct nodeNew], newTail, newHistLen, newCost, ...
				newSt, newIp, 1, 1);
			
			newOpen	= cat(1, newOpen, ...
				[nodeNew newHIndx (newCost + heur(nodeNew))]);				% Set of newly opened (node, hist) pairs
			nHExp	= nHExp + 1;
		end
		% Intentionally split these two ifs: program executes faster
	end
	
	nodeExp = cat(1, nodeExp, nodeAct);
	if (VCell(nodeAct,3) >= dMax)
		drawCdV05(VCell, fineAxes, 'r', 3, 0, nodeAct, 'y')
	else
		drawCdV05(VCell, fineAxes, 'k', 2, 0, nodeAct, 'y')
	end
	nodeData(nodeAct).mk(histAct) = 2;										% Mark this pair closed
	openL(1, :)	= [];														% Remove (nodeAct, histAct) from OPEN
	nOpen		= nOpen - 1;
	
	tmpOpen = binSort(openL(1:nOpen, :), newOpen, 3);						% Add newly opened pairs to sorted open list
	if numel(tmpOpen) == 0
		nOpen	= 0;
		openL	= [];
	else
		nOpen	= size(tmpOpen, 1);
		openL(1:nOpen, :)	= tmpOpen;
	end
% 	disp(openL(1:nOpen,	(openL(1:nOpen,3) < Inf)))
end
nData = [nIter nStExp nHExp];
% tConn/nStExp

% ------- Trace Optimal Path -------
nodeG	= nodeData(N+1).hist(1, 1);
optHist = nodeData(N+1).hist(1, 2);
optCost = nodeData(nodeG).d(optHist);

if (numel(optCost) == 0) || (optCost == Inf)
	fprintf('No path found.\n');
	optPath = [];
	return;
end

tail	= nodeData(nodeG).ot(optHist);
head	= nodeG;
hIndx	= optHist;
hLen	= nodeData(nodeG).ohL(optHist);
optPath.nodes	= nodeData(nodeG).ohst(hIndx, 1:hLen);
optPath.ip		= [];
optPath.indx	= [];
optPath.hLen	= [];
while tail ~= nodeS
	tail	= nodeData(head).ot(hIndx);										% Get tail	
	nextH	= [tail nodeData(head).ohst(hIndx, 1:(hLen-1))];
	
	optPath.nodes	= [tail optPath.nodes];
	optPath.ip		= [nodeData(head).ipSeq(hIndx).ip optPath.ip];
	optPath.indx	= [hIndx optPath.indx];
	optPath.hLen	= [hLen optPath.hLen];
	head	= nodeData(head).ohst(hIndx, (hLen-1));							% Step back one
		
	[~, hIndx]	= ismember(nextH, nodeData(head).ohst, 'rows');				% The next index is the one that corresponds to the history
																			% formed by tail and first H nodes of current history
	
	if hIndx
		hLen		= nodeData(head).ohL(hIndx);
	end
end
%--------------------------------------------------------------------------
%**************************************************************************
function UL = sizeConstrHistoryV02(G, VCell, dMax, L, nodeV, UH)
% Specialized for finding histories of a single node

UL	= [];
lnV = getAllHistV04(G, nodeV, L+1, 0, nodeV, []);
if ~numel(lnV), return; end

for n = 1:L
	lnV((VCell(lnV(:,n),3) >= dMax), :) = [];								% Remove nodes corresponding to big cells
	if ~numel(lnV), break; end
end
if ~numel(lnV), return; end

for n = 1:L																	% Remove big-to-small transitions
	dDiff = VCell(lnV(:, n+1), 3) - VCell(lnV(:, n), 3);
	lnV((dDiff < 0), :) = [];
	if ~numel(lnV), break; end
end
if ~numel(lnV), return; end

if numel(UH)
	projRows		= ismember(lnV, UH(:, 1:L+1), 'rows');					% Remove projections
	lnV(projRows,:)	= [];
end

UL = lnV;
%--------------------------------------------------------------------------
%**************************************************************************
function nodeData = updateNodeData(nodeData, nodeNew, newHIndx, ...
	newHist, newTail, newHistLen, newCost, newSt, newIp, newMk, chOrg)

H	=	numel(nodeData(nodeNew).hist(1,:)) - 1;
nodeData(nodeNew).hist(newHIndx, 1:(H+1))		= zeros(1, H+1);
nodeData(nodeNew).hist(newHIndx, (1:newHistLen))= newHist;
nodeData(nodeNew).h(newHIndx)					= newTail;
nodeData(nodeNew).hLen(newHIndx, :)				= newHistLen;
nodeData(nodeNew).d(newHIndx)					= newCost;
nodeData(nodeNew).st(:, newHIndx)				= newSt;
nodeData(nodeNew).ipSeq(newHIndx).ip			= newIp;
nodeData(nodeNew).mk(newHIndx)					= newMk;

if chOrg
	nodeData(nodeNew).ohL(newHIndx)		= newHistLen;
	nodeData(nodeNew).ohst(newHIndx, 1:newHistLen) = newHist;
	nodeData(nodeNew).ot(newHIndx)		= newTail;
end
%--------------------------------------------------------------------------