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

Description: A* algorithm;Handles mmultiple goals when cost-to-come for
each goal is required 
%}
function nodeData = aStarV04(G, nodeS, nodeG, heur)
% G	: Trasition cost matrix
% nodeS		: Start node
% nodeG		: Goal node set
% heur	: Heuristic (N x 1 vector, where N is number of nodes)

%	nodeData(n).mk	= marker, 0 = NEW, 1 = OPEN, 2 = CLOSED
%	nodeData(n).d	= cost to come
%	nodeData(n).b	= backpointer

nNodes		= size(G, 1);
nodeStruct	= struct('mk', 0, 'd', Inf, 'b', []);
nodeData	= repmat(nodeStruct, 1, nNodes);
nodeData(nodeS).mk	= 1;	nodeData(nodeS).d	= 0;

nOpen	= 1;
openL	= [nodeS heur(nodeS)];
goalCl	= 0;

nIter	= 0;
while (nOpen ~= 0) && (~goalCl)
	nIter	= nIter + 1;
	nodeAct	= openL(1, 1);													% Get node from top of (sorted) open stack
	nodeData(nodeAct).mk = 2;												% Mark that node as dead
	
	nOpen		= nOpen - 1;
	openL(1, :) = [];
	
	nhbrs	= find(G(nodeAct,:));
	for nodeNew = nhbrs														% For all neighbors		
		newCost	= G(nodeAct,nodeNew);										% Cost to go from act to new
		
		if nodeData(nodeNew).mk == 0										% Unvisited
			nodeData(nodeNew).mk	= 1;									% Mark open
			nodeData(nodeNew).d		= nodeData(nodeAct).d + newCost;		% Update c2come of newly visited state
			nodeData(nodeNew).b		= nodeAct;
			
			tmpOpen = binSort(openL(1:nOpen, :), [nodeNew nodeData(nodeNew).d + heur(nodeNew)], 2);
			if numel(tmpOpen) == 0
				nOpen	= 0;
				openL	= [];
			else
				nOpen	= size(tmpOpen, 1);
				openL(1:nOpen, :)	= tmpOpen;								% Add [nodeNew cost] to sorted open list
			end			
		elseif nodeData(nodeNew).mk == 1									% Already open, update c2come if necessary
			if nodeData(nodeNew).d > nodeData(nodeAct).d + newCost
				nodeData(nodeNew).d	= nodeData(nodeAct).d + newCost;
				nodeData(nodeNew).b	= nodeAct;
				
				[~, loc] = ismember(nodeNew, openL(1:nOpen, 1));
				openL(loc, :)= [];		nOpen = nOpen - 1;
				
				tmpOpen = binSort(openL(1:nOpen, :), [nodeNew nodeData(nodeNew).d + heur(nodeNew)], 2);
				if numel(tmpOpen) == 0
					nOpen	= 0;
					openL	= [];
				else
					nOpen	= size(tmpOpen, 1);
					openL(1:nOpen, :)	= tmpOpen;							% Add [nodeNew cost] to sorted open list
				end
			end
		end
	end
	
	goalCl = 1;
	for k = 1:numel(nodeG)
		if nodeData(nodeG(k)).mk ~= 2, goalCl = 0; break; end
	end
end

% GoalC2Come = [];
% for k = 1:numel(G)
% 	GoalC2Come = cat(1, GoalC2Come, nodeData(G(k)).d);
% end
% [LCost, LGoal] = min(GoalC2Come);
% 
% if LCost < Inf
% 	SPath = G(LGoal);
% 	while SPath(1) ~= S
% 		SPath = cat(2, nodeData(SPath(1)).b, SPath);
% 	end
% else
% 	SPath = [];
% end