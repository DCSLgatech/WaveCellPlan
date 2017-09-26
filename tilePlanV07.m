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

Description: Aircraft dynamics tile planner
%}

function [newCost, newSt, newIp] = tilePlanV07(hist, stAct, VCell, cellScale)

global vCruise g

% if (nodeNew == size(VHMR, 1)) || (nodeAct == (size(VHMR, 1) - 1))	
% 	newCost = GHMR(nodeAct, nodeNew);
% 	newSt	= stAct;
% 	newIp	= [];
% 	return;
% end

tileData.nodes = hist;
% fprintf('\tCurrent tile \t\t:');  disp(tileData.nodes);
% fprintf('\tCurrent state \t\t: (%f, %f, %f, %f)\n\n', stAct(1), ...
% 	stAct(2), stAct(3)*180/pi, stAct(4));
[feas, newIp, ~, ~, newSt] = pdTilePlanV12(stAct, tileData, VCell, cellScale);
if feas
% 	newCost = sum(newIp(3, :));
% 	newCost = 1e-3;
	m		= hist(2);
% 	newCost = ((256 - VCell(m,end))*VCell(m,3) + 0.1*(VCell(m,3))^2)/10 ...
% 		+ sum(newIp(3, :))/3e3;												% Copied from adj2cost06
	newCost = ((256 - VCell(m,end))*VCell(m,3) + 0.1*(VCell(m,3))^2)/10;	% Copied from adj2cost06
else
	newCost	= Inf;
	newSt	= Inf(4,1);
	newIp	= [];
end
% fprintf('\b\tFeasible? \t\t\t: %i\n\n', feas);
% if feas
% 	fprintf('\tNew state \t\t: (%f, %f, %f, %f)\n\n', newSt(1), ...
% 		newSt(2), newSt(3)*180/pi, newSt(4));
% end
