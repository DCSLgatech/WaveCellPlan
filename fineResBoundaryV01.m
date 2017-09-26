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

Description: Find boundary nodes of a certain size of the fine resolution
"window" of a given multi-resolution cell dec 
%}
function bdData = fineResBoundaryV01(G, VCell, dMax)

% bdData struct:
% 	bdData.node	= node that forms the "large-cell" boundary of the fine res
% 	window 
% 	bdData.nhbrs= "small cell", interior nhbrs of that node
%	bdData.cost = cost-to-go from this boundary node to goal
%	bdData.optP = shortest path in G from this node to goal

dMaxNodes	= find(VCell(:, 3) == dMax);
bdStruct	= struct('node', [], 'nhbrs', [], 'cost', 0, 'optP', []);
bdData		= repmat(bdStruct, 1, numel(dMaxNodes));

nBdNodes	= 0;
for m = dMaxNodes'
	if any(VCell((G(:,m) ~= 0), 3) < dMax)
		nBdNodes				= nBdNodes + 1;
		bdData(nBdNodes).node	= m;
		bdData(nBdNodes).nhbrs	= find(VCell((G(:,m) ~= 0), 3) < dMax);
	end
end

bdData((nBdNodes+1):end) = [];