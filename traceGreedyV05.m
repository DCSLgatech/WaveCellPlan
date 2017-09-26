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

Description: 
%}
function optP= traceGreedyV05(nodeS, nodeG, nodeData, knownNodes)

if nargin == 3
	nodeAct = nodeG;
	optP	= nodeG;

	while (nodeAct ~= nodeS)
		nodeAct = nodeData(nodeAct).b;
		optP	= cat(2, nodeAct, optP);
	end
	return
end

[hasGoal, idAct] = ismember(nodeG, knownNodes);
if ~hasGoal, error('knownNodes does not contain goal!'); end

optP = nodeG;
while (idAct ~= 1)
	idAct	= nodeData(idAct).b;
	optP	= cat(2, knownNodes(idAct), optP);
end
