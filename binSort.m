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

Description: Brute force binary sort
%}

%**************************************************************************
function A = binSort(A, B, c)
%--------------------------------------------------------------------------
% The rows of B are inserted into A, while sorting (ascending) according to
% column c. Both A and B have the same number of columns. A is assumed to
% sorted ascending.

[rA, cA] = size(A);
[rB, cB] = size(B);

if numel(A) == 0, A = B; return; end
if numel(B) == 0, return; end
if cB ~= cA, error('A and B must have same number of columns!\n'); end

for count = 1:rB
	thisIns		= B(count, :);
	thisCost	= thisIns(1, c);
	
	if ismember(thisIns, A, 'rows'), 
% 		fprintf('This one came back!\t\t'); disp(thisIns)
% 		redn = redn + 1;
		continue;
	end

	if A(rA, c) <= thisCost													% Insert after last row
		A	= cat(1, A, thisIns);
		rA	= rA + 1;
		continue;
	elseif A(1, c) >= thisCost												% Insert before first row
		A	= cat(1, thisIns, A);
		rA	= rA + 1;
		continue;
	end
	
	nCand	= rA;															% Number of candidate rows in A that can have greater cost
	testRow	= 0;
	dirn	= 1;	
	while nCand > 1
		p		= floor(nCand/2);
		testRow = testRow + dirn*p;
		dirn	= 2*(A(testRow, c) < thisCost) - 1;
		nCand	= nCand - p;
	end	

	insRow = testRow + (dirn + 1)/2;										% Insert before this row in A
	A	= cat(1, A(1:(insRow - 1), :), thisIns, A(insRow:rA, :));
	rA	= rA + 1;
end
%**************************************************************************