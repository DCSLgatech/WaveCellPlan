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
function [G, V] = findAdj03(inNewnotOldV, inOldnotNewV, remIndx, VOld, GOld);

GOld(remIndx, :) = [];
GOld(:, remIndx) = [];
VOld(remIndx,:)	= [];
V	= [VOld; inNewnotOldV];

if ~numel(inNewnotOldV)
	G = GOld;
	return;
end

nNew	= size(inNewnotOldV,1);
nOld	= size(VOld, 1);
newRec	= [];

for m1 = 1:nNew
	for m2 = (m1 + 1):nNew
		v1 = inNewnotOldV(m1, 1:3);
		v2 = inNewnotOldV(m2, 1:3);
		
		if v1(3) > v2(3)
			temp = v1;
			v1 = v2;
			v2 = temp;
		end
		indx = (v1(1:2) - v2(1:2))./v1(3);	
		
		
		if ((indx(1) <= v2(3)/v1(3)) && (indx(2) < v2(3)/v1(3))...
				&& (indx(1) >= -1) && (indx(2) > -1)) || ...
			((indx(2) <= v2(3)/v1(3)) && (indx(1) < v2(3)/v1(3))...
				&& (indx(2) >= -1) && (indx(1) > -1))						% 4-connectivity
			newRec = cat(1, newRec, [m1 m2 1; m2 m1 1]);
		end
		
	end
end
GNew = sparse(newRec(:,1), newRec(:,2), newRec(:,3), nNew, nNew);

newRec	= [];
for m1 = 1:nNew
	for m2 = 1:nOld
		v1 = inNewnotOldV(m1, 1:3);
		v2 = VOld(m2, 1:3);
		
		if v1(3) > v2(3)
			temp = v1;
			v1 = v2;
			v2 = temp;
		end
		indx = (v1(1:2) - v2(1:2))./v1(3);	
		
		if ((indx(1) <= v2(3)/v1(3)) && (indx(2) < v2(3)/v1(3))...
				&& (indx(1) >= -1) && (indx(2) > -1)) || ...
			((indx(2) <= v2(3)/v1(3)) && (indx(1) < v2(3)/v1(3))...
				&& (indx(2) >= -1) && (indx(1) > -1))						% 4-connectivity
			newRec = cat(1, newRec, [m2 m1 1]);
		end
		
	end
end
GNO = sparse(newRec(:,1), newRec(:,2), newRec(:,3), nOld, nNew);

G	= [GOld GNO; GNO' GNew];
