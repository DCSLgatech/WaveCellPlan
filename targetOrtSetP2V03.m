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

Description: arget set computation, Problem 2 (across adjacent edges)
%}
function [alfa, thtaPSmp, thtaMSmp, xFeasSmp] = ...
	targetOrtSetP2V03(w, d, r, xSmp, btaSmp)
%--------------------------------------------------------------------------

% -------- Intersect with segment with non-Inf bta
[~, btaNotInfIndx] = remInf(btaSmp(1,:));
if ~numel(btaNotInfIndx)
	alfa	= [-Inf; Inf];
	thtaPSmp= [];			thtaMSmp = [];
	xFeasSmp= [];
	return;
end

% -------- Portion of exit segment reachable by C arcs alone
% Rightmost point directly reachable by C- assumed to be C
% Leftmost point directly reachable by C+ assumed to be D 
% Rightmost point directly reachable by C+ assumed to be C
n1L		= r - sqrt(r^2 - w^2);												% Leftmost point directly reachable by C-
smpN1L	= min(max((findSample(xSmp, n1L)), btaNotInfIndx(1)), ...
	btaNotInfIndx(end));
smpN2U	= btaNotInfIndx(end);

if smpN2U <= smpN1L
	alfa	= [-Inf; Inf];
	thtaPSmp= [];			thtaMSmp = [];
	xFeasSmp= [];
	return;
end

% [smpN1L xSmp(smpN1L) n1L]
% -------- Calculate angles g+ and g-
gamPSmp = zeros(1, numel(btaNotInfIndx));									% Lowest angle possible at X, by C+
gamMSmp = zeros(1, numel(btaNotInfIndx));									% Highest angle possible at X, by C- 

for m = btaNotInfIndx(1):(smpN1L-1)
	x	= xSmp(m);
	C0	= x^2 + w^2;
	C1	= sqrt(4*(r^2)*C0 - C0^2);
	C2	= (2*r*x + C0);
	gamP= 2*atan((-2*r*w + C1)/C2);											% Angle at X of C+ arc from W to X	
	gamM= acos(1 - x/r);
	
	gamPSmp(m - btaNotInfIndx(1) + 1) = gamP;
	gamMSmp(m - btaNotInfIndx(1) + 1) = gamM;
end

for m = smpN1L:smpN2U
	x	= xSmp(m);
	C0	= x^2 + w^2;
	C1	= sqrt(4*(r^2)*C0 - C0^2);
	C2	= (2*r*x + C0);
	gamP= 2*atan((-2*r*w + C1)/C2);											% Angle at X of C+ arc from W to X
	gamM= pi2pi(pi + 2*atan((-2*r*w -C1)/C2));								% Angle at X of C- arc from W to X
	
	gamPSmp(m - btaNotInfIndx(1) + 1) = gamP;
	gamMSmp(m - btaNotInfIndx(1) + 1) = gamM;
end
% gamMSmp*180/pi

% -------- Find critical points (works because of continuity)
diffPU	= gamPSmp - btaSmp(1,btaNotInfIndx);	diffPL	= gamPSmp - btaSmp(2,btaNotInfIndx);
diffML	= gamMSmp - btaSmp(2,btaNotInfIndx);	diffMU	= gamMSmp - btaSmp(1,btaNotInfIndx);

xPU		= findZeros(diffPU);
xML		= findZeros(diffML);
xC		= [0 sort([xPU xML]) numel(btaNotInfIndx)];

% -------- Find type of interval
	% Note: because of continuity, one sample of each interval suffices to
	% figure out which "category" the whole interval belongs to
intMark = [];
for m = 1:(numel(xC)-1)
	if (diffPU(xC(m)+1) > 0) || (diffML(xC(m)+1) < 0)						% g+ > bta_U or g- < bta_L; no solution
		intMark = cat(2, intMark, 0);
	else
		intMark = cat(2, intMark, 1);
	end
end

% -------- Find largest interval in N1-N2
if intMark(1), intsGood = [xC(1)+1 xC(2)]; else intsGood = []; end
for m = 2:(numel(xC)-1)
	if intMark(m) && intMark(m-1)
		intsGood(end,2) = xC(m+1);
	elseif intMark(m) && ~intMark(m-1)
		intsGood		= cat(1, intsGood, [xC(m)+1 xC(m+1)]);
	end
end
if ~numel(intsGood)
	alfa = [-Inf; Inf];
	thtaPSmp = []; thtaMSmp = [];
	xFeasSmp = [];
	return;
end

intsGood = cat(2, intsGood, intsGood(:,2) - intsGood(:,1));
intsGood = sortrows(intsGood, 3);
xIntMaxSmp1 = intsGood(1, 1);	xIntMaxSmp2 = intsGood(1, 2);
xFeasSmp =  (xIntMaxSmp1:xIntMaxSmp2) + btaNotInfIndx(1) - 1;

% --------- Max and min of Ux' and Lx' on interval found above
intSmp	= xIntMaxSmp1:xIntMaxSmp2;

if ~numel(intSmp)
	alfa	= [-Inf; Inf];
	thtaPSmp= [];			thtaMSmp = [];
	xFeasSmp= [];
	return;
end

thtaPSmp= zeros(numel(intSmp),1);
thtaMSmp= zeros(numel(intSmp),1);
% intSmp

n = 0;
for m = intSmp
	n = n + 1;
	x = xSmp(m + btaNotInfIndx(1) - 1);
	if diffPL(m) < 0
		bta	= btaSmp(2,m + btaNotInfIndx(1) - 1);
		C0	= x^2 + w^2;
		A	= -cos(bta) - x/r;
		B	= sin(bta) - w/r;
		C	= 1 - (x*cos(bta) - w*sin(bta))/r - C0/(2*r^2);
		thtaP	= 2*atan((A + sqrt(A^2 + B^2 - C^2))/(B+C));				% Tangent angle of C+C- path at W
	else
		thtaP	= gamMSmp(m) - pi/2;										% Tangent angle of C+ path at W
	end
	thtaPSmp(n) = thtaP;
	
	if diffMU(m) > 0		
		bta = btaSmp(1,m + btaNotInfIndx(1) - 1);		
		C0	= x^2 + w^2;
		A	= -cos(bta) + x/r;
		B	= sin(bta) + w/r;
		C	= 1 + (x*cos(bta) - w*sin(bta))/r - C0/(2*r^2);
		thtaM	= max(-pi/2, 2*atan((A - sqrt(A^2 + B^2 - C^2))/(B+C)));	% Tangent angle of C-C+ path at W
	else
		if (m + btaNotInfIndx(1) - 1) >= smpN1L
			thtaM = gamPSmp(m) - pi/2;										% Tangent angle of C- path at W
		else
			thtaM = -pi/2;
		end
	end
	thtaMSmp(n) = thtaM;
end

alfa(1) = max(thtaPSmp);
alfa(2) = min(thtaMSmp);
%--------------------------------------------------------------------------