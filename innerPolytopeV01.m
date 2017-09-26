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

Description: Inner convex polytope approximation to target set
%}

function [btaTerm, thtaf] = innerPolytopeV01(targetSetData, tileData)
% function btaTerm = innerPolytopeV01(targetSetData, tileData)

xSmp	= targetSetData.xSmp(1,:);
alfaSmp = targetSetData.btaSmp(1:2, :);

% -------- Transform target set orientation to inertial coordinates
if tileData.pType(1) == 2
	alfaSmp = alfaSmp - pi/2;
end

% -------- Least squares fit with constraint of 'above' or 'below'
[lfitd1, notInfIndx]	= remInf(alfaSmp(1,:));
lfitd2					= alfaSmp(2, notInfIndx);
lfitC					= [(xSmp(notInfIndx))' ones(numel(notInfIndx),1)];

if ~numel(notInfIndx)
	btaTerm = [-1 0 -1; 1 0 1];												% Infeasible constraints
	thtaf	= Inf;
	return;
end

lsqOpt		= optimset('Display', 'off');
alfaLfit1	= lsqlin(lfitC, lfitd1', lfitC, lfitd1', ...
	[], [], [], [], [], lsqOpt);											% Upper bound
alfaLfit2	= lsqlin(lfitC, lfitd2', -lfitC, -lfitd2', ...
	[], [], [], [], [], lsqOpt);											% Lower bound
% 
% alfaLfit1	= lsqlin(lfitC, lfitd1', [], [], ...
% 	[], [], [], [], [], lsqOpt);											% Upper bound
% alfaLfit2	= lsqlin(lfitC, lfitd2', [], [], ...
% 	[], [], [], [], [], lsqOpt);											% Lower bound

% alfaL		= [lfitC*alfaLfit1 lfitC*alfaLfit2];
% figure;  grid on; hold on;
% plot((xSmp(notInfIndx))', alfaL(:,1)*180/pi, 'Color', [0 0 0.5], 'LineWidth', 2, 'LineStyle', '--');
% plot(xSmp, alfaSmp(1,:)*180/pi, 'Color', [0 0 0.5], 'LineWidth', 2);
% plot((xSmp(notInfIndx))', alfaL(:,2)*180/pi, 'Color', [0 0.5 0], 'LineWidth', 2, 'LineStyle', '--');
% plot(xSmp, alfaSmp(2,:)*180/pi, 'Color', [0 0.5 0], 'LineWidth', 2);

btaTerm	= [-alfaLfit1(1) 1 -alfaLfit1(2); alfaLfit2(1) -1 alfaLfit2(2); ...
	-1 0 xSmp(notInfIndx(1)); 1 0 -xSmp(notInfIndx(end))];
sc		= (xSmp(notInfIndx(1)) + xSmp(notInfIndx(end)))/2;
thtaf	= sc*(alfaLfit1(1) + alfaLfit2(1))/2  + (alfaLfit1(2) + alfaLfit2(2))/2;