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

Description: Aircraft point mass model (reduced planar motion)
%}
function zdot = mpcPlant09(~,z,u,modelParams)

% z = (x, y, psi, v)
% u = (T, phi)
% x,y	: inertial position
% psi	: heading (same as orientation, zero sideslip assumed)

global g																	% gravitational acceleration, m/s2

%----- Aircraft and flight parameters
mass	= modelParams.mass;													% mass, kg
CD0		= modelParams.CD0;													% skin friction drag coefficient, dim-less
K		= modelParams.K;													% induced drag coefficient, dim-less
S		= modelParams.S;													% wing area, m2
rhoD	= modelParams.rhoD;													% air density, kg/m3

%----- States in readable form
psi	= z(3);	v = z(4);

%----- Control inputs in readable form
T = u(1);	phi = u(2);

%----- Equations of motion
zdot(1,1) = v*cos(psi);
zdot(2,1) = v*sin(psi);
zdot(3,1) = -g*tan(phi)/v;
zdot(4,1) = (T - 0.5*rhoD*(v^2)*S*CD0 - ...
	K*((mass*g)^2)/(0.5*rhoD*(v^2)*S*(cos(phi)^2)))/mass;