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

Description: MPC Tile Planner for aircraft model
%}
function [feas, uMPC, targetSetData, tileData, zf] = ...
	pdTilePlanV12(z0, tileData, VCell, cellScale)

% 1. tileData struct:
%		tileData.nodes	= nodes that make up the tile, ordered first to
%			last (given)
%		tileData.rad	= min radius of turn within each cell in the tile
%			(to be set in this function)
%		tileData.cVert	= coordinates of A,B,C,D vertices of second cell of
%			the tile, through which traj is to be planned (returned by
%			'targetSetAnalysisV04')
%		tileData.pType	= whether second cell has traversal across parallel
%			(1) or adjacent edges (2), (returned by
%			'targetSetAnalysisV04')
% 2. targetSetData struct:
% 		targetSetData.alfaSmp	= reachable sets for cells 2 through H+1
% 		targetSetData.btaSmp	= target sets for cells 2 through H+1; for
%			cell C = H+1, target set is ([-pi/2, pi/2] x entire exit segment)
% 		targetSetData.nSol		= number of points on each entry segment
%			for which reachable orientation set is computed
% 		targetSetData.nSmp		= number of grid points on each entry and
%			exit segments
% 		targetSetData.wL		= lower boundary of entry segments
% 		targetSetData.wU		= upper boundary of entry segments
% 		targetSetData.xSmp		= grid points on exit segments
% 		targetSetData.wSmp		= grid points on entry segments
%		targetSetData evaluated by "targetSetAnalysisV02"
% 3. VCell = location (1,2), dimensions (3) of cells, friction circle
%		data (ft = 4, fr = 5) speed limit data (vmin = 6, vmax = 7), and
%		elevation (8)
% 4. z0 = initial condition for tile plan

global g																	% gravitational acceleration, m/s2; % gravitational acceleration, m/s2
global vCruise

%----- Aircraft and flight parameters
% Altitude arbitrary, rest from NASA atmosphere model
altitude	= 1e3;															% m
temperature = 15.04 - 0.00649*altitude + 273.1;								% Kelvin
pressure	= 101.29*((temperature/288.08)^5.256);							% kPa
stdrhoD		= pressure/(0.2869*temperature);								% kg/m3

% Example 2.11 from B. N. Pamadi
modelParams.mass	= 5e4/g;												% mass, kg
modelParams.CD0		= 0.02;													% skin friction drag coefficient, dim-less
modelParams.K		= 0.04;													% induced drag coefficient, dim-less
modelParams.S		= 30;													% wing area, m2
modelParams.rhoD	= stdrhoD;												% air density, kg/m3
modelParams.CLmax	= 1.2;													% max lift coefficient, dim-less

mass	= modelParams.mass;
CD0		= modelParams.CD0;
K		= modelParams.K;
SAr		= modelParams.S;
rhoD	= modelParams.rhoD;
CLmax	= modelParams.CLmax;

% vCruise	= 100;																% Cruise speed
T0		= 0.5*rhoD*(vCruise^2)*SAr*CD0 + ...
	K*((mass*g)^2)/(0.5*rhoD*(vCruise^2)*SAr);								% Constant thrust

%----- Input constraints
phimin	= -45*pi/180;													% Min bank, rad
% phimin	= -20*pi/180;														% Min bank, rad
phimax	= min(20*pi/180, acos((mass*g)/(0.5*rhoD*(z0(4)^2)*SAr*CLmax)));	% Max bank, rad
phirad	= min(abs(phimin), abs(phimax));

%-------- Obtain target sets
% cellScale	= 1e3;															% Size of each fine res. cell
VCell(:, 1:3) = VCell(:, 1:3).*(cellScale);

C		= (numel(tileData.nodes) - 2);	tileData.rad = zeros(1, C+2);
tileData.rad(2:C+1) = (vCruise^2)/g/abs(tan(phirad));


% disp(tileData.nodes)
% 
% disp(tileData.rad)

tileData.rad	= tileData.rad/cellScale;
VCell(:, 1:3)	= VCell(:, 1:3)./(cellScale);
% disp(tileData.rad)

%-------- Find target configurations
[targetSetData, tileData] = targetSetAnalysisV03(tileData, VCell);

%-------- Scale target set numbers
targetSetData.wL	= (targetSetData.wL).*cellScale;
targetSetData.wU	= (targetSetData.wU).*cellScale;
targetSetData.xSmp	= (targetSetData.xSmp).*cellScale;
targetSetData.wSmp	= (targetSetData.wSmp).*cellScale;
targetSetData.wSol	= (targetSetData.wSol).*cellScale;

%-------- Inner polytope approximation to target set of first cell
[btaTermConf, thtaf]= innerPolytopeV01(targetSetData, tileData);			% Orientation constraint as function of length parameter
nbtaTerm			= size(btaTermConf, 1);

VCell(:, 1:3)		= VCell(:, 1:3).*(cellScale);

%-------- Determine terminal speed constraint
xMax= VCell(tileData.nodes(2), 3);	xMin= -(1e-1)*xMax;
yMax= VCell(tileData.nodes(2), 3);	yMin= 0;		

vfMax = VCell(tileData.nodes(2), 7); vfMin = VCell(tileData.nodes(2), 6);

% -------- Transform z0 to cell coordinates
Rtrans(:, :, 1)	= [0 -1; 1 0];												% Rotate -90
Rtrans(:, :, 2)	= [0 1; -1 0];												% Rotate 90
Rtrans(:, :, 3)	= [1 0; 0 -1];												% Flip H
Rtrans(:, :, 4)	= [-1 0; 0 1];												% Flip V
pos(1:2, 1) = z0(1:2, 1) - (VCell(tileData.nodes(2), 1:2))' - ...
	(VCell(tileData.nodes(2), 3))*[0.5;0.5];								% From inertial axes to history cell axes
ort			= z0(3);
S			= eye(2);
rTran		= tileData.rTran(1,:);
for m = 1:2
	if rTran(1, m) ~= 0
		if (m == 2) && (rTran(1, m) <= 2) && (rTran(1, m-1) > 2)
			S	= (Rtrans(:, :, rTran(1, m)))'*S;
		else
			S	= Rtrans(:, :, rTran(1, m))*S;
		end
	end	
	if m == 1
		switch rTran(1, m)
			case {1, 2}
				ort = ort - sign(rTran(1,m) - 1.5)*pi/2;
			case 3
				ort = -ort;
			case 4
				ort	= pi - ort;			
		end
	end
	if m == 2
		switch rTran(1, m)
			case {1, 2}
				if rTran(1, m-1) > 2
					ort = ort + sign(rTran(1,m) - 1.5)*pi/2;
				else
					ort = ort - sign(rTran(1,m) - 1.5)*pi/2;
				end					
			case 3
				ort = -ort;
			case 4
				ort	= pi - ort;			
		end
	end
end
pos = S*pos + (VCell(tileData.nodes(2), 3))*[0.5; 0.5];						% At this point, [pos ort] is the configuration in the axes
ort	= pi2pi(ort);															% system attached to the second cell of the given history
z0	= [pos; ort; z0(4)];

if (z0(1) < xMin) || (z0(1) > xMax)
	feas = 0;	uMPC = [];	zf = [];
	return;
end

% -------- Check if initial orientation within (first) reachable set
% w0Smp	= findSample(targetSetData.wSmp(1,:), z0(2));
% a0U		= targetSetData.alfaSmp(1, w0Smp);
% a0L		= targetSetData.alfaSmp(2, w0Smp);
% if (z0(3) < (a0L - 0.05*abs(a0L))) || ...
% 		(z0(3) > (a0U + 0.05*abs(a0U)))
% 	feas = 0;	uMPC = [];	zf = [];
% 	return;
% end

% -------- Initial MPC Parameters
mpcParams.nSt = 4;	mpcParams.lIn	= 2;	mpcParams.mOp = 4;				% Dimensions of state, input, and output spaces: particle dynamical model
dlta2	= (7.5e-2);
tS		= 1;																% Sampling time chosen arbitrarily

Qblk	= 0.0001*eye(mpcParams.mOp);	Qblk(3,3)= 1; Qblk(4, 4) = 0;
% Qblk	= 15*eye(mpcParams.mOp);	Qblk(3,3)= 3; Qblk(4, 4) = 0.015;
% Qblk	= 10*eye(mpcParams.mOp);	Qblk(3,3)= 5; Qblk(4, 4) = 0.010;
mpcParams.Hu0	= 5;	mpcParams.Hw0	= 1;								% Prediction and control horizons
mpcParams.tS	= tS;	mpcParams.HpEst = 1;
mpcParams.lsqOpt= optimset('Display', 'off');

% -------- State constraints, hard coded for nSt = 4 (particle dynamics)
mpcConstraints.vLim = VCell(tileData.nodes(2), 6:7);						% Speed constraints in cell
mpcConstraints.vfMax= vfMax;	mpcConstraints.vfMin= vfMin;
mpcConstraints.Gu	= [[-1; 1] zeros(2,3) [xMin; -xMax]; ...				% x >= xMin, x <= xMax 
	[0; 0] [-1; 1] zeros(2) [yMin; -yMax]; ...								% y >= yMin, y <= yMax
	zeros(2,3) [-1; 1] [mpcConstraints.vLim(1); -mpcConstraints.vLim(2)]];	% v >= vMin, v <= vMax

% -------- Appropriate target for MPC (not very crucial, decides "turns")
turnLR = zeros(C-1, 1); 
for m = 2:C	
	if tileData.pType(m) == 1, continue; end
	turnLR(m-1) = 1;
	for l = 1:2
		if tileData.rTran(m, l) > 2
			turnLR(m-1) = -turnLR(m-1);
		end
	end
end
for l = 1:2
	if tileData.rTran(1, l) > 2
		turnLR = -turnLR;
	end
end

initPosnDef = 0; turnMult = 2; %stdDef = (min(0.75, (C/4)^2))*xMax/2;
stdDef = (min(0.75, (C/4)^2))*xMax/2;
for m = 2:C
	turnMult = turnMult/2;
	initPosnDef = initPosnDef - turnLR(m-1)*turnMult*stdDef;
end

switch tileData.pType(1)
	case 1
		mpcParams.zf= [xMax; xMax/2 + initPosnDef; thtaf; vCruise];			% Target config.
		btaTerm		= [zeros(nbtaTerm,1) btaTermConf(:,1:2) ...
			zeros(nbtaTerm, (mpcParams.nSt - 3)) btaTermConf(:,3)];
% 		GTerm		= [[-1; 1] zeros(2,3) [xMax*(1 - 2*dlta2); -xMax]; ...	% xf >= xMax(1 - dlta2), xf <= xMax
% 			[dlta2; dlta2] [-1; 1] zeros(2) [yMin; -yMax]; ...				% y >= yMin, y <= yMax
% 			btaTerm];														% target set constraints
		GTerm		= [[-1; 1] zeros(2,3) [xMax*(1 - 2*dlta2); -xMax]; ...	% xf >= xMax(1 - dlta2), xf <= xMax
			[0; 0] [-1; 1] zeros(2) [yMin; -yMax]; ...				% y >= yMin, y <= yMax
			btaTerm];														% target set constraints
	case 2
		mpcParams.zf= [xMax/2 + initPosnDef; 0; thtaf; vCruise];
		btaTerm		= [btaTermConf(:,1) zeros(nbtaTerm,1) btaTermConf(:,2) ...
			zeros(nbtaTerm, (mpcParams.nSt - 3)) btaTermConf(:,3)];
% 		GTerm		= [[-1; 1] zeros(2,3) [xMin; -xMax]; ...				% xf >= xMin, xf <= xMax
% 			[0; 0] [-1; 1] zeros(2) [yMin; -yMin*(1 + 2*dlta2)];... 		% y >= yMin, y <= yMin + dlta2
% 			btaTerm];														% target set constraints
		GTerm		= [[-1; 1] zeros(2,3) [xMin; -xMax]; ...				% xf >= xMin, xf <= xMax
			[0; 0] [-1; 1] zeros(2) [yMin; -(yMin + 2*xMax*dlta2)];... 		% y >= yMin, y <= yMin + dlta2
			btaTerm];														% target set constraints
		mpcParams.Hu0	= 5;
end
mpcParams.Qblk	= Qblk; mpcParams.Rblk	= 0;								% Weights for quadratic cost
mpcConstraints.GTerm= GTerm;

%-------- MPC algorithm
u0	= [T0; 0]; uMPC = [u0; 0];	feas= 0;	k = 1;	z0k	= z0;	zf = [];
r	= (vCruise^2)/g/tan(phirad);
while (1)
	u0k		= [T0; 0];
	ukPrev	= uMPC(1:2, k) - u0k;
	
	mpcModel= linDiscModel(z0k, u0k, tS, mpcParams.nSt, mpcParams.lIn, modelParams);
		
% 	fprintf('MPC Iteration k = %i\n', k);
% 	fprintf('\tVehicle state z0k \t= (%f, %f, %f, %f)\n', z0k(1), z0k(2), z0k(3)*180/pi, z0k(4));
% 	fprintf('\tPrevious input \t\t= (%f, %f)\n', ukPrev(1), ukPrev(2));
	
	switch tileData.pType(1)
		case 1
			lTmp	= sqrt((xMax - z0k(1))^2 + max(((yMin - z0k(2))^2), ...
				((yMax - z0k(2))^2)));										% For estimating maximum Hp
			lTmp2	= abs(xMax - z0k(1) - dlta2*xMax);						% For estimating minimum Hp in no accel case
			yBar	= z0k(2) + (xMax - z0k(1))*tan(z0k(3));
			
			xyBS	= findSample(targetSetData.xSmp(1,:), yBar);
			btaUBS	= targetSetData.btaSmp(1,xyBS);
			btaLBS	= targetSetData.btaSmp(2,xyBS);
		case 2
			lTmp	= sqrt((yMin - z0k(2))^2 + max(((xMin - z0k(1))^2), ...
				((xMax - z0k(1))^2)));
			lTmp2	= abs(z0k(2) - dlta2*xMax);
			xBar	= z0k(1) - z0k(2)/tan(z0k(3));
			
			xyBS	= findSample(targetSetData.xSmp(1,:), xBar);
			btaUBS	= targetSetData.btaSmp(1,xyBS) - pi/2;
			btaLBS	= targetSetData.btaSmp(2,xyBS) - pi/2;
	end
	lMax			= 2*r*atan(lTmp/2/r); 
	mpcParams.HpMax = max(1, ceil(lMax/vfMin/tS));
% 	mpcParams.HpMax = 40;
	

	% ----- Run MPC with u1 = 0 (only bank angle control)
% 	fprintf('\tNo acceleration...\n');
	flipC = 1;
	for m = 1:2
		if (tileData.rTran(1,m) == 3) || (tileData.rTran(1,m) == 4)
			flipC = -flipC;
		end
	end
	if flipC == -1
		mpcConstraints.Fu	= [1 phimin; -1 -phimax];
	else
		mpcConstraints.Fu	= [1 -phimax; -1 phimin];						% Constraint on u2 with u1 = 0;
% 		mpcConstraints.Fu	= [1/phimax -1; -1/phimin 1];						% Constraint on u2 with u1 = 0;
% 		mpcConstraints.Fu	= [0 -1; 0 -1];
	end
	mpcParams.HpEst		= max([1 (floor(lTmp2/z0k(4)/tS) - 2) ...
		(mpcParams.HpEst - 1)]);

	[feas, mpcData, ~, mpcParams]= oneControlMPC(z0k, ...
		ukPrev(2), mpcConstraints, mpcModel, mpcParams, 2);
	
	if ~feas,
		break;
	end
	
	phik	= (mpcData.ukSer(1) + ukPrev(2));
	Tk		= 0.5*rhoD*(vCruise^2)*SAr*CD0 + ...
		K*((mass*g)^2)/(0.5*rhoD*(vCruise^2)*SAr*(cos(phik))^2);								% Constant thrust
	uMPC(:, k+1)= [Tk phik tS]';
	
	% ----- Step forward	
	uk	= uMPC(1:2, k+1);
	tSk = uMPC(3, k+1);
	for m = 1:2
		if (tileData.rTran(1,m) == 3) || (tileData.rTran(1,m) == 4)
			uMPC(2, k+1) = -uMPC(2, k+1);									% ***** THIS (AND ASSOCIATED CONSTRAINT BEFORE)
		end																	% CAN SCREW UP ASYMMETRIC CONSTRAINTS ******
	end
	k	= k+1;

	% ----- Simulate nonlinear model with MPC input
	[~, zkSer]	= ode45(@(t,z) mpcPlant09(t, z, uk, modelParams), [0 tSk], z0k);	
	z0k			= zkSer(end, :)';

	if (tileData.pType(1) == 1) && (abs(z0k(1) - xMax) <= 1.005*dlta2*xMax), break; end		% ... or HERE
	if (tileData.pType(1) == 2) && (abs(z0k(2)) <= 1.005*dlta2*xMax), break; end				% ... or HERE
end

%-------- Post processing
if feas
% 	fprintf('\tVehicle state z0k \t= (%f, %f, %f, %f)\n', z0k(1), z0k(2), z0k(3)*180/pi, z0k(4));
	
	%----- Last time step: move from within error band to boundary of cell
% 	uLast = uMPC(1:2,end);
	uLast = [T0; 0];
	[tLastSer, zLastSer] = ode45(@(t,z) mpcPlant09(t, z, uLast, ...
		modelParams), linspace(0, 3*tS, 25000), z0k);	
	switch tileData.pType(1)
		case 1
			[~,lastRow] = min(abs(zLastSer(((zLastSer(:,2) > 0) & ...
				(zLastSer(:,2) < yMax)), 1) - xMax));
		case 2
			[~,lastRow] = min(abs(zLastSer(((zLastSer(:,1) > 0) & ...
				(zLastSer(:,1) < xMax)) ,2)));
	end
% 	zf	= zLastSer(lastRow, :)';
	
% 	while ((tileData.pType(1) == 1) && ((zf(1) - xMax) < -(1e-3)*xMax)) ...
% 			|| ((tileData.pType(1) == 2) && (zf(2) > (1e-3)*xMax))
% 		[tLastSer, zLastSer] = ode45(@(t,z) mpcPlant09(t, z, uLast, ...
% 			modelParams), linspace(0,tS,10000), zf);
% 		switch tileData.pType(1)
% 			case 1
% 				[~,lastRow] = min(abs(zLastSer(:,1) - xMax));
% 			case 2
% 				[~,lastRow] = min(abs(zLastSer(:,2)));
% 		end
% 		zf	= zLastSer(lastRow, :)';
% 	end
	tSLast	= tLastSer(lastRow);
	
	if tSLast > 1e-9
		uMPC(:, k+1)= [uLast' tSLast]';		
	end
	uMPC(:, 1)		= [];
	zf	= zLastSer(lastRow, :)';
	
	if ((tileData.pType(1) == 1) && (abs(zf(1) - xMax) > (1e-3)*xMax)) ...
			|| ((tileData.pType(1) == 2) && (abs(zf(2)) > (1e-3)*xMax) )
		feas= 0;
		zf	= [];
		uMPC= [];
		return;
	end
	
	% -------- Transform zf back to inertial coordinates
	posF= zf(1:2, :) - (VCell(tileData.nodes(2), 3))/2;
	ortF= zf(3);  S	= eye(2);
	for m = 2:-1:1
		if rTran(1, m) ~= 0
			if (m == 2) && (rTran(1, m) <= 2) && (rTran(1, m-1) > 2)
				S	= Rtrans(:, :, rTran(1, m))*S;
			else
				S	= (Rtrans(:, :, rTran(1, m)))'*S;		
			end
		end
		if m == 2 
			switch rTran(1, m)
				case {1, 2}
					ortF	= ortF + sign(rTran(1, m) - 1.5)*pi/2;
				case 3
					ortF	= -ortF;
				case 4
					ortF	= pi - ortF;
			end
		end
		if m == 1
			switch rTran(1, m)			
				case {1, 2}
					ortF	= ortF + sign(rTran(1, m) - 1.5)*pi/2;
				case 3
					ortF	= -ortF;										% Technically, the same thing as below should be here,
																			% but in no transformation does FH arise first
				case 4
					if (rTran(1, m+1) == 1) || (rTran(1, m+1) == 2)			% If previously rotated, then the flip is about the
						ortF	= -ortF;									% intermediate axes, hence FV becomes FH (intermediate axes)
					else
						ortF	= pi - ortF;
					end
			end
		end
	end
	posF= S*posF + (VCell(tileData.nodes(2), 1:2))' + (VCell(tileData.nodes(2), 3))/2;
	zf	= [posF; pi2pi(ortF); zf(4)];
else
	uMPC = [];
end
% -------------------------------------------------------------------------
%**************************************************************************
function mpcModel = linDiscModel(z0k, u0k, tS, nSt, lIn, modelParams)

global g																	% gravitational acceleration, m/s2

mass= modelParams.mass;														% mass, kg
CD0 = modelParams.CD0;														% skin friction drag coefficient, dim-less
K	= modelParams.K;														% induced drag coefficient, dim-less
SAr	= modelParams.S;														% wing area, m2
rhoD= modelParams.rhoD;														% air density, kg/m3

T0	= u0k(1);
phi0= u0k(2);

v0	= z0k(4);
psi0= z0k(3);

Ac	= [0 0 -v0*sin(psi0) cos(psi0); ...										% Linearized particle dynamics
	0 0 v0*cos(psi0) sin(psi0); zeros(1,4);
	zeros(1,3) -(rhoD*v0*SAr*CD0 - 4*K*((mass*g)^2)/(rhoD*(v0^3)*SAr))/mass];
B1c	= [zeros(1,3) 1/mass; zeros(1,2) -g/v0 0]';
B2c = [v0*cos(psi0); v0*sin(psi0); -g*tan(phi0)/v0; (T0 - ...
	0.5*rhoD*(v0^2)*SAr*CD0 - ...
	K*((mass*g)^2)/(0.5*rhoD*(v0^2)*SAr*(cos(phi0)^2)))/mass] - ...
	Ac*z0k - B1c*u0k;
Cc	= eye(4);

tmpS= expm([[Ac B1c B2c]*tS; zeros((lIn+1), (nSt+lIn+1))]);
Ad	= tmpS(1:nSt, 1:nSt);
Bd	= tmpS(1:nSt, (nSt+1):(nSt+lIn+1));

B1d = Bd(:,1:end-1);
B2d = Bd(:,end);

mpcModel.Ad	= Ad;	mpcModel.B1d= B1d; 
mpcModel.B2d= B2d;	mpcModel.Cd	= Cc;
mpcModel.Ac	= Ac;	mpcModel.B1c= B1c; 
mpcModel.B2c= B2c;	mpcModel.Cc	= Cc;
% -------------------------------------------------------------------------
%**************************************************************************
function [Psi,Lam1,Lam2,Theta] = linPredictors(mpcModel, nSt, lIn, Hp, Hu)

Ad	= mpcModel.Ad;	B1d = mpcModel.B1d;
B2d = mpcModel.B2d;	Cd	= mpcModel.Cd;

Cdiag	= eye(nSt*Hp); % This is for Cd = eye(4); in general, Cdiag	= kron(eye(Hp), Cd);
Atil	= zeros(nSt*Hp, nSt);
B1til	= zeros(nSt*Hp, lIn);	B2til	= zeros(nSt*Hp, 1);
AiB1	= zeros(nSt,lIn);		AiB2	= zeros(nSt,1);
for m = 1:Hp
	AiB1 = AiB1 + (Ad^(m-1))*B1d;
	AiB2 = AiB2 + (Ad^(m-1))*B2d;
	
	Atil(((m-1)*nSt + 1):(m*nSt), :)	= Ad^m;
	B1til(((m-1)*nSt + 1):(m*nSt), :)	= AiB1;
	B2til(((m-1)*nSt + 1):(m*nSt), :)	= AiB2;	
end	
Psi		= Cdiag*Atil;														% nSt*Hp x nSt
Lam1	= Cdiag*B1til;														% nSt*Hp x lIn
Lam2	= Cdiag*B2til;														% nSt*Hp x 1

Ttil = zeros(nSt*Hp, lIn*Hp);
for m = 1:Hp
	Ttil(:, ((m-1)*lIn + 1):(m*lIn)) = ...
		[zeros((m-1)*nSt,lIn); B1til((1:(Hp-m+1)*nSt),:)];
end
Ttil	= Ttil(:, 1:lIn*Hu);
Theta	= Cdiag*Ttil;														% nSt*Hp x lIn*Hu
% -------------------------------------------------------------------------
%**************************************************************************
function [feas, mpcData, zPredf, mpcParams] = oneControlMPC(z0k, ukPrev, ...
	mpcConstraints, mpcModel, mpcParams, uID)

HpEst	= mpcParams.HpEst;	Hp	= HpEst;	foundBestHp = 0;
nSt		= mpcParams.nSt;	lIn = 1;		mOp = mpcParams.mOp;
mpcModel.B1d = mpcModel.B1d(:,uID);

while (~foundBestHp) && (Hp <= mpcParams.HpMax)
	Hu	= min(mpcParams.Hu0, Hp);		Hw	= min(mpcParams.Hw0, Hp);		% Set control horizon and cost window
	mpcData.Hu = Hu;					mpcData.Hw = Hw;
	
% 	fprintf('\t\tHp = %i\n', Hp);
	
	
% 	Q	= zeros((Hp-Hw)*mOp);	Q = blkdiag(Q, mpcParams.Qblk);
% 	Q	= [zeros((Hp-Hw)*mOp) zeros((Hp-Hw)*mOp, mOp); ...
% 		zeros(mOp, (Hp-Hw)*mOp) mpcParams.Qblk];
	Q	= kron(eye(Hp-Hw)*mOp, mpcParams.Qblk); Q = blkdiag(Q, mpcParams.Qblk);
	R	= kron(eye(Hu), mpcParams.Rblk);
	S_Q	= sqrt(Q);				S_R	= sqrt(R);
	
	T	= (kron(mpcParams.zf',ones(Hp-Hw+1,1)))';	T	= T(:);
	
	% --- Set up constraint multiplier matrices	
	F	= [kron(eye(Hu), mpcConstraints.Fu(:,1:lIn)) ...
		kron(ones(Hu,1), mpcConstraints.Fu(:,end))];
	
	Fsc	= [];																% Input constraint
	for m = 1:Hu															% ***** REPLACE WITH kron *****
		Fscm	= zeros(size(F,1), lIn);		
		for l = m:Hu
			Fscm = Fscm + F(:,(l-1)*lIn+1:l*lIn);			
		end
		Fsc	= cat(2, Fsc, Fscm);		
		if m == 1, Fsc1 = Fscm; end
	end
	fLast	= F(:, end);	

	G	= [kron(eye(Hp-Hw), mpcConstraints.Gu(:,1:mOp))...
		kron([zeros(Hp-Hw,mOp) ones(Hp-Hw,1)], mpcConstraints.Gu(:,end))];	% State constraint	
	G	= cat(1, G, [kron([zeros(1,(Hp-Hw)) 1], ...
		mpcConstraints.GTerm(:,1:mOp)) mpcConstraints.GTerm(:,end)]);
	Gam	= G(:, 1:end-1);	gLast = G(:,end);		

	[mpcData.Psi, mpcData.Lam1, mpcData.Lam2, mpcData.Theta] ...
		= linPredictors(mpcModel, nSt, lIn, Hp, Hu);
	
% 	[~, LamU1, ~, ~] = linPredictors(mpcModel2, nSt, lIn, Hp, Hu);
	Ups	= mpcData.Psi*z0k + mpcData.Lam1*ukPrev + mpcData.Lam2;% + LamU1*uk1Prev;
	mpcData.Ups = Ups;

	E	= (T - Ups);
	for k2 = 3:nSt:size(E,1)
		E(k2) = pi2pi(E(k2));
	end

	% --- Run LSq minimization with terminal constraints
	OmegaL	= [Fsc; Gam*mpcData.Theta];
	omegaR	= [-Fsc1*ukPrev - fLast; -Gam*Ups - gLast];

	[mpcData.ukSer, ~, ~, exitFlag] = lsqlin([S_Q*mpcData.Theta; S_R], ...
		[S_Q*E; zeros(lIn*Hu, 1)], OmegaL, omegaR, ...
		[], [], [], [], [], mpcParams.lsqOpt);

	switch exitFlag
		case 1																% Feasible with terminal constraints
% 			fprintf('\t\t*** Feasible with TC, \tHp = %i\n', Hp);
			foundBestHp = 1;			
			feas	= 1;				
		case -2																% Infeasible with terminal constraints
			Hp	= Hp + 1;
			feas= 0;
		otherwise			
			error('Unknown exitFlag (with terminal constraint).');
	end
end
if feas
	zPred	= mpcData.Ups + mpcData.Theta*mpcData.ukSer;					% Predicted state trajectory with MPC input
	zPred	= reshape(zPred, 4, Hp - mpcData.Hw + 1);
	zPredf	= zPred(:,end);
	mpcParams.HpEst = Hp;
else
	zPredf = [];
end
% -------------------------------------------------------------------------