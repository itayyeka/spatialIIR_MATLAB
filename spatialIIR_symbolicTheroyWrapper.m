function [] = spatialIIR_symbolicTheroyWrapper(cfgIn,scriptFlags)
	%% configure	
	try
		scriptFlags.crlbSimpleSim;
	catch
		scriptFlags.crlbSimpleSim 	= 1;
	end
	try
		nSensors = cfgIn.nSensors;
	catch
		nSensors = 2;
	end
	
	%% symbolics
	alpha 				= sym('alpha',[nSensors 1]);
	alphaT 				= transpose(alpha);
	alphaH				= conj(alphaT);
	beta 				= sym('beta',[nSensors 1]);
	betaT 				= transpose(beta);
	betaH				= conj(betaT);
	syms 				omega;
	assume(				omega, 'real');
	syms 				tau;
	assume(				tau, 'real');
	syms 				c;
	assume(				c, 'real');
	syms 				D;
	assume(				D, 'real');
	syms 				theta;
	assume(				theta, 'real');
	tauTheta 			= omega*D*cos(theta)/c;
	sensorID_zeroBased 	= 0:(nSensors-1));
	sensorID_oneBased 	= 1:nSensors;
	d 					= exp(-1i*tauTheta*();
	dT 					= trnspose(d);
	dH 					= conj(dT);
	g 					= betaT*d*exp(-1i*omega*tau);
	A 					= -1i*omega*D*sin(theta)*diag(sensorID_zeroBased)/c;
	AT 					= transpose(A);
	AH 					= conj(AT); 
	B 					= d*dT*A-A*d*dT;
	BT 					= transpose(B);
	BH 					= conj(BT);
	
	crlbVariables.nSensors 				= nSensors;
	crlbVariables.alpha 				= alpha;
	crlbVariables.alphaT 				= alphaT;
	crlbVariables.alphaH 				= alphaH;
	crlbVariables.beta 					= beta;
	crlbVariables.betaT 				= betaT;
	crlbVariables.betaH 				= betaH;
	crlbVariables.omega 				= beta;
	crlbVariables.tau 					= tau;
	crlbVariables.c 					= c;
	crlbVariables.D 					= D;
	crlbVariables.theta 				= beta;
	crlbVariables.tauTheta 				= tauTheta;
	crlbVariables.sensorID_zeroBased 	= sensorID_zeroBased;
	crlbVariables.sensorID_oneBased 	= sensorID_oneBased;
	crlbVariables.d 					= d;
	crlbVariables.dT 					= dT;
	crlbVariables.dH 					= dH;
	crlbVariables.g 					= g;
	crlbVariables.A 					= A;
	crlbVariables.AT 					= AT;
	crlbVariables.AH 					= AH;
	crlbVariables.B 					= B;
	crlbVariables.BT 					= BT;
	crlbVariables.BH 					= BH;	
	
	%% sub scripts
	if scriptFlags.crlbSimpleSim
		crlbSimpleSim_output = crlbSimpleSim(crlbVariables,cfgIn);
	end
end