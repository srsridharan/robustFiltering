
classdef (ConstructOnLoad=false) codFreeFilterClass


	properties(SetAccess = private)
	%% the values below are simply initalizations.. they are never used since they are replaced in the constructor 
	% function	below 

		limitOnNumberOfQuadratics = 20; 
		QuadraticStored = [];
		systemSpecificParam = {};	 % next ver:  can change this to be an object of the system class
		stateEstimate = [];
		covarEstimate = [];
		filterWeightOnMeasNoise =  [];  % this is the Qw parameter
		filterWeightOnDynNoise =  []; % this is the R parameter
		filterWeightOnState = [];
		currentNumberOfQuadratics = 0;
		currentApproxToOutputPositive=[];
		currentApproxToOutputNegative=[];
		currentApproxToOutputSquared=[];
		outputFnToFit=[];
	end % of properties
	methods
	
	function obj =	codFreeFilterClass(systemSpecificParam,filterWeightOnMeasNoise,filterWeightOnDynNoise, filterWeightOnState,...
											initialStateEstimate,initialStateCovarEstimate, initialQuadratic,outputFnToFit)
		%% This is the constructor  for the codFreeFilterClass. It performs argument checking (number and nonempty) and initializes some of 
		% the class members
		
		
		errorIDforThisFunction = 'robustFilteringCODfree:constructor';	
%         keyboard;
		% argument number check
        try
            error(nargchk(8, 8, nargin,'struct'))
        catch excep1
            exception2 = MException(errorIDforThisFunction, 'incorrect number of inputs to constructor' );
			excep2 = addCause(excep1,exception2); 
			throw(excep2);
        end
			
		
		% check for empty arguments	
			argList = {'systemSpecificParam','filterWeightOnMeasNoise','filterWeightOnDynNoise',' filterWeightOnState',...
											'initialStateEstimate','initialStateCovarEstimate',' initialQuadratic','outputFnToFit'};	
			emptyChecker = cellfun(@(x)isempty(x),{systemSpecificParam,filterWeightOnMeasNoise,filterWeightOnDynNoise, filterWeightOnState,...
											initialStateEstimate,initialStateCovarEstimate, initialQuadratic,outputFnToFit},'UniformOutput',1);
			assert(~sum(emptyChecker),['empty argument(s): ', argList(logical(emptyChecker))]);
				
			

		% now initialize the arguments
				obj.systemSpecificParam = systemSpecificParam; % these contain the sys dynamics.
				obj.filterWeightOnMeasNoise =  filterWeightOnMeasNoise;
				obj.filterWeightOnDynNoise = filterWeightOnDynNoise;
				obj.filterWeightOnState =  filterWeightOnState;
				obj.stateEstimate = initialStateEstimate; 
				obj.covarEstimate = initialStateCovarEstimate;
				obj.QuadraticStored{1} = initialQuadratic;
				obj.currentNumberOfQuadratics = 1;
				obj.outputFnToFit = outputFnToFit;
				obj.currentApproxToOutputPositive = nan; 
				obj.currentApproxToOutputNegative = nan;
				obj.currentApproxToOutputSquared = nan;

  end % of constructor codFreeFilterClass

  
	%% Now we have the functions that are used to set parameter values
    function obj = set.systemSpecificParam(obj,value)
				errorIDforThisFunction = 'robustFilteringCODfree:systemSpecificParam';
				assert(~isempty(value),[errorIDforThisFunction,': value empty in assignment']);
				assert(iscell(value),[errorIDforThisFunction,': assigned value not a cell']);
        
				obj.systemSpecificParam = value;
    end
    function obj = set.filterWeightOnMeasNoise(obj,value)
      	errorIDforThisFunction = 'robustFilteringCODfree:filterWeightOnMeasNoise';
				assert(~isempty(value),[errorIDforThisFunction,': value empty in assignment']);
			  obj.filterWeightOnMeasNoise = value;
    end
    function obj = set.filterWeightOnDynNoise(obj,value)
        errorIDforThisFunction = 'robustFilteringCODfree:filterWeightOnDynNoise';
				assert(~isempty(value),[errorIDforThisFunction,': value empty in assignment']);
				obj.filterWeightOnDynNoise = value;
    end
    function obj = set.filterWeightOnState(obj,value)
        errorIDforThisFunction = 'robustFilteringCODfree:filterWeightOnState';
				assert(~isempty(value),[errorIDforThisFunction,': value empty in assignment']);
				obj.filterWeightOnState = value;
    end
    function obj = set.stateEstimate(obj,value)
        errorIDforThisFunction = 'robustFilteringCODfree:stateEstimate';
				assert(~isempty(value),[errorIDforThisFunction,': value empty in assignment']);
				obj.stateEstimate = value;
    end
    function obj = set.covarEstimate(obj,value)
				errorIDforThisFunction = 'robustFilteringCODfree:covarEstimate';
				assert(~isempty(value),[errorIDforThisFunction,': value empty in assignment']);
        obj.covarEstimate = value;
    end
    function obj = set.QuadraticStored(obj,value)
				errorIDforThisFunction = 'robustFilteringCODfree:QuadraticStored';
				assert(~isempty(value),[errorIDforThisFunction,': value empty in assignment']);
        obj.QuadraticStored = value;
    end
    function obj = set.currentNumberOfQuadratics(obj,value)
				errorIDforThisFunction = 'robustFilteringCODfree:currentNumberOfQuadratics';
				assert(~isempty(value),[errorIDforThisFunction,': value empty in assignment']);
        obj.currentNumberOfQuadratics = value;
    end
    
    function obj = set.currentApproxToOutputPositive(obj,value)
				errorIDforThisFunction = 'robustFilteringCODfree:currentApproxToOutputPositive';
				assert(~isempty(value),[errorIDforThisFunction,': value empty in assignment']);
        obj.currentApproxToOutputPositive = value;
    end
    
    function obj = set.currentApproxToOutputNegative(obj,value)
				errorIDforThisFunction = 'robustFilteringCODfree:currentApproxToOutputNegative';
				assert(~isempty(value),[errorIDforThisFunction,': value empty in assignment']);
        obj.currentApproxToOutputNegative = value;
    end
    
    function obj = set.currentApproxToOutputSquared(obj,value)
				errorIDforThisFunction = 'robustFilteringCODfree:currentApproxToOutputSquared';
				assert(~isempty(value),[errorIDforThisFunction,': value empty in assignment']);
        obj.currentApproxToOutputSquared = value;
    end
    
    function obj = set.outputFnToFit(obj,value)
				errorIDforThisFunction = 'robustFilteringCODfree:filterWeightOnMeasNoise';
				assert(~isempty(value),[errorIDforThisFunction,': value empty in assignment']);
        obj.outputFnToFit = value;
    end
   
	% end of  functions that are used to set parameter values
    
	
	function newSetOfValFnQuadratics = generateQuadraticsInNextTimeStep(obj,yval)
			%% this function creates the quadratics at the next time step from the quadratics at the current time step

			errorIDforThisFunction = 'robustFilteringCODfree:generateQuadraticsInNextTimeStep:';	
			
			counterToBuildQuadratics = 1;
			newSetOfValFnQuadratics = {}; % this will contain the quadratic terms which will be generated from
						% the current value fn quadratics and the output quadratics (ref paper)

			Atilde = obj.systemSpecificParam{1};
			Btilde = obj.systemSpecificParam{2};


			try
				% now for each of the quadratics (available in the current time step) in the approxmation of the value function 				
				% obtain the new/updated quadratic
				
				for(counterValFnQuadratic = 1:obj.currentNumberOfQuadratics)
						% decompose the valueFnQuadratic in the current loop iteration into its components	
							[N,L,Phi] = helperNLPfromQuadratic(obj.QuadraticStored{counterValFnQuadratic});		
							
						% now find the optimal value of the noise  (derived in the paper)
								wcommonmultFactor =  -inv(Btilde'*N'*Btilde + obj.filterWeightOnMeasNoise');
							% first the constant part	 of the optimal measurement noise estimate
								s = wcommonmultFactor * (Btilde'*N'*Atilde);
							
							% now the linear part of the measurement noise	
								z = wcommonmultFactor * (Btilde' * L');


							A1 = Atilde + Btilde *s;
								
							Ntilde = A1'*N*A1 + s'*obj.filterWeightOnMeasNoise*s; % the new covariance estimate (purely quadratic term)
							Ltilde = z'*Btilde'*N*A1  + L*A1 + z'*obj.filterWeightOnMeasNoise*s; % the linear term 
							Phitilde = Phi + yval'*obj.filterWeightOnDynNoise*yval + 2*L*Btilde*z + z'*Btilde'*N*Btilde*z + z'*obj.filterWeightOnMeasNoise*z;
								%the constant term
							% now use the above to form the quadratic kernal
							Qtilde = [Ntilde Ltilde';Ltilde  Phitilde];
							
							% first generate the output function approximation for the current value of the output
								qapproxforOutput = obj.generateQuadOutputApprox(yval);

							% now run through a for loop to  generate an approximation for Qtilde

							for(counteroutputQuadratic = 1: length(qapproxforOutput))
								newSetOfValFnQuadratics{counterToBuildQuadratics} = qapproxforOutput{counteroutputQuadratic}...
																									+ Qtilde;
								counterToBuildQuadratics = counterToBuildQuadratics + 1;	
							end
						
				end % of counterValFnQuadratic loop

			catch excep1
				if(isempty(obj.QuadraticStored{counterValFnQuadratic}))
					exception2 = MException(errorIDforThisFunction, 'emptyQuadraticAttemptingToIterateToNextStep' );
					excep2 = addcause(exception2,excep1); 
					throw(excep2);
				else
					try 
						getInv = inv(Btilde'*N'*Btilde + obj.filterWeightOnMeasNoise');
					catch inverr
						exception2 = addcause(inverror,excep1);
						throw(exception2);

					end
					rethrow(excep1);

				end % of if 

			end % of the try catch block

	end  % end of fn: generateQuadraticsInNextTimeStep

	function qapproxForOutput = generateQuadOutputApprox(obj,yval)
		% Given a output scalar yval, this function generates the corresponding (appropriately chosen) set of
		% quadratics that approximate the function  - y C(x) +  (1/2)C^2(x). For more details on the theory, ref to the paper

			qapproxForOutput ={};

			ctrforqgrowth = 1;
			if(yval>=0)
				qcoutputapproxtouse = obj.currentApproxToOutputPositive;
			else
				qcoutputapproxtouse = obj.currentApproxToOutputNegative;
			end
			
			qcsqapproxtouse = obj.currentApproxToOutputSquared;
				
			for(qcapproxctr = 1:length(qcoutputapproxtouse))
				for(qcsqapproxctr = 1:length(qcsqapproxtouse))
					qapproxForOutput{ctrforqgrowth} =  obj.filterWeightOnDynNoise*qcsqapproxtouse{qcsqapproxctr} + ...
											2*abs(yval)*obj.filterWeightOnDynNoise*qcoutputapproxtouse{qcapproxctr};
					ctrforqgrowth = ctrforqgrowth + 1;
				end % of for qcsqapproxctr
			end % of for qcapproxctr

	end  % end of genQuadOutputApprox function



	function obj = pruneQuadraticsAndGenerateNewEst(obj,newSetOfValFnQuadratics)
			% this function prunes the existing quadratics in order to obtain the needed number of quadratics. During this
			% process it also generates estimates for the state and covarince matrices.
			errorIDforThisFunction = 'robustFilteringCODfree:pruneQuadraticsAndGenerateNewEst:';	
	

			% check if there are a new set of Value function quadratics
				assert(~isempty(newSetOfValFnQuadratics),[errorIDforThisFunction,': the new set of quadratics is empty']);
				

				[minval,stateest] = cellfun(@(q) helpergenerateMinimaAndCorrespondingStateEst(q),...
															newSetOfValFnQuadratics,'UniformOutput',false);
				CoordinatesOfPointsToCluster = cell2mat(stateest)'; 

			% set the number of quadratics to retain. 
				numberOfQuadraticsToBeKept = min([obj.limitOnNumberOfQuadratics,length(newSetOfValFnQuadratics)]);
				[clusteridx, clusterCentroids] = kmeans(CoordinatesOfPointsToCluster,numberOfQuadraticsToBeKept, 'emptyaction','singleton');
			 
			% now for each cluster find and retain the quadratic with the lowest value
				for(klusterindxCounter = 1:numberOfQuadraticsToBeKept)
						rowselect= find(clusteridx == klusterindxCounter);
						[val,indxsort] = sort(cell2mat(minval(rowselect)));
						indxOfQuadraticsToKeep(klusterindxCounter) = rowselect(indxsort(1));
				end

	 		% Now sort and choose the intra cluster minima (from each cluster) and store them
				[val,indxsort ] = sort(cell2mat(minval(indxOfQuadraticsToKeep)));
				obj.stateEstimate  = stateest{indxOfQuadraticsToKeep(indxsort(1))};
				obj.currentNumberOfQuadratics = obj.limitOnNumberOfQuadratics; 
				obj.QuadraticStored = newSetOfValFnQuadratics(indxOfQuadraticsToKeep);

	end %of function pruneQuadratics



	function obj= generate_quadapprox_script(obj,xhat) 
			%% This function returns the quadratics which approximate a desired output fn and its square at the current state (ref paper)
			%% TO DO: the interface to this function can be further cleaned up. It can simply call an
			% external fn 
			fntofit =obj.outputFnToFit;  
			numberOfQuadratics = 7;
			numOfPtsAroundEachQuadratic = 10;
			rangeAroundXhat = 2;
			obj.currentApproxToOutputNegative = fitForFn(fntofit, xhat, numberOfQuadratics, numOfPtsAroundEachQuadratic, rangeAroundXhat );


			fntofit =@(x)(-1.*obj.outputFnToFit(x)); %@(x)  -x.^3/40;
			obj.currentApproxToOutputPositive = fitForFn(fntofit, xhat, numberOfQuadratics, numOfPtsAroundEachQuadratic, rangeAroundXhat );

			fntofit =@(x) obj.outputFnToFit(x).^2 ; %  (x.^3/40).^2;
			obj.currentApproxToOutputSquared = fitForFnSq(fntofit, xhat, numberOfQuadratics, numOfPtsAroundEachQuadratic, rangeAroundXhat );

	end % of function generate_quadapprox_script
		
	end % of methods 
		
	
end % end of class def



%% Helper functions

function [minval,stateest]= helpergenerateMinimaAndCorrespondingStateEst(Q)
		% this function generates the minima of the quadratic (having a kernal Q) and the corresponding state 
		% xhatmin. The value of the quadratic form for any state xhatmin is [xhatmin;1]'*Q*[xhatmin;1]
		
			q11 = Q(1:end-1,1:end-1);
			q21=Q(end,1:end-1);
			q12 = Q(1:end-1,end);
			q22 = Q(end,end);
			xhat0 = -inv(q11+q11')*(q12+ q21');
			stateest = xhat0;
 			minval = 0.5*[stateest;1]'*Q*[stateest;1];
			return;

end % of function generateMinimaAndCorrespondingStateEst


function [N,L,Phi] = helperNLPfromQuadratic(Q)
	%% this function decomposes  the kernal of a quadratic form  into its constituent blocks:
	% N  : the state covariance  (the purely 2nd order terms)
	% L  : the linear term
	% Phi : the constant term
	N = Q(1:end-1,1:end-1);
	L = Q(end,1:end-1);
	Phi = Q(end,end);

end
