%% this defines the codfree class

classdef (ConstructOnLoad=false) codfreeobj


	properties(SetAccess = private)
		% private datamembers with default values where needed.
		delta=0.1;
		numofsysqubits=2;
		numTsteps = 2;
		numToAllow = [11,20]; %number of quadratics to allow at  each stage
		numAtPresent; % keeps a running count of the number of quadratics 
		penalty=20; 
		cellofAvailableControlMatrices;
		bigListOfControlSeq ;
		costOperator;
		bigListofControlIndexes;
		controlcostindexlist;
		cValList;
		pathAndMatFileName;
		pathToFigFiles;
		pathToGmatrices;
		pruningZoneThreshold;

	end

	
	
	methods
		%constructor
		function obj = codfreeobj(delta,numofsysqubits,numTsteps,numToAllow,penalty,...
			cellofAvailableControlMatrices,controlcostindexlist,pathAndMatFileName,...
			pathToFigFiles,pathToGmatrices,pruningZoneThreshold)
			if(nargin < 11)
					error('codfree:constructor:insufficientArgumentsPassed');
			end
				
			obj.numofsysqubits = numofsysqubits;
			obj.numTsteps = numTsteps;
			obj.numToAllow = numToAllow;
			obj.penalty = penalty;
			obj.cellofAvailableControlMatrices = cellofAvailableControlMatrices; 
			obj.controlcostindexlist = controlcostindexlist;
			obj.numAtPresent = 1;
			obj.bigListOfControlSeq = cell(int32(sum(numToAllow)),1);
			obj.bigListofControlIndexes= cell(int32(sum(numToAllow)),1);
			obj.cValList=zeros(length(obj.bigListOfControlSeq),1);
			obj.pathAndMatFileName = pathAndMatFileName;
			obj.pathToFigFiles= pathToFigFiles;
			obj.pathToGmatrices = pathToGmatrices;
			obj.pruningZoneThreshold = pruningZoneThreshold;
			%now define the costOperator	
				dimOfMatrices=2^(numofsysqubits);% 4 for 2 qubits
				if dimOfMatrices~=size(cellofAvailableControlMatrices{1},1)
					error('codfree:constructor:sizeofcontrolmatrixDoesNotMatchThatForRequiredNumQubits');
				end
				Mtest=eye(dimOfMatrices);
				Mtest=Mtest(:)';
			obj.costOperator= (-2*penalty*[ Mtest zeros(1,dimOfMatrices*dimOfMatrices)]);
	
		end % of constructor function


		function obj = pruningProcess(obj)
				obj.bigListOfControlSeq{1}=eye(2^(obj.numofsysqubits));
				obj.bigListofControlIndexes{1}=0; % the control which does nothing.
				obj.cValList(1)=0;
				tmpLencellofAvailableControlMatrices = length(obj.cellofAvailableControlMatrices);
				obj.bigListOfControlSeq(2:tmpLencellofAvailableControlMatrices+1)= obj.cellofAvailableControlMatrices;
				obj.bigListofControlIndexes(2:tmpLencellofAvailableControlMatrices+1)=...
									num2cell([1:tmpLencellofAvailableControlMatrices]');
				cellOfIndicesOfAvailableControls = num2cell([1:tmpLencellofAvailableControlMatrices]');

				obj.cValList(2:tmpLencellofAvailableControlMatrices+1)= obj.controlcostindexlist;
			
				obj.numAtPresent= tmpLencellofAvailableControlMatrices+1; 
				% the above is the num of controls (unpruned) in the bigListOfControlSeq cell.
			
				controlSeqAfterLastStep = obj.cellofAvailableControlMatrices;
				controlindexListAfterLastStep = num2cell([1:tmpLencellofAvailableControlMatrices]');
				numberPruneableinEachIteration = zeros(length(obj.numToAllow),1);



				for(counterCurrentTimeStep=2:length(obj.numToAllow))
						if(~mod(counterCurrentTimeStep,5))
								 msg='five more done OOP_codfreecode';
								 sendemail(msg);	
						end 
						if(isempty(controlSeqAfterLastStep))
								break
						end
								
						[newBStepBeforePruning,newCindexListBeforePruning] = ...
												createCellOfNewPotentialControls(controlSeqAfterLastStep,...
												controlindexListAfterLastStep,obj.cellofAvailableControlMatrices,...
												cellOfIndicesOfAvailableControls)	;  
						costBStepInPruning = zeros(length(newBStepBeforePruning),1);
						newClambdabar = cellfun(@(x)sum(x),newCindexListBeforePruning);

						% now to obtain the costs for each of the new B elements
						for(counterCurrentQuadratic=1:length(costBStepInPruning))
								% subtract P_i - P_j to get elements of \bar{J} i.e P_\lambda
								% the notation P, lambda etc correspond to those in the conference/journal papers
								
								fixedP= obj.costOperator*convertOperator(newBStepBeforePruning{counterCurrentQuadratic});
								tmpvarForcfixed= newCindexListBeforePruning{counterCurrentQuadratic};
								fixedc = sum(obj.controlcostindexlist(int32(tmpvarForcfixed)));
								PLambda = cellfun(@(x)((obj.costOperator)*convertOperator(x)-fixedP),...
													[obj.bigListOfControlSeq(1:obj.numAtPresent);newBStepBeforePruning(...
													setdiff(1:end,counterCurrentQuadratic))],'UniformOutput',false);
								cLambda = cellfun(@(x)(x-fixedc),num2cell([(obj.cValList(1:obj.numAtPresent));...
													newClambdabar(setdiff(1:end,counterCurrentQuadratic))]),'UniformOutput',false);
								t1 = pruneAction(PLambda,cLambda,obj.pruningZoneThreshold,...
													obj.numofsysqubits,obj.pathToGmatrices,obj.delta);
					
								costBStepInPruning(counterCurrentQuadratic) = t1; % take the least margin of contribution 
															
						end % of counterCurrentQuadratic 
											
						% now to sort the costs of the new B elements to prune them
						[pruneabilityMetric indx1] = sort(costBStepInPruning,'descend');
						[indexListOfPruneableSeq]=find(pruneabilityMetric<=0); % find the pruneable negative values
						
						numberPruneableinEachIteration(counterCurrentTimeStep) = length(indexListOfPruneableSeq);
						%% now pick the ones to retain in the list of controls
						
						% pick the index of the top NumToKeep of them
						% choose those elements of BPossible and return those
					
						if(isempty(min(indexListOfPruneableSeq))) % i.e if there are no pruneable parameters
								obj.bigListOfControlSeq(obj.numAtPresent+1:obj.numAtPresent+obj.numToAllow(counterCurrentTimeStep)) = ...
												newBStepBeforePruning(indx1(1:int32(obj.numToAllow(counterCurrentTimeStep)))); 
								obj.bigListofControlIndexes(obj.numAtPresent+1:obj.numAtPresent+obj.numToAllow(counterCurrentTimeStep))=...
												newCindexListBeforePruning(int32(indx1(1:int32(obj.numToAllow(counterCurrentTimeStep))))) ;
								tt11tvmp= cell2mat(newCindexListBeforePruning(indx1(1:int32(obj.numToAllow(counterCurrentTimeStep)))));
								obj.cValList(obj.numAtPresent+1:obj.numAtPresent+obj.numToAllow(counterCurrentTimeStep)) =...
												sum(obj.controlcostindexlist(int32(tt11tvmp)),2);
								obj.numAtPresent = obj.numAtPresent+obj.numToAllow(counterCurrentTimeStep);
								controlSeqAfterLastStep = newBStepBeforePruning(indx1(1:obj.numToAllow(counterCurrentTimeStep)));
								controlindexListAfterLastStep= newCindexListBeforePruning(indx1(1:obj.numToAllow(counterCurrentTimeStep)));
						 else 
								% some of them are pruneable and indexListOfPruneableSeq has the locations of the pruneable ones
								tm1 =  min(indexListOfPruneableSeq)-1;
								obj.bigListOfControlSeq(obj.numAtPresent+1:obj.numAtPresent+min([obj.numToAllow(...
									 counterCurrentTimeStep),tm1])) = ...
													newBStepBeforePruning(indx1(1:int32(min([obj.numToAllow(counterCurrentTimeStep),tm1]))));
								obj.bigListofControlIndexes(obj.numAtPresent+1:obj.numAtPresent+...
																			min([obj.numToAllow(counterCurrentTimeStep),tm1]))=...
													newCindexListBeforePruning(indx1(1:int32(min([obj.numToAllow(counterCurrentTimeStep),tm1]))));
								tt11tvmp= cell2mat(newCindexListBeforePruning(indx1(1:int32(...
													min([obj.numToAllow(counterCurrentTimeStep),tm1])))));
								obj.cValList(obj.numAtPresent+1:obj.numAtPresent+int32(min([obj.numToAllow(...
															counterCurrentTimeStep),tm1]))) =sum(obj.controlcostindexlist(int32(tt11tvmp)),2);
								obj.numAtPresent = obj.numAtPresent+int32(min([obj.numToAllow(counterCurrentTimeStep),tm1]));
								controlSeqAfterLastStep = newBStepBeforePruning(indx1(1:int32(min(...
															[obj.numToAllow(counterCurrentTimeStep),tm1]))));
								controlindexListAfterLastStep= newCindexListBeforePruning(...
													indx1(1:int32(min([obj.numToAllow(counterCurrentTimeStep),tm1]))));

								if(isequal(int32(tm1),0))
									sendemail('no controls left after pruning');
									keyboard;
									break;
								end
						
						end % of if there are no pruneable params

					save([obj.pathAndMatFileName,'.mat'],'obj');
				end % of for loop with counterCurrentTimeStep
			
		end % end of fn pruningProcess i.e  the main file for pruning





	end
	%end of methods in class



end
% end of classDefinition


function ret = convertOperator(B)
% this converts any operation from SU(2^2) space to a 32 x 32 matrix which can operate on the 
% 32 element vector representation of SU(2^2)

		R= real(B);
		S= imag(B);

		ret = [blkdiag(R,R,R,R), blkdiag(-S,-S,-S,-S);blkdiag(S,S,S,S),blkdiag(R,R,R,R)]; 
		 
end % of function convertOperator

function sendemail(topic,attachmentPath,from,to)
			global  SMTPname FromEmailID ToEmailID; 
			if nargin > 4
					error('codfree:sendemail:TooManyInputs', ...
							'requires at most 4 inputs');
			end

			setpref('Internet','SMTP_Server',SMTPname);
			switch nargin %  what to do depending on number of inputs to the mail routine 
					case 1
							from= FromEmailID ;
							to = ToEmailID;
							setpref('Internet','E_mail',from);
							sendmail(to,topic);
					case 2
							from= FromEmailID ;
							to = ToEmailID; 
							setpref('Internet','E_mail',from);
							sendmail(to,topic);
					case 3
							error('codfree:sendemail:NeedTocreateFromToOptions:tbd');

			end %  what to do depending on number of inputs to the mail routine 

end % of function sendemail



function [ret,retc]=createCellOfNewPotentialControls(controlSeqAfterLastStep,controlindexListAfterLastStep,cellofAvailableControlMatrices,cellOfIndicesOfAvailableControls)
% cellofAvailableControlMatrices is the control signal cell array
% i.e controlSeqAfterLastStep is larger so better to use that as the fixed matrix
		ret1=cellfun(@(x)multbycellofAvailableControlMatrices(x,controlSeqAfterLastStep),...
								cellofAvailableControlMatrices,'UniformOutput',false);
		ret2=([ret1{:}]);
		ret=ret2(:);
		ret1=cellfun(@(x)appendtocellOfIndicesOfAvailableControls(x,...
								controlindexListAfterLastStep),cellOfIndicesOfAvailableControls,'UniformOutput',false);
		ret2=([ret1{:}]);
		retc=ret2(:);
end


function re=multbycellofAvailableControlMatrices(x,B)
	re =(mat2cell(x*[B{:}],[length(x)],length(x)*ones(1,length(B))))';
end
function re=appendtocellOfIndicesOfAvailableControls(x,iB)
	re =(mat2cell([cell2mat(iB),ones(length(iB),1)*x],ones(1,length(iB)),[length(iB{end})+1]))';
end





function ret = pruneAction(PLambda,cLambda,objpruningZoneThreshold,numofsysqubits,PathToGmatrices,delta)
		%% this function prunes the max-plus basis directions (i.e. a max-plus Projection) to reduce
		%the growth due to the control set

		% For consistency the notation in this function c0,A0,B0,M, G etc are similar to those used in the paper 

		% In this function we also use persistent variables to speed up repeated calls to this function (especially through
		% avoiding the loading of the constraint set each time)
	
		persistent var_NumberOfDimensions;
    persistent var_fasterPruning; 
    persistent G1cell G2cell;
    persistent A0 B0 c0 M0 cPruning MPruning bPruning pruningZoneThreshold;
	  newConstraints = [];

		if isempty(var_fasterPruning)
			% This uses the persistent var to avoid loading the file repeatedly
				var_NumberOfDimensions= 2^(2*numofsysqubits+1);
				var_fasterPruning = 0;

				% load the matrices which encode the special unitary group constraints
				load(PathToGmatrices);
				A0 = 0*speye(1+var_NumberOfDimensions); % takes into account  y =(x,\xi) 
				B0 = [0*speye(var_NumberOfDimensions,1);1];
				c0 = 0;
				M0 = [c0 B0'; B0 A0];
				sizeofm= 2^(numofsysqubits);
				cPruning = 2*sizeofm;
				bPruning = -[reshape(eye(sizeofm),[sizeofm^2,1]); zeros(sizeofm^2,1);0]; 
				pruningZoneThreshold=  objpruningZoneThreshold; % for 3 steps .. even though only 2 are used.
				%2.4263; % for 8 steps 
				%1.3973;%  3.6776; %for 6 IyIy.
				% 3.0271;% for 9 IyIy;  %7.4341;% 8-2*real(trace(expm(-1i*.1*5*IyIy)))
				MPruning = [cPruning-pruningZoneThreshold,bPruning';bPruning A0];
		end % end of isempty var_fasterPruning
	%begin of cvx code    

	cvx_begin
	cvx_solver sdpt3
	cvx_precision high
	cvx_quiet(true)
	variable Y(2+var_NumberOfDimensions,2+var_NumberOfDimensions) symmetric

	maximize(trace(M0*Y))
	subject to
		Y(1,1)==1
		lambda_min(Y)>=0
		%Y in positivesemidefinite
	for(constraintCtr = 1:length(G2cell))
			A = blkdiag(G1cell{constraintCtr},0);
			M = blkdiag(G2cell{constraintCtr},A);
			trace(M*Y)==0;
	end
	for(ctrprune = 1: length(cLambda))
			A = 0*speye(1+var_NumberOfDimensions);
			b = [PLambda{ctrprune}'; -1];
			c = delta*cLambda{ctrprune};
			M = [c,b';b A];
			trace(M*Y)>=0;
	end
				

	%% Now we constrain the region over which we perform the comparison for pruning

	% Tr[ (U-I) * \dag{(U-I)} ]  <= Threshold
	% this leads to   cPruning = 8,  
	% bPruning = -[reshape(eye(4),[16,1]); zeros(16,1]

	trace(MPruning*Y)<=0;
		
	% 	trace(MPruning*Y)>=0;   
			
	cvx_end
	%disp('test');
	% keyboard;
	% end of cvx code
	ret = trace(M0*double(Y));

end
