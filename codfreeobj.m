%% this defines the codfree class
% the from and to email id in the code must be filled in as required

classdef (ConstructOnLoad=false) codfreeobj


	properties(SetAccess = private)
		delta=0.1;
		numofsysqubits=2;
		numTsteps = 2;
		numToAllow = [11,20];
		numAtPresent;
		penalty=20; 
		Am;
		BList ;
		costOperator;
		controlindexBigList;
		controlcostindexlist;
		cValList;
		pathAndMatFileName;
		pathToFigFiles;
		pathToGmatrices;
		ThresholdPruning;

	%	costMatrixForControls	
	end

	
	
	methods
		%constructor
		function obj = codfreeobj(delta,numofsysqubits,numTsteps,numToAllow,penalty,Am,controlcostindexlist,pathAndMatFileName,pathToFigFiles,pathToGmatrices,ThresholdPruning)
		if(nargin < 11)
				error('codfree:constructor:insufficientArgumentsPassed');
			end
			
			obj.numofsysqubits = numofsysqubits;
			obj.numTsteps = numTsteps;
			obj.numToAllow = numToAllow;
			obj.penalty = penalty;
			obj.Am = Am; 
			obj.controlcostindexlist = controlcostindexlist;
			obj.numAtPresent = 1;
			obj.BList = cell(int32(sum(numToAllow)),1);
			obj.controlindexBigList= cell(int32(sum(numToAllow)),1);
			obj.cValList=zeros(length(obj.BList),1);
			obj.pathAndMatFileName = pathAndMatFileName;
			obj.pathToFigFiles= pathToFigFiles;
			obj.pathToGmatrices = pathToGmatrices;
			obj.ThresholdPruning = ThresholdPruning;
		  %now define costOperator	
			dimOfMatrices=2^(numofsysqubits);% 4 for 2 qubits
			if dimOfMatrices~=size(Am{1},1)
				error('codfree:constructor:sizeofcontrolmatrixDoesNotMatchThatForRequiredNumQubits');
			end
			Mtest=eye(dimOfMatrices);
			Mtest=Mtest(:)';
			obj.costOperator= (-2*penalty*[ Mtest zeros(1,dimOfMatrices*dimOfMatrices)]);
	
		end


		function obj = pruningProcess(obj)

				obj.BList{1}=eye(2^(obj.numofsysqubits));
				obj.controlindexBigList{1}=0; % the control which does nothing.
				obj.cValList(1)=0;
				%keyboard;
				tmpLenAm = length(obj.Am);
				obj.BList(2:tmpLenAm+1)= obj.Am;
				obj.controlindexBigList(2:tmpLenAm+1)= num2cell([1:tmpLenAm]');
				indexofAm = num2cell([1:tmpLenAm]');

				obj.cValList(2:tmpLenAm+1)= obj.controlcostindexlist;
				obj.numAtPresent= tmpLenAm+1;  % this is the num of controls (unpruned) in the BList cell.
				BAfterLastStep = obj.Am;
				controlindexListAfterLastStep = num2cell([1:tmpLenAm]');
				numberPruneableinEachIteration = zeros(length(obj.numToAllow),1);



			for(ctr1=2:length(obj.numToAllow))
					if(~mod(ctr1,5))
							 msg='five more done OOP_codfreecode';
							 sendemail(msg);	
					end 
					if(isempty(BAfterLastStep))
							break
					end
							
					[newBStepBeforePruning,newCindexListBeforePruning] = ...
											usingYalmipCreateBPossibleCell4(BAfterLastStep,controlindexListAfterLastStep,obj.Am,indexofAm)	;  %%!! verify that this function works!!
					costBStepInPruning = zeros(length(newBStepBeforePruning),1);
					newClambdabar = cellfun(@(x)sum(x),newCindexListBeforePruning);
				% now to obtain the costs for each of the new B elements
					for(ctr2=1:length(costBStepInPruning))
							%		new stuff from the use of maxplus and cvx:
							% subtract P_i - P_j to get elements of \bar{J} i.e P_\lambda
							fixedP= obj.costOperator*convertOperator(newBStepBeforePruning{ctr2});
							tmpvarForcfixed= newCindexListBeforePruning{ctr2};
							fixedc = sum(obj.controlcostindexlist(int32(tmpvarForcfixed)));
							PLambda = cellfun(@(x)((obj.costOperator)*convertOperator(x)-fixedP),[obj.BList(1:obj.numAtPresent);newBStepBeforePruning(setdiff(1:end,ctr2))],...
										'UniformOutput',false);
							cLambda = cellfun(@(x)(x-fixedc),num2cell([(obj.cValList(1:obj.numAtPresent));newClambdabar(setdiff(1:end,ctr2))]),'UniformOutput',false);
							t1 = pruneAction(PLambda,cLambda,obj.ThresholdPruning,obj.numofsysqubits,obj.pathToGmatrices,obj.delta);
				
							costBStepInPruning(ctr2) = t1; % take the least margin of contribution 
														
					end % of ctr2 
										
					% now to sort the costs of the new B elements to prune them
					[a1 indx1] = sort(costBStepInPruning,'descend');
					[indx2]=find(a1<=0); % find the pruneable negative values
					
					numberPruneableinEachIteration(ctr1) = length(indx2);
					%% now pick the ones to retain in the list of controls
					
					% pick the index of the top NumToKeep of them
					% choose those elements of BPossible and return those
				
					if(isempty(min(indx2))) % i.e if there are no pruneable parameters
							obj.BList(obj.numAtPresent+1:obj.numAtPresent+obj.numToAllow(ctr1)) = newBStepBeforePruning(indx1(1:int32(obj.numToAllow(ctr1)))); 
											%closer to the identity, the more its worth
							obj.controlindexBigList(obj.numAtPresent+1:obj.numAtPresent+obj.numToAllow(ctr1))=...
																										newCindexListBeforePruning(int32(indx1(1:int32(obj.numToAllow(ctr1))))) ;
							tt11tvmp= cell2mat(newCindexListBeforePruning(indx1(1:int32(obj.numToAllow(ctr1)))));
									obj.cValList(obj.numAtPresent+1:obj.numAtPresent+obj.numToAllow(ctr1)) = sum(obj.controlcostindexlist(int32(tt11tvmp)),2);
							obj.numAtPresent = obj.numAtPresent+obj.numToAllow(ctr1);
							BAfterLastStep = newBStepBeforePruning(indx1(1:obj.numToAllow(ctr1)));
							controlindexListAfterLastStep= newCindexListBeforePruning(indx1(1:obj.numToAllow(ctr1)));
				   else 
							% some of them are pruneable and indx2 has the locations of the pruneable ones
							tm1 =  min(indx2)-1;
							obj.BList(obj.numAtPresent+1:obj.numAtPresent+min([obj.numToAllow(ctr1),tm1])) = newBStepBeforePruning(indx1(1:int32(min([obj.numToAllow(ctr1),tm1]))));
							obj.controlindexBigList(obj.numAtPresent+1:obj.numAtPresent+min([obj.numToAllow(ctr1),tm1]))=...
																															newCindexListBeforePruning(indx1(1:int32(min([obj.numToAllow(ctr1),tm1]))));
							tt11tvmp= cell2mat(newCindexListBeforePruning(indx1(1:int32(min([obj.numToAllow(ctr1),tm1])))));
							obj.cValList(obj.numAtPresent+1:obj.numAtPresent+int32(min([obj.numToAllow(ctr1),tm1]))) = sum(obj.controlcostindexlist(int32(tt11tvmp)),2);
							obj.numAtPresent = obj.numAtPresent+int32(min([obj.numToAllow(ctr1),tm1]));
							BAfterLastStep = newBStepBeforePruning(indx1(1:int32(min([obj.numToAllow(ctr1),tm1]))));
							controlindexListAfterLastStep= newCindexListBeforePruning(indx1(1:int32(min([obj.numToAllow(ctr1),tm1]))));
							if(isequal(int32(tm1),0))
								sendemail('no controls left after pruning');
								keyboard;
								break;
							end
					
					end											% of if there are no pruneable params

				save([obj.pathAndMatFileName,'.mat'],'obj');
			end % of ctr 1	
			

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
 

end

function sendemail(topic,attachmentPath,from,to)
global  SMTPname FromEmailID ToEmailID; 
	if nargin > 4
			error('codfree:sendemail:TooManyInputs', ...
					'requires at most 4 inputs');
	end

	setpref('Internet','SMTP_Server',SMTPname);
	switch nargin
			case 1
					from= FromEmailID ;
					to = ToEmailID;
					setpref('Internet','E_mail',from);
					sendmail(to,topic);
				case 2
					from= ;
					to =;
					setpref('Internet','E_mail',from);
					sendmail(to,topic);
				case 3
					error('codfree:sendemail:NeedTocreateFromToOptions:tbd');

		end
				

end



function [ret,retc]=usingYalmipCreateBPossibleCell4(BAfterLastStep,controlindexListAfterLastStep,Am,indexofAm)
% Am is the control signal cell array
% i.e BAfterLastStep is larger so better to use that as the fixed matrix
ret1=cellfun(@(x)multbyAm(x,BAfterLastStep),Am,'UniformOutput',false);
ret2=([ret1{:}]);
ret=ret2(:);
ret1=cellfun(@(x)appendtoindexofAm(x,controlindexListAfterLastStep),indexofAm,'UniformOutput',false);
ret2=([ret1{:}]);
retc=ret2(:);
end


function re=multbyAm(x,B)
re =(mat2cell(x*[B{:}],[length(x)],length(x)*ones(1,length(B))))';
end
function re=appendtoindexofAm(x,iB)
re =(mat2cell([cell2mat(iB),ones(length(iB),1)*x],ones(1,length(iB)),[length(iB{end})+1]))';
end





function ret = pruneAction(PLambda,cLambda,objThresholdPruning,numofsysqubits,PathToGmatrices,delta)
    persistent vvar1;
    persistent speederpruning; 
    persistent G1cell G2cell;
    persistent A0 B0 c0 M0 cPruning MPruning bPruning ThresholdPruning;
	  newConstraints = [];

		if isempty(speederpruning)
			% This uses the persistent var to avoid loading the file repeatedly
				vvar1= 2^(2*numofsysqubits+1);
				speederpruning = 0;
				load(PathToGmatrices);
				A0 = 0*speye(1+vvar1); % takes into account  y =(x,\xi) 
				B0 = [0*speye(vvar1,1);1];
				c0 = 0;
				M0 = [c0 B0'; B0 A0];
				sizeofm= 2^(numofsysqubits);
				cPruning = 2*sizeofm;
				bPruning = -[reshape(eye(sizeofm),[sizeofm^2,1]); zeros(sizeofm^2,1);0]; 
				ThresholdPruning=  objThresholdPruning; % for 3 steps .. even though only 2 are used.
		%2.4263; % for 8 steps 

		%1.3973;%  3.6776; %for 6 IyIy.
		% 3.0271;% for 9 IyIy;  %7.4341;% 8-2*real(trace(expm(-1i*.1*5*IyIy)))
				MPruning = [cPruning-ThresholdPruning,bPruning';bPruning A0];
		end % end of isempty speederpruning
   
	%begin of cvx code    

	cvx_begin
	cvx_solver sdpt3
	cvx_precision high
	cvx_quiet(true)
	variable Y(2+vvar1,2+vvar1) symmetric

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
			
					A = 0*speye(1+vvar1);
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
