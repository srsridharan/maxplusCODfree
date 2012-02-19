function mainFileForControlGenSU4(costMatrixForControls,pathtoBaseDir,fileName,delta,numofsysqubits,numTsteps,numToAllow,penalty)
%% mainFileForControlGenSU4(costMatrixForControls,pathtoBaseDir,fileName,delta,...
%									numofsysqubits,numTsteps,numToAllow,penalty)
% This function is responsible for synthesizing the control sequences, pruning them using the max-plus ideas,
% and storing them.


		cellofAvailableControlMatrices = CreateControlSignalCellArray4(delta);
		controlcostindexlist =costMatrixForControls;
		pathToSaveMatFiles = fileName ;
		pathToFigFiles = [pathtoBaseDir,'/figfiles/'];

		% now set the path to the constraints that are imposed by the Lie group i.e. special
		% unitary group structure.(This is used in the pruning step)
		pathToGmatrices = [pathtoBaseDir,'/Gmatrices.mat'];
		
		% The following controls the zone around the identity element  that the pruning procedure 
		% studies the pruneability of various control sequences (This is used in the pruning step)
		pruningZoneThreshold =0.3573; 
		 



		Obj1=codfreeobj(delta,numofsysqubits,numTsteps,numToAllow,penalty,cellofAvailableControlMatrices,...
						controlcostindexlist,pathToSaveMatFiles,pathToFigFiles,pathToGmatrices,pruningZoneThreshold);
		Obj1.pruningProcess();

end







function cellofAvailableControlMatrices= CreateControlSignalCellArray4(delta)
	%% the following code creates the control signal cell array cellofAvailableControlMatrices ...
	%	(i.e B^0_k's i.e the possible controls)
	% the varargin may contain costs for the various control signals. The costs are modeled as various magnitudes on the 
	% controls a cost of 5 means that the amplitude is 1/5 as large.
	
	% now generate the fundamental Pauli operators and tensor products thereof. These are chosen such that
	%	they form	the control signals available. 
	Ix=sigmax;
	Iy=sigmay;
	Iz=sigmaz;
	I2=qo(eye(2));

	IxI= tensor(Ix,I2);
	IIx= tensor(I2,Ix);
	IIz=tensor(I2,Iz);
	IzI=tensor(Iz,I2);
	IxIx=tensor(Ix,Ix);


	cellofAvailableControlDirections={full(IxI(:,:)), full(-IxI(:,:)),...
										full(IzI(:,:)),full(-IzI(:,:)),full(IIx(:,:)),...
										full(-IIx(:,:)), full(IIz(:,:)),full(-IIz(:,:)),full( IxIx(:,:)),full(-IxIx(:,:)) };
	numControlsDim=length(cellofAvailableControlDirections);
	y=eye(numControlsDim); 
	A1=mat2cell(y,ones(length(y),1),[numControlsDim]);

	cellofAvailableControlMatricesat=cell(1,numControlsDim); % contains the...
										% matrices of the control directions of the Lie alg
	
	cellofAvailableControlMatrices=cellfun(@(x)makecontrol(x,cellofAvailableControlDirections,...
										numControlsDim,delta),A1,'UniformOutput',false); % a cell array of control matrices


end % of function cellofAvailableControlMatrices

function ret = makecontrol(x,cellofAvailableControlDirections,numControlsDim,delta)
%% makecontrol is a helper function that generates the control unitary (state transiion matrix)
%	from the control matrices
		atmp=zeros(length(cellofAvailableControlDirections{1}));
		for(ctr1=1:numControlsDim)
				atmp=atmp+cellofAvailableControlDirections{ctr1}*x(ctr1);
		end

		ret = expm(-i*delta*atmp);
end
 
    
    
    

