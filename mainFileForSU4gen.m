function mainFileForSU4gen(costMatrixForControls,pathtoBaseDir,fileName,delta,numofsysqubits,numTsteps,numToAllow,penalty)

Am = CreateControlSignalCellArray4(delta);
controlcostindexlist =costMatrixForControls;
pathToSaveMatFiles = fileName ;
pathToFigFiles = [pathtoBaseDir,'/figfiles/'];
pathToGmatrices = [pathtoBaseDir,'/Gmatrices.mat'];
ThresholdPruning =0.3573; 
 



Obj1=codfreeobj(delta,numofsysqubits,numTsteps,numToAllow,penalty,Am,controlcostindexlist,pathToSaveMatFiles,pathToFigFiles,pathToGmatrices,ThresholdPruning);
Obj1.pruningProcess();

end







function Am= CreateControlSignalCellArray4(delta)
	%% the following code creates the control signal cell array Am (i.e B^0_k's i.e the possible controls)
	% the varargin may contain costs for the various control signals. The costs are modeled as various magnitudes on the 
	% controls a cost of 5 means that the amplitude is 1/5 as large.

	Ix=sigmax;
	Iy=sigmay;
	Iz=sigmaz;
	I2=qo(eye(2));

	IxI= tensor(Ix,I2);
	IIx= tensor(I2,Ix);
	IIz=tensor(I2,Iz);
	IzI=tensor(Iz,I2);
	IxIx=tensor(Ix,Ix);


	Amatrix={full(IxI(:,:)), full(-IxI(:,:)), full(IzI(:,:)),full(-IzI(:,:)),full(IIx(:,:)),...
										full(-IIx(:,:)), full(IIz(:,:)),full(-IIz(:,:)),full( IxIx(:,:)),full(-IxIx(:,:)) };
	numControlsDim=length(Amatrix);
	y=eye(numControlsDim); 
	A1=mat2cell(y,ones(length(y),1),[numControlsDim]);
	Amat=cell(1,numControlsDim); % contains the matrices of the control directions of the Lie alg
	Am=cellfun(@(x)makecontrol(x,Amatrix,numControlsDim,delta),A1,'UniformOutput',false); % a 243*1 cell array of control matrices


end

function ret = makecontrol(x,Amatrix,numControlsDim,delta)
		atmp=zeros(length(Amatrix{1}));
		for(ctr1=1:numControlsDim)
				atmp=atmp+Amatrix{ctr1}*x(ctr1);
		end

		ret = expm(-i*delta*atmp);
end
 
    
    
    

