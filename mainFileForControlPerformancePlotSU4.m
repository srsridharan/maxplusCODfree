function mainFileForControlPerformancePlotSU4(costMatrixForControls,pathtoBaseDir,...
						fileName,delta,numofsysqubits,numTsteps,numToAllow,penalty)
%% mainFileForControlPerformancePlotSU4(costMatrixForControls,pathtoBaseDir,...
%						fileName,delta,numofsysqubits,numTsteps,numToAllow,penalty)
% This function is responsible for performing the error analysis for the curse of dimensionality 
% free approach. The value function is observed over a 2d plane (in the much larger space) using the
% controls (quadratics) chosen via pruning. Various metrics for the error are plotted and the results are emailed
% back to the user as requested



%% Now for the plots:
					Iz=sigmaz;
					I2=qo(eye(2));
					Ix = sigmax;
					Iy = sigmay;
					IxI= tensor(Ix,I2);
					IIx= tensor(I2,Ix);

					IIz=tensor(I2,Iz);
					IzI=tensor(Iz,I2);
					IxIx=tensor(Ix,Ix);
					


					IyI=tensor(Iy,I2);
					IIy=tensor(I2,Iy);
					IyIy=tensor(Iy,Iy);
					IzIz=tensor(Iz,Iz);


					Ix = full(Ix(:,:));
					Iy = full(Iy(:,:));
					Iz = full(Iz(:,:));
					IxI = full(IxI(:,:));
					IIx = full(IIx(:,:));
					IIy = full(IIy(:,:));
					IyI = full(IyI(:,:));
					IIz = full(IIz(:,:));
					IzI = full(IzI(:,:));
					IzIz = full(IzIz(:,:));
					IyIy = full(IyIy(:,:));
					IxIx = full(IxIx(:,:));


	pathToSaveMatFiles = fileName ;%[pathtoBaseDir,'/matfiles/',fileName];
	pathToFigFiles = [pathtoBaseDir,'/figfiles/'];
	
	Obj2 = load(pathToSaveMatFiles);
	Obj1 = Obj2.obj;	

	% now initialize the plane over which the plot is to be generated
	APlotPlane={IxI,IzI};

	% mention the filename desired for saving the results and the plot titles etc
	toSaveAs=[pathToFigFiles,'IxIvsIzI'];
	xLab='I_x\otimes I';
	yLab ='I_z\otimes I';
	plotTitle='Plot of I_x\otimes I vs I_z\otimes I';

	% set the plot axes limits and the number of grid points
	plotThreshold = 0.2;
	numberOfGridPts = 10;

	% create a plotcodfreeobj object and invoke the plot member function
	cplotobj = plotcodfreeobj(APlotPlane, toSaveAs, xLab,yLab,plotTitle,plotThreshold,numberOfGridPts);
	cplotobj.plotAll(Obj1);




end
% end of mainFileForControlPerformancePlotSU4
