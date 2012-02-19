classdef (ConstructOnLoad=false) plotcodfreeobj < handle
% This file defines the plotcodfreeobj class to construct plotting options and then generate plots 
% via the use of specific member functions thereof. This allows seperating the plotting from the 
% control synthesis. Due to the nature of how the plots are to be treated (i.e all passing of the object
% to a function is intented to be via reference to the same object), we derive this class from the  handle class. 
	properties(SetAccess = private)
		% set the private data members with default values where possible
		planeToPlotOn;
		toSaveAs;
		xAxisLabel;
		yAxisLabel;
		tickLabVal;
		tickLab;
		plotTitle;
		plotThreshold;
		numberOfGridPts = 10;
		azimuthSet =-41.5000;
    el=74;
		fontsize = 16;
		alphaVal = .9; % controls transparency
		numTicks = 5;
		indxOfStartingTimeStep = 0;
		numOfTimeStepsTotalTk  = 2; 
		costmat = [];
		indexplot =[];
		hjberrplot = [];
		quadraticindexplot = [];
		X = [];
		Y = [];
  end % of properties setting

	methods
		function obj = plotcodfreeobj(planeToPlotOn, toSaveAs, xAxisLabel,yAxisLabel,plotTitle,plotThreshold,numberOfGridPts)
			obj.planeToPlotOn = planeToPlotOn;
			obj.toSaveAs=toSaveAs;
			obj.xAxisLabel= xAxisLabel;
			obj.yAxisLabel = yAxisLabel;
			obj.plotTitle = plotTitle;
			obj.plotThreshold = plotThreshold;
			obj.numberOfGridPts = numberOfGridPts;
		end
		% end fo fn plotcodfreeobj (constructor)
		function plotAll(obj,codfreeObj)
		%% this function plots all the required plots given a cod free object. 
				mainPlotterFn(obj,codfreeObj);
				obj.plotofControlChoice();
				obj.plotofhjberror();
				obj.plotofcostmat();

		end
	
		function mainPlotterFn(obj,codfreeObj)

			d=linspace(-obj.plotThreshold,obj.plotThreshold,obj.numberOfGridPts);
			d=union(d,0);
			[obj.X,obj.Y] = meshgrid(d,d);  % A cube

			gridstepsizeforhjb = range(d)/obj.numberOfGridPts;
			obj.tickLab= linspace((min(d)),(max(d)),obj.numTicks);
			obj.tickLabVal = cellfun(@(x)num2str(x),num2cell(obj.tickLab),'UniformOutput',false);

			cMat = arrayfun(@(x,y)obj.findcost(x,y,codfreeObj),obj.X,obj.Y,'UniformOutput', false);


			indexplot = zeros(size(obj.X)); % control index used
			hjberrplot = zeros(size(obj.X)); % hjb error
			gradientsqerrorplot = zeros(size(obj.X));
			costmat = zeros(size(obj.X));
			quadraticindexplot = zeros(size(obj.X));
			hjberrNormalizedRatioplot = zeros(size(obj.X));

			for(k1=1:length(cMat))
					for(k2 = 1:length(cMat))
							costmat(k1,k2) = cMat{k1,k2}{1};
							indexplot(k1,k2) = cMat{k1,k2}{2}; % costu0
							hjberrplot(k1,k2) = cMat{k1,k2}{3};
							gradientneterrorplot(k1,k2) = cMat{k1,k2}{6};
							gradienterr(k1,k2) = cMat{k1,k2}{5};
							quadraticindexplot(k1,k2)= cMat{k1,k2}{7};
							hjberrNormalizedRatioplot(k1,k2)= cMat{k1,k2}{8};

					end    
			end
			obj.costmat = costmat;
			obj.indexplot = indexplot;
			obj.hjberrplot = hjberrplot;
			obj.quadraticindexplot = quadraticindexplot;

			save([obj.toSaveAs,'.mat'],'obj');

		end %end of mainPlotterFn


		function plotofControlChoice(obj)
			plotTitle = 'Plot of control choice';
			figure;surf(obj.X,obj.Y,obj.indexplot);
			title(plotTitle); view(obj.azimuthSet,obj.el);
			set(gca,'Ylim',[-obj.plotThreshold,obj.plotThreshold]);
			set(gca,'Xlim',[-obj.plotThreshold,obj.plotThreshold]);
			set(gca,'fontsize',obj.fontsize);
			xlabel(obj.xAxisLabel,'fontsize',obj.fontsize); ylabel(obj.yAxisLabel,'fontsize',obj.fontsize);
			set(gca,'ZtickLabel',[]);
			set(gca,'XTick',obj.tickLab);
			set(gca,'YTick',obj.tickLab);
			set(gca,'XtickLabel',obj.tickLabVal);
			set(gca,'YtickLabel',obj.tickLabVal);
			colorbar;
			saveas(gcf,[obj.toSaveAs,'plotofcontrolchoice.jpg']);
			close(gcf);
		end % end of plot of controlchoice 
 	
     
		function plotofquadraticchoice(obj)

			plotTitle = 'Plot of quadratic  chosen';
			figure;surf(obj.X,obj.Y,obj.quadraticindexplot);
			title(plotTitle); view(obj.azimuthSet,obj.el);
			set(gca,'Ylim',[-obj.plotThreshold,obj.plotThreshold]);
			set(gca,'Xlim',[-obj.plotThreshold,obj.plotThreshold]);
	    set(gca,'fontsize',16);
	   	alpha(obj.alphaVal);
			xlabel(obj.xAxisLabel,'fontsize',obj.fontsize); ylabel(obj.yAxisLabel,'fontsize',obj.fontsize);
			set(gca,'ZtickLabel',[]);
			set(gca,'XTick',obj.tickLab);
			set(gca,'YTick',obj.tickLab);
			set(gca,'XtickLabel',obj.tickLabVal);
			set(gca,'YtickLabel',obj.tickLabVal);
			colorbar;
			saveas(gcf,[obj.toSaveAs,'plotofquadraticchoice.jpg']);
			close(gcf);
		end % end of plot of quadratic choice  
 	
		function plotofhjberror(obj)

				plotTitle = 'Plot of hjb error';
				figure;surf(obj.X,obj.Y,obj.hjberrplot);
				title(plotTitle); view(obj.azimuthSet,obj.el);
				set(gca,'Ylim',[-obj.plotThreshold,obj.plotThreshold]);
				set(gca,'Xlim',[-obj.plotThreshold,obj.plotThreshold]);
				set(gca,'fontsize',16);
				alpha(obj.alphaVal);
				xlabel(obj.xAxisLabel,'fontsize',obj.fontsize); ylabel(obj.yAxisLabel,'fontsize',obj.fontsize);
				set(gca,'ZtickLabel',[]);
				set(gca,'XTick',obj.tickLab);
				set(gca,'YTick',obj.tickLab);
				set(gca,'XtickLabel',obj.tickLabVal);
				set(gca,'YtickLabel',obj.tickLabVal);
				colorbar;
				saveas(gcf,[obj.toSaveAs,'plotofhjberror.jpg']);
				close(gcf);
		end % end of plot of  hjb error
		

		function plotofcostmat(obj)
					plotTitle = 'Plot of cost funct ion';
					figure;surf(obj.X,obj.Y,obj.costmat);
					title(plotTitle); view(obj.azimuthSet,obj.el);
					set(gca,'Ylim',[-obj.plotThreshold,obj.plotThreshold]);
					set(gca,'Xlim',[-obj.plotThreshold,obj.plotThreshold]);
					set(gca,'fontsize',16);
					alpha(obj.alphaVal);
					xlabel(obj.xAxisLabel,'fontsize',obj.fontsize); ylabel(obj.yAxisLabel,'fontsize',obj.fontsize);
					set(gca,'ZtickLabel',[]);
					set(gca,'XTick',obj.tickLab);
					set(gca,'YTick',obj.tickLab);
					set(gca,'XtickLabel',obj.tickLabVal);
					set(gca,'YtickLabel',obj.tickLabVal);
					colorbar;
					saveas(gcf,[obj.toSaveAs,'costfnplot.jpg']);
					close(gcf);
		end % end of plot of cost function
			

					
		function resultvector = findcost(obj,x,y,codfreeObj)

						persistent indexusable_tk0 indexusable_tk0plus1 Ix Iy Iz IxI IIx IzI IIz IyI IIy IxIx IyIy IzIz;
					 	

							 
						if(isempty(indexusable_tk0))
											%% begin: initialize ix iy etc...
												Ix=sigmax;
												Iy=sigmay;
												Iz=sigmazimuthSet;
												I2=qo(eye(2));

												IxI= tensor(Ix,I2);
												IIx= tensor(I2,Ix);
												IIz=tensor(I2,Iz);
												IzI=tensor(Iz,I2);
												IxIx=tensor(Ix,Ix);

												IxIy=tensor(Ix,Iy);
												IxIz=tensor(Ix,Iz);

												IyIx=tensor(Iy,Ix);
												IyIz=tensor(Iy,Iz);
												IzIx=tensor(Iz,Ix);
												IzIy=tensor(Iz,Iy);  

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

												IxIy=full(IxIy(:,:));
												IxIz=full(IxIz(:,:));
												IyIx=full(IyIx(:,:));
												IyIz=full(IyIz(:,:));
												IzIx=full(IzIx(:,:));
												IzIy=full(IzIy(:,:));


												IzIz = full(IzIz(:,:));
												IyIy = full(IyIy(:,:));
												IxIx = full(IxIx(:,:));
												%% end initialize ix, iy etc

								%% now find out the index of usable controls for t0
								% the possible lengths of controls for obj.indxOfStartingTimeStep is found by searching for all the controls that have length atmost 
								% that corresponding to obj.numOfTimeStepsTotalTk-obj.indxOfStartingTimeStep
								indexusable_tk0 = find(cellfun(@(x)length(x),codfreeObj.controlindexBigList(1:codfreeObj.numAtPresent))...
													<=obj.numOfTimeStepsTotalTk-obj.indxOfStartingTimeStep,1,'last'); 
							
								%% now find the index of controls for obj.indxOfStartingTimeStep+1
								% the idea is the same as in the code above
								indexusable_tk0plus1 = find(cellfun(@(x)length(x),codfreeObj.controlindexBigList(1:codfreeObj.numAtPresent))...
												<=obj.numOfTimeStepsTotalTk-obj.indxOfStartingTimeStep-1,1,'last'); 
						end

						% now generate the unitary element from the x,y values given
						U0= expm(-1i*(obj.planeToPlotOn{1}*x+obj.planeToPlotOn{2}*y));
						% Now find V_{obj.indxOfStartingTimeStep}(u0) and V_{tk+1}(u0), V_{obj.indxOfStartingTimeStep}(utilde) 

						% first V_{obj.indxOfStartingTimeStep}(u0) 
						cost1= cellfun(@(M,N)(N*codfreeObj.delta + 2*codfreeObj.penalty*...
											abs(real(trace(eye(4)-M*U0)))),codfreeObj.BList(1:indexusable_tk0),...
											num2cell(codfreeObj.cValList(1:indexusable_tk0)));
						[costu0tk0,indxcontrolsigu0tk0]  = min(cost1);
						quadraticIndex_tk0 = indxcontrolsigu0tk0;
						indxcontrolsigu0tk0 = codfreeObj.controlindexBigList{indxcontrolsigu0tk0}(1);

						% next V_{obj.indxOfStartingTimeStep+1}(u0) 
						cost1= cellfun( @(M,N)(N*codfreeObj.delta + 2*codfreeObj.penalty*...
											abs(real(trace(eye(4)-M*U0)))),codfreeObj.BList(1:indexusable_tk0plus1),...
											num2cell(codfreeObj.cValList(1:indexusable_tk0plus1)));
						[costu0tk0plus1,indxcontrolsigu0tk0plus1]  = min(cost1);
						indxcontrolsigu0tk0plus1 = codfreeObj.controlindexBigList{indxcontrolsigu0tk0plus1}(1);

											
						pDV_pDt =  (costu0tk0plus1 - costu0tk0)/codfreeObj.delta;
						
						%now find V_{obj.indxOfStartingTimeStep}(exp()\cdot u(obj.indxOfStartingTimeStep)) 
						%i.e. V_{obj.indxOfStartingTimeStep}(utilde)
						if(indxcontrolsigu0tk0) %non zero optimal control
								 tildeU = codfreeObj.Am{indxcontrolsigu0tk0}*U0; 
						else
								 tildeU = U0;
						end
						cost1 = cellfun( @(M,N)(N*codfreeObj.delta + 2*codfreeObj.penalty*...
												abs(real(trace(eye(4)-M*tildeU)))),codfreeObj.BList(1:indexusable_tk0),...
												num2cell(codfreeObj.cValList(1:indexusable_tk0)));
						[costutildetk0,b]  = min(cost1);% this is V_{tk0}(utilde)
						% next V_{obj.indxOfStartingTimeStep+1}(utilde) 
						cost1= cellfun(@(M,N)(N*codfreeObj.delta + 2*codfreeObj.penalty* abs(real(trace(eye(4)-M*tildeU)))),...
											codfreeObj.BList(1:indexusable_tk0plus1),num2cell(codfreeObj.cValList(1:indexusable_tk0plus1)));
						[costutildetk0plus1,indxcontrolsigutildetk0plus1]  = min(cost1);
						indxcontrolsigutildetk0plus1=codfreeObj.controlindexBigList{indxcontrolsigutildetk0plus1}(1);
									
						DVDT = (costutildetk0plus1-costu0tk0 )/codfreeObj.delta;
									 
						pDV_pDx_times_Dx_Dt = (costutildetk0 - costu0tk0)/codfreeObj.delta;
									 
									 
						if(indxcontrolsigu0tk0) %if non zero optimal control
							 grad_err = (costutildetk0-costu0tk0)/codfreeObj.delta;
							 hjberrorNonSteadyState = codfreeObj.controlcostindexlist(indxcontrolsigu0tk0) + DVDT ;%
							 grad_sensitive = abs(abs(grad_err)- ( codfreeObj.controlcostindexlist(indxcontrolsigu0tk0)/norm(IxIx)));
							 ratioofhjbnorm = hjberrorNonSteadyState/( codfreeObj.controlcostindexlist(indxcontrolsigu0tk0)+ ...
																							abs(pDV_pDt)+ abs(pDV_pDx_times_Dx_Dt) );
								
						else % if the optimal control turns out to be 0 i.e. to do nothing
							 grad_err = 0;
							 grad_sensitive = 0;
							 hjberrorNonSteadyState=DVDT ; %pDV_pDt;
							 ratioofhjbnorm = 1;
						end %if non zero optimal control

					 
						resultvector = {costu0tk0,indxcontrolsigu0tk0,(hjberrorNonSteadyState),...
						norm(hjberrorNonSteadyState),grad_err, grad_sensitive,quadraticIndex_tk0,ratioofhjbnorm};
								 
		end	% end of fn findcost 
						
  
	
	end % of methods setting




end 
% end of classdef
