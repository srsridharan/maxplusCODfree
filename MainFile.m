function MainFile
%% MainFile  performs the control design using the curse of dimensionality free appraoch on the SU4 group
%% and then develops the the performance plots

			global SMTPname FromEmailID ToEmailID;
			% these global variables are to set the details of the email routines to send
			% reports/results/figures  once the simulation is complete
			SMTPname = ;
			FromEmailID = ;
			ToEmailID =  ;
			try 
							% initialize the parameters for the control design and run the control design
							% and performance plot studies
							ctmain = 3;
							% main program for SU4 tests
							% IxI IIx IIz IzI IxIx
							global ctmainVariable SMTPname FromEmailID ToEmailID; 

							%% now set the parameters required for the simulation
							% all directions except the directions that use both qubits are weighed by (1/ctmain)
							costMatrixForControls= [ones(1,8)/ctmain, 1,1];

							delta = 0.2;
							numofsysqubits = 2;
							numTsteps = 2;
							numToAllow = [11,20];
							penalty= 20;


							try % run the control design and plot the performance
								ctmainVariable = ctmain;
								pathtoBaseDir= pwd;

								% create  file name including the cost details and number of steps in the name
								pathtoMatFiletoSave = [pwd,'/matfiles/sdplessprunedfullcontrol',...
														num2str(ctmain),'cost',num2str(ctmain),num2str(numTsteps),'steps'];
								
								% now run the main object creation and dimensionality reduction approach on the 
								% SU4 unitary group
								mainFileForControlGenSU4(costMatrixForControls,pathtoBaseDir,...
														fileName,delta,numofsysqubits,numTsteps,numToAllow,penalty);
								% now plot the results of the control design
								mainFileForControlPerformancePlotSU4(costMatrixForControls,pathtoBaseDir,...
														pathtoMatFiletoSave,delta,numofsysqubits,numTsteps,numToAllow,penalty);
							catch % run the control design and plot the performance
								errs=lasterror;
								setpref('Internet','SMTP_Server',SMTPname);
								setpref('Internet','E_mail',FromEmailID); 
								sendmail(ToEmailID,'error while running the control design and performance plot');
								keyboard;

							end  % run the control design and plot the performance

							
							try % try saving and emailing the results
								zip('./matfiles/backup',{'./matfiles/sdplessprunedfullcontrol*.mat'});
								zip('./matfiles/backup',{'./figfiles/*.jpg'});
								setpref('Internet','SMTP_Server',SMTPname);
								setpref('Internet','E_mail',FromEmailID);
								sendmail(ToEmailID,...
												 'Test subject','Test message',...
												 {['./matfiles/backup.zip']});
							catch % try saving and emailing the results
								setpref('Internet','SMTP_Server',SMTPname);
								setpref('Internet','E_mail',FromEmailID); 
								sendmail(ToEmailID,...
												 'error in sending the zipped result files');
							end % try saving and emailing the results

			catch % initialize the parameters for the control design and run the control design and performance plot studies
					errs=lasterror;
					setpref('Internet','SMTP_Server',SMTPname);
					setpref('Internet','E_mail',FromEmailID);
					sendmail(ToEmailID,...
								 'err_macmachinemainfile',errs.message);
			end % initialize the parameters for the control design and run the control design and performance plot studies

end % end of function MainFile
