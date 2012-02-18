function usingYalmipMainSu4Iterations(ctmain)
% main program for SU4 tests
% IxI IIx IIz IzI IxIx
global ctmainVariable SMTPname FromEmailID ToEmailID; 

costMatrixForControls= [ones(1,8)/ctmain, 1,1];
delta = 0.2;
numofsysqubits = 2;
numTsteps = 2;
numToAllow = [11,20];
penalty= 20;


 	try

        ctmainVariable = ctmain;
				pathtoBaseDir= pwd;
				pathtoMatFiletoSave = [pwd,'/matfiles/sdplessprunedfullcontrol', num2str(ctmain),'cost',num2str(ctmain),num2str(numTsteps),'steps'];
			mainFileForSU4gen(costMatrixForControls,pathtoBaseDir,fileName,delta,numofsysqubits,numTsteps,numToAllow,penalty);
  mainFileForSU4plot(costMatrixForControls,pathtoBaseDir,pathtoMatFiletoSave,delta,numofsysqubits,numTsteps,numToAllow,penalty);




 	catch

 		errs=lasterror;
		setpref('Internet','SMTP_Server',SMTPname);
		setpref('Internet','E_mail',FromEmailID); 
		 sendmail(ToEmailID,...
		         'save error2');

		keyboard;
	%	save(['./matfiles/', pathtoMatFiletoSave,'errcondition.mat'],'errs');
 	end
	
	

try
zip('./matfiles/backup',{'./matfiles/sdplessprunedfullcontrol*.mat'});
zip('./matfiles/backup',{'./figfiles/*.jpg'});
setpref('Internet','SMTP_Server',SMTPname);
setpref('Internet','E_mail',FromEmailID);
sendmail(ToEmailID,...
         'Test subject','Test message',...
         {['./matfiles/backup.zip']});


			 catch


		setpref('Internet','SMTP_Server',SMTPname);
		setpref('Internet','E_mail',FromEmailID); 
		 sendmail(ToEmailID,...
		         'error in sending the zipped file');

			 end

end
