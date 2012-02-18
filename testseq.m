global SMTPname FromEmailID ToEmailID;
SMTPname = ;
FromEmailID = ;
ToEmailID =  ;
try
 		usingYalmipMainSu4Iterations(3);
	catch
		errs=lasterror;
	    setpref('Internet','SMTP_Server',SMTPname);
	    setpref('Internet','E_mail',FromEmailID);
        sendmail(ToEmailID,...
	         'err_macmachinemainfile',errs.message);
end

