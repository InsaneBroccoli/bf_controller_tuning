function autoclean(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% M.E.Peter, 05.09.2011 
% autoclean('path'): angegebenen Verzeichniss
% autoclean        : aktuelles Verzeichniss
% Es wird gelöscht: .asv
%                   .mexw32
%                   .mexw64
%                   .autosave
%                   XPC-GPA-files    
%                   slprj              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fprintf('+++ autoclean +++\n');
	
	if(nargin == 0)
		directory = dir();
		dirName = pwd();
	elseif(nargin == 1)
		dirName = varargin{1};
		if(strcmp(dirName(end-1:end),'\') ||...
			strcmp(dirName(end-1:end),'/'))
			dirName = dirName(1:end-1);
		end
		
		directory = dir(varargin{1});
	else
		error('Falsche Anzahl von Argumenten übergeben');
	end

	for i=1:length(directory)
		name = directory(i).name;
		if (length(name) >= 4 && strcmpi(name(end-3:end),'.asv'))
                delete([dirName,'\', name]);
                fprintf('   %s\n',name);
            elseif (length(name) >= 7 && strcmpi(name(end-6:end),'.mexw32'))
                delete([dirName,'\', name]);
                fprintf('   %s\n',name);
            elseif (length(name) >= 7 && strcmpi(name(end-6:end),'.mexw64'))
                delete([dirName,'\', name]);
                fprintf('   %s\n',name);
            elseif (length(name) >= 5 && strcmpi(name(end-4:end),'slprj'))
                rmdir([dirName,'\', name],'s');
                fprintf('   %s\n',name);
            elseif (length(name) >= 8 && strcmpi(name(end-7:end),'_xpc_rtw'))
                rmdir([dirName,'\', name],'s');
                fprintf('   %s\n',name);
            elseif (length(name) >= 9 && strcmpi(name(end-8:end),'.autosave'))
                delete([dirName,'\', name]);
                fprintf('   %s\n',name);
                
        end
        if (length(name) >= 4 && strcmpi(name(end-3:end),'.mdl'))
            for k=1:length(directory)
                XPCname = directory(k).name;
                if length(name) <= length(XPCname) 
                    if k~=i && strcmpi(name(1:length(name)-4),XPCname(1:length(name)-4))
                        if (length(XPCname) >= 4 && strcmpi(XPCname(end-3:end),'.dlm')) ||...
                           (length(XPCname) >= 4 && strcmpi(XPCname(end-3:end),'.xml')) ||...
                           (length(XPCname) >= 5 && strcmpi(XPCname(end-4:end),'bio.m')) ||...
                           (length(XPCname) >= 4 && strcmpi(XPCname(end-3:end),'pt.m')) ||...
                           (length(XPCname) >= 5 && strcmpi(XPCname(end-4:end),'ref.m'))
                            delete([dirName,'\', XPCname]);
                            fprintf('   %s\n',XPCname);
                        end
                    end
                end
            end
        end
        if (length(name) >= 4 && strcmpi(name(end-3:end),'.slx'))
            for k=1:length(directory)
                XPCname = directory(k).name;
                if length(name) <= length(XPCname) 
                    if k~=i && strcmpi(name(1:length(name)-4),XPCname(1:length(name)-4))
                        if (length(XPCname) >= 4 && strcmpi(XPCname(end-3:end),'.dlm')) ||...
                           (length(XPCname) >= 4 && strcmpi(XPCname(end-3:end),'.xml')) ||...
                           (length(XPCname) >= 5 && strcmpi(XPCname(end-4:end),'bio.m')) ||...
                           (length(XPCname) >= 4 && strcmpi(XPCname(end-3:end),'pt.m')) ||...
                           (length(XPCname) >= 5 && strcmpi(XPCname(end-4:end),'ref.m'))
                            delete([dirName,'\', XPCname]);
                            fprintf('   %s\n',XPCname);
                        end
                    end
                end
            end
        end
	end