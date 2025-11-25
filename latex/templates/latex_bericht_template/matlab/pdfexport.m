function pdfexport(pdfname,pdfsize,fontsize,pdfpath)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Erzeugung eine pdf-Dokumentes aus einer Matlab-Grafik
% M.E.Peter, 25.01.2012
%
% pdfexport(pdfname,pdfsize,fontsize,pdfpath)
% -------------------------------------------
%
% pdfname  := String mit Namen, default: 'pdfplot'
% pdfsize  := Grösse pdf in cm, default: [15 10]
% fontsize := Schriftgrösse, default: 11
% pdfpath  := String mit Speicherpfad, default: pwd (aktueller Pfad)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('pdfname','var');
    pdfname = 'pdfplot';
end
if ~exist('pdfsize','var');
%     pdfsize = [14 9]; % Bode
    pdfsize = [15 10]; % Zeitsignale
end
if ~exist('pdfpath','var');
    pdfpath = 'E:\pmic\temp\latex_bericht_template\images';
%     pdfpath = pwd;
end
if ~exist('fontsize','var');
    fontsize = 11;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

colorVec = [0.1 0.1 0.1];
lineobj = findobj('type','line');
set(lineobj, 'Linewidth',1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = gcf;
hchild = get(h,'children');
for k = 1:length(hchild)
    set(hchild(k),'linewidth',0.5,'fontsize',fontsize,'fontname','times',...
        'xcolor',colorVec,'ycolor',colorVec,'zcolor',colorVec);                 
    set(get(hchild(k),'xlabel'),'fontname','times',...
        'fontSize',fontsize,'color',colorVec);
    set(get(hchild(k),'ylabel'),'fontname', 'times',...
        'fontSize',fontsize,'color',colorVec);
    set(get(hchild(k),'title'),'fontname', 'times',...
        'fontSize',fontsize,'color',colorVec);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(h,'PaperUnits','centimeters');
set(gca,'LooseInset',get(gca,'TightInset'));
set(h,'PaperPosition',[0 0 pdfsize(1) pdfsize(2)]); 
set(h,'PaperSize',pdfsize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

current_path = pwd;
print(h,'-dpdf',pdfname);   
cd(pdfpath);                              
print(h,'-dpdf',pdfname);              
cd(current_path);                         
fprintf('+++ pdfexport +++\n');
display([pdfname,'.pdf successfully save to: ',pwd]);
display([pdfname,'.pdf successfully save to: ',pdfpath]);

return

