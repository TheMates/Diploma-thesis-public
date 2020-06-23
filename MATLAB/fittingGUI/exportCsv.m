function exportCsv(fname,Bm,Am,FIR)

separator = csvsep;
dlmwrite([fname '_Am.csv'],Am','delimiter',separator,'precision',15);
dlmwrite([fname '_Bm.csv'],Bm','delimiter',separator,'precision',15);
dlmwrite([fname '_FIR.csv'],FIR,'delimiter',separator,'precision',15);
end