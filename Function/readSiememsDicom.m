function [fid,h] = readSiememsDicom(file2read)

h = simpledicom(file2read);

frequency = h.f.SF/1e3;
 
fid=h.fid;