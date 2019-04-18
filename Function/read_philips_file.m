function [signal,step,frequency]= read_philips_file(file2read)

fin	= fopen([(file2read) '.spar'],'rt');
%modification Andy Devos - 11/10/2004: facility to read .SPAR/.SDAT files
% modification D. Sima - 04/2010: save position information form MRSI data
%                      - 01/2013: positopn corrected for 3D MRSI
if fin==-1
  fin	= fopen([(file2read) '.SPAR'],'rt')
end;

tmsrf	= 0; nrows=1; ncols=1; 
while 1
    [dummy,count] = fscanf(fin,'%s',1);
    if count < 1
        break;
    end
    if strcmp(dummy,'nucleus') == 1
        dummy = fscanf(fin,'%s',1);
        dummy = fscanf(fin,'%s',1);
        if strcmp(dummy,'"1H"') | strcmp(dummy,'"1h"') | strcmp(dummy,'1H')
            nuc = 'h';
        elseif strcmp(dummy,'"31P"') | strcmp(dummy,'"31p"') | strcmp(dummy,'31P')
            nuc = 'p';
        elseif strcmp(dummy,'"13C"') | strcmp(dummy,'"13c"') | strcmp(dummy,'13C')
            nuc = 'c';
        elseif strcmp(dummy,'"19F"') | strcmp(dummy,'"19f"') | strcmp(dummy,'19F')
            nuc = 'f';
        else
            nuc = dummy(2:length(dummy)-1);
        end
    elseif strcmp(dummy,'synthesizer_frequency') == 1
        dummy    = fscanf(fin,'%s',1);
        frequency	= fscanf(fin,'%f',1)/1000;
    elseif strcmp(dummy,'offset_frequency') == 1
        dummy    = fscanf(fin,'%s',1);
        woff	= fscanf(fin,'%f',1);
    elseif strcmp(dummy,'sample_frequency') == 1
        dummy    = fscanf(fin,'%s',1);
        sw	= fscanf(fin,'%f',1);
    elseif strcmp(dummy,'samples') == 1
        dummy    = fscanf(fin,'%s',1);
        ndp	= fscanf(fin,'%d',1);
    elseif strcmp(dummy,'averages') == 1
        dummy    = fscanf(fin,'%s',1);
        scans	= fscanf(fin,'%d',1);
    elseif strcmp(dummy,'ap_size') == 1
        dummy = fscanf(fin, '%s',1);
        ap_size= fscanf(fin, '%d',1);
    elseif strcmp(dummy,'lr_size') == 1
        dummy = fscanf(fin, '%s',1);
        lr_size= fscanf(fin, '%d',1);
    elseif strcmp(dummy,'cc_size') == 1
        dummy = fscanf(fin, '%s',1);
        cc_size= fscanf(fin, '%d',1);
    elseif strcmp(dummy,'nr_time_series_spectra') == 1
        dummy    = fscanf(fin,'%s',1);
        nfids	= fscanf(fin,'%d',1);
        if nfids < 1; nfids	= 1; end
        tmsrf	= 1; 
    end
    if tmsrf == 0
        if strcmp(dummy,'spec_num_row') == 1
            dummy    = fscanf(fin,'%s',1);
            nfids	= fscanf(fin,'%d',1);
            nrows = nfids;
            ncols = 1;
        end
        if strcmp(dummy,'nr_of_slices_for_multislice') == 1
            dummy    = fscanf(fin,'%s',1);
            slices	= fscanf(fin,'%d',1);
        end        
        if strcmp(dummy,'dim2_pnts') == 1
            dummy    = fscanf(fin,'%s',1);
            nrows	= fscanf(fin,'%d',1);
        end
        if strcmp(dummy,'dim3_pnts') == 1
            dummy    = fscanf(fin,'%s',1);
            ncols	= fscanf(fin,'%d',1);
        end
    end
end
fclose(fin);

if (nrows*ncols*slices > 1)
        pos1 =[]; pos2 = []; pos3 = [];
        for l = 1:slices
            for i = 1:nrows
                pos1 = [pos1; [1:ncols]'];
                pos2 = [pos2; i*ones(1,ncols)'];
            end
            pos3 = [pos3; l*ones(nrows*ncols,1)];
        end
        position = [[1:nfids]' pos1 pos2 pos3];
    else
        position = [1 1 1 1];
end
    

%Rene zegt: bij philips begin keihard op 0 zetten!

begin=0;


%gedeeld door 1000 omdat in de .spar file de samplefreq gegeven is in Hz, terwijl
%MRUI en zo met kHz werkt.
step	= 1000.0 / sw;

pmsdat=[(file2read) '.sdat'];
pmspar=[(file2read) '.spar'];

bpar=[1000 1 0];


pmsfid       = zeros(nfids,ndp);
%fin = fopen(pmsdat,'r','d');
fin = fopen(pmsdat,'r','ieee-le'); 
%modification Andy Devos - 11/10/2004: facility to read .SPAR/.SDAT
%files
if fin==-1
  pmsdat=[(file2read) '.SDAT'];
  %fin = fopen(pmsdat,'r','d'); %  works not after Matlab2008b
  fin = fopen(pmsdat,'r','ieee-le'); 
end;

i =  sqrt(-1);
for scnt = 1:nfids
    [sig1]       = freadVAXD(fin,[2 ndp],'float32');
    pmsfid(scnt,1:ndp) = bpar(2) * (sig1(1,1:ndp) + i *sig1(2,1:ndp));
end
fclose(fin);
ftit         = 'Philips data file, converted format (*.sdat)';

signal=pmsfid; 

volume=ap_size*lr_size*cc_size/1000;  %volume in ml!!!!!
