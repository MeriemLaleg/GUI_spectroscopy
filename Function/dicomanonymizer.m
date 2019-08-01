%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Christian Labadie -- 2-Nov-2004, 11-January-2008
% DICOM Siemens Trio anonymizer

function dicomanonymizer(fdcm)


if nargin==0, fdcm=uigetfile('*.*','Pick a DICOM file'); end

f = dir( fdcm ); 

fp = fopen( fdcm ,'r');
hdr = fread(fp,f.bytes,'uint8=>char')';
fclose(fp);
 

% ------------------ DICOM header ------------------


% Date, Time, Instance number
k = strfind(hdr,[char(hex2dec(['08';'00';'22';'00']))' 'DA' ]); % date
disp(sprintf('Datum    : %s', hdr(k+8:k+15)))

k = strfind(hdr,[char(hex2dec(['08';'00';'32';'00']))' 'TM' ]); % time
disp(sprintf('Time     : %s', hdr(k+8:k+15)))

k = strfind(hdr,[char(hex2dec(['20';'00';'10';'00']))' 'SH' char(hex2dec(['02';'00']))']);
disp(sprintf('StudyNr  : %s', deblank(hdr(k+8:k+10))))

k = strfind(hdr,[char(hex2dec(['20';'00';'11';'00']))' 'IS' char(hex2dec(['02';'00']))']);
disp(sprintf('SeriesNr : %s', deblank(hdr(k+8:k+10))))

k = strfind(hdr,[char(hex2dec(['20';'00';'12';'00']))' 'IS' char(hex2dec(['02';'00']))']);
disp(sprintf('AcqNr    : %s', deblank(hdr(k+8:k+10))))

k = strfind(hdr,[char(hex2dec(['20';'00';'13';'00']))' 'IS' char(hex2dec(['02';'00']))']);
disp(sprintf('ImaNr    : %s',deblank(hdr(k+8:k+10))));

k = strfind(hdr, [ char(hex2dec(['10';'00';'10';'00']))' ]);
disp('following fields are erased:')
while strcmp(hdr(k:k+1), char(hex2dec(['10';'00']))' )
    n = double(hdr(k+6)) + double(hdr(k+7))*256;
    disp(sprintf('%.2x,%.2x;%.2x,%.2x (%2d) : %s',double(hdr(k:k+3)), n, hdr(k+8:k+7+n)))
    for i=k+8:k+7+n, hdr(i)=' '; end
    k = k+8+n;
end

fp = fopen( [ fdcm '-anonymous' ],'w');
fwrite(fp,hdr,'uint8');
fclose(fp);

disp(['saved in ' fdcm '-anonymous']);