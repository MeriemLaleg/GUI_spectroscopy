%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Christian Labadie -- 2-Nov-2004, 11-January-2008
% DICOM Siemens Trio: reads IMA header, svs and mri
% Usage: 
%   h = simpledicom(uigetfile('Pick a DICOM');
%   if h.isimage
%       imagesc(h.image), axis image
%   else
%     plot(h.ppm, abs(h.spec)), set(gca,'XDir','reverse');
%   end

function h=simpledicom(fdcm)

f = dir( fdcm ); 

SF = 0; DW = 0;
fp = fopen( fdcm ,'r');
hdr = fread(fp,f.bytes,'uint8=>char')';

% ------------------ DICOM header ------------------

% Date, Time, Instance number
k = strfind(hdr,[char(hex2dec(['08';'00';'22';'00']))' 'DA' ]); % date
Datum = hdr(k+8:k+15);

k = strfind(hdr,[char(hex2dec(['08';'00';'32';'00']))' 'TM' ]); % time
Time = hdr(k+8:k+15);

k = strfind(hdr,[char(hex2dec(['20';'00';'10';'00']))' 'SH' char(hex2dec(['02';'00']))']);
StudyNr = deblank(hdr(k+8:k+10));   %strtrim

k = strfind(hdr,[char(hex2dec(['20';'00';'11';'00']))' 'IS' char(hex2dec(['02';'00']))']);
SeriesNr = deblank(hdr(k+8:k+10));  %strtrim

k = strfind(hdr,[char(hex2dec(['20';'00';'12';'00']))' 'IS' char(hex2dec(['02';'00']))']);
AcqNr = deblank(hdr(k+8:k+10));     %strtrim

k = strfind(hdr,[char(hex2dec(['20';'00';'13';'00']))' 'IS' char(hex2dec(['02';'00']))']);
ImaNr = deblank(hdr(k+8:k+10));     %strtrim

h.f.name = f.name;
h.f.bytes = f.bytes;

% extract MrProtocol from the DICOM header
ibeg = strfind(hdr, ['### ASCCONV BEGIN ###'] );
iend = strfind(hdr,'### ASCCONV END ###');
if(length(ibeg)>1 || length(iend)>1)
    count_beg = 0;
    count_end = 0;
    found = 0;
    while(~found)
        if(ibeg(count_beg)>iend(count_end))
            count_end = count_end+1;
        elseif(ibeg(count_beg+1)<iend(count_end))
            count_beg = count_beg+1;
        else
            ibeg = ibeg(count_beg);
            iend = iend(count_end);
            found = 1;
        end
    end
end
pro = hdr(ibeg:iend+19);
str{1} = fileparts(fdcm);
str{2} = f.name;
str{3} = [num2str(f.bytes) ' bytes' ];
if ~isempty(pro)
    h.f.study = str2num(StudyNr);
    h.f.series = str2num(SeriesNr);
    h.f.acq = str2num(AcqNr);
    h.f.ima = str2num(ImaNr);
%     h.f.date = [ Datum(1:4) '.' Datum(5:6) '.' Datum(7:8) ];
%     h.f.time = [ Time(1:2) ':' Time(3:4) ':' Time(5:end) ];
%     str{4} = [ 'Study:' StudyNr '   Series:' SeriesNr ...
%             '   Acquisition:' AcqNr '   Instance/Image:' ImaNr ];
%     str{5} = ['Date: ' h.f.date  '  Time: ' h.f.time ];
end
idx = strfind(pro,char(10)); istart=1;
for k=1:length(idx)
    str{k+6} = pro(istart:idx(k)-1);
    istart = idx(k)+1;
    ieq = strfind(str{k}, ' = ');
    if length(ieq) == 1
        key = deblank(str{k}(1:ieq));   %strtrim
        value = str{k}(ieq+3:end);
        str{k} = [ key ' = ' value ];
        if strcmp(key,'sTXSPEC.asNucleusInfo[0].lFrequency'), SF = str2num(value); end
        if strcmp(key,'sRXSPEC.alDwellTime[0]'), DW = str2num(value); end
        if strcmp(key,'sSpecPara.dDeltaFrequency'), DeltaFreq = str2num(value); end
    end
end
h.f.hdr = str';


% DICOM data
fid = strfind(hdr,[char(hex2dec(['e1';'7f';'10';'10']))' 'OB' char(hex2dec(['00';'00']))']);
pic = strfind(hdr,[char(hex2dec(['e0';'7f';'10';'00']))' 'OW' char(hex2dec(['00';'00']))']);

h.isimage = false;
if length(fid) == 1
    h = mrsvs(h, fp, fid, DW, SF);
elseif length(pic) == 1
    row = strfind(hdr,[char(hex2dec(['28';'00';'10';'00']))' 'US' char(hex2dec(['02';'00']))']);
    col = strfind(hdr,[char(hex2dec(['28';'00';'11';'00']))' 'US' char(hex2dec(['02';'00']))']);
    h = mrimage(h, fp, pic, row(1), col(1));
else
    h.fid = 'unknown format';
end
fclose(fp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MR single voxel spectroscopy
function h = mrsvs(h, fp, k, DW, SF)

fseek(fp, k+7, 'bof');  fidsize = fread(fp,1,'int32','ieee-le') / 4; % 4 bytes for float32
fseek(fp, k+11, 'bof'); fid = fread(fp, fidsize,'float32','ieee-le');
fid=reshape(fid,2,length(fid)/2);
h.fid = complex(fid(1,:) ,fid(2,:));
NP = length(fid);

h.spec = fftshift(fft(h.fid, length(h.fid) ));  % fft(X,N) N-point FFT, padded with zeros

h.ppm = 1:-1/(NP-1):0;
if SF > 0 && DW > 0
    DW = DW / 1000000000; % nanoseconds
    SW = 1 / 2 / DW;
    h.ppm = h.ppm * SW / SF * 1000000;
end

h.f.NP = NP;
h.f.SF = SF;
h.f.DW = DW;
%h.f.SW = SW;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MR image
function h=mrimage(h, fp, k, row, col)
fseek(fp, row + 7, 'bof');  row = fread(fp,1,'int16','ieee-le');
fseek(fp, col + 7, 'bof');  col = fread(fp,1,'int16','ieee-le');
fseek(fp, k+7, 'bof');  picsize = fread(fp,1,'int32','ieee-le') / 2;       % 2 bytes for int16
fseek(fp, k+11, 'bof'); 
h.image = reshape(fread(fp,picsize,'int16','ieee-le'), col, row);
h.isimage = true;
