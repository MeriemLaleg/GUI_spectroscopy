% SVSVIEWER for Siemens spectroscopy VA25 or VB15 DICOM files
% Christian Labadie -- 9-13 Juli 2007, updated 9 April 2008
% Max Planck Institute for Human Cognitive and Brain Sciences, Leipzig
%
% Tipps: - the current DICOM data is accessible in the workspace with:
%          >> disp(h.f)
%        - [copy] button makes a copy of the MrProtocol in the clipboard
%        - [save as...] exports in MRUI text format or anonymized DICOM

function h = svsviewer(func, hg, arg)

% your prefered DICOM folders can be listed here:
myLocations = {  '/scr/elfi1/Daten/Hepta' ...
                 '/scr/elfi1/Daten/TimTrio' '/scr/elfi1/Daten/Trio' ...
                 'example/unix/folder' '~/mydicom/folder' ...
                 'C:\example\PC\folder' };

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin==2 && strcmp(func,'folder') && exist(hg,'dir')

    evalin('base','if exist(''h'',''var''), assignin(''caller'',''h'',h); end');
    h.curdir = hg;
    addcd = true;
    func='';

elseif nargin < 2

    if nargin==1, myLocations = { func myLocations{:} }; end
    h = startgui(myLocations);
    ima = folders(h);
    h = svsview(h, fullfile(h.curdir, ima));
    if nargout == 0, assignin('base', 'h', h); clear; end
    return

elseif nargin==3 && strcmp(func,'cd') && exist(arg,'dir')

    disp(arg)
    h = get(hg, 'UserData');
    h.curdir = arg;
    addcd=true;
    func='';

elseif nargin==2 && strcmp(func,'cpclip')

    h = get( get(hg,'Parent'), 'UserData' );
    copyhdr(h);
    return

elseif nargin==2 && strcmp(func,'saveas')

    h = get( get(hg,'Parent'), 'UserData' );
    dcmsaveas(h);
    return

else

    h = get( get(hg,'Parent'), 'UserData' );
    str = get(hg, 'String'); str = str( get(hg, 'Value') );
    addcd=false;

end

if strcmp(func,'cd') && ~strcmp(h.curdir, str{1})

    cd( h.curdir )
    d = dir(str{1});
    if ~isempty(d)
        cd( str{1} )
        h.curdir = pwd;
        cd( h.mdir )
        addcd = true;
        uicontrol(h.gui.cf)  % this could be a problem for MATLAB 6
    else
        cd( h.mdir )
    end

elseif strcmp(func, 'cr')

    directoryname = uigetdir(h.curdir, 'Choose a Directory');
    if directoryname ~= 0
        cd( directoryname )
        h.curdir = pwd;
        cd(h.mdir)
        addcd=true;
    end

elseif strcmp(func, 'cl')

    directoryname = menu('Choose a Directory', h.folders);
    if directoryname ~= 0
        cd( h.folders{directoryname} )
        h.curdir = pwd;
        cd(h.mdir)
        addcd=true;
    end

elseif strcmp(func, 'cf')

    if ~isempty(str), h = svsview(h, fullfile(h.curdir,parseImaListing(str{1}))) ; end

end

if addcd
    ima = folders(h);
    h = svsview(h, [ h.curdir filesep ima ]);
    % uicontrol(h.gui.cf)
    if h.f.isdicom
        n = length(h.folders);
        ipos=n+1;
        for i=1:n, if strcmp(h.folders{i},h.curdir), ipos=0; end; end
        if ipos>0, h.folders{ipos}=h.curdir; end
    end
end

set(h.gui.fig,'UserData', h);

if nargout == 0
    assignin('base', 'h', h);
    clear
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reads IMA header, displays the FID and FFT
function h = svsview(h, IMAname)

f = dir( IMAname );

if length(f) ~= 1 || f.isdir
    guiclear(h);
    set(h.gui.pro,'String','');
    set(h.gui.cs,'String','');
    return
end
   
[dSt,fSt,eSt] = fileparts(IMAname);

if f.bytes>1024*1024*1024, fInfo = sprintf(' (%.2f GB)', f.bytes/1024/1024/1024);
elseif f.bytes>1024*1024, fInfo = sprintf(' (%.2f MB)', f.bytes/1024/1024);
elseif f.bytes>1024, fInfo = sprintf(' (%.2f KB)', f.bytes/1024);
else fInfo = sprintf(' (%.d bytes)', f.bytes);
end

set(h.gui.cs,'String', [fSt eSt fInfo]);
if strcmpi(eSt,'.zip'), guiclear(h); return; end

fp = fopen( IMAname ,'r');
hdr = fread(fp,f.bytes,'uint8=>char')';
h.f = dicominfo(f,hdr);  % DICOM header
h.f.hdr{1} = h.curdir;
h.f.hdr{2} = f.name;
h.f.hdr{3} = [ 'File:' num2str( get(h.gui.cf,'Value') ) ' ' fInfo];

if h.f.isdicom
    h.f.hdr{4} = [ 'Study:' h.f.study '  Series:' h.f.series ...
        '  Acquisition:' h.f.acq '  Instance:' h.f.ima ...
        '  Date: ' h.f.date  '  Time: ' h.f.time ];
    h.f.hdr{5} = ['Sequence: ' h.f.sequence '  Description: ' h.f.description ];
    set(h.gui.cs, 'String', [ '[' h.f.series '] ' max(h.f.acq,h.f.ima) ', ' h.f.description ' (' h.f.date ' ' h.f.time ')']);
end

isel = get(h.gui.pro, 'Value');
if isel > length(h.f.hdr), isel=1; end
set(h.gui.pro, 'String', h.f.hdr, 'Value', isel);

% DICOM data
fid = strfind(hdr,[char(hex2dec(['e1';'7f';'10';'10']))' 'OB' char(hex2dec(['00';'00']))']);
pic = strfind(hdr,[char(hex2dec(['e0';'7f';'10';'00']))' 'OW' char(hex2dec(['00';'00']))']);

%% multivoxel

if length(fid) == 1
    h = mrsvs(h, fp, fid);
    fclose(fp);
    h.f.isdicom = true;
elseif length(pic) == 1
    row = strfind(hdr,[char(hex2dec(['28';'00';'10';'00']))' 'US' char(hex2dec(['02';'00']))']);
    col = strfind(hdr,[char(hex2dec(['28';'00';'11';'00']))' 'US' char(hex2dec(['02';'00']))']);
    mrimage(h, fp, pic, row(1), col(1));
    fclose(fp);
    h.f.isdicom = true;
else
    fclose(fp);
    if strcmpi(eSt,'.bmp') || strcmpi(eSt,'.jpg') || strcmpi(eSt,'.png') || strcmpi(eSt,'.gif')
        axes(h.gui.fft), cla, set(h.gui.fft, 'Visible', 'off')
        axes(h.gui.fid), cla, set(h.gui.fid, 'Visible', 'off')
        axes(h.gui.mag), set(h.gui.mag, 'Visible', 'on')
        image(imread(IMAname));
        axis(h.gui.mag,'off'), axis(h.gui.mag,'image');
    else
        guiclear(h)
    end
    h.f.isdicom = false;
end

if h.f.isdicom
    set(h.gui.td, 'Visible', 'on')
else
    set(h.gui.td, 'Visible', 'off')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MR single voxel spectroscopy
function h = mrsvs(h, fp, k)

axes(h.gui.mag), cla, set(h.gui.mag, 'Visible', 'on')
axes(h.gui.fft), cla, set(h.gui.fft, 'Visible', 'on')
axes(h.gui.fid), cla, set(h.gui.fid, 'Visible', 'on')

fseek(fp, k+7, 'bof');  fidsize = fread(fp,1,'int32','ieee-le') / 4; % 4 bytes for float32
fseek(fp, k+11, 'bof'); fid = fread(fp, fidsize,'float32','ieee-le');
fid=reshape(fid,2,length(fid)/2);
fid=complex(fid(1,:) ,fid(2,:));
NP = length(fid);

xr = 1:NP; xi = xr + NP;
axes(h.gui.fid), plot(xr, real(fid), xi, imag(fid), 'g')
v = axis;
axis( [ 0 2*NP+1 v(3) v(4) ] )
grid on

spec = fftshift(fft(fid, length(fid) ));  % fft(X,N) N-point FFT, padded with zeros
[ phi0 phased ] = simplephase(spec);
axes(h.gui.fft), plot(xr, real(phased), xi, imag(phased), 'g');
v = axis;
axis( [ 0 2*length(fid)+1 v(3) v(4) ] )
grid on

axes(h.gui.mag)
ppm = 1:-1/(NP-1):0;
if h.f.SF > 0 && h.f.DW > 0
    DW = h.f.DW / 1e+9; % DwellTime is given in nanoseconds
    SW = 1 / (2 * DW);
    ppm = ppm * SW / h.f.SF * 1e+6;
end
plot(ppm, abs(spec));
v = axis;
axis( [ ppm(end) ppm(1) v(3) v(4) ] )
set(gca,'XDir','reverse')
xlabel('ppm')
title('magnitude')

h.f.NP = NP;
h.f.DW = DW;
h.f.SW = SW;
h.f.ppm = ppm;
h.f.fid = fid;
h.f.spec = spec;
h.f.phi0 = phi0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simple SVS zero order phase
function [ phi0 , phased ] = simplephase(spec)
best = 0; dphi = pi / 1.5; phi0 = 0;
for step = 1:4
    for phi = phi0-1.5*dphi:dphi / 10:phi0+1.5*dphi
        dum = exp(1i * phi) * spec;
        signal = sum (real( dum ) );
        if signal > best, best = signal; phi0 = phi; phased = dum; end
    end
    dphi = dphi / 10;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MR image
function mrimage(h, fp, k, row, col)

axes(h.gui.fft), cla, set(h.gui.fft, 'Visible', 'off')
axes(h.gui.fid), cla, set(h.gui.fid, 'Visible', 'off')
axes(h.gui.mag), set(h.gui.mag, 'Visible', 'on')

fseek(fp, row + 7, 'bof');  row = fread(fp,1,'int16','ieee-le');
fseek(fp, col + 7, 'bof');  col = fread(fp,1,'int16','ieee-le');
fseek(fp, k+7, 'bof');  picsize = fread(fp,1,'int32','ieee-le') / 2;       % 2 bytes for int16
fseek(fp, k+11, 'bof');  imagesc(reshape(fread(fp,picsize,'int16','ieee-le'), col, row));
axis(h.gui.mag,'off'), axis(h.gui.mag,'image'), colormap gray;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start GUI and create layout
function h = startgui(myLocations)

h.mdir = fileparts(mfilename('fullpath'));
h.mfile = mfilename;

% check prefered locations and determine PC and Unix paths
nf=1;
for i=1:length(myLocations)
    if exist([myLocations{i}],'dir'), p{nf}=myLocations{i}; nf=nf+1; disp(myLocations{i}); end
end
p{nf} = pwd;
h.curdir = p{1};

if strcmp(computer,'PCWIN')  % PC standard path
    [s,w] = dos('echo %USERPROFILE%'); w = [deblank(w) '\Desktop'];
    if s==0 && exist(w,'dir') && ~strcmp(w,p{1}), p{length(p)+1} = w; end
    tmp='CDEFGHIJKLMNOPQRSTUVWXYZAB';
    for i=1:length(tmp)
        w=[tmp(i) ':\']; if exist(w,'dir'), p{length(p)+1} = w; end
    end
else  % unix standard path
    p{nf+1} = '~/'; p{nf+2} = '/';
end
h.folders = p;

% determine screen size
p = get(0,'ScreenSize');
h.gui.fig = figure('Name', 'SVS viewer', 'Position',[ 50 p(4)-670-55 925 640 ],'Toolbar','Figure');

h.gui.mag = axes('Parent',h.gui.fig, 'Position',[ 0.038 0.580 0.440 0.375 ],'Box','on');
h.gui.fft = axes('Parent',h.gui.fig, 'Position',[ 0.534 0.737 0.439 0.222 ],'Box','on'); set(gca,'XTickLabel',[])
h.gui.fid = axes('Parent',h.gui.fig, 'Position',[ 0.534 0.449 0.439 0.222 ],'Box','on');

% current file selection
h.gui.cs = uicontrol('Parent', h.gui.fig, 'Units', 'normalized', 'Position', [0.431+0.055 0.356 0.542-0.055 0.021], ...
    'Style', 'Text', 'HorizontalAlignment', 'left','Backgroundcolor',[0.8 0.8 0.8]);

% current directory
h.gui.ct = uicontrol('Parent', h.gui.fig, 'Units', 'normalized', 'Position', [0.431 0.378   0.542 0.021], ...
    'String', h.curdir, 'Style', 'Text', 'HorizontalAlignment', 'left','Backgroundcolor',[0.8 0.8 0.8]);

% protocol viewer
h.gui.pro = uicontrol('Parent', h.gui.fig, 'Units', 'normalized', 'Position', [0.431 0.019 0.542 0.336 ], ...
    'String', {}, 'Style','listbox', 'Backgroundcolor','w', 'Min',0,'Max',2);

% root selection
h.gui.cr = uicontrol('Parent', h.gui.fig, 'Units', 'normalized', 'Position', [0.157 0.454 0.041 0.043], ...
    'String', 'cd...', 'Style', 'pushbutton', 'Callback', [mfilename '(''cr'', gcbo);']);

% root selection
h.gui.cl = uicontrol('Parent', h.gui.fig, 'Units', 'normalized', 'Position', [0.205 0.454 0.066 0.042], ...
    'String', 'folders', 'Style', 'pushbutton', 'Callback', [mfilename '(''cl'', gcbo);']);

% root selection
h.gui.td = uicontrol('Parent', h.gui.fig, 'Units', 'normalized', 'Position', [0.278 0.454 0.066 0.042], ...
    'String', 'save as...', 'Style', 'pushbutton', 'Callback', [ mfilename '(''saveas'', gcbo);'], ...
    'Visible','off');

% file selection
h.gui.cf = uicontrol('Parent', h.gui.fig,'Units', 'normalized', 'Position',[0.157 0.019 0.246 0.419], ...
    'String', {}, 'Style','listbox', 'Callback',[ mfilename '(''cf'', gcbo);']);

% directory selection
h.gui.cd = uicontrol('Parent', h.gui.fig, 'Units', 'normalized', 'Position', [0.017 0.019 0.125 0.48], ...
    'String', {}, 'Style','listbox', 'Callback',[mfilename '(''cd'', gcbo);']);

% copytoclipboard
h.gui.clip = uicontrol('Parent', h.gui.fig, 'Units', 'normalized', 'Position', [0.431 0.356 0.05 0.021], ...
    'Style', 'PushButton', 'String','copy', 'Callback',[mfilename '(''cpclip'', gcbo);']);

h.f.isdicom = false;

set(h.gui.fig, 'UserData', h);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create GUI layout
function guiclear(h)
axes(h.gui.mag), cla, set(h.gui.mag, 'Visible', 'off')
axes(h.gui.fft), cla, set(h.gui.fft, 'Visible', 'off')
axes(h.gui.fid), cla, set(h.gui.fid, 'Visible', 'off')
set(h.gui.td, 'Visible', 'off')
% set(h.gui.pro,'String','');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create file selection listboxes
function ima = folders(h)

skey=[];
str={ h.curdir }; str{2} = '..';
fstr = {}; fima = {};
idir = 3; ifile = 1; ima = 1;
dum = dir( h.curdir );
n=0;
for i=1:length(dum), if dum(i).isdir && dum(i).name(1) ~= '.', n = n+1; end; end
dumord=1:length(dum); if n>21, dumord=length(dum):-1:1; end
for k=dumord

    if dum(k).isdir && dum(k).name(1) ~= '.'
        str{idir} = dum(k).name;
        idir = idir + 1;

    elseif dum(k).name(1) ~= '.'

        idx = strfind(dum(k).name,'.');
        if length(idx) > 8
            skey(ima) = 1000 * str2double(dum(k).name(idx(3)+1:idx(4)-1)) + ...
                str2double(dum(k).name(idx(4)+1:idx(5)-1));
            fima{ima} = dicomheader(fullfile(h.curdir, dum(k).name));
            ima = ima + 1;
        else
            fstr{ifile} = dum(k).name;
            ifile = ifile + 1;
        end

    end

end

set(h.gui.cd, 'String', str, 'Value', 1);

% sort the IMA files according to their Series and Instance numbers
ima = '';
if ~isempty(skey)
    [skey,idx] = sort(skey);
    fstr = {fima{ idx } fstr{:}};
end

if ~isempty(fstr), ima = parseImaListing(fstr{1}); end
set(h.gui.cf, 'String', fstr, 'Value', 1);

st = h.curdir;
if length(st) > 80, st=[ '... ' st(end-80:end) ]; end
st=sprintf('%s  (%d files)', st, length(fstr));
set(h.gui.ct, 'String', st);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse a Siemens MrProtocol string for display
function st=strprot(st)
if strcmp(st(1:2),'""'), st=st(3:end); end
if strcmp(st(1),'"'), st=st(2:end); end
if strcmp(st(end-1:end),'""'), st=st(1:end-2); end
if strcmp(st(end:end),'"'), st=st(1:end-1); end
st=regexprep(st, '+AF8-', '_');
st=regexprep(st, '+AEA-', '@');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% partial check of a DICOM header
function [listing, dcm]=dicomheader(IMAname)
f = dir( IMAname );
b = min(8*1024,f.bytes);  % 8KB should be enough to read some DICOM info
fp = fopen( IMAname ,'r');
hdr = fread(fp,b,'uint8=>char')';
dcm = dicominfo(f,hdr);
fclose(fp);
listing=['[' dcm.series '] ' num2str(max(str2double(dcm.acq),str2double(dcm.ima))) ', ' dcm.description ' | ' f.name];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% quick parse of a DICOM header
function dcm=dicominfo(f,hdr)
dcm.name = f.name;
dcm.bytes = f.bytes;
dcm.isdicom = false;
dcm.instrument = '';
dcm.description = '';
dcm.date = '';
dcm.time = '';
dcm.study = '';
dcm.series = '';
dcm.acq = '';
dcm.ima = '';
dcm.hdr = {};
dcm.sequence = '';
dcm.protocol = '';
dcm.B0 = '';
dcm.nucleus = '';
dcm.SF = 1;
dcm.DW = 1;

if length(hdr)<132 || ~strcmp(hdr(129:132),'DICM'), return; end
dcm.isdicom = true;

% parse some DICOM header fields
% date
k = strfind(hdr,[char(hex2dec(['08';'00';'22';'00']))' 'DA' ]);
if isempty(k), return; end
k = hdr(k+8:k+15);
dcm.date = [ k(1:4) '.' k(5:6) '.' k(7:8) ];

% time
k = strfind(hdr,[char(hex2dec(['08';'00';'32';'00']))' 'TM' ]);
if isempty(k), return; end
k = hdr(k+8:k+15);
dcm.time = [ k(1:2) ':' k(3:4) ':' k(5:end) ];

% manufacturer model's name
k = strfind(hdr,char(hex2dec(['08';'00';'90';'10']))');
if isempty(k), return; end
n = double(hdr(k+6)) + double(hdr(k+7))*256;
dcm.instrument = deblank(hdr(k+8:k+7+n));

% series description
k = strfind(hdr,char(hex2dec(['08';'00';'3E';'10']))');
if isempty(k), return; end
n = double(hdr(k+6)) + double(hdr(k+7))*256;
dcm.description = deblank(hdr(k+8:k+7+n));

% study, series, instance, acquisition
k = strfind(hdr,[char(hex2dec(['20';'00';'10';'00']))' 'SH' char(hex2dec(['02';'00']))']);
if isempty(k), return; end
dcm.study = deblank(hdr(k+8:k+10));

k = strfind(hdr,[char(hex2dec(['20';'00';'11';'00']))' 'IS' char(hex2dec(['02';'00']))']);
if isempty(k), return; end
dcm.series = deblank(hdr(k+8:k+10));

k = strfind(hdr,[char(hex2dec(['20';'00';'12';'00']))' 'IS' char(hex2dec(['02';'00']))']);
if isempty(k), return; end
dcm.acq = deblank(hdr(k+8:k+10));

k = strfind(hdr,[char(hex2dec(['20';'00';'13';'00']))' 'IS' char(hex2dec(['02';'00']))']);
if isempty(k), return; end
dcm.ima = deblank(hdr(k+8:k+10));

% extract MrProtocol from the DICOM header
ibeg = strfind(hdr, ['### ASCCONV BEGIN ###' char(10) 'ulVersion'] );
if isempty(ibeg), return; end  % MrProtocol not found
iend = strfind(hdr,'### ASCCONV END ###');
if isempty(iend), iend=length(hdr); else iend=iend+18; end
hdr = hdr(ibeg(1):iend(1));
hdr(length(hdr)+1) = char(10);
idx = strfind(hdr,char(10));
istart = 1;
for k=1:length(idx)
    dcm.hdr{k+6} = hdr(istart:idx(k)-1);
    istart = idx(k)+1;
    ieq = strfind(dcm.hdr{k+6}, ' = ');
    if length(ieq) == 1
        key = deblank(dcm.hdr{k+6}(1:ieq));   %strtrim
        value = dcm.hdr{k+6}(ieq+3:end);
        dcm.hdr{k+6}=[key ' = ' value];
        if strcmp(key,'tSequenceFileName'), dcm.sequence = strprot(value); end
        if strcmp(key,'tProtocolName'), dcm.protocol = strprot(value); end
        if strcmp(key,'sTXSPEC.asNucleusInfo[0].lFrequency'), dcm.SF = str2num(value); end
        if strcmp(key,'sTXSPEC.asNucleusInfo[0].tNucleus'), dcm.nucleus = value(3:end-2); end
        if strcmp(key,'sRXSPEC.alDwellTime[0]'), dcm.DW = str2num(value); end
        if strcmp(key,'sSpecPara.dDeltaFrequency'), DeltaFreq = value; end
        if strcmp(key,'sProtConsistencyInfo.flNominalB0'), dcm.B0 = value; end
        if strcmp(key,'sSpecPara.lVectorSize'), dcm.VectorSize = str2num(value); end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filename in GUI listbox
function f=parseImaListing(f)
[a,b]=strtok(f,'|');
if length(b)>2, f=b(3:end); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% copy header to the clipboard
function copyhdr(h)
nl=sprintf('\n');
if strcmp(computer,'PCWIN'), nl=sprintf('\r\n'); end
data = [h.f.hdr{1} nl];
for i=2:length(h.f.hdr)
    data = [ data h.f.hdr{i} nl];
end
if strcmp(computer,'PCWIN')  % PC standard path
    clipboard('copy',data)
else
    cpobj = com.mathworks.page.utils.ClipboardHandler;
    cpobj.copy(data);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% saveas user interface
function dcmsaveas(h)
if ~h.f.isdicom, return; end

fp = [ h.f.series '_' h.f.description '_' max(h.f.acq,h.f.ima) ];
[fp, d] = uiputfile({'*.txt', 'jMRUI Text files (*.txt)'; ...
    '*.dcm', 'DICOM anonymized (*.dcm)'}, 'Save as', fullfile(h.curdir,fp));
if ~fp, return; end
[a,b,c]=fileparts(fp);
if strcmpi(c,'.dcm')
    saveasdicomanon(fullfile(d,fp), h)
elseif strcmpi(c,'.txt')
    saveasmrui(fullfile(d,fp), h.f)
else
    disp(['unsupported format : ' c]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% saveas jMRUI MRSI text file
function saveasmrui(fp, f)
fp=fopen(fp,'w');
fprintf(fp,'jMRUI Data Textfile\n\nFilename: %s\n\n',f.name);
fprintf(fp,'PointsInDataset: %d\n',f.VectorSize); % updated for CSI
fprintf(fp,'DatasetsInFile: %d\n',f.NP/f.VectorSize); % updated for CSI
fprintf(fp,'SamplingInterval: %f\n',f.DW*1000*2); % 2*DW converted to milliseconds
fprintf(fp,'ZeroOrderPhase: 0E0\n');   % f.phi0/pi*180
fprintf(fp,'BeginTime: 0E0\n');
fprintf(fp,'TransmitterFrequency: %d\n',f.SF);
fprintf(fp,'MagneticField: %s\n', f.B0);
fprintf(fp,'TypeOfNucleus: %s\n', f.nucleus);
fprintf(fp,'NameOfPatient: anonymous\n');
fprintf(fp,'DateOfExperiment: %s %s\n',f.date,f.time);
fprintf(fp,'Spectrometer: %s\n',f.instrument);

fprintf(fp,'DicomStudy: %s\n',f.study);
fprintf(fp,'DicomSeries: %s\n',f.series);
fprintf(fp,'DicomAcquisition: %s\n',f.acq);
fprintf(fp,'DicomInstance: %s\n',f.ima);
fprintf(fp,'PulseSequence: %s\n',f.sequence);
fprintf(fp,'ProtocolName: %s\n',f.protocol);
fprintf(fp,'AdditionalInfo: Siemens MrProtocol follows\n');
for i=7:length(f.hdr)
    fprintf(fp,'%s\n',f.hdr{i});
end

fprintf(fp,'\nSignal\nsig(real)\tsig(imag)\n')
% save the complex conjugate of the Siemens FID
for i = 1:f.NP/f.VectorSize
    fprintf(fp,'Signal %d out of %d in file\n',i,f.NP/f.VectorSize);
    for k=1:f.VectorSize
        fprintf(fp,'%e\t%e\n',real(f.fid((i-1)*f.VectorSize+k)),-imag(f.fid((i-1)*f.VectorSize+k))); % complex conjugate
    end
end
fclose(fp);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% saveas DICOM anonymous
function saveasdicomanon(fp,h)

fn = fullfile(h.curdir,h.f.name);
f = dir(fp); 
fn = fopen(fn,'r');
hdr = fread(fn,f.bytes,'uint8=>char')';
fclose(fn);

k = strfind(hdr, [ char(hex2dec(['10';'00';'10';'00']))' ]);
nk=0;
while strcmp(hdr(k:k+1), char(hex2dec(['10';'00']))' )
    nk=nk+1; patient(nk).k = k;
    switch hdr(k+2:k+3)
        case char(hex2dec(['10';'00']))', key='Name';
        case char(hex2dec(['20';'00']))', key='ID';
        case char(hex2dec(['30';'00']))', key='Birthdate';
        case char(hex2dec(['40';'00']))', key='Sex';
        case char(hex2dec(['10';'10']))', key='Age';
        case char(hex2dec(['30';'10']))', key='Weight';
        case char(hex2dec(['00';'40']))', key='Comments';
        otherwise, key='unkown';
    end
    patient(nk).n = double(hdr(k+6)) + double(hdr(k+7))*256;
    prompt{nk}=sprintf('%s (%.2x,%.2x %2d bytes) : %s', ...
        key, hdr(k+2:k+3), patient(nk).n, hdr(k+8:k+7+patient(nk).n) );
    def{nk}='';
    k = k+8+patient(nk).n;
end

value = inputdlg(prompt,'Patient (10,00)',1,def);

if isempty(value), return; end

for i=1:nk
    k=patient(i).k;
    n=patient(i).n;
    arg = sprintf(sprintf('%%-%ds',n),value{i});
    hdr(k+8:k+7+n) = arg;
end

fp = fopen(fp,'w');
fwrite(fp,hdr,'uint8');
fclose(fp);
