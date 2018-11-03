
function GUI_start()
addpath ./Function
addpath ./Function/SCSA
% This function is used for starting the water suppression GUI

% Author: Sourav Bhaduri
% University of Ghent
% email: Sourav.Bhaduri@UGent.be
% Feb 2018; Last revision: Feb 2018

clear all
clear global all
close all
clc

global mainGUI_handles;
global details_para;
global globalpara;
global data_value;
global show_plot;show_plot=0;

globalpara=[];

Pix_SS = get(0,'screensize');
x_max = Pix_SS(:,3)./2;
y_max = 90;

data_value.scsa_flag=0;
data_value.sim_flag=0;
details_para.slice_no=1;
details_para.image_flag=0;
data_value.save_flag=0;
data_value.sense_flag=0;
details_para.fit_flag=0;

% main figure creation
mainGUI_handles.fig = figure('Visible','off','Position',[0,(Pix_SS(:,4)-120),x_max,y_max],'Color',[0.68,0.82,.68],'Name','Water Suppression 1.0','NumberTitle', 'off');
set(mainGUI_handles.fig,'MenuBar','none'); 

% menu set-up
% set(mainGUI_handles.fig,'MenuBar','none'); 
% mainGUI_handles.project = uimenu(mainGUI_handles.fig,'Label','Project');
% mainGUI_handles.new_project = uimenu(mainGUI_handles.project,'Label','New Project','Callback',@new_project_callback);
% mainGUI_handles.load_project = uimenu(mainGUI_handles.project,'Label','Load Project','Callback',@load_project_callback);
% mainGUI_handles.save_project = uimenu(mainGUI_handles.project,'Label','Save Project','Callback',@save_project_callback);
% mainGUI_handles.save_figure = uimenu(mainGUI_handles.project,'Label','Save Figure','Callback',@save_main_fig_callback);
% 
% mainGUI_handles.load = uimenu(mainGUI_handles.fig,'Label','Load');
% mainGUI_handles.load_data = uimenu(mainGUI_handles.load,'Label','MRS Data','Callback',@load_data_callback);
% mainGUI_handles.load_mri = uimenu(mainGUI_handles.load,'Label','MRI Image','Callback',@MRI_load_callback);
% 
% mainGUI_handles.settings = uimenu(mainGUI_handles.fig,'Label','Settings');
% mainGUI_handles.save_setting = uimenu(mainGUI_handles.settings,'Label','Create and Save Settings','Callback',@save_settings_callback);
% mainGUI_handles.load_setting = uimenu(mainGUI_handles.settings,'Label','Load Settings');%,'Callback',@load_setting_callback);
% mainGUI_handles.load_default_setting = uimenu(mainGUI_handles.load_setting,'Label','Default','Callback',@load_default_setting_callback);
% mainGUI_handles.load_user_setting = uimenu(mainGUI_handles.load_setting,'Label','User Defined','Callback',@load_user_setting_callback);
% mainGUI_handles.Method = uimenu(mainGUI_handles.fig,'Label','Method');
% mainGUI_handles.automatic = uimenu(mainGUI_handles.Method ,'Label','Automatic','Callback',@auto_callback);
% mainGUI_handles.manual= uimenu(mainGUI_handles.Method ,'Label','Manual','Callback',@manual_callback);
% 
% mainGUI_handles.mathematics = uimenu(mainGUI_handles.fig,'Label','Mathematics');
% mainGUI_handles.save_setting = uimenu(mainGUI_handles.mathematics,'Label','Sum','Callback',@sum_callback);
% mainGUI_handles.load_setting = uimenu(mainGUI_handles.mathematics,'Label','Average','Callback',@average_callback);
% mainGUI_handles.load_default_setting = uimenu(mainGUI_handles.mathematics,'Label','Subtract','Callback',@subtract_callback);
% 
% mainGUI_handles.segmentation = uimenu(mainGUI_handles.fig,'Label','Segmentation');
% mainGUI_handles.segment = uimenu(mainGUI_handles.segmentation,'Label','Segment','Callback',@segment_callback);
% mainGUI_handles.segment_display = uimenu(mainGUI_handles.segmentation,'Label','Display','Callback',@segment_display_callback);
% mainGUI_handles.segment_bias = uimenu(mainGUI_handles.segmentation,'Label','Bias correct','Callback',@segment_bias_callback);


%axis and voxel list setup
mainGUI_handles.slice_selection_panel = uipanel('Title','','Units','normalized','Position',[.04 .05 .9 .8],'BackgroundColor',[0.68,0.72,1],...
    'FontSize',9,'FontWeight','Bold','TitlePosition','centertop');
% mainGUI_handles.slice_selection_list = uicontrol('Parent',mainGUI_handles.slice_selection_panel,'Style','listbox','Units','normalized','Position',...
%     [.05 .027 .9 .97],'Callback', @slice_selection_callback,'Interruptible','off','BusyAction','cancel');
mainGUI_handles.slice_load = uicontrol(mainGUI_handles.slice_selection_panel,'Style', 'pushbutton','Units','normalized','Position', [.02 .3 .2 .5],...
'String', 'Load','FontSize',11,'FontWeight','Bold', 'Interruptible','off','BusyAction','cancel','Callback',@slice_load_Callback);
mainGUI_handles.slice_preprocess = uicontrol(mainGUI_handles.slice_selection_panel,'Style', 'pushbutton','Units','normalized','Position', [.25 .3 .2 .5],...
'String', 'Preprocess','FontSize',11,'FontWeight','Bold', 'Interruptible','off','BusyAction','cancel','Callback',@slice_preprocess_Callback);
mainGUI_handles.slice_fit = uicontrol(mainGUI_handles.slice_selection_panel,'Style', 'pushbutton','Units','normalized','Position', [.48 .3 .2 .5],...
'String', 'Fit','FontSize',9,'FontWeight','Bold', 'Interruptible','off','BusyAction','cancel','Callback',@slice_fit_Callback);
% mainGUI_handles.slice_map = uicontrol(mainGUI_handles.slice_selection_panel,'Style', 'pushbutton','Units','normalized','Position', [.09 .05 .8 .2],...
% 'String', 'Metabolite map','FontSize',9,'FontWeight','Bold', 'Interruptible','off','BusyAction','cancel','Callback',@slice_map_Callback);

mainGUI_handles.close = uicontrol(mainGUI_handles.fig,'Style', 'pushbutton','Units','normalized','Position',[670/x_max 5/y_max 45/x_max 25/y_max],...
'String', 'Close','FontSize',9,'FontWeight','Bold', 'Interruptible','off','BusyAction','cancel','Callback',@close_Callback);
mainGUI_handles.Help = uicontrol(mainGUI_handles.fig,'Style', 'pushbutton','Units','normalized','Position',[550/x_max 5/y_max 45/x_max 25/y_max],...
'String', 'Help','FontSize',9,'FontWeight','Bold', 'Interruptible','off','BusyAction','cancel','Callback',@help_Callback);
mainGUI_handles.save = uicontrol(mainGUI_handles.fig,'Style', 'pushbutton','Units','normalized','Position',[610/x_max 5/y_max 45/x_max 25/y_max],...
'String', 'Save','FontSize',9,'FontWeight','Bold', 'Interruptible','off','BusyAction','cancel','Callback',@save_Callback);

set(mainGUI_handles.fig,'Visible','on');
set(mainGUI_handles.slice_load, 'Enable', 'on');
set(mainGUI_handles.slice_preprocess, 'Enable', 'off');
set(mainGUI_handles.slice_fit, 'Enable', 'off');
% set(mainGUI_handles.slice_map, 'Enable', 'off');

set(mainGUI_handles.save, 'Enable', 'off');
set(mainGUI_handles.fig,'CloseRequestFcn',@close_Callback);
% undecorateFig(mainGUI_handles.fig);

waitfor(mainGUI_handles.fig);

%------------------------------------------------------------------------------------------------------------------------------------------------------------------------
function slice_load_Callback(hObject, eventdata)

global details_para;
global mainGUI_handles;
global data_value;
global preprocess_handles;
global preprocess_sel_vox;
global FID_FD_copy;
global met_tick;
global fitting_handles;
global preprocess_display;
global fitting_display;
global cur_vox;
global ind_allvox;
global zerofill_value;
global default_parameters;
global pha_cor_data;
global fit_para;
global apod_para;
global pha_para;
global ws_para;
global ds_para;
global met_para;
global model_para;
global cur_data;

% set(mainGUI_handles.slice_preprocess, 'Enable', 'on');
set(mainGUI_handles.slice_load, 'Enable', 'off');
set(mainGUI_handles.slice_fit, 'Enable', 'off');
% set(mainGUI_handles.slice_map, 'Enable', 'off');

try
    if (data_value.save_flag==0)
        if (data_value.sense_flag==0)
            info = 'Please select binary file to read';
            [fname,pathname]=uigetfile([pwd,'\Input_data\*.mat'],info);
            file_exten_ind=strfind(fname,'.');
            file_exten=fname(file_exten_ind:end);
            if (strcmp(file_exten,'.mat'))
                data_value.file_flag=1;
                data_value.temp_filename=[pathname fname];
                data_value.spect_name=fname;
            end
            try
                data_load_suppress_init();
            catch
                return;
            end
            if details_para.error==1
                return;
            end
        else
            rapid_sense_gui();
        end
    end
catch
    set(mainGUI_handles.slice_preprocess, 'Enable', 'off');
    return;
end
try
if (details_para.error==0)
    set(mainGUI_handles.slice_preprocess, 'Enable', 'on');
else
    set(mainGUI_handles.slice_preprocess, 'Enable', 'off');
    set(mainGUI_handles.slice_fit, 'Enable', 'off');
    set(mainGUI_handles.slice_map, 'Enable', 'off');
end
if (data_value.save_flag==1)
    
        load(data_value.temp_filename,'details_para','data_value','preprocess_handles','FID_FD_copy','met_tick','preprocess_sel_vox','fitting_handles','apod_para','fit_para','ind_allvox','pha_cor_data','default_parameters','zerofill_value','cur_vox','fitting_display','preprocess_display','pha_para','ws_para','ds_para','met_para','model_para','cur_data');
        set(mainGUI_handles.slice_preprocess, 'Enable', 'on');
end
set(mainGUI_handles.slice_load, 'Enable', 'on');
catch
    return;
end


%------------------------------------------------------------------------------------------------------------------------------------------------------

function slice_preprocess_Callback(hObject, eventdata)

global details_para;
global mainGUI_handles;
global data_value;
global preprocess_handles;

% set(mainGUI_handles.fig,'Visible','off');
set(mainGUI_handles.slice_fit, 'Enable', 'off');
% set(mainGUI_handles.slice_map, 'Enable', 'off');
set(mainGUI_handles.slice_preprocess, 'Enable', 'off');
set(mainGUI_handles.slice_load, 'Enable', 'off');

if (isempty(preprocess_handles))
     data_load_suppress();
else
    set(mainGUI_handles.slice_preprocess, 'Enable', 'off');
     data_load_suppress();
end
try
    set(mainGUI_handles.slice_preprocess, 'Enable', 'on');
catch
    return;
end   
if details_para.preprocess_done==1;
    set(mainGUI_handles.slice_fit, 'Enable', 'on');
end
% set(mainGUI_handles.fig,'Visible','on');
% set(mainGUI_handles.load_image, 'Enable', 'on');
set(mainGUI_handles.slice_load, 'Enable', 'on');

%------------------------------------------------------------------------------------------------------------------------------------------------------

function slice_fit_Callback(hObject, eventdata)

global details_para;
global mainGUI_handles;
global data_value;
global fitting_handles;

% set(mainGUI_handles.fig,'Visible','off');
set(mainGUI_handles.slice_fit, 'Enable', 'off');
set(mainGUI_handles.slice_preprocess, 'Enable', 'off');
% set(mainGUI_handles.slice_map, 'Enable', 'off');
set(mainGUI_handles.slice_load, 'Enable', 'off');

if (isempty(fitting_handles))
    Quantify_HSVD_SCSA();
else
     set(mainGUI_handles.slice_fit, 'Enable', 'off');
     Quantify_HSVD_SCSA();
end
try
    set(mainGUI_handles.slice_fit, 'Enable', 'on');
catch
    return;
end
set(mainGUI_handles.slice_preprocess, 'Enable', 'on');
% set(mainGUI_handles.fig,'Visible','on');
% set(mainGUI_handles.slice_map, 'Enable', 'on');
if details_para.fit_flag==1;
    set(mainGUI_handles.save, 'Enable', 'on');
else
    set(mainGUI_handles.save, 'Enable', 'off');
end
set(mainGUI_handles.slice_load, 'Enable', 'on');
%------------------------------------------------------------------------------------------------------------------------------------------------------




function close_Callback(hObject, eventdata)
%Callback for 'close' button

try
    close all;
catch 
    return;
end



function help_Callback(hObject, eventdata)
%Callback for 'close' button

file=strcat(pwd,'\README.docx');

try
    open(file);
catch 
    return;
end

function save_Callback(hObject, eventdata)
%Callback for 'close' button
global details_para;
global mainGUI_handles;
global data_value;
global preprocess_handles;
global preprocess_sel_vox;
global FID_FD_copy;
global met_tick;
global fitting_handles;
global preprocess_display;
global fitting_display;
global cur_vox;
global ind_allvox;
global zerofill_value;
global default_parameters;
global pha_cor_data;
global fit_para;
global apod_para;
global pha_para;
global ws_para;
global ds_para;
global met_para;
global model_para;
global cur_data;

set(mainGUI_handles.save, 'Enable', 'off');

file_exten_ind=strfind(data_value.spect_name,'.');
file_exten=data_value.spect_name(1:1:file_exten_ind-1);
file=strcat(pwd,'\Output_data','\',file_exten,'_output','.mat');
file_xl=strcat(pwd,'\Output_data','\',file_exten,'_output','.csv');
data_value.save_flag=1;
cur_data(:,1)=data_value.peak_area;
cur_data(:,2)=data_value.snr;
cur_data(:,3)=data_value.NAA_Cr;
cur_data(:,4)=data_value.Cho_Cr;
cur_data(:,5)=details_para.se;
col_name = {'Area','SNR','NAA/Cr','Cho/Cr','Error'};
row_name = cell(1,details_para.num_vox);
for i=1:details_para.num_vox
    row_name{i} = strcat('voxel ', num2str(i));
end
h = waitbar(0,'Please wait...');
try
    save(file,'details_para','mainGUI_handles','data_value','preprocess_handles','FID_FD_copy','met_tick','preprocess_sel_vox','fitting_handles','apod_para','fit_para','ind_allvox','pha_cor_data','default_parameters','zerofill_value','cur_vox','fitting_display','preprocess_display','pha_para','ws_para','ds_para','met_para','model_para','cur_data');
    dlmwrite(file_xl,cur_data,'precision', 9);
catch 
    set(mainGUI_handles.save, 'Enable', 'on');
    return;
end
waitbar(1);
close(h);
set(mainGUI_handles.save, 'Enable', 'on');




%--------------------------------------------------------------------------------------------------------------------------------------------------------