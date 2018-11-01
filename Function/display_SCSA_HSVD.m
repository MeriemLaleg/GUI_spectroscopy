function[]= display_SCSA_HSVD()

% This function is used for pre-processing the data: Phase correction, Apodisation

% Author: Sourav Bhaduri
% University of Ghent
% email: Sourav.Bhaduri@UGent.be
% Feb 2018; Last revision: Feb 2018


global data_value;
global display_handles;
global details_para;
global apod_para;
global pha_para;
global ws_para;
global ds_para;
global ind_allvox;
global pha_cor_data;
global default_parameters;
global zerofill_value;
global info_struct;
global cur_vox;
global preprocess_display;
global FID_FD_copy;
global met_tick;
global mainGUI_handles;
global hsvd_para;
global scsa_para;


Pix_SS = get(0,'screensize');
x_max = Pix_SS(:,3)./2;
y_max = Pix_SS(:,4).*.86;
selected_vox = details_para.selected_voxels;
vox_no = cell(1,length(selected_vox));
for i = 1:length(selected_vox)
vox_no{i} =  num2str(selected_vox(1,i));
end

% Cerate GUI
display_handles.fig = figure('Visible','off','Position',[(x_max+2),35,x_max,y_max./1],'Color',[0.68,0.92,.8],'Name','Display window','NumberTitle', 'off');
% set(display_handles.fig,'MenuBar','none'); 

display_handles.axis = axes('Units','Pixels','Units','normalized','Position',[50/x_max,60/y_max,700/x_max,620/y_max]);
% display_handles.axis_text = uicontrol(display_handles.fig,'Style','text','Units','normalized','Position',[220/x_max 460/y_max 130/x_max 15/y_max],...
%     'String','Signal Spectrum','BackgroundColor',[0.68,0.92,1],'FontSize',9,'FontWeight','Bold');


display_handles.file_panel = uipanel(display_handles.fig,'Title','','Position',[.3 .93 .43 .07],'BackgroundColor','yellow','FontSize',5,'FontWeight','Bold');
display_handles.main_text = uicontrol(display_handles.file_panel,'Style','text','Units','normalized','Position',[0.005,0.25,0.1,0.1*2.5],'String',...
    'File name','FontSize',7,'FontWeight','Bold','BackgroundColor','yellow');
display_handles.main_val = uicontrol(display_handles.file_panel,'Style','edit','Units','normalized','Position',[0.15,0.15,0.95,0.35*2],'String',data_value.file_name{1},...
    'Callback',{@file_main_Callback},'Interruptible','off','BusyAction','cancel'); 

display_handles.close = uicontrol(display_handles.fig,'Style', 'pushbutton','Units','normalized','Position',[660/x_max 10/y_max 100/x_max 25/y_max],...
'String', 'OK/Close','FontSize',9,'FontWeight','Bold', 'Interruptible','off','BusyAction','cancel','Callback',@close_Callback);
% display_handles.met_map = uicontrol(display_handles.fig,'Style', 'pushbutton','Units','normalized','Position',[560/x_max 10/y_max 100/x_max 25/y_max],...
% 'String', 'metabolite map','FontSize',9,'FontWeight','Bold', 'Interruptible','off','BusyAction','cancel','Callback',@metabolite_map_Callback);

set(display_handles.fig,'Visible','on');
set(display_handles.main_val, 'Enable', 'off')
set(display_handles.fig,'CloseRequestFcn',@close_Callback);


 
preprocess_display = 1;

% Display the selected voxel 

function file_main_Callback(hObject, eventdata)
%Callback for 'first-order' edit box 
global pha_para;
global cur_vox;
global ind_allvox;
global display_handles;
global details_para;
 

function close_Callback(hObject, eventdata)
%Callback for 'close' button
global display_handles;
global met_tick;
global mainGUI_handles;
global details_para;
global preprocess_handles;
global fitting_handles;

if details_para.preprocess_going==1;
    details_para.preprocess_going=0;
    try
        close(preprocess_handles.fig);
    catch 
        return;
    end
end
if details_para.fitting_going==1;
    details_para.fitting_going=0;
    try
        close(fitting_handles.fig);
    catch 
        return;
    end
end
if (met_tick==1 && ~isempty(display_handles.met_fig))
    close(display_handles.fig)
    close(display_handles.met_fig);
else
     try
        close(display_handles.fig);
     catch 
        return;
     end  
  
end
%-----------------------------------------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------------------------------------
