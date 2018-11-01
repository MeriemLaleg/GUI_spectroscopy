function[]= Quantify_HSVD_SCSA()

% This function is used for fitting the data: HLSVD based 

% Author: Sourav Bhaduri
% University of Ghent
% email: Sourav.Bhaduri@UGent.be
% Feb 2018; Last revision: Feb 2018

global data_value;
global fitting_handles;
global details_para;
global apod_para;
global fit_para;
global ind_allvox_fit;
global pha_cor_data;
global default_parameters;
global zerofill_value;
global info_struct;
global cur_vox;
global fitting_display;
global FID_FD_copy;
global met_tick;
global model_para;
global met_para;
global mainGUI_handles;


FID_FD_copy = data_value.FID_FD;
details_para.fitting_going=1;

if details_para.fit_done==0
    default_parameters.pha_0=1.9;
    default_parameters.pha_1=2.1;
    default_parameters.pha_0_NAA=1.9;
    default_parameters.pha_1_NAA=2.3;
    default_parameters.pha_0_Cr=2.8;
    default_parameters.pha_1_Cr=3.2;
    default_parameters.pha_0_Cho=3.3;
    default_parameters.pha_1_Cho=3.5;
    default_parameters.pha_0_Lac=1.1;
    default_parameters.pha_1_Lac=1.4;
    default_parameters.pha_0_mi=3.4;
    default_parameters.pha_1_mi=3.65;
    default_parameters.apod_para_lor=0;
    default_parameters.apod_para_gaus=0;
    default_parameters.apod_para_sigm=0;
    default_parameters.apod_para_type=1;
    ind_allvox_fit=1;
    met_tick=0;
    info_struct = [];
    details_para.selected_voxels=1:size(data_value.FID_FD,2);
    details_para.num_vox=size(data_value.FID_FD,2);
    details_para.N=size(data_value.FID_FD,1);
    cur_vox = details_para.selected_voxels(1,1);
    details_para.disp_ppm_seg=[-2+2.01 4+2.01];
    details_para.disp_freq_seg(1) = ((details_para.disp_ppm_seg(1) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
    details_para.disp_seg(1) = ceil((details_para.disp_freq_seg(1)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
    details_para.disp_freq_seg(2) = ((details_para.disp_ppm_seg(2) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
    details_para.disp_seg(2) = ceil((details_para.disp_freq_seg(2)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
    details_para.disp_fit(1)=default_parameters.pha_0;
    details_para.disp_fit(2)=default_parameters.pha_1;
    details_para.fit_method=0;
    for i = 1:length(details_para.selected_voxels)
        vox = details_para.selected_voxels(1,i);
        apod_para.type=1*ones(details_para.num_vox,1);
        fit_para.pha_0 = default_parameters.pha_0*ones(details_para.num_vox,1);
        fit_para.pha_1 = default_parameters.pha_1*ones(details_para.num_vox,1);
        fit_para.pha_0_NAA = default_parameters.pha_0_NAA*ones(details_para.num_vox,1);
        fit_para.pha_1_NAA = default_parameters.pha_1_NAA*ones(details_para.num_vox,1);
        fit_para.pha_0_Cr = default_parameters.pha_0_Cr*ones(details_para.num_vox,1);
        fit_para.pha_1_Cr = default_parameters.pha_1_Cr*ones(details_para.num_vox,1);
        fit_para.pha_0_Cho = default_parameters.pha_0_Cho*ones(details_para.num_vox,1);
        fit_para.pha_1_Cho = default_parameters.pha_1_Cho*ones(details_para.num_vox,1);
        fit_para.pha_0_Lac = default_parameters.pha_0_Lac*ones(details_para.num_vox,1);
        fit_para.pha_1_Lac = default_parameters.pha_1_Lac*ones(details_para.num_vox,1);
        fit_para.pha_0_mi = default_parameters.pha_0_mi*ones(details_para.num_vox,1);
        fit_para.pha_1_mi = default_parameters.pha_1_mi*ones(details_para.num_vox,1);
        fit_para.order=15*ones(details_para.num_vox,1);
        fit_para.max_fun=50*ones(details_para.num_vox,1);
        fit_para.max_iter=50*ones(details_para.num_vox,1);
        model_para.type=1*ones(details_para.num_vox,1);
        met_para.type=1*ones(details_para.num_vox,1);
    end
    details_para.fit_freq_NAA(1) = ((fit_para.pha_0_NAA(cur_vox) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
    details_para.fit_freq_NAA(2) = ((fit_para.pha_1_NAA(cur_vox) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
    details_para.fit_freq_Cr(1) = ((fit_para.pha_0_Cr(cur_vox) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
    details_para.fit_freq_Cr(2) = ((fit_para.pha_1_Cr(cur_vox) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
    details_para.fit_freq_Cho(1) = ((fit_para.pha_0_Cho(cur_vox) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
    details_para.fit_freq_Cho(2) = ((fit_para.pha_1_Cho(cur_vox) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
    details_para.fit_freq_Lac(1) = ((fit_para.pha_0_Lac(cur_vox) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
    details_para.fit_freq_Lac(2) = ((fit_para.pha_1_Lac(cur_vox) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
    details_para.fit_freq_mi(1) = ((fit_para.pha_0_mi(cur_vox) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
    details_para.fit_freq_mi(2) = ((fit_para.pha_1_mi(cur_vox) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
    details_para.fit_disp_seg_NAA(1) = floor((details_para.fit_freq_NAA(1)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
    details_para.fit_disp_seg_NAA(2) = ceil((details_para.fit_freq_NAA(2)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
    details_para.fit_disp_seg_Cr(1) = floor((details_para.fit_freq_Cr(1)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
    details_para.fit_disp_seg_Cr(2) = ceil((details_para.fit_freq_Cr(2)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
    details_para.fit_disp_seg_Cho(1) = floor((details_para.fit_freq_Cho(1)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
    details_para.fit_disp_seg_Cho(2) = ceil((details_para.fit_freq_Cho(2)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
    details_para.fit_disp_seg_Lac(1) = floor((details_para.fit_freq_Lac(1)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
    details_para.fit_disp_seg_Lac(2) = ceil((details_para.fit_freq_Lac(2)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
    details_para.fit_disp_seg_mi(1) = floor((details_para.fit_freq_mi(1)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
    details_para.fit_disp_seg_mi(2) = ceil((details_para.fit_freq_mi(2)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
end
Pix_SS = get(0,'screensize');
x_max = Pix_SS(:,3)./2;
y_max = Pix_SS(:,4).*.78;
selected_vox = details_para.selected_voxels;
vox_no = cell(1,length(selected_vox));
for i = 1:length(selected_vox)
vox_no{i} =  num2str(selected_vox(1,i));
end


% Cerate GUI
fitting_handles.fig = figure('Visible','off','Position',[0,35,x_max,y_max./1],'Color',[0.68,0.92,.8],'Name','Fitting window','NumberTitle', 'off');
set(fitting_handles.fig,'MenuBar','none'); 
% undecorateFig(fitting_handles.fig);
% fitting_handles.save_figure = uimenu(fitting_handles.fig,'Label','Save Figure','Callback',@save_preprocess_fig_callback);

fitting_handles.axis1 = axes('Units','Pixels','Units','normalized','Position',[120/x_max 2/y_max 140/x_max 80/y_max]);
% display_handles.axis_text = uicontrol(display_handles.fig,'Style','text','Units','normalized','Position',[220/x_max 460/y_max 130/x_max 15/y_max],...
%     'String','Signal Spectrum','BackgroundColor',[0.68,0.92,1],'FontSize',9,'FontWeight','Bold');
fitting_handles.axis2 = axes('Units','Pixels','Units','normalized','Position',[250/x_max 2/y_max 150/x_max 80/y_max]);

myImage=imread([pwd,'\Function\ugent_logo.png']);
myImage= imresize(myImage,1);
imshow(myImage,'Parent',fitting_handles.axis1);
myImage=imread([pwd,'\Function\kaust_logo.jpg']);
myImage= imresize(myImage,1);
imshow(myImage,'Parent',fitting_handles.axis2);


fitting_handles.h1 = uipanel(fitting_handles.fig,'Title','Voxel No','Units','normalized','Position',[0.91 .58 .08 .4],...
    'BackgroundColor','yellow','FontSize',8, 'FontWeight','Bold');
% fitting_handles.voxel_list = uicontrol(fitting_handles.h1,'Style', 'listbox','Units','normalized','Position', [.05 .027 .9 .97],'String', vox_no,...
%     'Max', details_para.num_vox, 'Min', 1, 'Callback', @voxel_list_callback, 'Interruptible','off','BusyAction','cancel'); 
fitting_handles.voxel_list = uicontrol(fitting_handles.h1,'Style', 'listbox','Units','normalized','Position', [.05 .027 .9 .97],'String', vox_no,...
    'value',1, 'Callback', @voxel_list_callback, 'Interruptible','off','BusyAction','cancel'); 


fitting_handles.display_panel = uipanel(fitting_handles.fig,'Title','Display Properties','Position',[.03 .58 .87 .4],'BackgroundColor','yellow','FontSize',9,'FontWeight','Bold');
fitting_handles.seg_min_text = uicontrol(fitting_handles.display_panel,'Style','text','Units','normalized','Position',[0.4,0.4,0.2,0.1*2.5],'String',...
    'Min','FontSize',8,'FontWeight','Bold','BackgroundColor','yellow');
fitting_handles.seg_max_text = uicontrol(fitting_handles.display_panel,'Style','text','Units','normalized','Position',[0.7,0.4,0.2,0.1*2.5],'String',...
    'Max','FontSize',8,'FontWeight','Bold','BackgroundColor','yellow');
fitting_handles.display_seg_mim = uicontrol(fitting_handles.display_panel,'Style','edit','Units','normalized','Position',[0.45,0.2,0.15,0.15*2],'String','0',...
    'Callback',{@display_seg_mim_Callback},'Interruptible','off','BusyAction','cancel'); 
fitting_handles.display_seg_max = uicontrol(fitting_handles.display_panel,'Style','edit','Units','normalized','Position',[0.75,0.2,0.15,0.15*2],'String','0',...
    'Callback',{@display_seg_max_Callback},'Interruptible','off','BusyAction','cancel'); 
% preprocess_handles.zerofill_text = uicontrol(preprocess_handles.fig,'Style','text','Units','normalized','Position',[557/x_max 460/y_max 70/x_max 15/y_max],...
%     'String','Zerofill','BackgroundColor',[0.68,0.92,1],'FontSize',9,'FontWeight','Bold');
% preprocess_handles.zerofill_options = uicontrol(preprocess_handles.fig,'Style','popup','String','No Fill|1 times|2 times|3 times','Units','normalized',...
%     'Position',[555/x_max 405/y_max 80/x_max 15/y_max],'Value',zerofill_value,'Interruptible','off','BusyAction','cancel','Callback', @zerofill_options_callback); 

fitting_handles.TF_FD_display = uibuttongroup(fitting_handles.display_panel,'Units','normalized','Position',[0.04 0.76 .7 0.16],...
    'FontSize',9,'BackgroundColor',[0.68,0.92,.8],'BorderType','none','SelectionChangeFcn',@TD_FD_display_callback);
fitting_handles.FD_disp_real = uicontrol(fitting_handles.TF_FD_display,'Style','radiobutton','Units','normalized','Position',[.01 .1 .3 .9],'String',...
    'FD signal real','FontWeight','Bold','BackgroundColor',[0.68,0.92,.8]);
fitting_handles.FD_disp_abs = uicontrol(fitting_handles.TF_FD_display,'Style','radiobutton','Units','normalized','Position',[.35 .1 .3 .9],'String',...
    'FD signal abs','FontWeight','Bold','BackgroundColor',[0.68,0.92,.8]);
fitting_handles.TD_disp = uicontrol(fitting_handles.TF_FD_display,'Style','radiobutton','Units','normalized','String','TD signal','Position',...
    [.71 .1 .3 .9],'FontWeight','Bold','BackgroundColor',[0.68,0.92,.8]);
% apodization panel
fitting_handles.h2 = uipanel(fitting_handles.fig,'Title','Metabolite','Units','normalized','Position',[20/x_max 305/y_max 170/x_max 80/y_max],...
    'BackgroundColor','yellow','FontSize',9, 'FontWeight','Bold');
fitting_handles.metabolite = uicontrol(fitting_handles.h2,'Style','text','Units','normalized','String','Metabolite Type','Units','normalized',...
    'Position',[0.1 0.65 0.8 0.3],'BackgroundColor','yellow','FontSize',8,'FontWeight','Bold');
fitting_handles.metabolite_type = uicontrol(fitting_handles.h2,'Style','popup','String','NAA|Lactate|Creatine|Choline|Myo','Units','normalized','Position',...
    [0.1 0.19 0.8 0.3],'Value',met_para.type(cur_vox), 'Callback',{@metabolite_type_Callback});

% model panel
fitting_handles.h_model = uipanel(fitting_handles.fig,'Title','','Units','normalized','Position',[190/x_max 305/y_max 170/x_max 80/y_max],...
    'BackgroundColor','yellow','FontSize',9, 'FontWeight','Bold');
fitting_handles.model_type = uicontrol(fitting_handles.h_model,'Style','text','Units','normalized','String','Fitting Model','Units','normalized',...
    'Position',[0.1 0.52 0.8 0.3],'BackgroundColor','yellow','FontSize',8,'FontWeight','Bold');
fitting_handles.model_select = uicontrol(fitting_handles.h_model,'Style','popup','String','Lorentz|Gauss|Voigt','Units','normalized','Position',...
    [0.1 0.19 0.8 0.3],'Value',model_para.type(cur_vox), 'Callback',{@model_select_Callback});


% phase correction panel
fitting_handles.h3 = uipanel(fitting_handles.fig,'Title','Fitting range','Units','normalized','Position',[20/x_max 145/y_max 620/x_max 150/y_max],...
    'BackgroundColor','yellow','FontSize',9, 'FontWeight','Bold');
fitting_handles.pha_corr0_text = uicontrol(fitting_handles.h3,'Style','text','Units','normalized','Position',[0.1 0.5 0.2 .2],'String',...
    'start ','FontSize',8,'FontWeight','Bold','BackgroundColor','yellow');
fitting_handles.pha_corr1_text = uicontrol(fitting_handles.h3,'Style','text','Units','normalized','Position',[0.1 0.2 0.2 0.2],'String',...
    'end','FontSize',8,'FontWeight','Bold','BackgroundColor','yellow');
fitting_handles.pha_corr_value_text = uicontrol(fitting_handles.h3,'Style','text','Units','normalized','Position',[0.31 0.71 0.2 .2],'String',...
    'Value (PPM)','FontSize',8,'FontWeight','Bold','BackgroundColor','yellow');
fitting_handles.pha_corr0_value = uicontrol(fitting_handles.h3,'Style','edit','Units','normalized','Position',[0.34 0.5 .15 .2],...
    'String',num2str(fit_para.pha_0(cur_vox)),'Callback',{@pha0_val_Callback});
fitting_handles.pha_corr1_value = uicontrol(fitting_handles.h3,'Style','edit','Units','normalized','Position',[0.34 0.2 .15 .2],...
    'String',num2str(fit_para.pha_1(cur_vox)),'Callback',{@pha1_val_Callback});
fitting_handles.order_text = uicontrol(fitting_handles.h3,'Style','text','Units','normalized','Position',[0.71 0.7 0.2 .2],'String',...
    'order: ','FontSize',8,'FontWeight','Bold','BackgroundColor','yellow');
fitting_handles.order_value = uicontrol(fitting_handles.h3,'Style','edit','Units','normalized','Position',[0.85 0.75 .1 .2],...
    'String',num2str(fit_para.order(cur_vox)),'Callback',{@order_val_Callback});

fitting_handles.max_fun_text = uicontrol(fitting_handles.h3,'Style','text','Units','normalized','Position',[0.75 0.4 0.1 .2],'String',...
    'max_fun: ','FontSize',8,'FontWeight','Bold','BackgroundColor','yellow');
fitting_handles.max_fun_value = uicontrol(fitting_handles.h3,'Style','edit','Units','normalized','Position',[0.85 0.45 .1 .2],...
    'String',num2str(fit_para.max_fun(cur_vox)),'Callback',{@max_fun_val_Callback});
fitting_handles.max_iter_text = uicontrol(fitting_handles.h3,'Style','text','Units','normalized','Position',[0.74 0.1 0.11 .2],'String',...
    'iterations: ','FontSize',8,'FontWeight','Bold','BackgroundColor','yellow');
fitting_handles.max_iter_value = uicontrol(fitting_handles.h3,'Style','edit','Units','normalized','Position',[0.85 0.15 .1 .2],...
    'String',num2str(fit_para.max_iter(cur_vox)),'Callback',{@max_iter_val_Callback});


fitting_handles.buttongroup = uibuttongroup(fitting_handles.fig,'Units','normalized','Position',[32/x_max .6 130/x_max .1],...
    'BackgroundColor','yellow','FontSize',9, 'FontWeight','Bold');%,'SelectionChangeFcn',@ind_allvox_fit_value);
fitting_handles.all_vox = uicontrol('parent',fitting_handles.buttongroup,'Style','radiobutton','String','All Voxel',...
    'Units','normalized','Position',[0.1 0.2 0.7 0.3],'BackgroundColor','yellow');
fitting_handles.ind_vox = uicontrol('parent',fitting_handles.buttongroup,'Style','radiobutton','String','Individual Voxel',...
    'Units','normalized','Position',[0.1 0.6 0.8 0.3],'BackgroundColor','yellow');

fitting_handles.developer_display = uipanel(fitting_handles.fig,'Title','Developers:','Units','normalized','Position',[20/x_max 85/y_max 620/x_max 60/y_max],...
    'BackgroundColor','yellow','FontSize',9, 'FontWeight','Bold');
fitting_handles.developer_text = uicontrol(fitting_handles.developer_display,'Style','text','Units','normalized','Position',[0.01 0.01 0.7 0.7],'String',...
    'Sourav Bhaduri, Abderrazak Chahid','FontSize',8,'FontWeight','Bold','BackgroundColor','yellow');


fitting_handles.close = uicontrol(fitting_handles.fig,'Style', 'pushbutton','Units','normalized','Position',[660/x_max 10/y_max 100/x_max 25/y_max],...
'String', 'OK/Close','FontSize',9,'FontWeight','Bold', 'Interruptible','off','BusyAction','cancel','Callback',@close_Callback);


fitting_handles.fit = uicontrol(fitting_handles.fig,'Style', 'pushbutton','Units','normalized','Position',[350/x_max 220/y_max 100/x_max 25/y_max],...
'String', 'Fit','FontSize',9,'FontWeight','Bold', 'Interruptible','off','BusyAction','cancel','Callback',@fitting_Callback);

% fitting_handles.wavelet = uicontrol(fitting_handles.fig,'Style', 'pushbutton','Units','normalized','Position',[240/x_max 140/y_max 100/x_max 25/y_max],...
% 'String', 'Denoise','FontSize',9,'FontWeight','Bold', 'Interruptible','off','BusyAction','cancel','Callback',@denoise_Callback);

set(fitting_handles.buttongroup,'SelectionChangeFcn',@ind_allvox_fit_value);



if ind_allvox_fit == 1
    set(fitting_handles.buttongroup,'SelectedObject',fitting_handles.ind_vox);
else
    set(fitting_handles.buttongroup,'SelectedObject',fitting_handles.all_vox);
end

set(fitting_handles.fig,'Visible','on');
set(fitting_handles.display_seg_mim, 'Enable', 'on','string',num2str(details_para.disp_ppm_seg(1)))
set(fitting_handles.display_seg_max, 'Enable', 'on','string',num2str(details_para.disp_ppm_seg(2)))
set(fitting_handles.pha_corr0_value, 'Enable', 'on','string',num2str(details_para.disp_fit(1)))
set(fitting_handles.pha_corr1_value, 'Enable', 'on','string',num2str(details_para.disp_fit(2)))


% set(fitting_handles.main_val, 'Enable', 'off')
% set(fitting_handles.ref_val, 'Enable', 'off')
set(fitting_handles.model_select,'value',1);
 
fitting_display = 1;

set(fitting_handles.fig,'CloseRequestFcn',@close_Callback);
% Display the selected voxel 
display_SCSA_HSVD();
display_data()
waitfor(fitting_handles.fig);

function voxel_list_callback(hObject, eventdata)
% Callback for 'voxel no.' get the selected voxel and update the plot according to the voxel selected
global details_para;
global fitting_handles;
global apod_para;
global fit_para;
global cur_vox;
global model_para;
global met_para;

% get the selected voxel
preprocess_sel_vox = get(hObject,'Value');

if(isempty(preprocess_sel_vox))
    preprocess_sel_vox = 1;
    set(fitting_handles.voxel_list,'value',1);
end

selected_vox = details_para.selected_voxels;
cur_vox = selected_vox(1,preprocess_sel_vox);
display_data()


% Update the GUI based on the selected voxel
temp = apod_para.type(cur_vox);
set(fitting_handles.metabolite_type,'value',temp);

set(fitting_handles.pha_corr0_value,'string',num2str(fit_para.pha_0(cur_vox)));
set(fitting_handles.pha_corr1_value,'string',num2str(fit_para.pha_1(cur_vox)));
set(fitting_handles.model_select,'value',model_para.type(cur_vox));
set(fitting_handles.metabolite_type,'value',met_para.type(cur_vox));
set(fitting_handles.order_value,'string',num2str(fit_para.order(cur_vox)));
set(fitting_handles.max_fun_value,'string',num2str(fit_para.max_fun(cur_vox)));
set(fitting_handles.max_iter_value,'string',num2str(fit_para.max_iter(cur_vox)));
set(fitting_handles.pha_corr0_value, 'Enable', 'on','string',num2str(fit_para.pha_0(cur_vox)));
set(fitting_handles.pha_corr1_value, 'Enable', 'on','string',num2str(fit_para.pha_1(cur_vox)));
% set(fitting_handles.display_seg_mim, 'Enable', 'on','string',num2str(details_para.disp_ppm_seg(1)))
% set(fitting_handles.display_seg_max, 'Enable', 'on','string',num2str(details_para.disp_ppm_seg(2)))



function metabolite_type_Callback(hObject, eventdata)
%Callback for 'Apodisation type' 
global fitting_handles;
global apod_para;
global ind_allvox_fit;
global details_para;
global cur_vox;
global fit_para;
global met_para;
% global data_value;

% get the selected apodization type
temp = get(hObject,'Value');
if(ind_allvox_fit == 1)
    met_para.type(cur_vox) = temp;
else
    met_para.type(details_para.selected_voxels(1,:)) = temp;
end

% enable/disable based on the selected type
switch(temp)
  case 1 %NAA
        phi_0 = 1.9;
         if(ind_allvox_fit == 1)
            fit_para.pha_0(cur_vox) = phi_0; 
         else
            fit_para.pha_0(details_para.selected_voxels(1,:)) = phi_0; 
         end
        phi_1 = 2.1;
         if(ind_allvox_fit == 1)
            fit_para.pha_1(cur_vox) = phi_1; 
         else
            fit_para.pha_1(details_para.selected_voxels(1,:)) = phi_1; 
         end  
        set(fitting_handles.pha_corr0_value, 'Enable', 'on','string',num2str(phi_0));
        set(fitting_handles.pha_corr1_value, 'Enable', 'on','string',num2str(phi_1));
  case 2 %Lactate
         phi_0 = 1.1;
          if(ind_allvox_fit == 1)
            fit_para.pha_0(cur_vox) = phi_0; 
         else
            fit_para.pha_0(details_para.selected_voxels(1,:)) = phi_0; 
         end   
         phi_1 = 1.4;
          if(ind_allvox_fit == 1)
            fit_para.pha_1(cur_vox) = phi_1; 
         else
            fit_para.pha_1(details_para.selected_voxels(1,:)) = phi_1; 
         end  
        set(fitting_handles.pha_corr0_value, 'Enable', 'on','string',num2str(phi_0));
        set(fitting_handles.pha_corr1_value, 'Enable', 'on','string',num2str(phi_1));
  case 3 %Creatine
         phi_0 = 2.8;
         if(ind_allvox_fit == 1)
            fit_para.pha_0(cur_vox) = phi_0; 
         else
            fit_para.pha_0(details_para.selected_voxels(1,:)) = phi_0; 
         end
         phi_1 = 3.1;
         if(ind_allvox_fit == 1)
            fit_para.pha_1(cur_vox) = phi_1; 
         else
            fit_para.pha_1(details_para.selected_voxels(1,:)) = phi_1;   
         end
        set(fitting_handles.pha_corr0_value, 'Enable', 'on','string',num2str(phi_0));
        set(fitting_handles.pha_corr1_value, 'Enable', 'on','string',num2str(phi_1));
  case 4 %Choline
         phi_0 = 3.3;
          if(ind_allvox_fit == 1)
            fit_para.pha_0(cur_vox) = phi_0; 
         else
            fit_para.pha_0(details_para.selected_voxels(1,:)) = phi_0; 
         end    
         phi_1 = 3.4;
          if(ind_allvox_fit == 1)
            fit_para.pha_1(cur_vox) = phi_1; 
         else
            fit_para.pha_1(details_para.selected_voxels(1,:)) = phi_1; 
         end  
         set(fitting_handles.pha_corr0_value, 'Enable', 'on','string',num2str(phi_0));
         set(fitting_handles.pha_corr1_value, 'Enable', 'on','string',num2str(phi_1));
  case 5 %Myo inositol
         phi_0 = 3.4;
          if(ind_allvox_fit == 1)
            fit_para.pha_0(cur_vox) = phi_0; 
         else
            fit_para.pha_0(details_para.selected_voxels(1,:)) = phi_0; 
         end   
         phi_1 = 3.65;
          if(ind_allvox_fit == 1)
            fit_para.pha_1(cur_vox) = phi_1; 
         else
            fit_para.pha_1(details_para.selected_voxels(1,:)) = phi_1; 
         end  
         set(fitting_handles.pha_corr0_value, 'Enable', 'on','string',num2str(phi_0));
         set(fitting_handles.pha_corr1_value, 'Enable', 'on','string',num2str(phi_1));
end
%------------------------------------------------------------------------------------------------------------
function model_select_Callback(hObject, eventdata)
%Callback for 'Apodisation type' 
global fitting_handles;
global apod_para;
global ind_allvox_fit;
global details_para;
global cur_vox;
global fit_para;
global model_para;
% global data_value;

% get the selected apodization type
temp = get(hObject,'Value');
if(ind_allvox_fit == 1)
    model_para.type(cur_vox) = temp;
else
    model_para.type(details_para.selected_voxels(1,:)) = temp;
end
%-----------------------------------------------------------------------------------------------------------
function ind_allvox_fit_value(source,eventdata)
% call back for individual voxel and all voxel button group
global ind_allvox_fit;

% select whether to apply pre-processing on all the voxels or in the selected voxel 
temp = get(eventdata.NewValue,'string');
if(strcmp(temp,'Individual Voxel'))
    ind_allvox_fit = 1;
else
    ind_allvox_fit = 2;
end
%-------------------------------------------------------------------------------------------------------- 
% 
function pha0_val_Callback(hObject, eventdata)
%Callback for 'zero-order' edit box 
global fit_para;
global cur_vox;
global ind_allvox_fit;
global fitting_handles;
global details_para;

% read and check the zero-order phase value from the edit box 
phi_0 = str2double(get(hObject,'string'));
if((phi_0 < -10)|| (phi_0 > 10) || (isnan(phi_0)))
    errordlg('Invalid value');
    set(hObject,'string',num2str(fit_para.pha_0(cur_vox)))
    return;
end
fit_para.pha_0(details_para.selected_voxels(1,:)) = phi_0;   
set(fitting_handles.pha_corr0_value,'string',num2str(fit_para.pha_0(cur_vox)));

% 
function pha1_val_Callback(hObject, eventdata)
%Callback for 'first-order' edit box 
global fit_para;
global cur_vox;
global ind_allvox_fit;
global fitting_handles;
global details_para;

% read and check the first-order phase value from the edit box 
phi_1 = str2double(get(hObject,'string'));
if((phi_1 < -10)|| (phi_1 > 10) || (isnan(phi_1)))
    errordlg('Invalid value');
    set(hObject,'string',num2str(fit_para.pha_1(cur_vox)))
    return;
end
fit_para.pha_1(details_para.selected_voxels(1,:)) = phi_1; 
set(fitting_handles.pha_corr1_value,'string',num2str(fit_para.pha_1(cur_vox)));

function order_val_Callback(hObject, eventdata)
%Callback for 'first-order' edit box 
global fit_para;
global cur_vox;
global ind_allvox_fit;
global fitting_handles;
global details_para;

% read and check the first-order phase value from the edit box 
ord_1 = str2double(get(hObject,'string'));
if((ord_1 < 0)|| (ord_1 > 100) || (isnan(ord_1)))
    errordlg('Invalid value');
    set(hObject,'string',num2str(fit_para.order(cur_vox)))
    return;
end
if(ind_allvox_fit == 1)
    fit_para.order(cur_vox) = ord_1; 
else
    fit_para.order(details_para.selected_voxels(1,:)) = ord_1; 
end
set(fitting_handles.order_value,'string',num2str(fit_para.order(cur_vox)));

function max_fun_val_Callback(hObject, eventdata)
%Callback for 'first-order' edit box 
global fit_para;
global cur_vox;
global ind_allvox_fit;
global fitting_handles;
global details_para;

% read and check the first-order phase value from the edit box 
ord_1 = str2double(get(hObject,'string'));
if((ord_1 < 0)|| (ord_1 > 10000) || (isnan(ord_1)))
    errordlg('Invalid value');
    set(hObject,'string',num2str(fit_para.max_fun(cur_vox)))
    return;
end
if(ind_allvox_fit == 1)
    fit_para.max_fun(cur_vox) = ord_1; 
else
    fit_para.max_fun(details_para.selected_voxels(1,:)) = ord_1; 
end
set(fitting_handles.max_fun_value,'string',num2str(fit_para.max_fun(cur_vox)));

function max_iter_val_Callback(hObject, eventdata)
%Callback for 'first-order' edit box 
global fit_para;
global cur_vox;
global ind_allvox_fit;
global fitting_handles;
global details_para;

% read and check the first-order phase value from the edit box 
ord_1 = str2double(get(hObject,'string'));
if((ord_1 < 0)|| (ord_1 > 10000) || (isnan(ord_1)))
    errordlg('Invalid value');
    set(hObject,'string',num2str(fit_para.max_iter(cur_vox)))
    return;
end
if(ind_allvox_fit == 1)
    fit_para.max_iter(cur_vox) = ord_1; 
else
    fit_para.max_iter(details_para.selected_voxels(1,:)) = ord_1; 
end
set(fitting_handles.max_iter_value,'string',num2str(fit_para.max_iter(cur_vox)));



function TD_FD_display_callback(hObject, eventdata)
global fitting_display;

temp = get(get(hObject,'SelectedObject'),'String');

if strcmp(temp,'FD signal real')
    fitting_display = 1;
elseif strcmp(temp,'FD signal abs')
    fitting_display = 2;
else
    fitting_display = 0;
end
display_data()

function display_data()
global fitting_display;
global fitting_handles;
global cur_vox;
global details_para;
global data_value;
global apod_para;
global cur_data;
global display_handles;

seg = details_para.disp_seg;
ppm = details_para.ppm_referenced;

if (details_para.fit_flag(:,cur_vox)==0)
    if (fitting_display == 1)
        cla(display_handles.axis);
        plot_var = real(data_value.FID_FD(seg(1):seg(2),cur_vox));
        plot(display_handles.axis,ppm(seg(1):seg(2))',plot_var);
        max_p = max(plot_var);
        min_p = min(plot_var);
        axis(display_handles.axis,[ppm(seg(1)), ppm(seg(2)),min_p - 0.01*max_p, max_p + 0.01*max_p])
        set(display_handles.axis,'XDir','reverse');
        legend(display_handles.axis,'Water suppressed and denoised');
        xlabel(display_handles.axis,'PPM');
        ylabel(display_handles.axis,'Amplitude');
    elseif (fitting_display == 2)
        cla(display_handles.axis);
        plot_var = abs(data_value.FID_FD(seg(1):seg(2),cur_vox));
        plot(display_handles.axis,ppm(seg(1):seg(2))',plot_var);
        max_p = max(plot_var);
        min_p = min(plot_var);
        axis(display_handles.axis,[ppm(seg(1)), ppm(seg(2)),min_p - 0.01*max_p, max_p + 0.01*max_p])
        set(display_handles.axis,'XDir','reverse');
        legend(display_handles.axis,'Water suppressed and denoised');
        xlabel(display_handles.axis,'PPM');
        ylabel(display_handles.axis,'Amplitude');
    else
        cla(display_handles.axis);
        plot_var = real(data_value.FID_TD(:,cur_vox));
    %     plot_var= ifftshift(real(data_value.FID_TD(:,cur_vox))+flipud(real(data_value.FID_TD(:,cur_vox))));
    %     filter = apod_para.function_weight(:,cur_vox);
        max_p = max(abs(plot_var));
        min_p = min(abs(plot_var));
    %     if (max(filter)> max_p)
    %         max_p = max(filter);
    %         min_p = -max_p;
    %     else
    %         max_p = max(plot_var);
    %         min_p = -max_p;
    %     end 
        set(display_handles.axis,'Xdir','Normal');
        axis(display_handles.axis,[details_para.t(1), details_para.t(end),min_p - 0.01*max_p, max_p + 0.01*max_p]);
        plot(display_handles.axis,details_para.t',plot_var);
        legend(display_handles.axis,'Water suppressed and denoised');
        xlabel(display_handles.axis,'time');
        ylabel(display_handles.axis,'Amplitude');
    %     hold on;
    %     plot(display_handles.axis,details_para.t',filter,'r');
    %     axis(display_handles.axis,[details_para.t(1), details_para.t(details_para.N),min_p , max_p])
    %     hold off
    end
else
    Pix_SS = get(0,'screensize');
    x_max = Pix_SS(:,3)./2;
    y_max = Pix_SS(:,4).*.75;
    col_name = {'Area','SNR','NAA/Cr','Cho/Cr','Error'};
    row_name = cell(1,details_para.num_vox);
    for i=1:details_para.num_vox
        row_name{i} = strcat('voxel ', num2str(i));
    end
    diplay_handles.para_table = uitable(fitting_handles.fig,'Units','normalized','position',[370/x_max,290/y_max,320/x_max,130/y_max], 'data',cur_data,...
    'RowName', row_name, 'ColumnName', col_name,'ColumnWidth',{85});
     if (fitting_display == 1)
        cla(display_handles.axis);
        plot_var = real(data_value.FID_FD(seg(1):seg(2),cur_vox));
        plot_fit = real(data_value.fit_FD(seg(1):seg(2),cur_vox));
        plot(display_handles.axis,ppm(seg(1):seg(2))',plot_var);
        hold (display_handles.axis,'on'); 
        plot(display_handles.axis,ppm(seg(1):seg(2))',plot_fit);
        max_p = max(plot_var);
        min_p = min(plot_var);
        axis(display_handles.axis,[ppm(seg(1)), ppm(seg(2)),min_p - 0.01*max_p, max_p + 0.01*max_p])
        set(display_handles.axis,'XDir','reverse');
        legend(display_handles.axis,'Water suppressed and denoised','Fitted+baseline');
        xlabel(display_handles.axis,'PPM');
        ylabel(display_handles.axis,'Amplitude');
    elseif (fitting_display == 2)
        cla(display_handles.axis);
        plot_var = abs(data_value.FID_FD(seg(1):seg(2),cur_vox));
        plot_fit = abs(data_value.fit_FD(seg(1):seg(2),cur_vox));
        plot(display_handles.axis,ppm(seg(1):seg(2))',plot_var);
        hold (display_handles.axis,'on'); 
        plot(display_handles.axis,ppm(seg(1):seg(2))',plot_fit);
        max_p = max(plot_var);
        min_p = min(plot_var);
        axis(display_handles.axis,[ppm(seg(1)), ppm(seg(2)),min_p - 0.01*max_p, max_p + 0.01*max_p])
        set(display_handles.axis,'XDir','reverse');
        legend(display_handles.axis,'Water suppressed and denoised','Fitted+baseline');
        xlabel(display_handles.axis,'PPM');
        ylabel(display_handles.axis,'Amplitude');
    else
        cla(display_handles.axis);
        plot_var = real(data_value.FID_TD(:,cur_vox));
        plot_fit = real(data_value.fit_TD(:,cur_vox));
    %     plot_var= ifftshift(real(data_value.FID_TD(:,cur_vox))+flipud(real(data_value.FID_TD(:,cur_vox))));
    %     filter = apod_para.function_weight(:,cur_vox);
        max_p = max(abs(plot_var));
        min_p = min(abs(plot_var));
    %     if (max(filter)> max_p)
    %         max_p = max(filter);
    %         min_p = -max_p;
    %     else
    %         max_p = max(plot_var);
    %         min_p = -max_p;
    %     end  
        set(display_handles.axis,'Xdir','Normal');
        axis(display_handles.axis,[details_para.t(1), details_para.t(end),min_p - 0.5*max_p, max_p + 0.1*max_p])
        plot(display_handles.axis,details_para.t',plot_var);
        hold (display_handles.axis,'on'); 
        plot(display_handles.axis,details_para.t',plot_fit);
        legend(display_handles.axis,'Water suppressed and denoised','Fitted+baseline');
        xlabel(display_handles.axis,'time');
        ylabel(display_handles.axis,'Amplitude');
    %     hold on;
    %     plot(display_handles.axis,details_para.t',filter,'r');
    %     axis(display_handles.axis,[details_para.t(1), details_para.t(details_para.N),min_p , max_p])
    %     hold off
    end
end

function save_preprocess_fig_callback(hObject, eventdata)
%Callback for 'save figure'
global fitting_handles;
global log_file;
ini_path = log_file.saving_path;
[temp,s_path] = uiputfile(ini_path);
if(temp == 0)
    return;
else
    dot_pos = strfind(temp, '.');
    s_fln = temp(1:dot_pos-1);
    file_name = strcat(s_path,s_fln);
    saveas(fitting_handles.fig,file_name,'png');
    log_file.saving_path = s_path;
end



function close_Callback(hObject, eventdata)
%Callback for 'close' button
global fitting_handles;
global met_tick;
global data_value;
global details_para;
global display_handles;

if (met_tick==1 && ~isempty(fitting_handles.met_fig))
    try
        close(display_handles.fig);
    catch 
        return;
    end
    close(fitting_handles.fig)
    close(fitting_handles.met_fig);
%     close(display_handles.fig);
else
    try
        close(display_handles.fig);
        close(fitting_handles.fig);
     catch 
        return;
     end
%      close(fitting_handles.fig);
%      close(display_handles.fig);
end
details_para.fit_done=1;
details_para.fitting_going=0;
    
  
%-----------------------------------------------------------------------------------------------------------------------------------
function fitting_Callback(hObject, eventdata)

global data_value;
global details_para;
global fitting_handles;
global preprocess_sel_vox;
global FID_FD_copy;
global fit_para;
global cur_vox;
global ind_allvox_fit;
global model_para;
global cur_data;


details_para.fit_freq(1) = ((fit_para.pha_0(cur_vox) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
details_para.fit_freq(2) = ((fit_para.pha_1(cur_vox) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
details_para.fit_freq_NAA(1) = ((fit_para.pha_0_NAA(cur_vox) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
details_para.fit_freq_NAA(2) = ((fit_para.pha_1_NAA(cur_vox) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
details_para.fit_freq_Cr(1) = ((fit_para.pha_0_Cr(cur_vox) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
details_para.fit_freq_Cr(2) = ((fit_para.pha_1_Cr(cur_vox) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
details_para.fit_freq_Cho(1) = ((fit_para.pha_0_Cho(cur_vox) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
details_para.fit_freq_Cho(2) = ((fit_para.pha_1_Cho(cur_vox) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
details_para.fit_freq_Lac(1) = ((fit_para.pha_0_Lac(cur_vox) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
details_para.fit_freq_Lac(2) = ((fit_para.pha_1_Lac(cur_vox) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
details_para.fit_freq_mi(1) = ((fit_para.pha_0_mi(cur_vox) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
details_para.fit_freq_mi(2) = ((fit_para.pha_1_mi(cur_vox) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
details_para.fit_disp_seg(1) = floor((details_para.fit_freq(1)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
details_para.fit_disp_seg(2) = ceil((details_para.fit_freq(2)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
details_para.fit_disp_seg_NAA(1) = floor((details_para.fit_freq_NAA(1)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
details_para.fit_disp_seg_NAA(2) = ceil((details_para.fit_freq_NAA(2)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
details_para.fit_disp_seg_Cr(1) = floor((details_para.fit_freq_Cr(1)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
details_para.fit_disp_seg_Cr(2) = ceil((details_para.fit_freq_Cr(2)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
details_para.fit_disp_seg_Cho(1) = floor((details_para.fit_freq_Cho(1)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
details_para.fit_disp_seg_Cho(2) = ceil((details_para.fit_freq_Cho(2)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
details_para.fit_disp_seg_Lac(1) = floor((details_para.fit_freq_Lac(1)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
details_para.fit_disp_seg_Lac(2) = ceil((details_para.fit_freq_Lac(2)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
details_para.fit_disp_seg_mi(1) = floor((details_para.fit_freq_mi(1)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
details_para.fit_disp_seg_mi(2) = ceil((details_para.fit_freq_mi(2)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples

% set(fitting_handles.fig,'Visible','off');
h = waitbar(0,'Performing fitting operation....');
cnt=1;
try
    if(ind_allvox_fit == 1)
        if details_para.fit_method==0
            
                cnt=cur_vox;

                re_bl = ssa(real(data_value.FID_FD(:,cur_vox)),40); % estimate the baseline for real part
        %         sp2 = csaps(ppm,imag(data_value.FID_FD(:,i)),30);
        %         im_bl = fnval(sp2,ppm);
                im_bl = ssa(imag(data_value.FID_FD(:,cur_vox)),40); % estimate the baseline for imaginary part
                data_value.FID_FD_base(:,cur_vox)=re_bl + 1i*im_bl;
                data_value.FID_FD_b(:,cur_vox)=data_value.FID_FD(:,cur_vox)-data_value.FID_FD_base(:,cur_vox);
                data_value.FID_TD_b(:,cur_vox)=ifft(ifftshift(data_value.FID_FD_b(:,cur_vox)));
                baseline(:,1)=data_value.FID_FD_base(:,cur_vox);
                baseline_t(:,1)=ifft(ifftshift(data_value.FID_FD_base(:,cur_vox)));
                [ fit_sig,amp,Ind_sig,freq,Z_sig,FWHM] = HSVDPK_new(data_value.FID_TD_b(:,cur_vox),details_para.fres,details_para.Fs,fit_para.order(cur_vox));
                Initx = [(amp) (FWHM/1)' freq'];
                lb= [((amp) - abs(amp./4)) ((FWHM/1)- 10)' (freq-10)'];
                ub= [((amp) + abs(amp./4)) ((FWHM/1) + 10)' (freq+10)'];
                options = optimset('MaxFunEvals', fit_para.max_fun(cur_vox), 'MaxIter', fit_para.max_iter(cur_vox), 'Algorithm','trust-region-reflective', 'Display', 'off'); % optimization parameters 
                [fit_param1,~,residual,~,~,~,J] = lsqnonlin(@(x)model_estimate_equation2(x,details_para.t,(data_value.FID_TD(:,cur_vox)),model_para.type(cur_vox)),Initx,lb,ub,options); % non-linear least squares algorithm
                [sig_sum,C,Z_sig] = model_estimate(fit_param1,details_para.t,model_para.type(cur_vox));
                fit_sig=sig_sum;
                Ind_sig=C;
                signal_fitted(:,1)=fftshift(fft(fit_sig));
                count=2;
                while(1)
        %                 sp1 = csaps(ppm,real(data_value.FID_FD(:,i)-signal_fitted(:,count-1)),30);
        %                 re_bl =fnval(sp1,ppm);
                        re_bl = ssa(real(data_value.FID_FD(:,cur_vox)-signal_fitted(:,count-1)),30); % estimate the baseline for real part
        %                 sp2 = csaps(ppm,imag(data_value.FID_FD(:,i)-signal_fitted(:,count-1)),30);
        %                 im_bl =fnval(sp2,ppm);
                        im_bl = ssa(imag(data_value.FID_FD(:,cur_vox)-signal_fitted(:,count-1)),30); % estimate the baseline for imaginary part
                        baseline(:,count)=re_bl + 1i*im_bl;
                        baseline_t(:,count)=ifft(ifftshift(baseline(:,count)));
                        [ fit_sig,amp,Ind_sig,freq,Z_sig,FWHM] = HSVDPK_new(data_value.FID_TD(:,cur_vox)-baseline_t(:,count),details_para.fres,details_para.Fs,fit_para.order(cur_vox));
                        Initx = [amp (FWHM/1)' freq'];
                        lb= [(amp - abs(abs(amp./4))) ((FWHM/1)- 10)' (freq-10)'];
                        ub= [(amp + abs(abs(amp./4))) ((FWHM/1) + 10)' (freq+10)'];
                        options = optimset('MaxFunEvals', fit_para.max_fun(cur_vox), 'MaxIter', fit_para.max_iter(cur_vox), 'Algorithm','trust-region-reflective', 'Display', 'off'); % optimization parameters 
                        [fit_param1,~,residual,~,~,~,J] = lsqnonlin(@(x)model_estimate_equation2(x,details_para.t,(data_value.FID_TD(:,cur_vox)),model_para.type(cur_vox)),Initx,lb,ub,options); % non-linear least squares algorithm
                        [sig_sum,C,Z_sig] = model_estimate(fit_param1,details_para.t,model_para.type(cur_vox));
                        fit_sig=sig_sum;
                        Ind_sig=C;
                        signal_fitted(:,count)=fftshift(fft(fit_sig));
                        if (std(abs(baseline(:,count)-baseline(:,count-1)))<=1e-5 || count==70)
                            break;
                        end
                        count=count+1;
                end
                se = calculate_SE(fit_param1,residual,J) ;
                se=se(1:3:end);
        %         [ fit_sig,amp,Ind_sig,freq,Z_sig,FWHM] = HSVDPK_new(data_value.FID_TD_b(:,i),details_para.fres,details_para.Fs,fit_para.order(cur_vox));
                 model_select=1;
                Initx = [amp (FWHM)' freq'];
                freq=fit_param1(:,3);
            %     lb= [(amp-(amp/4)) ((FWHM)-50)' (freq-20)'];
            %     ub= [(amp+(amp/4)) ((FWHM)+50)' (freq+20)'];
            %      options = optimset('MaxFunEvals', 500, 'MaxIter', 20, 'Algorithm','trust-region-reflective', 'Display', 'off'); % optimization parameters 
            %     [fit_param11,~,residual,~,~,~,J] = lsqnonlin(@(x)model_estimate_equation2(x,details_para.t,data_value.FID_TD(:,i),model_select),Initx,[],[],options); % non-linear least squares algorithm
            %     [fit_sig,Ind_sig,Z_sig] = model_estimate(fit_param11,details_para.t,model_select);
                ind=find((freq>=details_para.fit_freq(1)-20 & freq<=details_para.fit_freq(2)+20));
                ind_NAA=find((freq>=details_para.fit_freq_NAA(1)-20 & freq<=details_para.fit_freq_NAA(2)+20));
                ind_Cr=find((freq>=details_para.fit_freq_Cr(1)-20 & freq<=details_para.fit_freq_Cr(2)+20));
                ind_Cho=find((freq>=details_para.fit_freq_Cho(1)-20 & freq<=details_para.fit_freq_Cho(2)+20));
                ind_Myo=find((freq>=details_para.fit_freq_Cr(1)-20 & freq<=details_para.fit_freq_Cr(2)+20));
                ind_Lac=find((freq>=details_para.fit_freq_Cho(1)-20 & freq<=details_para.fit_freq_Cho(2)+20));
                half_sig2=sum(Ind_sig(:,ind),2);
                data_value.fit_TD(:,cur_vox)=half_sig2;
                data_value.fit_FD_NAA(:,cur_vox)=fftshift(fft(sum(Ind_sig(:,ind_NAA),2)));
                data_value.fit_FD_Cr(:,cur_vox)=fftshift(fft(sum(Ind_sig(:,ind_Cr),2)));
                data_value.fit_FD_Cho(:,cur_vox)=fftshift(fft(sum(Ind_sig(:,ind_Cho),2)));
                data_value.fit_FD_Myo(:,cur_vox)=fftshift(fft(sum(Ind_sig(:,ind_Myo),2)));
                data_value.fit_FD_Lac(:,cur_vox)=fftshift(fft(sum(Ind_sig(:,ind_Lac),2)));
                data_value.fit_FD(:,cur_vox)=fftshift(fft(data_value.fit_TD(:,cur_vox)));
                data_value.fit_FD_temp(:,cur_vox)=fftshift(fft(data_value.fit_TD(:,cur_vox)));
                details_para.se(:,cur_vox)=sum(se(ind));
        %         min_ref=min(real(data_value.fit_FD(1:30,i)));
        %         data_value.fit_FD(:,i)=data_value.fit_FD(:,i);
        %         data_value.peak_area(:,i)=trapz(details_para.ppm,real(data_value.fit_FD(:,i)));
                data_value.peak_area(:,cur_vox)=trapz(details_para.ppm(details_para.fit_disp_seg(1):details_para.fit_disp_seg(2)),abs(real(data_value.fit_FD(details_para.fit_disp_seg(1):details_para.fit_disp_seg(2),cur_vox))));
        %         base_temp=real(data_value.FID_FD(:,i)) - real(fftshift(fft(baseline_t(:,count))));
        %         data_value.peak_area(:,i)=trapz(details_para.ppm(details_para.fit_disp_seg(1):details_para.fit_disp_seg(2)),abs(real(base_temp(details_para.fit_disp_seg(1):details_para.fit_disp_seg(2)))));
        %         data_value.peak_area(:,i)=trapz(details_para.ppm(details_para.fit_disp_seg(1):details_para.fit_disp_seg(2)),abs(data_value.FID_FD(details_para.fit_disp_seg(1):details_para.fit_disp_seg(2),i))-abs(baseline(details_para.fit_disp_seg(1):details_para.fit_disp_seg(2),count-1)));%details_para.fit_disp_seg(1)
                data_value.peak_area_NAA(:,cur_vox)=trapz(details_para.ppm,abs(real(data_value.fit_FD_NAA(:,cur_vox))));
                data_value.peak_area_Cr(:,cur_vox)=trapz(details_para.ppm,abs(real(data_value.fit_FD_Cr(:,cur_vox))));
                data_value.peak_area_Cho(:,cur_vox)=trapz(details_para.ppm,abs(real(data_value.fit_FD_Cho(:,cur_vox))));
        %         data_value.peak_area_NAA(:,i)=trapz(details_para.ppm_referenced(details_para.fit_disp_seg_NAA(1):details_para.fit_disp_seg_NAA(2)),abs(real(data_value.fit_FD_NAA(details_para.fit_disp_seg_NAA(1):details_para.fit_disp_seg_NAA(2),i))));
        %         data_value.peak_area_NAA(:,i)=max(real(data_value.fit_FD_NAA(:,i)));
        %         data_value.peak_area_NAA(:,i)=max(real(data_value.fit_FD(details_para.fit_disp_seg_NAA(1):details_para.fit_disp_seg_NAA(2),i)));
        %         data_value.peak_area_Cr(:,i)=trapz(details_para.ppm_referenced(details_para.fit_disp_seg(1):details_para.fit_disp_seg(2)),real(data_value.fit_FD_Cr(details_para.fit_disp_seg(1):details_para.fit_disp_seg(2),i)));
        %         data_value.peak_area_Cho(:,i)=trapz(details_para.ppm_referenced(details_para.fit_disp_seg(1):details_para.fit_disp_seg(2)),real(data_value.fit_FD_Cho(details_para.fit_disp_seg(1):details_para.fit_disp_seg(2),i)));
                data_value.peak_area_Myo(:,cur_vox)=trapz(details_para.ppm_referenced(details_para.fit_disp_seg(1):details_para.fit_disp_seg(2)),real(data_value.fit_FD_Myo(details_para.fit_disp_seg(1):details_para.fit_disp_seg(2),cur_vox)));
                data_value.peak_area_Lac(:,cur_vox)=trapz(details_para.ppm_referenced(details_para.fit_disp_seg(1):details_para.fit_disp_seg(2)),real(data_value.fit_FD_Lac(details_para.fit_disp_seg(1):details_para.fit_disp_seg(2),cur_vox)));
                data_value.corr_baseline(:,cur_vox)=real(data_value.FID_FD(:,cur_vox)) - real(fftshift(fft(baseline_t(:,count))));
        %         data_value.peak_area_NAA(:,i)=max(data_value.corr_baseline(:,i));
                data_value.noise(:,cur_vox)=std(data_value.corr_baseline(end-50:end-10,cur_vox));%std(real(data_value.FID_FD(:,i))-real(fftshift(fft(fit_sig))));
        %         data_value.noise(:,i)=std(real(data_value.FID_FD(end-30:end,i)));
                data_value.snr(:,cur_vox)=data_value.peak_area_NAA(:,cur_vox)./data_value.noise(:,cur_vox);
                data_value.NAA_Cr(:,cur_vox)=max(abs(real(data_value.FID_FD(details_para.fit_disp_seg_NAA(1):details_para.fit_disp_seg_NAA(2),cur_vox))))./max(abs(real(data_value.FID_FD(details_para.fit_disp_seg_Cr(1):details_para.fit_disp_seg_Cr(2),cur_vox))));%max(abs(real(data_value.fit_FD_NAA(:,cur_vox))))./max(abs(real(data_value.fit_FD_Cr(:,cur_vox))));
                data_value.fit_temp(:,cur_vox)=data_value.fit_TD(:,cur_vox);
                data_value.Cho_Cr(:,cur_vox)=max(abs(real(data_value.FID_FD(details_para.fit_disp_seg_Cho(1):details_para.fit_disp_seg_Cho(2),cur_vox))))./max(abs(real(data_value.FID_FD(details_para.fit_disp_seg_Cr(1):details_para.fit_disp_seg_Cr(2),cur_vox))));%max(abs(real(data_value.fit_FD_Cho(:,cur_vox))))./max(abs(real(data_value.fit_FD_Cr(:,cur_vox))));
                data_value.Myo_Cr(:,cur_vox)=data_value.peak_area_Myo(:,cur_vox)./data_value.peak_area_Cr(:,cur_vox);
                data_value.Lac_Cr(:,cur_vox)=data_value.peak_area_Lac(:,cur_vox)./data_value.peak_area_Cr(:,cur_vox);
                data_value.fit_TD(:,cur_vox)=1*data_value.fit_TD(:,cur_vox) + 1*baseline_t(:,count);%baseline_t(:,count-1);%data_value.fit_TD(:,i) + ifft(ifftshift(data_value.FID_FD_base(:,i)));
                data_value.fit_FD(:,cur_vox)=fftshift(fft(data_value.fit_TD(:,cur_vox)));
                data_value.fit_baseline(:,cur_vox)= data_value.fit_FD(:,cur_vox)-fftshift(fft(data_value.fit_temp(:,cur_vox)));
                waitbar(1);
                details_para.fit_flag(:,cur_vox)=1;
                cur_data(cnt,1)=data_value.peak_area(:,cnt);
                details_para.mean_peak_area=mean(cur_data(cnt,1));
                details_para.std_peak_area=std(cur_data(cnt,1));
                cur_data(cnt,2)=data_value.snr(:,cnt);
                cur_data(cnt,3)=data_value.NAA_Cr(:,cnt);
                cur_data(cnt,4)=data_value.Cho_Cr(:,cnt);
        elseif (details_para.fit_method==1)
            for i = 1:details_para.num_vox
                details_para.fla=0;
                details_para.seg = [1 length(data_value.FID_TD(:,i))];
        %         if i==1
        %              [ fit_sig,amp,Ind_sig,freq,Z_sig,FWHM] = HSVDPK_new(data_value.FID_TD(:,i),details_para.fres,details_para.Fs,fit_para.order(cur_vox));
                    [No_pks,pk_loc,FWHM,signal_power_each] = find_peaks_new(data_value.FID_TD(:,i),details_para.fres);
        %         end
                freq=details_para.fres(pk_loc);
                FWHM=FWHM/1;
                [fit_sig,~,amp_est] = model_varpro(data_value.FID_TD(:,i),details_para.fres(pk_loc), FWHM, ones(1,length(pk_loc)),details_para.t);
                amp_est=amp_est;
                Initx = [amp_est' (FWHM/1)' freq'];
                lb= [(amp_est- data_value.FID_TD(1,i))' ((FWHM/1)- 10)' (freq-10)'];
                ub= [(amp_est + data_value.FID_TD(1,i))' ((FWHM/1) + 10)' (freq+10)'];
                options = optimset('MaxFunEvals', 5000, 'MaxIter', 5000, 'Algorithm','trust-region-reflective', 'Display', 'off'); % optimization parameters 
                [fit_param1,~,residual,~,~,~,J] = lsqnonlin(@(x)model_estimate_equation2(x,details_para.t,(data_value.FID_TD(:,i)),1),Initx,lb,ub,options); % non-linear least squares algorithm
                [sig_sum,C,Z_sig] = model_estimate(fit_param1,details_para.t,1);
                Ind_sig=C;
                ind=find((freq>=details_para.fit_freq(1)-20 & freq<=details_para.fit_freq(2)+20));
                ind_NAA=find((freq>=details_para.fit_freq_NAA(1)-20 & freq<=details_para.fit_freq_NAA(2)+20));
                ind_Cr=find((freq>=details_para.fit_freq_Cr(1)-20 & freq<=details_para.fit_freq_Cr(2)+20));
                ind_Cho=find((freq>=details_para.fit_freq_Cho(1)-20 & freq<=details_para.fit_freq_Cho(2)+20));
                half_sig2=sum(Ind_sig(:,ind),2);
                data_value.fit_TD(:,i)=half_sig2;
                data_value.fit_FD_NAA(:,i)=fftshift(fft(sum(Ind_sig(:,ind_NAA),2)));
                data_value.fit_FD_Cr(:,i)=fftshift(fft(sum(Ind_sig(:,ind_Cr),2)));
                data_value.fit_FD_Cho(:,i)=fftshift(fft(sum(Ind_sig(:,ind_Cho),2)));
                data_value.fit_FD(:,i)=fftshift(fft(data_value.fit_TD(:,i)));
                data_value.peak_area(:,i)=trapz(details_para.ppm_referenced,abs(data_value.fit_FD(:,i)));
                data_value.peak_area_NAA(:,i)=trapz(details_para.ppm_referenced,abs(data_value.fit_FD_NAA(:,i)));
                data_value.peak_area_Cr(:,i)=trapz(details_para.ppm_referenced,abs(data_value.fit_FD_Cr(:,i)));
                data_value.peak_area_Cho(:,i)=trapz(details_para.ppm_referenced,abs(data_value.fit_FD_Cho(:,i)));
                data_value.noise(:,i)=std(data_value.FID_FD(:,i)-fftshift(fft(fit_sig)));
                data_value.snr(:,i)=data_value.peak_area_NAA(:,i)./data_value.noise(:,i);
                data_value.NAA_Cr(:,i)=data_value.peak_area_NAA(:,i)./data_value.peak_area_Cr(:,i);
                data_value.Cho_Cr(:,i)=data_value.peak_area_Cho(:,i)./data_value.peak_area_Cr(:,i);
                waitbar(i/details_para.num_vox);
            end
        else
            for i = 1:details_para.num_vox
                details_para.fla=0;
                data_value.FID_FD_base(:,i)=ssa(data_value.FID_FD(:,i),100);
                data_value.FID_FD_b(:,i)=data_value.FID_FD(:,i)-data_value.FID_FD_base(:,i);
                data_value.FID_TD_b(:,i)=ifft(ifftshift(data_value.FID_FD_b(:,i)));
                details_para.seg = [1 length(data_value.FID_TD(:,i))];
        %         if i==1
                    [No_pks,pk_loc,FWHM(1,:),signal_power(1,:)] = find_peaks_new(data_value.FID_TD_b(:,i),details_para.fres);
        %         end
        %         [ fit_sig,amp,Ind_sig,freq,Z_sig,FWHM] = HSVDPK_new(data_value.FID_TD(:,i),details_para.fres,details_para.Fs,fit_para.order(cur_vox));
                freq=details_para.fres(pk_loc);
                FWHM=1*FWHM*(1E6/details_para.Tf)./2;
                FWHM(1:end)=FWHM(1);
                for kk=1:length(FWHM)%length(pk_loc)
                    Data2bFit=(fftshift(fft(data_value.FID_TD_b(:,i))))/1;
                    PEAK_ppm=4.7 + freq(kk)*(1E6/details_para.Tf);%details_para.ppm_referenced(pk_loc(kk));
                    z=abs(details_para.ppm_referenced-(PEAK_ppm+FWHM(kk)));
                    ub=find(min(z)==z);
                    z=abs(details_para.ppm_referenced-(PEAK_ppm-FWHM(kk)));
                    lb=find(min(z)==z); 
                    freqrange = details_para.ppm_referenced(lb:ub);
                    SumSpec_half = (Data2bFit(lb:ub,:));
                    max_peak=max(SumSpec_half);
                    Init = [max_peak./2 FWHM(kk) PEAK_ppm 0 0 0 ];
                    dataParams_tot(kk,:) = FitPeaks(freqrange, SumSpec_half, Init);
                    dataParams_tot(kk,4:6)=0;
        %             dataParams_tot(kk,2)=0;
                    fitted_signal_tot(kk,:)=LorentzModel(dataParams_tot(kk,:), details_para.ppm_referenced);
                    Ind_sig(:,kk)=fitted_signal_tot(kk,:);
                    Ind_sig(:,kk)=ifft(ifftshift(Ind_sig(:,kk)));
                end
                fit_sig=sum(Ind_sig,2);
        %         fit_sig=ifft(ifftshift(fit_sig));
                ind=find((freq>=details_para.fit_freq(1)-20 & freq<=details_para.fit_freq(2)+20));
                ind_NAA=find((freq>=details_para.fit_freq_NAA(1)-20 & freq<=details_para.fit_freq_NAA(2)+20));
                ind_Cr=find((freq>=details_para.fit_freq_Cr(1)-20 & freq<=details_para.fit_freq_Cr(2)+20));
                ind_Cho=find((freq>=details_para.fit_freq_Cho(1)-20 & freq<=details_para.fit_freq_Cho(2)+20));
                half_sig2=sum(Ind_sig(:,ind),2);
                data_value.fit_TD(:,i)=half_sig2;
                data_value.fit_FD_NAA(:,i)=fftshift(fft(sum(Ind_sig(:,ind_NAA),2)));
                data_value.fit_FD_Cr(:,i)=fftshift(fft(sum(Ind_sig(:,ind_Cr),2)));
                data_value.fit_FD_Cho(:,i)=fftshift(fft(sum(Ind_sig(:,ind_Cho),2)));
                data_value.fit_FD(:,i)=fftshift(fft(data_value.fit_TD(:,i)));
                data_value.peak_area(:,i)=trapz(details_para.ppm_referenced,abs(data_value.fit_FD(:,i)));
                data_value.peak_area_NAA(:,i)=trapz(details_para.ppm_referenced,abs(data_value.fit_FD_NAA(:,i)));
                data_value.peak_area_Cr(:,i)=trapz(details_para.ppm_referenced,abs(data_value.fit_FD_Cr(:,i)));
                data_value.peak_area_Cho(:,i)=trapz(details_para.ppm_referenced,abs(data_value.fit_FD_Cho(:,i)));
                data_value.noise(:,i)=std(data_value.FID_FD(:,i)-fftshift(fft(fit_sig)));
                data_value.snr(:,i)=data_value.peak_area_NAA(:,i)./data_value.noise(:,i);
                data_value.NAA_Cr(:,i)=data_value.peak_area_NAA(:,i)./data_value.peak_area_Cr(:,i);
                data_value.Cho_Cr(:,i)=data_value.peak_area_Cho(:,i)./data_value.peak_area_Cr(:,i);
                data_value.fit_TD(:,i)=data_value.fit_TD(:,i) + ifft(ifftshift(data_value.FID_FD_base(:,i)));
                data_value.fit_FD(:,i)=fftshift(fft(data_value.fit_TD(:,i)));
                waitbar(i/details_para.num_vox);
            end
        end
        
        
    else
        
        if details_para.fit_method==0
            for i = 1:details_para.num_vox
          
                cnt=i;
                re_bl = ssa(real(data_value.FID_FD(:,i)),40); % estimate the baseline for real part
        %         sp2 = csaps(ppm,imag(data_value.FID_FD(:,i)),30);
        %         im_bl = fnval(sp2,ppm);
                im_bl = ssa(imag(data_value.FID_FD(:,i)),40); % estimate the baseline for imaginary part
                data_value.FID_FD_base(:,i)=re_bl + 1i*im_bl;
                data_value.FID_FD_b(:,i)=data_value.FID_FD(:,i)-data_value.FID_FD_base(:,i);
                data_value.FID_TD_b(:,i)=ifft(ifftshift(data_value.FID_FD_b(:,i)));
                baseline(:,1)=data_value.FID_FD_base(:,i);
                baseline_t(:,1)=ifft(ifftshift(data_value.FID_FD_base(:,i)));
                [ fit_sig,amp,Ind_sig,freq,Z_sig,FWHM] = HSVDPK_new(data_value.FID_TD_b(:,i),details_para.fres,details_para.Fs,fit_para.order(cur_vox));
                Initx = [(amp) (FWHM/1)' freq'];
                lb= [((amp) - abs(amp./4)) ((FWHM/1)- 10)' (freq-10)'];
                ub= [((amp) + abs(amp./4)) ((FWHM/1) + 10)' (freq+10)'];
                options = optimset('MaxFunEvals', fit_para.max_fun(cur_vox), 'MaxIter', fit_para.max_iter(cur_vox), 'Algorithm','trust-region-reflective', 'Display', 'off'); % optimization parameters 
                [fit_param1,~,residual,~,~,~,J] = lsqnonlin(@(x)model_estimate_equation2(x,details_para.t,(data_value.FID_TD(:,i)),model_para.type(cur_vox)),Initx,lb,ub,options); % non-linear least squares algorithm
                [sig_sum,C,Z_sig] = model_estimate(fit_param1,details_para.t,model_para.type(cur_vox));
                fit_sig=sig_sum;
                Ind_sig=C;
                signal_fitted(:,1)=fftshift(fft(fit_sig));
                count=2;
                while(1)
        %                 sp1 = csaps(ppm,real(data_value.FID_FD(:,i)-signal_fitted(:,count-1)),30);
        %                 re_bl =fnval(sp1,ppm);
                        re_bl = ssa(real(data_value.FID_FD(:,i)-signal_fitted(:,count-1)),30); % estimate the baseline for real part
        %                 sp2 = csaps(ppm,imag(data_value.FID_FD(:,i)-signal_fitted(:,count-1)),30);
        %                 im_bl =fnval(sp2,ppm);
                        im_bl = ssa(imag(data_value.FID_FD(:,i)-signal_fitted(:,count-1)),30); % estimate the baseline for imaginary part
                        baseline(:,count)=re_bl + 1i*im_bl;
                        baseline_t(:,count)=ifft(ifftshift(baseline(:,count)));
                        [ fit_sig,amp,Ind_sig,freq,Z_sig,FWHM] = HSVDPK_new(data_value.FID_TD(:,i)-baseline_t(:,count),details_para.fres,details_para.Fs,fit_para.order(cur_vox));
                        Initx = [amp (FWHM/1)' freq'];
                        lb= [(amp - abs(abs(amp./4))) ((FWHM/1)- 10)' (freq-10)'];
                        ub= [(amp + abs(abs(amp./4))) ((FWHM/1) + 10)' (freq+10)'];
                        options = optimset('MaxFunEvals', fit_para.max_fun(cur_vox), 'MaxIter', fit_para.max_iter(cur_vox), 'Algorithm','trust-region-reflective', 'Display', 'off'); % optimization parameters 
                        [fit_param1,~,residual,~,~,~,J] = lsqnonlin(@(x)model_estimate_equation2(x,details_para.t,(data_value.FID_TD(:,i)),model_para.type(cur_vox)),Initx,lb,ub,options); % non-linear least squares algorithm
                        [sig_sum,C,Z_sig] = model_estimate(fit_param1,details_para.t,model_para.type(cur_vox));
                        fit_sig=sig_sum;
                        Ind_sig=C;
                        signal_fitted(:,count)=fftshift(fft(fit_sig));
                        if (std(abs(baseline(:,count)-baseline(:,count-1)))<=1e-5 || count==70)
                            break;
                        end
                        count=count+1;
                end
                se = calculate_SE(fit_param1,residual,J) ;
                se=se(1:3:end);
        %         [ fit_sig,amp,Ind_sig,freq,Z_sig,FWHM] = HSVDPK_new(data_value.FID_TD_b(:,i),details_para.fres,details_para.Fs,fit_para.order(cur_vox));
                 model_select=1;
                Initx = [amp (FWHM)' freq'];
                freq=fit_param1(:,3);
            %     lb= [(amp-(amp/4)) ((FWHM)-50)' (freq-20)'];
            %     ub= [(amp+(amp/4)) ((FWHM)+50)' (freq+20)'];
            %      options = optimset('MaxFunEvals', 500, 'MaxIter', 20, 'Algorithm','trust-region-reflective', 'Display', 'off'); % optimization parameters 
            %     [fit_param11,~,residual,~,~,~,J] = lsqnonlin(@(x)model_estimate_equation2(x,details_para.t,data_value.FID_TD(:,i),model_select),Initx,[],[],options); % non-linear least squares algorithm
            %     [fit_sig,Ind_sig,Z_sig] = model_estimate(fit_param11,details_para.t,model_select);
                ind=find((freq>=details_para.fit_freq(1)-20 & freq<=details_para.fit_freq(2)+20));
                ind_NAA=find((freq>=details_para.fit_freq_NAA(1)-20 & freq<=details_para.fit_freq_NAA(2)+20));
                ind_Cr=find((freq>=details_para.fit_freq_Cr(1)-20 & freq<=details_para.fit_freq_Cr(2)+20));
                ind_Cho=find((freq>=details_para.fit_freq_Cho(1)-20 & freq<=details_para.fit_freq_Cho(2)+20));
                ind_Myo=find((freq>=details_para.fit_freq_Cr(1)-20 & freq<=details_para.fit_freq_Cr(2)+20));
                ind_Lac=find((freq>=details_para.fit_freq_Cho(1)-20 & freq<=details_para.fit_freq_Cho(2)+20));
                half_sig2=sum(Ind_sig(:,ind),2);
                data_value.fit_TD(:,i)=half_sig2;
                data_value.fit_FD_NAA(:,i)=fftshift(fft(sum(Ind_sig(:,ind_NAA),2)));
                data_value.fit_FD_Cr(:,i)=fftshift(fft(sum(Ind_sig(:,ind_Cr),2)));
                data_value.fit_FD_Cho(:,i)=fftshift(fft(sum(Ind_sig(:,ind_Cho),2)));
                data_value.fit_FD_Myo(:,i)=fftshift(fft(sum(Ind_sig(:,ind_Myo),2)));
                data_value.fit_FD_Lac(:,i)=fftshift(fft(sum(Ind_sig(:,ind_Lac),2)));
                data_value.fit_FD(:,i)=fftshift(fft(data_value.fit_TD(:,i)));
                data_value.fit_FD_temp(:,i)=fftshift(fft(data_value.fit_TD(:,i)));
                details_para.se(:,i)=sum(se(ind));
        %         min_ref=min(real(data_value.fit_FD(1:30,i)));
        %         data_value.fit_FD(:,i)=data_value.fit_FD(:,i);
        %         data_value.peak_area(:,i)=trapz(details_para.ppm,real(data_value.fit_FD(:,i)));
                data_value.peak_area(:,i)=trapz(details_para.ppm(details_para.fit_disp_seg(1):details_para.fit_disp_seg(2)),abs(real(data_value.fit_FD(details_para.fit_disp_seg(1):details_para.fit_disp_seg(2),i))));
        %         base_temp=real(data_value.FID_FD(:,i)) - real(fftshift(fft(baseline_t(:,count))));
        %         data_value.peak_area(:,i)=trapz(details_para.ppm(details_para.fit_disp_seg(1):details_para.fit_disp_seg(2)),abs(real(base_temp(details_para.fit_disp_seg(1):details_para.fit_disp_seg(2)))));
        %         data_value.peak_area(:,i)=trapz(details_para.ppm(details_para.fit_disp_seg(1):details_para.fit_disp_seg(2)),abs(data_value.FID_FD(details_para.fit_disp_seg(1):details_para.fit_disp_seg(2),i))-abs(baseline(details_para.fit_disp_seg(1):details_para.fit_disp_seg(2),count-1)));%details_para.fit_disp_seg(1)
                data_value.peak_area_NAA(:,i)=trapz(details_para.ppm,abs(real(data_value.fit_FD_NAA(:,i))));
                data_value.peak_area_Cr(:,i)=trapz(details_para.ppm,abs(real(data_value.fit_FD_Cr(:,i))));
                data_value.peak_area_Cho(:,i)=trapz(details_para.ppm,abs(real(data_value.fit_FD_Cho(:,i))));
        %         data_value.peak_area_NAA(:,i)=trapz(details_para.ppm_referenced(details_para.fit_disp_seg_NAA(1):details_para.fit_disp_seg_NAA(2)),abs(real(data_value.fit_FD_NAA(details_para.fit_disp_seg_NAA(1):details_para.fit_disp_seg_NAA(2),i))));
        %         data_value.peak_area_NAA(:,i)=max(real(data_value.fit_FD_NAA(:,i)));
        %         data_value.peak_area_NAA(:,i)=max(real(data_value.fit_FD(details_para.fit_disp_seg_NAA(1):details_para.fit_disp_seg_NAA(2),i)));
        %         data_value.peak_area_Cr(:,i)=trapz(details_para.ppm_referenced(details_para.fit_disp_seg(1):details_para.fit_disp_seg(2)),real(data_value.fit_FD_Cr(details_para.fit_disp_seg(1):details_para.fit_disp_seg(2),i)));
        %         data_value.peak_area_Cho(:,i)=trapz(details_para.ppm_referenced(details_para.fit_disp_seg(1):details_para.fit_disp_seg(2)),real(data_value.fit_FD_Cho(details_para.fit_disp_seg(1):details_para.fit_disp_seg(2),i)));
                data_value.peak_area_Myo(:,i)=trapz(details_para.ppm_referenced(details_para.fit_disp_seg(1):details_para.fit_disp_seg(2)),real(data_value.fit_FD_Myo(details_para.fit_disp_seg(1):details_para.fit_disp_seg(2),i)));
                data_value.peak_area_Lac(:,i)=trapz(details_para.ppm_referenced(details_para.fit_disp_seg(1):details_para.fit_disp_seg(2)),real(data_value.fit_FD_Lac(details_para.fit_disp_seg(1):details_para.fit_disp_seg(2),i)));
                data_value.corr_baseline(:,i)=real(data_value.FID_FD(:,i)) - real(fftshift(fft(baseline_t(:,count))));
        %         data_value.peak_area_NAA(:,i)=max(data_value.corr_baseline(:,i));
                data_value.noise(:,i)=std(data_value.corr_baseline(end-50:end-10,i));%std(real(data_value.FID_FD(:,i))-real(fftshift(fft(fit_sig))));
        %         data_value.noise(:,i)=std(real(data_value.FID_FD(end-30:end,i)));
                data_value.snr(:,i)=data_value.peak_area_NAA(:,i)./data_value.noise(:,i);
                data_value.NAA_Cr(:,i)=max(abs(real(data_value.FID_FD(details_para.fit_disp_seg_NAA(1):details_para.fit_disp_seg_NAA(2),i))))./max(abs(real(data_value.FID_FD(details_para.fit_disp_seg_Cr(1):details_para.fit_disp_seg_Cr(2),i))));%max(abs(real(data_value.fit_FD_NAA(:,i))))./max(abs(real(data_value.fit_FD_Cr(:,i))));
                data_value.fit_temp(:,i)=data_value.fit_TD(:,i);
                data_value.Cho_Cr(:,i)=max(abs(real(data_value.FID_FD(details_para.fit_disp_seg_Cho(1):details_para.fit_disp_seg_Cho(2),i))))./max(abs(real(data_value.FID_FD(details_para.fit_disp_seg_Cr(1):details_para.fit_disp_seg_Cr(2),i))));%max(abs(real(data_value.fit_FD_Cho(:,i))))./max(abs(real(data_value.fit_FD_Cr(:,i))));
                data_value.Myo_Cr(:,i)=data_value.peak_area_Myo(:,i)./data_value.peak_area_Cr(:,i);
                data_value.Lac_Cr(:,i)=data_value.peak_area_Lac(:,i)./data_value.peak_area_Cr(:,i);
                data_value.fit_TD(:,i)=1*data_value.fit_TD(:,i) + 1*baseline_t(:,count);%baseline_t(:,count-1);%data_value.fit_TD(:,i) + ifft(ifftshift(data_value.FID_FD_base(:,i)));
                data_value.fit_FD(:,i)=fftshift(fft(data_value.fit_TD(:,i)));
                data_value.fit_baseline(:,i)= data_value.fit_FD(:,i)-fftshift(fft(data_value.fit_temp(:,i)));
                waitbar(i/details_para.num_vox);
                details_para.fit_flag(:,i)=1;
                cur_data(cnt,1)=data_value.peak_area(:,cnt);
                details_para.mean_peak_area=mean(cur_data(cnt,1));
                details_para.std_peak_area=std(cur_data(cnt,1));
                cur_data(cnt,2)=data_value.snr(:,cnt);
                cur_data(cnt,3)=data_value.NAA_Cr(:,cnt);
                cur_data(cnt,4)=data_value.Cho_Cr(:,cnt);
            end
        elseif (details_para.fit_method==1)
            for i = 1:details_para.num_vox
                details_para.fla=0;
                details_para.seg = [1 length(data_value.FID_TD(:,i))];
        %         if i==1
        %              [ fit_sig,amp,Ind_sig,freq,Z_sig,FWHM] = HSVDPK_new(data_value.FID_TD(:,i),details_para.fres,details_para.Fs,fit_para.order(cur_vox));
                    [No_pks,pk_loc,FWHM,signal_power_each] = find_peaks_new(data_value.FID_TD(:,i),details_para.fres);
        %         end
                freq=details_para.fres(pk_loc);
                FWHM=FWHM/1;
                [fit_sig,~,amp_est] = model_varpro(data_value.FID_TD(:,i),details_para.fres(pk_loc), FWHM, ones(1,length(pk_loc)),details_para.t);
                amp_est=amp_est;
                Initx = [amp_est' (FWHM/1)' freq'];
                lb= [(amp_est- data_value.FID_TD(1,i))' ((FWHM/1)- 10)' (freq-10)'];
                ub= [(amp_est + data_value.FID_TD(1,i))' ((FWHM/1) + 10)' (freq+10)'];
                options = optimset('MaxFunEvals', 5000, 'MaxIter', 5000, 'Algorithm','trust-region-reflective', 'Display', 'off'); % optimization parameters 
                [fit_param1,~,residual,~,~,~,J] = lsqnonlin(@(x)model_estimate_equation2(x,details_para.t,(data_value.FID_TD(:,i)),1),Initx,lb,ub,options); % non-linear least squares algorithm
                [sig_sum,C,Z_sig] = model_estimate(fit_param1,details_para.t,1);
                Ind_sig=C;
                ind=find((freq>=details_para.fit_freq(1)-20 & freq<=details_para.fit_freq(2)+20));
                ind_NAA=find((freq>=details_para.fit_freq_NAA(1)-20 & freq<=details_para.fit_freq_NAA(2)+20));
                ind_Cr=find((freq>=details_para.fit_freq_Cr(1)-20 & freq<=details_para.fit_freq_Cr(2)+20));
                ind_Cho=find((freq>=details_para.fit_freq_Cho(1)-20 & freq<=details_para.fit_freq_Cho(2)+20));
                half_sig2=sum(Ind_sig(:,ind),2);
                data_value.fit_TD(:,i)=half_sig2;
                data_value.fit_FD_NAA(:,i)=fftshift(fft(sum(Ind_sig(:,ind_NAA),2)));
                data_value.fit_FD_Cr(:,i)=fftshift(fft(sum(Ind_sig(:,ind_Cr),2)));
                data_value.fit_FD_Cho(:,i)=fftshift(fft(sum(Ind_sig(:,ind_Cho),2)));
                data_value.fit_FD(:,i)=fftshift(fft(data_value.fit_TD(:,i)));
                data_value.peak_area(:,i)=trapz(details_para.ppm_referenced,abs(data_value.fit_FD(:,i)));
                data_value.peak_area_NAA(:,i)=trapz(details_para.ppm_referenced,abs(data_value.fit_FD_NAA(:,i)));
                data_value.peak_area_Cr(:,i)=trapz(details_para.ppm_referenced,abs(data_value.fit_FD_Cr(:,i)));
                data_value.peak_area_Cho(:,i)=trapz(details_para.ppm_referenced,abs(data_value.fit_FD_Cho(:,i)));
                data_value.noise(:,i)=std(data_value.FID_FD(:,i)-fftshift(fft(fit_sig)));
                data_value.snr(:,i)=data_value.peak_area_NAA(:,i)./data_value.noise(:,i);
                data_value.NAA_Cr(:,i)=data_value.peak_area_NAA(:,i)./data_value.peak_area_Cr(:,i);
                data_value.Cho_Cr(:,i)=data_value.peak_area_Cho(:,i)./data_value.peak_area_Cr(:,i);
                waitbar(i/details_para.num_vox);
                details_para.fit_flag(:,i)=1;
            end
        else
            for i = 1:details_para.num_vox
                details_para.fla=0;
                data_value.FID_FD_base(:,i)=ssa(data_value.FID_FD(:,i),100);
                data_value.FID_FD_b(:,i)=data_value.FID_FD(:,i)-data_value.FID_FD_base(:,i);
                data_value.FID_TD_b(:,i)=ifft(ifftshift(data_value.FID_FD_b(:,i)));
                details_para.seg = [1 length(data_value.FID_TD(:,i))];
        %         if i==1
                    [No_pks,pk_loc,FWHM(1,:),signal_power(1,:)] = find_peaks_new(data_value.FID_TD_b(:,i),details_para.fres);
        %         end
        %         [ fit_sig,amp,Ind_sig,freq,Z_sig,FWHM] = HSVDPK_new(data_value.FID_TD(:,i),details_para.fres,details_para.Fs,fit_para.order(cur_vox));
                freq=details_para.fres(pk_loc);
                FWHM=1*FWHM*(1E6/details_para.Tf)./2;
                FWHM(1:end)=FWHM(1);
                for kk=1:length(FWHM)%length(pk_loc)
                    Data2bFit=(fftshift(fft(data_value.FID_TD_b(:,i))))/1;
                    PEAK_ppm=4.7 + freq(kk)*(1E6/details_para.Tf);%details_para.ppm_referenced(pk_loc(kk));
                    z=abs(details_para.ppm_referenced-(PEAK_ppm+FWHM(kk)));
                    ub=find(min(z)==z);
                    z=abs(details_para.ppm_referenced-(PEAK_ppm-FWHM(kk)));
                    lb=find(min(z)==z); 
                    freqrange = details_para.ppm_referenced(lb:ub);
                    SumSpec_half = (Data2bFit(lb:ub,:));
                    max_peak=max(SumSpec_half);
                    Init = [max_peak./2 FWHM(kk) PEAK_ppm 0 0 0 ];
                    dataParams_tot(kk,:) = FitPeaks(freqrange, SumSpec_half, Init);
                    dataParams_tot(kk,4:6)=0;
        %             dataParams_tot(kk,2)=0;
                    fitted_signal_tot(kk,:)=LorentzModel(dataParams_tot(kk,:), details_para.ppm_referenced);
                    Ind_sig(:,kk)=fitted_signal_tot(kk,:);
                    Ind_sig(:,kk)=ifft(ifftshift(Ind_sig(:,kk)));
                end
                fit_sig=sum(Ind_sig,2);
        %         fit_sig=ifft(ifftshift(fit_sig));
                ind=find((freq>=details_para.fit_freq(1)-20 & freq<=details_para.fit_freq(2)+20));
                ind_NAA=find((freq>=details_para.fit_freq_NAA(1)-20 & freq<=details_para.fit_freq_NAA(2)+20));
                ind_Cr=find((freq>=details_para.fit_freq_Cr(1)-20 & freq<=details_para.fit_freq_Cr(2)+20));
                ind_Cho=find((freq>=details_para.fit_freq_Cho(1)-20 & freq<=details_para.fit_freq_Cho(2)+20));
                half_sig2=sum(Ind_sig(:,ind),2);
                data_value.fit_TD(:,i)=half_sig2;
                data_value.fit_FD_NAA(:,i)=fftshift(fft(sum(Ind_sig(:,ind_NAA),2)));
                data_value.fit_FD_Cr(:,i)=fftshift(fft(sum(Ind_sig(:,ind_Cr),2)));
                data_value.fit_FD_Cho(:,i)=fftshift(fft(sum(Ind_sig(:,ind_Cho),2)));
                data_value.fit_FD(:,i)=fftshift(fft(data_value.fit_TD(:,i)));
                data_value.peak_area(:,i)=trapz(details_para.ppm_referenced,abs(data_value.fit_FD(:,i)));
                data_value.peak_area_NAA(:,i)=trapz(details_para.ppm_referenced,abs(data_value.fit_FD_NAA(:,i)));
                data_value.peak_area_Cr(:,i)=trapz(details_para.ppm_referenced,abs(data_value.fit_FD_Cr(:,i)));
                data_value.peak_area_Cho(:,i)=trapz(details_para.ppm_referenced,abs(data_value.fit_FD_Cho(:,i)));
                data_value.noise(:,i)=std(data_value.FID_FD(:,i)-fftshift(fft(fit_sig)));
                data_value.snr(:,i)=data_value.peak_area_NAA(:,i)./data_value.noise(:,i);
                data_value.NAA_Cr(:,i)=data_value.peak_area_NAA(:,i)./data_value.peak_area_Cr(:,i);
                data_value.Cho_Cr(:,i)=data_value.peak_area_Cho(:,i)./data_value.peak_area_Cr(:,i);
                data_value.fit_TD(:,i)=data_value.fit_TD(:,i) + ifft(ifftshift(data_value.FID_FD_base(:,i)));
                data_value.fit_FD(:,i)=fftshift(fft(data_value.fit_TD(:,i)));
                waitbar(i/details_para.num_vox);
                details_para.fit_flag(:,i)=1;
            end
        end
    end
catch
        close(h);
        details_para.fit_flag(:,cnt)=0;
        data_value.peak_area(cnt)=0;
        data_value.snr(cnt)=0;
        data_value.NAA_Cr(cnt)=0;
        data_value.Cho_Cr(cnt)=0;
        details_para.se(cnt)=0;
        cur_data(cnt,1)=data_value.peak_area(cnt);
        details_para.mean_peak_area=mean(cur_data(cnt,1));
        details_para.std_peak_area=std(cur_data(cnt,1));
        cur_data(cnt,2)=data_value.snr(cnt);
        cur_data(cnt,3)=data_value.NAA_Cr(cnt);
        cur_data(cnt,4)=data_value.Cho_Cr(cnt);
        cur_data(cnt,5)=details_para.se(cnt);
        set(fitting_handles.fig,'Visible','on');
        return;
end
close(h);
set(fitting_handles.fig,'Visible','on');
cur_data(:,1)=data_value.peak_area;
details_para.mean_peak_area=mean(cur_data(:,1));
details_para.std_peak_area=std(cur_data(:,1));
cur_data(:,2)=data_value.snr;
cur_data(:,3)=data_value.NAA_Cr;
cur_data(:,4)=data_value.Cho_Cr;
cur_data(:,5)=details_para.se;
FID_FD_copy=data_value.FID_FD;
Pix_SS = get(0,'screensize');
x_max = Pix_SS(:,3)./2;
y_max = Pix_SS(:,4).*.75;
col_name = {'Area','SNR','NAA/Cr','Cho/Cr','Error'};
row_name = cell(1,details_para.num_vox);
for i=1:details_para.num_vox
    row_name{i} = strcat('voxel ', num2str(i));
end
diplay_handles.para_table = uitable(fitting_handles.fig,'Units','normalized','position',[370/x_max,290/y_max,320/x_max,130/y_max], 'data',cur_data,...
    'RowName', row_name, 'ColumnName', col_name,'ColumnWidth',{85});
% seg = details_para.disp_seg;
% ppm = details_para.ppm_referenced;
% cla(fitting_handles.axis);
% plot_var = real(data_value.FID_FD(seg(1):seg(2)));
% plot(fitting_handles.axis,ppm(seg(1):seg(2))',plot_var);
% max_p = max(plot_var);
% mim_p = min(plot_var);
% axis(fitting_handles.axis,[ppm(seg(1)), ppm(seg(2)),mim_p - 0.1*max_p, max_p + 0.1*max_p])
% set(fitting_handles.axis,'XDir','reverse')
display_data()
% set(fitting_handles.display_grid , 'Enable', 'on');
%------------------------------------------------------------------------------------------------------------------------------
function display_seg_mim_Callback(hObject, eventdata)

global details_para;

temp = get(hObject,'string');
temp1 = str2double(temp);
if(temp1<details_para.ppm_referenced(1) || temp1>details_para.ppm_referenced(details_para.N) || (temp1>details_para.disp_ppm_seg(2))||(isnan(temp1)))
    errordlg('Invalid value');
    set(hObject,'string',num2str(details_para.disp_ppm_seg(1)))
else
    details_para.disp_ppm_seg(1) = temp1; %in ppm
    details_para.disp_freq_seg(1) = ((details_para.disp_ppm_seg(1) - details_para.ref)*details_para.Tf)/1E6; %ppm to frequency
    details_para.disp_seg(1) = floor((details_para.disp_freq_seg(1)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
    
    display_data();
end
%------------------------------------------------------------------------------------------------------------------------------
function display_seg_max_Callback(hObject, eventdata)
global details_para;

temp = get(hObject,'string');
temp1 = str2double(temp); 
if(temp1<details_para.ppm_referenced(1) || temp1>details_para.ppm_referenced(details_para.N) || (temp1<details_para.disp_ppm_seg(1))||(isnan(temp1)))
    errordlg('Invalid value');
    set(hObject,'string',num2str(details_para.disp_ppm_seg(2)))
else
    details_para.disp_ppm_seg(2) = temp1;%ppm
    details_para.disp_freq_seg(2) = ((details_para.disp_ppm_seg(2) - details_para.ref)*details_para.Tf)/1E6;% ppm to frequency
    details_para.disp_seg(2) = ceil((details_para.disp_freq_seg(2)+ details_para.Fs/2)*details_para.N/details_para.Fs); %frequency to samples
    
    display_data();
end
%------------------------------------------------------------------------------------------------------------------------------------------
