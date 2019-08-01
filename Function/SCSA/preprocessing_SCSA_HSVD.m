function[]= preprocessing_SCSA_HSVD()

% This function is used for pre-processing the data: Phase correction, Apodisation

% Author: Sourav Bhaduri
% University of Ghent
% email: Sourav.Bhaduri@UGent.be
% Feb 2018; Last revision: Feb 2018
addpath ./Function


global data_value;
global preprocess_handles;
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

FID_FD_copy = data_value.FID_FD;
details_para.preprocess_going=1;

if details_para.preprocess_done==0
    default_parameters.pha_0=0;
    default_parameters.pha_1=0;
    default_parameters.apod_para_lor=0;
    default_parameters.apod_para_gaus=0;
    default_parameters.apod_para_sigm=0;
    default_parameters.apod_para_type=1;
    default_parameters.ws_para_SCSA=17;
    default_parameters.ws_para_HLSVD=10;
    default_parameters.ws_para_type=1;
    default_parameters.ds_para_SCSA=1;
    default_parameters.ds_para_HLSVD=3;
    default_parameters.ds_para_type=1;
    ind_allvox=1;
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
    pha_para.temp_pha0=default_parameters.pha_0*ones(details_para.num_vox,1);
    pha_para.temp_pha1=default_parameters.pha_1*ones(details_para.num_vox,1);
    pha_para.pha_0=default_parameters.pha_0*ones(details_para.num_vox,1);
    pha_para.pha_1=default_parameters.pha_1*ones(details_para.num_vox,1);
    apod_para.type = default_parameters.apod_para_type*ones(details_para.num_vox,1);
    apod_para.lor = default_parameters.apod_para_lor*ones(details_para.num_vox,1);
    apod_para.gaus = default_parameters.apod_para_gaus*ones(details_para.num_vox,1);
    apod_para.sigm = default_parameters.apod_para_sigm*ones(details_para.num_vox,1);
    apod_para.temp_lor = default_parameters.apod_para_lor*ones(details_para.num_vox,1);
    apod_para.temp_gaus = default_parameters.apod_para_gaus*ones(details_para.num_vox,1);
    apod_para.temp_sigm = default_parameters.apod_para_sigm*ones(details_para.num_vox,1);
    ws_para.type = default_parameters.ws_para_type*ones(details_para.num_vox,1);
    ws_para.SCSA = default_parameters.ws_para_SCSA*ones(details_para.num_vox,1);
    ws_para.HLSVD = default_parameters.ws_para_HLSVD*ones(details_para.num_vox,1);
    ds_para.type = default_parameters.ds_para_type*ones(details_para.num_vox,1);
    ds_para.SCSA = default_parameters.ds_para_SCSA*ones(details_para.num_vox,1);
    ds_para.HLSVD = default_parameters.ds_para_HLSVD*ones(details_para.num_vox,1);
    pha_cor_data =  data_value.FID_TD;
    hsvd_para.order=10;
    scsa_para.order=17;
if(isempty(apod_para)) % if pre-processing is not done previously load default pre-processing parameters
    
    apod_para.type = default_parameters.apod_para_type*ones(details_para.num_vox,1);
    apod_para.lor = default_parameters.apod_para_lor*ones(details_para.num_vox,1);
    apod_para.gaus = default_parameters.apod_para_gaus*ones(details_para.num_vox,1);
    apod_para.sigm = default_parameters.apod_para_sigm*ones(details_para.num_vox,1);
    ws_para.type = default_parameters.ws_para_type*ones(details_para.num_vox,1);
    ws_para.SCSA = default_parameters.ws_para_SCSA*ones(details_para.num_vox,1);
    ws_para.HLSVD = default_parameters.ws_para_HLSVD*ones(details_para.num_vox,1);
    ds_para.type = default_parameters.ds_para_type*ones(details_para.num_vox,1);
    ds_para.SCSA = default_parameters.ds_para_SCSA*ones(details_para.num_vox,1);
    ds_para.HLSVD = default_parameters.ds_para_HLSVD*ones(details_para.num_vox,1);
    apod_para.function_weight = ones(size(data_value.FID_FD));
    pha_para.pha_0 = default_parameters.pha_0*ones(details_para.num_vox,1);
    pha_para.pha_1 = default_parameters.pha_1*ones(details_para.num_vox,1);
    pha_para.temp_pha0=default_parameters.pha_0*ones(details_para.num_vox,1);
    pha_para.temp_pha1=default_parameters.pha_1*ones(details_para.num_vox,1);
    pha_cor_data =  data_value.FID_TD;
else % if pre-processing is done use the current pre-processing parameters 
    % and perform pre-processing using those values 
    for i = 1:length(details_para.selected_voxels)
        vox = details_para.selected_voxels(1,i);
%         pha_cor_data(:,vox) = phase_correction(FID_FD_copy(:,vox),pha_para.pha_0(vox),pha_para.pha_1(vox),details_para.fres); %phase correction is done on frequency domain non-fftshifted data
%         [temp2, function_weight] = apodize(pha_cor_data(:,vox),apod_para.type(vox),[details_para.Fs,apod_para.lor(vox),apod_para.gaus(vox),apod_para.sigm(vox)]);
%         sht_sig = fftshift(fft(temp2),1);
%         apod_para.function_weight(:, vox) = function_weight;
%         data_value.FID_TD(:,vox) = temp2;
%         data_value.FID_FD(:,vox) = sht_sig;
        apod_para.type=1*ones(details_para.num_vox,1);
        ws_para.type=1*ones(details_para.num_vox,1);
        ds_para.type=1*ones(details_para.num_vox,1);
    end
end
end
Pix_SS = get(0,'screensize');
x_max = Pix_SS(:,3)./2;
y_max = Pix_SS(:,4).*.75;
selected_vox = details_para.selected_voxels;
vox_no = cell(1,length(selected_vox));
for i = 1:length(selected_vox)
vox_no{i} =  num2str(selected_vox(1,i));
end

% Cerate GUI
preprocess_handles.fig = figure('Visible','off','Position',[0,35,x_max,y_max./1],'Color',[0.68,0.92,.8],'Name','Preprocessing window','NumberTitle', 'off');
set(preprocess_handles.fig,'MenuBar','none'); 
% undecorateFig(preprocess_handles.fig);
preprocess_handles.reload = uimenu(preprocess_handles.fig,'Label','Reload','Callback',@reload_callback);
% preprocess_handles.ws_tool = uimenu(preprocess_handles.tool,'Label','Water suppression');
% preprocess_handles.denoise_tool = uimenu(preprocess_handles.tool,'Label','Denoise');
% preprocess_handles.scsa_ws_tool = uimenu(preprocess_handles.ws_tool,'Label','SCSA','Callback',@SCSA_ws_callback);
% preprocess_handles.hsvd_ws_tool = uimenu(preprocess_handles.ws_tool,'Label','HLSVD-PRO','Callback',@HSVD_ws_callback);
% preprocess_handles.scsa_ws_tool = uimenu(preprocess_handles.denoise_tool,'Label','SCSA','Callback',@SCSA_denoise_callback);
% preprocess_handles.hsvd_ws_tool = uimenu(preprocess_handles.denoise_tool,'Label','SVD','Callback',@HSVD_denoise_callback);
% preprocess_handles.save_figure = uimenu(preprocess_handles.fig,'Label','Save Figure','Callback',@save_preprocess_fig_callback);


preprocess_handles.axis1 = axes('Units','Pixels','Units','normalized','Position',[20/x_max 2/y_max 140/x_max 80/y_max]);
% display_handles.axis_text = uicontrol(display_handles.fig,'Style','text','Units','normalized','Position',[220/x_max 460/y_max 130/x_max 15/y_max],...
%     'String','Signal Spectrum','BackgroundColor',[0.68,0.92,1],'FontSize',9,'FontWeight','Bold');
preprocess_handles.axis2 = axes('Units','Pixels','Units','normalized','Position',[150/x_max 2/y_max 150/x_max 80/y_max]);
% 
myImage=imread([pwd,'\Function','\ugent_logo.png']);
myImage= imresize(myImage,1);
imshow(myImage,'Parent',preprocess_handles.axis1);
myImage=imread([pwd,'\Function','\kaust_logo.jpg']);
myImage= imresize(myImage,1);
imshow(myImage,'Parent',preprocess_handles.axis2);

preprocess_handles.h1 = uipanel(preprocess_handles.fig,'Title','Voxel No','Units','normalized','Position',[0.91 .58 .08 .4],...
    'BackgroundColor','yellow','FontSize',8, 'FontWeight','Bold');
% preprocess_handles.voxel_list = uicontrol(preprocess_handles.h1,'Style', 'listbox','Units','normalized','Position', [.05 .027 .9 .97],'String', vox_no,...
%     'Max', details_para.num_vox, 'Min', 1, 'Callback', @voxel_list_callback, 'Interruptible','off','BusyAction','cancel'); 
preprocess_handles.voxel_list = uicontrol(preprocess_handles.h1,'Style', 'listbox','Units','normalized','Position', [.05 .027 .9 .97],'String', vox_no,...
    'value',1, 'Callback', @voxel_list_callback, 'Interruptible','off','BusyAction','cancel'); 


preprocess_handles.display_panel = uipanel(preprocess_handles.fig,'Title','Display Properties','Position',[.03 .58 .87 .4],'BackgroundColor','yellow','FontSize',9,'FontWeight','Bold');
preprocess_handles.seg_min_text = uicontrol(preprocess_handles.display_panel,'Style','text','Units','normalized','Position',[0.4,0.4,0.2,0.1*2.5],'String',...
    'Min','FontSize',8,'FontWeight','Bold','BackgroundColor','yellow');
preprocess_handles.seg_max_text = uicontrol(preprocess_handles.display_panel,'Style','text','Units','normalized','Position',[0.7,0.4,0.2,0.1*2.5],'String',...
    'Max','FontSize',8,'FontWeight','Bold','BackgroundColor','yellow');
preprocess_handles.display_seg_mim = uicontrol(preprocess_handles.display_panel,'Style','edit','Units','normalized','Position',[0.45,0.2,0.15,0.15*2],'String','0',...
    'Callback',{@display_seg_mim_Callback},'Interruptible','off','BusyAction','cancel'); 
preprocess_handles.display_seg_max = uicontrol(preprocess_handles.display_panel,'Style','edit','Units','normalized','Position',[0.75,0.2,0.15,0.15*2],'String','0',...
    'Callback',{@display_seg_max_Callback},'Interruptible','off','BusyAction','cancel'); 
% preprocess_handles.zerofill_text = uicontrol(preprocess_handles.fig,'Style','text','Units','normalized','Position',[557/x_max 460/y_max 70/x_max 15/y_max],...
%     'String','Zerofill','BackgroundColor',[0.68,0.92,1],'FontSize',9,'FontWeight','Bold');
% preprocess_handles.zerofill_options = uicontrol(preprocess_handles.fig,'Style','popup','String','No Fill|1 times|2 times|3 times','Units','normalized',...
%     'Position',[555/x_max 405/y_max 80/x_max 15/y_max],'Value',zerofill_value,'Interruptible','off','BusyAction','cancel','Callback', @zerofill_options_callback); 

preprocess_handles.TF_FD_display = uibuttongroup(preprocess_handles.display_panel,'Units','normalized','Position',[0.04 0.76 .7 0.16],...
    'FontSize',9,'BackgroundColor',[0.68,0.92,.8],'BorderType','none','SelectionChangeFcn',@TD_FD_display_callback);
preprocess_handles.FD_disp_real = uicontrol(preprocess_handles.TF_FD_display,'Style','radiobutton','Units','normalized','Position',[.01 .1 .2 .9],'String',...
    'FD (real)','FontWeight','Bold','BackgroundColor',[0.68,0.92,.8]);
preprocess_handles.FD_disp_abs = uicontrol(preprocess_handles.TF_FD_display,'Style','radiobutton','Units','normalized','Position',[.3 .1 .2 .9],'String',...
    'FD (abs)','FontWeight','Bold','BackgroundColor',[0.68,0.92,.8]);
preprocess_handles.TD_disp = uicontrol(preprocess_handles.TF_FD_display,'Style','radiobutton','Units','normalized','String','FD (imag)','Position',...
    [.53 .1 .2 .9],'FontWeight','Bold','BackgroundColor',[0.68,0.92,.8]);
preprocess_handles.FD_disp_imag = uicontrol(preprocess_handles.TF_FD_display,'Style','radiobutton','Units','normalized','String','TD','Position',...
    [.75 .1 .2 .9],'FontWeight','Bold','BackgroundColor',[0.68,0.92,.8]);

% apodization panel
preprocess_handles.h2 = uipanel(preprocess_handles.fig,'Title','Apodization','Units','normalized','Position',[380/x_max 40/y_max (x_max-360)/(3.2*x_max) 320/y_max],...
    'BackgroundColor','yellow','FontSize',9, 'FontWeight','Bold');
preprocess_handles.Apodize = uicontrol(preprocess_handles.h2,'Style','text','Units','normalized','String','Apodization Type','Units','normalized',...
    'Position',[0.1 0.65 0.8 0.3],'BackgroundColor','yellow','FontSize',8,'FontWeight','Bold');
preprocess_handles.Apodize_type = uicontrol(preprocess_handles.h2,'Style','popup','String','Exponential|Gauss|Gauss-Exponential|sigmoid','Units','normalized','Position',...
    [0.1 0.59 0.8 0.3],'Value',apod_para.type(cur_vox), 'Callback',{@Apodize_type_Callback});

preprocess_handles.Apodize_Lor_text = uicontrol(preprocess_handles.h2, 'Style','text','Units','normalized','Position',[0.1 0.7 0.8 0.07],'String',...
    'order','FontSize',8,'FontWeight','Bold','BackgroundColor','yellow');
preprocess_handles.Apodize_Lor_edit = uicontrol(preprocess_handles.h2,'Style','edit','Units','normalized','Position',[0.1 0.65 0.8 0.07],'String',...
    num2str(apod_para.lor(cur_vox)),'Callback',{@lor_Callback}); 
preprocess_handles.Apodize_Lor_slide = uicontrol(preprocess_handles.h2,'Style', 'slider','Min',0,'Max',100,'Value',apod_para.lor(cur_vox),'Units','normalized',...
    'Position', [0.1 0.58 0.8 0.07],'SliderStep',[0.001,0.01],'Callback', {@slider_lor_Callback}); 
%     
preprocess_handles.Apodize_Gaus_text = uicontrol(preprocess_handles.h2,'Style','text','Units','normalized','Position',[0.1 0.45 0.8 0.07],'String',....
    'order','FontSize',8,'FontWeight','Bold','BackgroundColor','yellow');
preprocess_handles.Apodize_Gaus_edit = uicontrol(preprocess_handles.h2,'Style','edit','Units','normalized','Position',[0.1 0.4 0.8 0.07],...
    'String',num2str(apod_para.gaus(cur_vox)),'Callback',{@gaus_Callback});  
preprocess_handles.Apodize_Gaus_slide  = uicontrol(preprocess_handles.h2,'Style', 'slider','Min',0,'Max',1000,'Value',apod_para.gaus(cur_vox),'Units','normalized',...
    'Position', [0.1 0.33 0.8 0.07],'SliderStep',[0.0001,0.001],'Callback', {@slider_gaus_Callback});  

preprocess_handles.Apodize_Sig_text = uicontrol(preprocess_handles.h2,'Style','text','Units','normalized','Position',[0.1 0.2 0.8 0.07],'String',...
    'order','FontSize',8,'FontWeight','Bold','BackgroundColor','yellow');
preprocess_handles.Apodize_Sigm_edit = uicontrol(preprocess_handles.h2,'Style','edit','Units','normalized','Position',[0.1 0.15 0.8 0.07],...
    'String',num2str(apod_para.sigm(cur_vox)),'Callback',{@sigm_Callback});  
preprocess_handles.Apodize_Sigm_slide = uicontrol(preprocess_handles.h2,'Style', 'slider','Min',0,'Max',500,'Value',apod_para.sigm(cur_vox),'Units','normalized',...
    'Position', [0.1 0.08 0.8 0.07],'SliderStep',[0.001,0.01],'Callback', {@slider_sigm_Callback});

% ws panel
preprocess_handles.ws = uipanel(preprocess_handles.fig,'Title','Water suppression','Units','normalized','Position',[(380/x_max)+(x_max-360)/(3.2*x_max) 40/y_max (x_max-360)/(3.2*x_max) 320/y_max],...
    'BackgroundColor','yellow','FontSize',9, 'FontWeight','Bold');
preprocess_handles.ws_orient = uicontrol(preprocess_handles.ws,'Style','text','Units','normalized','String','WS Type','Units','normalized',...
    'Position',[0.1 0.65 0.8 0.3],'BackgroundColor','yellow','FontSize',8,'FontWeight','Bold');
preprocess_handles.ws_type = uicontrol(preprocess_handles.ws,'Style','popup','String','SCSA|HLSVD-PRO','Units','normalized','Position',...
    [0.1 0.59 0.8 0.3],'Value',ws_para.type(cur_vox), 'Callback',{@ws_type_Callback});

preprocess_handles.ws_SCSA_text = uicontrol(preprocess_handles.ws, 'Style','text','Units','normalized','Position',[0.1 0.7 0.8 0.07],'String',...
    'order','FontSize',8,'FontWeight','Bold','BackgroundColor','yellow');
preprocess_handles.ws_SCSA_edit = uicontrol(preprocess_handles.ws,'Style','edit','Units','normalized','Position',[0.1 0.65 0.8 0.07],'String',...
    num2str(ws_para.SCSA(cur_vox)),'Callback',{@SCSA_ws_Callback}); 
preprocess_handles.ws_SCSA_slide = uicontrol(preprocess_handles.ws,'Style', 'slider','Min',0,'Max',30,'Value',apod_para.lor(cur_vox),'Units','normalized',...
    'Position', [0.1 0.58 0.8 0.07],'SliderStep',[0.1,1],'Callback', {@slider_SCSA_Callback}); 
%     
preprocess_handles.ws_HLSVD_text = uicontrol(preprocess_handles.ws,'Style','text','Units','normalized','Position',[0.1 0.45 0.8 0.07],'String',....
    'order','FontSize',8,'FontWeight','Bold','BackgroundColor','yellow');
preprocess_handles.ws_HLSVD_edit = uicontrol(preprocess_handles.ws,'Style','edit','Units','normalized','Position',[0.1 0.4 0.8 0.07],...
    'String',num2str(ws_para.HLSVD(cur_vox)),'Callback',{@HLSVD_ws_Callback});  
preprocess_handles.ws_HLSVD_slide  = uicontrol(preprocess_handles.ws,'Style', 'slider','Min',0,'Max',30,'Value',apod_para.gaus(cur_vox),'Units','normalized',...
    'Position', [0.1 0.33 0.8 0.07],'SliderStep',[0.1,1],'Callback', {@slider_HLSVD_Callback}); 

% denoise panel
preprocess_handles.ds = uipanel(preprocess_handles.fig,'Title','Denoise','Units','normalized','Position',[(380/x_max)+(2*(x_max-360))/(3.2*x_max) 40/y_max (x_max-360)/(3.2*x_max) 320/y_max],...
    'BackgroundColor','yellow','FontSize',9, 'FontWeight','Bold');
preprocess_handles.ds_orient = uicontrol(preprocess_handles.ds,'Style','text','Units','normalized','String','Type','Units','normalized',...
    'Position',[0.1 0.65 0.8 0.3],'BackgroundColor','yellow','FontSize',8,'FontWeight','Bold');
preprocess_handles.ds_type = uicontrol(preprocess_handles.ds,'Style','popup','String','SCSA|SVD','Units','normalized','Position',...
    [0.1 0.59 0.8 0.3],'Value',ds_para.type(cur_vox), 'Callback',{@ds_type_Callback});

preprocess_handles.ds_SCSA_text = uicontrol(preprocess_handles.ds, 'Style','text','Units','normalized','Position',[0.1 0.7 0.8 0.07],'String',...
    'level','FontSize',8,'FontWeight','Bold','BackgroundColor','yellow');
preprocess_handles.ds_SCSA_edit = uicontrol(preprocess_handles.ds,'Style','edit','Units','normalized','Position',[0.1 0.65 0.8 0.07],'String',...
    num2str(ds_para.SCSA(cur_vox)),'Callback',{@SCSA_ds_Callback}); 
preprocess_handles.ds_SCSA_slide = uicontrol(preprocess_handles.ds,'Style', 'slider','Min',0.2,'Max',3,'Value',apod_para.lor(cur_vox),'Units','normalized',...
    'Position', [0.1 0.58 0.8 0.07],'SliderStep',[0.05,1],'Callback', {@slider_SCSA_ds_Callback}); 
%     
preprocess_handles.ds_HLSVD_text = uicontrol(preprocess_handles.ds,'Style','text','Units','normalized','Position',[0.1 0.45 0.8 0.07],'String',....
    'level','FontSize',8,'FontWeight','Bold','BackgroundColor','yellow');
preprocess_handles.ds_HLSVD_edit = uicontrol(preprocess_handles.ds,'Style','edit','Units','normalized','Position',[0.1 0.4 0.8 0.07],...
    'String',num2str(ds_para.HLSVD(cur_vox)),'Callback',{@HLSVD_ds_Callback});  
preprocess_handles.ds_HLSVD_slide  = uicontrol(preprocess_handles.ds,'Style', 'slider','Min',0,'Max',10,'Value',apod_para.gaus(cur_vox),'Units','normalized',...
    'Position', [0.1 0.33 0.8 0.07],'SliderStep',[0.1,1],'Callback', {@slider_HLSVD_ds_Callback});  

% phase correction panel
preprocess_handles.h3 = uipanel(preprocess_handles.fig,'Title','Phase correction','Units','normalized','Position',[20/x_max 160/y_max 360/x_max 200/y_max],...
    'BackgroundColor','yellow','FontSize',9, 'FontWeight','Bold');
preprocess_handles.pha_corr0_text = uicontrol(preprocess_handles.h3,'Style','text','Units','normalized','Position',[0.1 0.6 0.2 .2],'String',...
    'Zero-order phase correction ','FontSize',7.5,'FontWeight','Bold','BackgroundColor','yellow');
preprocess_handles.pha_corr1_text = uicontrol(preprocess_handles.h3,'Style','text','Units','normalized','Position',[0.1 0.3 0.2 0.2],'String',...
    'First-order phase correction','FontSize',7.5,'FontWeight','Bold','BackgroundColor','yellow');
preprocess_handles.pha_corr_value_text = uicontrol(preprocess_handles.h3,'Style','text','Units','normalized','Position',[0.31 0.61 0.2 .2],'String',...
    'Value','FontSize',8,'FontWeight','Bold','BackgroundColor','yellow');
preprocess_handles.pha_corr_slider_text = uicontrol(preprocess_handles.h3,'Style','text','Units','normalized','Position',[0.66 0.61 0.2 0.2],'String',...
    'Slider','FontSize',8,'FontWeight','Bold','BackgroundColor','yellow');
preprocess_handles.pha_corr0_value = uicontrol(preprocess_handles.h3,'Style','edit','Units','normalized','Position',[0.34 0.5 .15 .2],...
    'String',num2str(pha_para.pha_0(cur_vox)),'Callback',{@pha0_val_Callback});
preprocess_handles.pha_corr1_value = uicontrol(preprocess_handles.h3,'Style','edit','Units','normalized','Position',[0.34 0.2 .15 .2],...
    'String',num2str(pha_para.pha_1(cur_vox)),'Callback',{@pha1_val_Callback});
preprocess_handles.pha_corr0_slide = uicontrol(preprocess_handles.h3,'Style','slider','Min',-pi,'Max',pi,'Value',pha_para.pha_0(cur_vox),'Units','normalized','Position',...
    [0.55 0.5 .42 .2],'SliderStep',[0.001,0.01],'Callback', {@pha0_slider_Callback});  
preprocess_handles.pha_corr1_slide = uicontrol(preprocess_handles.h3,'Style','slider','Min',-10,'Max',10,'Value',pha_para.pha_1(cur_vox),'Units','normalized','Position',...
    [0.55 0.2 .42 .2],'SliderStep',[0.001,0.01],'Callback', {@pha1_slider_Callback});

preprocess_handles.developer_display = uipanel(preprocess_handles.fig,'Title','Developers:','Units','normalized','Position',[20/x_max 80/y_max 360/x_max 80/y_max],...
    'BackgroundColor','yellow','FontSize',9, 'FontWeight','Bold');
preprocess_handles.developer_text = uicontrol(preprocess_handles.developer_display,'Style','text','Units','normalized','Position',[0.01 0.01 0.8 0.8],'String',...
    'Sourav Bhaduri, Abderrazak Chahid','FontSize',8,'FontWeight','Bold','BackgroundColor','yellow');

preprocess_handles.close = uicontrol(preprocess_handles.fig,'Style', 'pushbutton','Units','normalized','Position',[((380/x_max)+(2*(x_max-360))/(3.0*x_max)) 10/y_max 100/x_max 25/y_max],...
'String', 'OK/Close','FontSize',9,'FontWeight','Bold', 'Interruptible','off','BusyAction','cancel','Callback',@close_Callback);
% preprocess_handles.met_map = uicontrol(preprocess_handles.fig,'Style', 'pushbutton','Units','normalized','Position',[560/x_max 10/y_max 100/x_max 25/y_max],...
% 'String', 'metabolite map','FontSize',9,'FontWeight','Bold', 'Interruptible','off','BusyAction','cancel','Callback',@metabolite_map_Callback);

preprocess_handles.buttongroup = uibuttongroup(preprocess_handles.fig,'Units','normalized','Position',[32/x_max .6 130/x_max .1],...
    'BackgroundColor','yellow','FontSize',9, 'FontWeight','Bold');%,'SelectionChangeFcn',@ind_allvox_value);

preprocess_handles.all_vox = uicontrol('parent',preprocess_handles.buttongroup,'Style','radiobutton','String','All Voxel',...
    'Units','normalized','Position',[0.1 0.2 0.7 0.3],'BackgroundColor','yellow');
preprocess_handles.ind_vox = uicontrol('parent',preprocess_handles.buttongroup,'Style','radiobutton','String','Individual Voxel',...
    'Units','normalized','Position',[0.1 0.6 0.8 0.3],'BackgroundColor','yellow');
% preprocess_handles.water_suppress = uicontrol(preprocess_handles.fig,'Style', 'pushbutton','Units','normalized','Position',[240/x_max 150/y_max 100/x_max 25/y_max],...
% 'String', 'Water suppress','FontSize',9,'FontWeight','Bold', 'Interruptible','off','BusyAction','cancel','Callback',@water_suppress_Callback);

preprocess_handles.water_suppress = uicontrol(preprocess_handles.fig,'Style', 'pushbutton','Units','normalized','Position',[((380/x_max)+(1*(x_max-360))/(3.0*x_max)) 70/y_max 80/x_max 25/y_max],...
'String', 'Water suppress','FontSize',7.5,'FontWeight','Bold', 'Interruptible','off','BusyAction','cancel','Callback',@water_suppress_Callback);
preprocess_handles.denoise = uicontrol(preprocess_handles.fig,'Style', 'pushbutton','Units','normalized','Position',[((380/x_max)+(2*(x_max-360))/(3.0*x_max)) 70/y_max 80/x_max 25/y_max],...
'String', 'Denoise','FontSize',7.5,'FontWeight','Bold', 'Interruptible','off','BusyAction','cancel','Callback',@denoise_Callback);
% preprocess_handles.wavelet = uicontrol(preprocess_handles.fig,'Style', 'pushbutton','Units','normalized','Position',[240/x_max 140/y_max 100/x_max 25/y_max],...
% 'String', 'Denoise','FontSize',9,'FontWeight','Bold', 'Interruptible','off','BusyAction','cancel','Callback',@denoise_Callback);
set(preprocess_handles.buttongroup,'SelectionChangeFcn',@ind_allvox_value);



if ind_allvox == 1
    set(preprocess_handles.buttongroup,'SelectedObject',preprocess_handles.ind_vox);
else
    set(preprocess_handles.buttongroup,'SelectedObject',preprocess_handles.all_vox);
end

set(preprocess_handles.ws_SCSA_slide,'value',ws_para.SCSA(cur_vox));
set(preprocess_handles.ws_HLSVD_slide,'value',ws_para.HLSVD(cur_vox));

set(preprocess_handles.ds_SCSA_slide,'value',ds_para.SCSA(cur_vox));
set(preprocess_handles.ds_HLSVD_slide,'value',ds_para.HLSVD(cur_vox));

if (apod_para.type(selected_vox(1,1)) == 1)
    set(preprocess_handles.Apodize_Lor_edit, 'Enable', 'on');
    set(preprocess_handles.Apodize_Lor_slide, 'Enable', 'on');
    set(preprocess_handles.Apodize_Sigm_edit, 'Enable', 'off');
    set(preprocess_handles.Apodize_Sigm_slide, 'Enable', 'off');
    set(preprocess_handles.Apodize_Gaus_edit, 'Enable', 'off');
    set(preprocess_handles.Apodize_Gaus_slide, 'Enable', 'off');
elseif (apod_para.type(selected_vox(1,1)) == 2)
    set(preprocess_handles.Apodize_Lor_edit, 'Enable', 'off');
    set(preprocess_handles.Apodize_Sigm_edit, 'Enable', 'off');
    set(preprocess_handles.Apodize_Lor_slide, 'Enable', 'off');
    set(preprocess_handles.Apodize_Sigm_slide, 'Enable', 'off');
    set(preprocess_handles.Apodize_Gaus_edit, 'Enable', 'on');
    set(preprocess_handles.Apodize_Gaus_slide, 'Enable', 'on');
elseif (apod_para.type(selected_vox(1,1)) == 3)
    set(preprocess_handles.Apodize_Lor_edit, 'Enable', 'on');
    set(preprocess_handles.Apodize_Sigm_edit, 'Enable', 'off');
    set(preprocess_handles.Apodize_Lor_slide, 'Enable', 'on');
    set(preprocess_handles.Apodize_Sigm_slide, 'Enable', 'off');
    set(preprocess_handles.Apodize_Gaus_edit, 'Enable', 'on');
    set(preprocess_handles.Apodize_Gaus_slide, 'Enable', 'on');
else
    set(preprocess_handles.Apodize_Lor_edit, 'Enable', 'off');
    set(preprocess_handles.Apodize_Sigm_edit, 'Enable', 'on');
    set(preprocess_handles.Apodize_Lor_slide, 'Enable', 'off');
    set(preprocess_handles.Apodize_Sigm_slide, 'Enable', 'on');
    set(preprocess_handles.Apodize_Gaus_edit, 'Enable', 'off');
    set(preprocess_handles.Apodize_Gaus_slide, 'Enable', 'off');
end

if (ws_para.type(selected_vox(1,1)) == 1)
    set(preprocess_handles.ws_SCSA_edit, 'Enable', 'on');
    set(preprocess_handles.ws_SCSA_slide, 'Enable', 'on');
    set(preprocess_handles.ws_HLSVD_edit, 'Enable', 'off');
    set(preprocess_handles.ws_HLSVD_slide, 'Enable', 'off');
elseif(ws_para.type(selected_vox(1,1)) == 2)
    set(preprocess_handles.ws_SCSA_edit, 'Enable', 'off');
    set(preprocess_handles.ws_SCSA_slide, 'Enable', 'off');
    set(preprocess_handles.ws_HLSVD_edit, 'Enable', 'on');
    set(preprocess_handles.ws_HLSVD_slide, 'Enable', 'on');
end

if (ds_para.type(selected_vox(1,1)) == 1)
    set(preprocess_handles.ds_SCSA_edit, 'Enable', 'on');
    set(preprocess_handles.ds_SCSA_slide, 'Enable', 'on');
    set(preprocess_handles.ds_HLSVD_edit, 'Enable', 'off');
    set(preprocess_handles.ds_HLSVD_slide, 'Enable', 'off');
elseif(ws_para.type(selected_vox(1,1)) == 2)
    set(preprocess_handles.ds_SCSA_edit, 'Enable', 'off');
    set(preprocess_handles.ds_SCSA_slide, 'Enable', 'off');
    set(preprocess_handles.ds_HLSVD_edit, 'Enable', 'on');
    set(preprocess_handles.ds_HLSVD_slide, 'Enable', 'on');
end
    

set(preprocess_handles.fig,'Visible','on');
set(preprocess_handles.display_seg_mim, 'Enable', 'on','string',num2str(details_para.disp_ppm_seg(1)))
set(preprocess_handles.display_seg_max, 'Enable', 'on','string',num2str(details_para.disp_ppm_seg(2)))

set(preprocess_handles.fig,'CloseRequestFcn',@close_Callback);
 
preprocess_display = 1;
% Display the selected voxel 
display_SCSA_HSVD();
display_data()
waitfor(preprocess_handles.fig);

function voxel_list_callback(hObject, eventdata)
% Callback for 'voxel no.' get the selected voxel and update the plot according to the voxel selected
global details_para;
global preprocess_handles;
global apod_para;
global pha_para;
global cur_vox;
global ws_para;
global ds_para;

% get the selected voxel
preprocess_sel_vox = get(hObject,'Value');

if(isempty(preprocess_sel_vox))
    preprocess_sel_vox = 1;
    set(preprocess_handles.voxel_list,'value',1);
end

selected_vox = details_para.selected_voxels;
cur_vox = selected_vox(1,preprocess_sel_vox);
display_data()


% Update the GUI based on the selected voxel
temp = apod_para.type(cur_vox);
temp1 = ws_para.type(cur_vox);
temp2 = ds_para.type(cur_vox);
set(preprocess_handles.Apodize_type,'value',temp);
set(preprocess_handles.ws_type,'value',temp1);
set(preprocess_handles.ds_type,'value',temp2);
switch(temp)
  case 1
        set(preprocess_handles.Apodize_Lor_edit, 'Enable', 'on','string',num2str(apod_para.lor(cur_vox)));
        set(preprocess_handles.Apodize_Gaus_edit, 'Enable', 'off');
        set(preprocess_handles.Apodize_Sigm_edit, 'Enable', 'off');
        set(preprocess_handles.Apodize_Lor_slide, 'Enable', 'on','value',apod_para.lor(cur_vox));
        set(preprocess_handles.Apodize_Sigm_slide, 'Enable', 'off');
        set(preprocess_handles.Apodize_Gaus_slide, 'Enable', 'off');
        
  case 2
        set(preprocess_handles.Apodize_Lor_edit, 'Enable', 'off');
        set(preprocess_handles.Apodize_Gaus_edit, 'Enable', 'on','string',num2str(apod_para.gaus(cur_vox)));
        set(preprocess_handles.Apodize_Sigm_edit, 'Enable', 'off');
        set(preprocess_handles.Apodize_Lor_slide, 'Enable', 'off');
        set(preprocess_handles.Apodize_Gaus_slide, 'Enable', 'on','value',apod_para.gaus(cur_vox));
        set(preprocess_handles.Apodize_Sigm_slide, 'Enable', 'off');
  case 3
        set(preprocess_handles.Apodize_Lor_edit, 'Enable', 'on','string',num2str(apod_para.lor(cur_vox)));
        set(preprocess_handles.Apodize_Gaus_edit, 'Enable', 'on','string',num2str(apod_para.gaus(cur_vox)));
        set(preprocess_handles.Apodize_Sigm_edit, 'Enable', 'off');
        set(preprocess_handles.Apodize_Lor_slide, 'Enable', 'on','value',apod_para.lor(cur_vox));
        set(preprocess_handles.Apodize_Gaus_slide, 'Enable', 'on','value',apod_para.gaus(cur_vox));
        set(preprocess_handles.Apodize_Sigm_slide, 'Enable', 'off');
  case 4
        set(preprocess_handles.Apodize_Lor_edit, 'Enable', 'off');
        set(preprocess_handles.Apodize_Gaus_edit, 'Enable', 'off');
        set(preprocess_handles.Apodize_Sigm_edit, 'Enable', 'on','string',num2str(apod_para.sigm(cur_vox)));
        set(preprocess_handles.Apodize_Lor_slide, 'Enable', 'off');
        set(preprocess_handles.Apodize_Gaus_slide, 'Enable', 'off');
        set(preprocess_handles.Apodize_Sigm_slide, 'Enable', 'on','value',apod_para.sigm(cur_vox));
end
switch(temp1)
   case 1 %SCSA
        set(preprocess_handles.ws_SCSA_edit, 'Enable', 'on');
        set(preprocess_handles.ws_HLSVD_edit, 'Enable', 'off');
        set(preprocess_handles.ws_SCSA_slide, 'Enable', 'on');
        set(preprocess_handles.ws_HLSVD_slide, 'Enable', 'off');
  case 2 %HLSVD
        set(preprocess_handles.ws_SCSA_edit, 'Enable', 'off');
        set(preprocess_handles.ws_HLSVD_edit, 'Enable', 'on');
        set(preprocess_handles.ws_SCSA_slide, 'Enable', 'off');
        set(preprocess_handles.ws_HLSVD_slide, 'Enable', 'on');
end
switch(temp2)
  case 1 %SCSA
        set(preprocess_handles.ds_SCSA_edit, 'Enable', 'on');
        set(preprocess_handles.ds_HLSVD_edit, 'Enable', 'off');
        set(preprocess_handles.ds_SCSA_slide, 'Enable', 'on');
        set(preprocess_handles.ds_HLSVD_slide, 'Enable', 'off');
  case 2 %HLSVD
        set(preprocess_handles.ds_SCSA_edit, 'Enable', 'off');
        set(preprocess_handles.ds_HLSVD_edit, 'Enable', 'on');
        set(preprocess_handles.ds_SCSA_slide, 'Enable', 'off');
        set(preprocess_handles.ds_HLSVD_slide, 'Enable', 'on');
end

set(preprocess_handles.pha_corr0_slide,'value',pha_para.pha_0(cur_vox));
set(preprocess_handles.pha_corr1_slide,'value',pha_para.pha_1(cur_vox));
set(preprocess_handles.pha_corr0_value,'string',num2str(pha_para.pha_0(cur_vox)));
set(preprocess_handles.pha_corr1_value,'string',num2str(pha_para.pha_1(cur_vox)));
set(preprocess_handles.ds_HLSVD_slide,'value',ds_para.HLSVD(cur_vox));
set(preprocess_handles.ds_SCSA_slide,'value',ds_para.SCSA(cur_vox));
set(preprocess_handles.ws_HLSVD_slide,'value',ws_para.HLSVD(cur_vox));
set(preprocess_handles.ws_SCSA_slide,'value',ws_para.SCSA(cur_vox));
set(preprocess_handles.ds_HLSVD_edit,'string',num2str(ds_para.HLSVD(cur_vox)));
set(preprocess_handles.ds_SCSA_edit,'string',num2str(ds_para.SCSA(cur_vox)));
set(preprocess_handles.ws_HLSVD_edit,'string',num2str(ws_para.HLSVD(cur_vox)));
set(preprocess_handles.ws_SCSA_edit,'string',num2str(ws_para.SCSA(cur_vox)));
% set(preprocess_handles.display_seg_mim, 'Enable', 'on','string',num2str(details_para.disp_ppm_seg(1)))
% set(preprocess_handles.display_seg_max, 'Enable', 'on','string',num2str(details_para.disp_ppm_seg(2)))


function zerofill_options_callback(hObject, eventdata)
%Callback for 'Zero-fill' drop down box
global zerofill_value;
global details_para;

% update the selected zero-fill value, But dont zero fill the signal.
% Zero-filling is done once the pre-processing GUI is closed
zerofill_value = get(hObject,'Value');
if(zerofill_value ~= 1)
    [a,b] = size(data_value.FID_TD);
    data_value.FID_TD = [data_value.FID_TD;zeros((zerofill_value-1)*a,b)]; %add zeros in the end
    data_value.FID_FD = fftshift(fft(data_value.FID_TD),1);
    
    %update parameters after zerofilling
    N = (zerofill_value)*(a); 
    n = 0:N-1;
    fres = (-details_para.Fs/2 + (details_para.Fs/N)*n);
    ref = details_para.ref;
    ppm = (fres)*(1E6/details_para.Tf);
    t = n/details_para.Fs;
    details_para.N = N;
    details_para.ppm = ppm;
    details_para.fres = fres;
    details_para.t = t;
    details_para.ppm_referenced = details_para.ref + ppm;
end
display_data();



function Apodize_type_Callback(hObject, eventdata)
%Callback for 'Apodisation type' 
global preprocess_handles;
global apod_para;
global ind_allvox;
global details_para;
global cur_vox;
% global data_value;

% get the selected apodization type
temp = get(hObject,'Value');
if(ind_allvox == 1)
    apod_para.type(cur_vox) = temp;
else
    apod_para.type(details_para.selected_voxels(1,:)) = temp;
end

% enable/disable based on the selected type
switch(temp)
  case 1 %Lorentz
        set(preprocess_handles.Apodize_Lor_edit, 'Enable', 'on');
        set(preprocess_handles.Apodize_Gaus_edit, 'Enable', 'off');
        set(preprocess_handles.Apodize_Sigm_edit, 'Enable', 'off');
        set(preprocess_handles.Apodize_Lor_slide, 'Enable', 'on');
        set(preprocess_handles.Apodize_Sigm_slide, 'Enable', 'off');
        set(preprocess_handles.Apodize_Gaus_slide, 'Enable', 'off');
  case 2 %Gaussian
        set(preprocess_handles.Apodize_Lor_edit, 'Enable', 'off');
        set(preprocess_handles.Apodize_Gaus_edit, 'Enable', 'on');
        set(preprocess_handles.Apodize_Sigm_edit, 'Enable', 'off');
        set(preprocess_handles.Apodize_Lor_slide, 'Enable', 'off');
        set(preprocess_handles.Apodize_Gaus_slide, 'Enable', 'on');
        set(preprocess_handles.Apodize_Sigm_slide, 'Enable', 'off');
  case 3 %Lor-Gauss
        set(preprocess_handles.Apodize_Lor_edit, 'Enable', 'on');
        set(preprocess_handles.Apodize_Gaus_edit, 'Enable', 'on');
        set(preprocess_handles.Apodize_Sigm_edit, 'Enable', 'off');
        set(preprocess_handles.Apodize_Lor_slide, 'Enable', 'on');
        set(preprocess_handles.Apodize_Gaus_slide, 'Enable', 'on');
        set(preprocess_handles.Apodize_Sigm_slide, 'Enable', 'off');
  case 4 %sigmoid
        set(preprocess_handles.Apodize_Lor_edit, 'Enable', 'off');
        set(preprocess_handles.Apodize_Gaus_edit, 'Enable', 'off');
        set(preprocess_handles.Apodize_Sigm_edit, 'Enable', 'on');
        set(preprocess_handles.Apodize_Lor_slide, 'Enable', 'off');
        set(preprocess_handles.Apodize_Gaus_slide, 'Enable', 'off');
        set(preprocess_handles.Apodize_Sigm_slide, 'Enable', 'on');
end
do_and_display_apod()% apodize and display the pre-processed signal
  
function lor_Callback(hObject, eventdata)
%Callback for 'Lorentz' edit box  
global preprocess_handles;
global apod_para;
global ind_allvox;
global details_para;
global cur_vox;

% read and check the Lorentz apod para from the edit box 
temp = str2double(get(hObject,'string'));
if((temp < 0) || (temp > 100) || (isnan(temp)))
    errordlg('Invalid value');
    set(hObject,'string',num2str(apod_para.lor(cur_vox)))
    return;
end
if(ind_allvox == 1)
    apod_para.lor(cur_vox) = temp;
    apod_para.temp_lor(cur_vox) = temp;
else
    apod_para.lor(details_para.selected_voxels(1,:)) = temp;
    apod_para.temp_lor(details_para.selected_voxels(1,:)) = temp;
end
set(preprocess_handles.Apodize_Lor_slide,'value',apod_para.lor(cur_vox));

do_and_display_apod(); % apodize and display the pre-processed signal

function slider_lor_Callback(hObject, eventdata)
%Callback for 'Lorentz' slider  
global preprocess_handles;
global apod_para;
global ind_allvox;
global details_para;
global cur_vox;

% read and ckeck the Lorentz apod para from the Lorentz slider 
temp = get(hObject,'Value');
if(ind_allvox == 1)
    apod_para.temp_lor(cur_vox)=temp - apod_para.lor(cur_vox);
    apod_para.lor(cur_vox) = temp;
else
    apod_para.temp_lor(details_para.selected_voxels(1,:))=temp- apod_para.lor(details_para.selected_voxels(1,:));
    apod_para.lor(details_para.selected_voxels(1,:)) = temp;  
end
set(preprocess_handles.Apodize_Lor_edit,'string',num2str(apod_para.lor(cur_vox)));
do_and_display_apod(); % apodize and display the pre-processed signal

function gaus_Callback(hObject, eventdata)
%Callback for 'Gaussian' edit box 
global preprocess_handles;
global apod_para;
global cur_vox;
global ind_allvox;
global details_para;

% read and check the Gaussian apod para from the edit box 
temp = str2double(get(hObject,'string'));
if((temp < 0) || (temp > 1000) || (isnan(temp)))
    errordlg('Invalid value');
    set(hObject,'string',num2str(apod_para.gaus(cur_vox)))
    return;
end
if(ind_allvox == 1)
    apod_para.gaus(cur_vox) = temp;
    apod_para.temp_gaus(cur_vox) = temp;
else
    apod_para.gaus(details_para.selected_voxels(1,:)) = temp;
    apod_para.temp_gaus(details_para.selected_voxels(1,:)) = temp;
end
set(preprocess_handles.Apodize_Gaus_slide,'value',apod_para.gaus(cur_vox));
do_and_display_apod();% apodize and display the pre-processed signal


function slider_gaus_Callback(hObject, eventdata)
%Callback for 'Gaussian' slider 
global preprocess_handles;
global apod_para;
global cur_vox;
global ind_allvox;
global details_para;

% read and check the Gaussian apod para from the slider 
temp = get(hObject,'Value');
if(ind_allvox == 1)
    apod_para.temp_gaus(cur_vox)=temp - apod_para.gaus(cur_vox);
    apod_para.gaus(cur_vox) = temp;
else
    apod_para.temp_gaus(details_para.selected_voxels(1,:))=temp- apod_para.gaus(details_para.selected_voxels(1,:));
    apod_para.gaus(details_para.selected_voxels(1,:)) = temp;  
end

set(preprocess_handles.Apodize_Gaus_edit,'string',num2str(apod_para.gaus(cur_vox)));
do_and_display_apod();% apodize and display the pre-processed signal

% 
function sigm_Callback(hObject, eventdata)
%Callback for 'Sigmoid' edit box 
global preprocess_handles;
global apod_para;
global cur_vox;
global ind_allvox;
global details_para;

% read and check the sigmoid apod para from the edit box 
temp = str2double(get(hObject,'string'));
if((temp < 0)|| (temp > 500) || (isnan(temp)))
    errordlg('Invalid value');
    set(hObject,'string',num2str(apod_para.sigm(cur_vox)))
    return;
end
if(ind_allvox == 1)
    apod_para.sigm(cur_vox) = temp;
    apod_para.temp_sigm(cur_vox) = temp;
else
    apod_para.sigm(details_para.selected_voxels(1,:)) = temp;
    apod_para.temp_sigm(details_para.selected_voxels(1,:)) = temp;
end
set(preprocess_handles.Apodize_Sigm_slide,'value',apod_para.sigm(cur_vox));
do_and_display_apod();% apodize and display the pre-processed signal


function slider_sigm_Callback(hObject, eventdata)
%Callback for 'Sigmoid' slider 
global preprocess_handles;
global apod_para;
global cur_vox;
global ind_allvox;
global details_para;

% read and check the sigmoid apod para from the slider
temp = get(hObject,'Value');
if(ind_allvox == 1)
    apod_para.temp_sigm(cur_vox)=temp - apod_para.sigm(cur_vox);
    apod_para.sigm(cur_vox) = temp;
else
    apod_para.temp_sigm(details_para.selected_voxels(1,:))=temp- apod_para.sigm(details_para.selected_voxels(1,:));
    apod_para.sigm(details_para.selected_voxels(1,:)) = temp;  
end
set(preprocess_handles.Apodize_Sigm_edit,'string',num2str(apod_para.sigm(cur_vox)));
do_and_display_apod();% apodize and display the pre-processed signal

% 
function pha0_val_Callback(hObject, eventdata)
%Callback for 'zero-order' edit box 
global pha_para;
global cur_vox;
global ind_allvox;
global preprocess_handles;
global details_para;

% read and check the zero-order phase value from the edit box 
phi_0 = str2double(get(hObject,'string'));
if((phi_0 < -pi)|| (phi_0 > pi) || (isnan(phi_0)))
    errordlg('Invalid value');
    set(hObject,'string',num2str(pha_para.pha_0(cur_vox)))
    return;
end

if(ind_allvox == 1)
    pha_para.pha_0(cur_vox) = phi_0;
    pha_para.temp_pha0(cur_vox)=phi_0;
else
    pha_para.pha_0(details_para.selected_voxels(1,:)) = phi_0;  
    pha_para.temp_pha0(details_para.selected_voxels(1,:))=phi_0;
end
set(preprocess_handles.pha_corr0_slide,'value',pha_para.pha_0(cur_vox));
do_and_display_phase() % Phase-correct and display the pre-processed signal



% 
function pha1_val_Callback(hObject, eventdata)
%Callback for 'first-order' edit box 
global pha_para;
global cur_vox;
global ind_allvox;
global preprocess_handles;
global details_para;

% read and check the first-order phase value from the edit box 
phi_1 = str2double(get(hObject,'string'));
if((phi_1 < -10)|| (phi_1 > 10) || (isnan(phi_1)))
    errordlg('Invalid value');
    set(hObject,'string',num2str(pha_para.pha_1(cur_vox)))
    return;
end
if(ind_allvox == 1)
    pha_para.pha_1(cur_vox) = phi_1;
    pha_para.temp_pha1(cur_vox)=phi_1;
else
    pha_para.pha_1(details_para.selected_voxels(1,:)) = phi_1; 
    pha_para.temp_pha1(details_para.selected_voxels(1,:))=phi_1;
end
set(preprocess_handles.pha_corr1_slide,'value',pha_para.pha_1(cur_vox));
do_and_display_phase() % Phase-correct and display the pre-processed signal


function pha0_slider_Callback(hObject, eventdata)
%Callback for 'zero-order' slider 
global preprocess_handles;
global pha_para;
global cur_vox;
global ind_allvox;
global details_para;

% read and check the zero-order phase value from the slider 
phi_0 = get(hObject,'Value');
if(ind_allvox == 1)
    pha_para.temp_pha0(cur_vox)=phi_0- pha_para.pha_0(cur_vox);
    pha_para.pha_0(cur_vox) = phi_0;
else
    pha_para.temp_pha0(details_para.selected_voxels(1,:))=phi_0- pha_para.pha_0(details_para.selected_voxels(1,:));
    pha_para.pha_0(details_para.selected_voxels(1,:)) = phi_0;  
end
set(preprocess_handles.pha_corr0_value,'string',num2str(pha_para.pha_0(cur_vox)));
do_and_display_phase();% Phase-correct and display the pre-processed signal



function pha1_slider_Callback(hObject, eventdata)
%Callback for 'first-order' slider 
global pha_para;
global cur_vox;
global ind_allvox;
global preprocess_handles;
global details_para;

% read and check the first-order phase value from the slider 
phi_1 = get(hObject,'Value');
if(ind_allvox == 1)
    pha_para.temp_pha1(cur_vox)=phi_1- pha_para.pha_1(cur_vox);
    pha_para.pha_1(cur_vox) = phi_1;
else
    pha_para.temp_pha1(details_para.selected_voxels(1,:))=phi_1- pha_para.pha_1(details_para.selected_voxels(1,:));
    pha_para.pha_1(details_para.selected_voxels(1,:)) = phi_1; 
end
set(preprocess_handles.pha_corr1_value,'string',num2str(pha_para.pha_1(cur_vox)));
do_and_display_phase()% Phase-correct and display the pre-processed signal




function ind_allvox_value(source,eventdata)
% call back for individual voxel and all voxel button group
global ind_allvox;

% select whether to apply pre-processing on all the voxels or in the selected voxel 
temp = get(eventdata.NewValue,'string');
if(strcmp(temp,'Individual Voxel'))
    ind_allvox = 1;
else
    ind_allvox = 2;
end

function do_and_display_phase()
% This function performs phase correction and displays the signal
global ind_allvox;
global cur_vox;
global pha_para;
global pha_cor_data;
global data_value;
global apod_para;
global details_para;
global preprocess_handles;
global FID_FD_copy;

if(ind_allvox == 1) % perform phase correction on the selected vovel
    %     do phase correction
    [td_sig] = phase_correction(data_value.FID_FD(:,cur_vox),pha_para.temp_pha0(cur_vox),pha_para.temp_pha1(cur_vox),details_para.fres);
    pha_cor_data(:,cur_vox) = td_sig;
%     pha_para.pha_0(cur_vox)=order;
%     do apodization
%     [temp2, function_weight] = apodize(td_sig,apod_para.type(cur_vox),[details_para.Fs,apod_para.lor(cur_vox),apod_para.gaus(cur_vox),apod_para.sigm(cur_vox)]);
    sht_sig = fftshift(fft(td_sig));
    data_value.FID_TD(:,cur_vox) =(td_sig);
    data_value.FID_FD(:,cur_vox) = (sht_sig);
%     apod_para.function_weight(:, cur_vox) = function_weight;
else % perform phase correction on the all vovel
    %     do phase correction
%     td_sig = phase_correction(data_value.FID_FD_org,pha_para.pha_0(preprocess_sel_vox),pha_para.pha_1(preprocess_sel_vox));
    for i=1:length(details_para.selected_voxels(1,:))
        td_sig = phase_correction(data_value.FID_FD(:,i),pha_para.temp_pha0(cur_vox),pha_para.temp_pha1(cur_vox),details_para.fres);
        pha_cor_data(:,i) = td_sig;
    %     do apodization
%     temp2=td_sig;
%         [temp2, function_weight]= apodize(td_sig,apod_para.type(cur_vox),[details_para.Fs,apod_para.lor(cur_vox),apod_para.gaus(cur_vox),apod_para.sigm(cur_vox)]);
        sht_sig = fftshift(fft(td_sig));
%     apod_para.function_weight(:, details_para.selected_voxels(1,:)) = function_weight;
        data_value.FID_TD(:,i) = td_sig;
        data_value.FID_FD(:,i) = sht_sig;  
    end
end
display_data()
set(preprocess_handles.pha_corr0_value,'string',num2str(pha_para.pha_0(cur_vox)));
set(preprocess_handles.Apodize_Sigm_edit,'string',num2str(apod_para.sigm(cur_vox)));

function do_and_display_apod()
% This function performs apodization and displays the signal
global details_para;
global preprocess_handles;
global data_value;
global cur_vox;
global apod_para;
global pha_cor_data;
global ind_allvox;


if(ind_allvox == 1) % perform apodization on the selected vovel
    temp1 = pha_cor_data(:,cur_vox);
        %     do apodization
    [temp2, function_weight] = apodize(temp1,apod_para.type(cur_vox),[details_para.Fs,apod_para.temp_lor(cur_vox),apod_para.temp_gaus(cur_vox),apod_para.temp_sigm(cur_vox)]);
    sht_sig = fftshift(fft(temp2),1);
    data_value.FID_TD(:,cur_vox) = temp2;
    data_value.FID_FD(:,cur_vox) = sht_sig;
    apod_para.function_weight(:, cur_vox) = function_weight;
    pha_cor_data(:,cur_vox) = temp2;
else %perform apodization on the all vovel
    temp1 = pha_cor_data(:,details_para.selected_voxels(1,:));
        %     do apodization
    [temp2, function_weight] = apodize(temp1,apod_para.type(cur_vox),[details_para.Fs,apod_para.temp_lor(cur_vox),apod_para.temp_gaus(cur_vox),apod_para.temp_sigm(cur_vox)]);
    sht_sig = fftshift(fft(temp2),1);
    apod_para.function_weight(:, details_para.selected_voxels(1,:)) = function_weight;
    data_value.FID_TD(:,details_para.selected_voxels(1,:)) = temp2;
    data_value.FID_FD(:,details_para.selected_voxels(1,:)) = sht_sig; 
    pha_cor_data(:,details_para.selected_voxels(1,:)) = temp2;
end

display_data()
set(preprocess_handles.Apodize_Sigm_edit,'string',num2str(apod_para.sigm(cur_vox)));


function TD_FD_display_callback(hObject, eventdata)
global preprocess_display;

temp = get(get(hObject,'SelectedObject'),'String');

if strcmp(temp,'FD (real)')
    preprocess_display = 1;
elseif strcmp(temp,'FD (abs)')
    preprocess_display = 2;
elseif strcmp(temp,'FD (imag)')
    preprocess_display = 3;
else
    preprocess_display = 0;
end
display_data()

function display_data()
global preprocess_display;
global preprocess_handles;
global cur_vox;
global details_para;
global data_value;
global apod_para;
global display_handles;


seg = details_para.disp_seg;
ppm = details_para.ppm_referenced;

if (details_para.ds_flag(:,cur_vox)==0)

    if (preprocess_display == 1)
        cla(display_handles.axis);
        plot_var = real(data_value.FID_FD(seg(1):seg(2),cur_vox));
        plot(display_handles.axis,ppm(seg(1):seg(2))',plot_var);
        max_p = max(plot_var);
        min_p = min(plot_var);
        axis(display_handles.axis,[ppm(seg(1)), ppm(seg(2)),min_p - 0.01*max_p, max_p + 0.01*max_p])
        set(display_handles.axis,'XDir','reverse');
        legend(display_handles.axis,'Original');
        xlabel(display_handles.axis,'PPM');
        ylabel(display_handles.axis,'Amplitude');
    elseif (preprocess_display == 2)
        cla(display_handles.axis);
        plot_var = abs(data_value.FID_FD(seg(1):seg(2),cur_vox));
        plot(display_handles.axis,ppm(seg(1):seg(2))',plot_var);
        max_p = max(plot_var);
        min_p = min(plot_var);
        axis(display_handles.axis,[ppm(seg(1)), ppm(seg(2)),min_p - 0.01*max_p, max_p + 0.01*max_p])
        set(display_handles.axis,'XDir','reverse');
        legend(display_handles.axis,'Original');
        xlabel(display_handles.axis,'PPM');
        ylabel(display_handles.axis,'Amplitude');
    elseif (preprocess_display == 3)
        cla(display_handles.axis);
        plot_var = imag(data_value.FID_FD(seg(1):seg(2),cur_vox));
        plot(display_handles.axis,ppm(seg(1):seg(2))',plot_var);
        max_p = max(plot_var);
        min_p = min(plot_var);
        axis(display_handles.axis,[ppm(seg(1)), ppm(seg(2)),min_p - 0.01*max_p, max_p + 0.01*max_p])
        set(display_handles.axis,'XDir','reverse');
        legend(display_handles.axis,'Original');
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
        legend(display_handles.axis,'Original');
        xlabel(display_handles.axis,'time');
        ylabel(display_handles.axis,'Amplitude');
    %     hold on;
    %     plot(display_handles.axis,details_para.t',filter,'r');
    %     axis(display_handles.axis,[details_para.t(1), details_para.t(details_para.N),min_p , max_p])
    %     hold off
    end
elseif(details_para.ds_flag(:,cur_vox)==1)
       if (preprocess_display == 1)
        cla(display_handles.axis);
        plot_var = real(data_value.FID_FD_org2(seg(1):seg(2),cur_vox));
        plot_ds = real(data_value.FID_FD(seg(1):seg(2),cur_vox));
        if details_para.ws_count==0
            plot_ws = real(data_value.FID_FD(seg(1):seg(2),cur_vox));
        else
            plot_ws = real(data_value.FID_FD_org1(seg(1):seg(2),cur_vox));
        end
        plot(display_handles.axis,ppm(seg(1):seg(2))',plot_var);
        hold (display_handles.axis,'on'); 
        if details_para.ws_count>0
            plot(display_handles.axis,ppm(seg(1):seg(2))',plot_ds);
            plot(display_handles.axis,ppm(seg(1):seg(2))',plot_ws,'g');
        else
            plot(display_handles.axis,ppm(seg(1):seg(2))',plot_ws,'g');
            plot(display_handles.axis,ppm(seg(1):seg(2))',plot_ds,'r');
        end
        max_p = max(plot_var);
        min_p = min(plot_ws);%- round((min(plot_ws).*2));
        if min_p<-20
                min_p=-30;
        end
        axis(display_handles.axis,[ppm(seg(1)), ppm(seg(2)),min_p - 0.01*max_p, max_p + 0.01*max_p])
        set(display_handles.axis,'XDir','reverse');
        if details_para.ws_count>0
            lgd=legend(display_handles.axis,'water suppressed','Denoised','Original');
        else
            lgd=legend(display_handles.axis,'water suppressed','Original','Denoised');
        end
        lgd.FontSize=8;
        xlabel(display_handles.axis,'PPM');
        ylabel(display_handles.axis,'Amplitude');
       elseif (preprocess_display == 2)
        cla(display_handles.axis);
        plot_var = abs(data_value.FID_FD_org2(seg(1):seg(2),cur_vox));
        plot_ds = abs(data_value.FID_FD(seg(1):seg(2),cur_vox));
        if details_para.ws_count==0
            plot_ws = abs(data_value.FID_FD(seg(1):seg(2),cur_vox));
        else
            plot_ws = abs(data_value.FID_FD_org1(seg(1):seg(2),cur_vox));
        end
        plot(display_handles.axis,ppm(seg(1):seg(2))',plot_var);
        hold (display_handles.axis,'on'); 
        if details_para.ws_count>0
            plot(display_handles.axis,ppm(seg(1):seg(2))',plot_ds);
            plot(display_handles.axis,ppm(seg(1):seg(2))',plot_ws,'g');
        else
            plot(display_handles.axis,ppm(seg(1):seg(2))',plot_ws,'g');
            plot(display_handles.axis,ppm(seg(1):seg(2))',plot_ds,'r');
        end
        max_p = max(plot_var);
        min_p = min(plot_ws);%-round((min(plot_ws).*2));
        if min_p<-20
                min_p=-30;
        end
        axis(display_handles.axis,[ppm(seg(1)), ppm(seg(2)),min_p - 0.01*max_p, max_p + 0.01*max_p])
        set(display_handles.axis,'XDir','reverse');
        if details_para.ws_count>0
            lgd=legend(display_handles.axis,'water suppressed','Denoised','Original');
        else
            lgd=legend(display_handles.axis,'water suppressed','Original','Denoised');
        end
        lgd.FontSize=8;
        xlabel(display_handles.axis,'PPM');
        ylabel(display_handles.axis,'Amplitude');
        elseif (preprocess_display == 3)
        cla(display_handles.axis);
        plot_var = imag(data_value.FID_FD_org2(seg(1):seg(2),cur_vox));
        plot_ds =imag(data_value.FID_FD(seg(1):seg(2),cur_vox));
        if details_para.ws_count==0
            plot_ws = imag(data_value.FID_FD(seg(1):seg(2),cur_vox));
        else
            plot_ws = imag(data_value.FID_FD_org1(seg(1):seg(2),cur_vox));
        end
        plot(display_handles.axis,ppm(seg(1):seg(2))',plot_var);
        hold (display_handles.axis,'on'); 
        if details_para.ws_count>0
            plot(display_handles.axis,ppm(seg(1):seg(2))',plot_ds);
            plot(display_handles.axis,ppm(seg(1):seg(2))',plot_ws,'g');
        else
            plot(display_handles.axis,ppm(seg(1):seg(2))',plot_ws,'g');
            plot(display_handles.axis,ppm(seg(1):seg(2))',plot_ds,'r');
        end
        max_p = max(plot_var);
        min_p = min(plot_ws);%-round((min(plot_ws).*2));
        if min_p<-20
                min_p=-30;
        end
        axis(display_handles.axis,[ppm(seg(1)), ppm(seg(2)),min_p - 0.01*max_p, max_p + 0.01*max_p])
        set(display_handles.axis,'XDir','reverse');
        if details_para.ws_count>0
            lgd=legend(display_handles.axis,'water suppressed','Denoised','Original');
        else
            lgd=legend(display_handles.axis,'water suppressed','Original','Denoised');
        end
        lgd.FontSize=8;
        xlabel(display_handles.axis,'PPM');
        ylabel(display_handles.axis,'Amplitude');
       else
        cla(display_handles.axis);
        plot_var = real(data_value.FID_TD_org2(:,cur_vox));
        plot_ds = real(data_value.FID_TD(:,cur_vox));
        if details_para.ws_count==0
            plot_ws = real(data_value.FID_TD(seg(1):seg(2),cur_vox));
        else
            plot_ws = real(data_value.FID_TD_org1(seg(1):seg(2),cur_vox));
        end
    %     plot_var= ifftshift(real(data_value.FID_TD(:,cur_vox))+flipud(real(data_value.FID_TD(:,cur_vox))));
    %     filter = apod_para.function_weight(:,cur_vox);
        max_p = max(abs(plot_var));
        min_p = min(abs(plot_ws));%-round((min(plot_ws).*2));
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
        hold (display_handles.axis,'on');  
        if details_para.ws_count>0
            plot(display_handles.axis,details_para.t',plot_ds);
            plot(display_handles.axis,details_para.t',plot_ws,'g');
        else
            plot(display_handles.axis,details_para.t',plot_ws,'g');
            plot(display_handles.axis,details_para.t',plot_ds,'r');
        end
        if details_para.ws_count>0
            lgd=legend(display_handles.axis,'water suppressed','Denoised','Original');
        else
            lgd=legend(display_handles.axis,'water suppressed','Original','Denoised');
        end
        lgd.FontSize=8;
        xlabel(display_handles.axis,'time');
        ylabel(display_handles.axis,'Amplitude');
    %     hold on;
    %     plot(display_handles.axis,details_para.t',filter,'r');
    %     axis(display_handles.axis,[details_para.t(1), details_para.t(details_para.N),min_p , max_p])
    %     hold off
       end
else
      if (preprocess_display == 1)
        cla(display_handles.axis);
        plot_var = real(data_value.FID_FD(seg(1):seg(2),cur_vox));
        plot_ws = real(data_value.FID_FD_org1(seg(1):seg(2),cur_vox));
        plot(display_handles.axis,ppm(seg(1):seg(2))',plot_var);
        hold (display_handles.axis,'on'); 
        plot(display_handles.axis,ppm(seg(1):seg(2))',plot_ws,'g');
        max_p = max(plot_var);
        min_p = min(plot_ws);%-round((min(plot_ws).*2));
        if min_p<-20
                min_p=-30;
        end
        axis(display_handles.axis,[ppm(seg(1)), ppm(seg(2)),min_p - 0.01*max_p, max_p + 0.01*max_p])
        set(display_handles.axis,'XDir','reverse');
        lgd=legend(display_handles.axis,'water suppressed','Original');
        lgd.FontSize=8;
        xlabel(display_handles.axis,'PPM');
        ylabel(display_handles.axis,'Amplitude');
       elseif (preprocess_display == 2)
        cla(display_handles.axis);
        plot_var = abs(data_value.FID_FD(seg(1):seg(2),cur_vox));
        plot_ws = abs(data_value.FID_FD_org1(seg(1):seg(2),cur_vox));
        plot(display_handles.axis,ppm(seg(1):seg(2))',plot_var);
        hold (display_handles.axis,'on'); 
        plot(display_handles.axis,ppm(seg(1):seg(2))',plot_ws,'g');
        max_p = max(plot_var);
        min_p = min(plot_ws);%-round((min(plot_ws).*2));
        if min_p<-20
                min_p=-30;
        end
        axis(display_handles.axis,[ppm(seg(1)), ppm(seg(2)),min_p - 0.01*max_p, max_p + 0.01*max_p])
        set(display_handles.axis,'XDir','reverse');
        lgd=legend(display_handles.axis,'water suppressed','Original');
        lgd.FontSize=8;
        xlabel(display_handles.axis,'PPM');
        ylabel(display_handles.axis,'Amplitude');
         elseif (preprocess_display == 3)
        cla(display_handles.axis);
        plot_var = imag(data_value.FID_FD(seg(1):seg(2),cur_vox));
        plot_ws = imag(data_value.FID_FD_org1(seg(1):seg(2),cur_vox));
        plot(display_handles.axis,ppm(seg(1):seg(2))',plot_var);
        hold (display_handles.axis,'on'); 
        plot(display_handles.axis,ppm(seg(1):seg(2))',plot_ws,'g');
        max_p = max(plot_var);
        min_p = min(plot_ws);%-round((min(plot_ws).*2));
        if min_p<-20
                min_p=-30;
        end
        axis(display_handles.axis,[ppm(seg(1)), ppm(seg(2)),min_p - 0.01*max_p, max_p + 0.01*max_p])
        set(display_handles.axis,'XDir','reverse');
        lgd=legend(display_handles.axis,'water suppressed','Original');
        lgd.FontSize=8;
        xlabel(display_handles.axis,'PPM');
        ylabel(display_handles.axis,'Amplitude');
       else
        cla(display_handles.axis);
        plot_var = real(data_value.FID_TD(:,cur_vox));
        plot_ws = real(data_value.FID_TD_org1(:,cur_vox));
    %     plot_var= ifftshift(real(data_value.FID_TD(:,cur_vox))+flipud(real(data_value.FID_TD(:,cur_vox))));
    %     filter = apod_para.function_weight(:,cur_vox);
        max_p = max(abs(plot_var));
        min_p = min(abs(plot_ws));%-round((min(plot_ws).*2));
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
        hold (display_handles.axis,'on'); 
        plot(display_handles.axis,details_para.t',plot_ws,'g');
        lgd=legend(display_handles.axis,'water suppressed','Original');
        lgd.FontSize=8;
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
global preprocess_handles;
global log_file;
ini_path = log_file.saving_path;
[temp,s_path] = uiputfile(ini_path);
if(temp == 0)
    return;
else
    dot_pos = strfind(temp, '.');
    s_fln = temp(1:dot_pos-1);
    file_name = strcat(s_path,s_fln);
    saveas(preprocess_handles.fig,file_name,'png');
    log_file.saving_path = s_path;
end



function close_Callback(hObject, eventdata)
%Callback for 'close' button
global preprocess_handles;
global met_tick;
global mainGUI_handles;
global details_para;
global display_handles;

if (met_tick==1 && ~isempty(preprocess_handles.met_fig))
    try
        close(display_handles.fig);
    catch 
        return;
    end
    close(preprocess_handles.fig)
    close(preprocess_handles.met_fig);
%     close(display_handles.fig);
else
     try
        close(display_handles.fig);
        close(preprocess_handles.fig);
     catch 
        return;
     end
%      close(preprocess_handles.fig);
%      close(display_handles.fig);
end
details_para.preprocess_done=1;
details_para.preprocess_going=0;
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
function metabolite_map_Callback(hObject, eventdata)

global data_value;
global details_para;
global preprocess_handles;
global preprocess_sel_vox;
global FID_FD_copy;
global met_tick;

for i = 1:details_para.num_vox
    [ fit_sig,amp,Ind_sig,freq] = HSVDPK_new(data_value.FID_TD(:,i),details_para.fres,details_para.Fs,10);
    ind=find((freq>=-80 & freq<=80));
    half_sig2=sum(Ind_sig(:,ind),2);
    data_value.FID_TD(:,i)=data_value.FID_TD(:,i)- half_sig2;
    data_value.FID_FD(:,i)=fftshift(fft(data_value.FID_TD(:,i)));
end
PE=details_para.PE;
RE=details_para.RE;
if round(details_para.Fs)==2000
    seg_NAA=160:180;
    seg_Lac=130:159;%210:220;
    seg_Cho=209:218;%205:210;%200:207;%205:210;%205:210;%205:210;%200:205;198:210;
    seg_Cr=196:205;%200:205;%200:205;198:210;
elseif round(details_para.Fs)==1600
    seg_NAA=135:165;
    seg_Lac=110:130;
    seg_Cho=185:200;
elseif round(details_para.Fs)==1000
    seg_NAA=60:95;
    seg_Lac=30:50;
    seg_Cho=145:165;
    seg_Cr=145:165;
elseif round(details_para.Fs)==2500
    seg_NAA=186:200;
    seg_Lac=165:185;
    seg_Cho=210:225;
end
for i=1:(PE/1)*RE
    temp1=abs(data_value.FID_FD(:,i));
    NAA_CSI(i)=max(temp1(seg_NAA))/1;
    Lac_CSI(i)=max(temp1(seg_Lac))/1;
    Cho_CSI(i)=max(temp1(seg_Cho))/1;
    Cr_CSI(i)=max(temp1(seg_Cr));

end
l=1;
for i=1:PE/1
    for j=1:RE
        NAA_map_csi(i,j)=NAA_CSI(l);
        Lac_map_csi(i,j)=Lac_CSI(l);
        Cho_map_csi(i,j)=Cho_CSI(l);
        Cr_map_csi(i,j)=Cr_CSI(l);
        l=l+1;
    end
end
met_tick=1;
preprocess_handles.met_fig=figure();
set(preprocess_handles.met_fig,'Name','RAPID 1.0','NumberTitle', 'off');
subplot(2,2,1)
imagesc(NAA_map_csi);%./Cr_map_csi);
title('NAA map');
colorbar
subplot(2,2,2)
imagesc(Lac_map_csi);%./Cr_map_csi);
title('Lac map');
colorbar
subplot(2,2,3)
imagesc(Cho_map_csi);%./Cr_map_csi);
title('Cho map');
colorbar
subplot(2,2,4)
imagesc(Cr_map_csi);%./Cr_map_csi);
title('Cr map');
colorbar

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
%------------------------------------------------------------------------------------------------------------------------------
function ws_type_Callback(hObject, eventdata)
%Callback for 'Apodisation type' 
global preprocess_handles;
global apod_para;
global ws_para;
global ind_allvox;
global details_para;
global cur_vox;
% global data_value;

% get the selected apodization type
temp = get(hObject,'Value');
if(ind_allvox == 1)
    ws_para.type(cur_vox) = temp;
else
    ws_para.type(details_para.selected_voxels(1,:)) = temp;
end

% enable/disable based on the selected type
switch(temp)
  case 1 %SCSA
        set(preprocess_handles.ws_SCSA_edit, 'Enable', 'on');
        set(preprocess_handles.ws_HLSVD_edit, 'Enable', 'off');
        set(preprocess_handles.ws_SCSA_slide, 'Enable', 'on');
        set(preprocess_handles.ws_HLSVD_slide, 'Enable', 'off');
  case 2 %HLSVD
        set(preprocess_handles.ws_SCSA_edit, 'Enable', 'off');
        set(preprocess_handles.ws_HLSVD_edit, 'Enable', 'on');
        set(preprocess_handles.ws_SCSA_slide, 'Enable', 'off');
        set(preprocess_handles.ws_HLSVD_slide, 'Enable', 'on');
end
%-----------------------------------------------------------------------------------------------------------------
function slider_SCSA_Callback(hObject, eventdata)
%Callback for 'Sigmoid' slider 
global preprocess_handles;
global apod_para;
global ws_para;
global cur_vox;
global ind_allvox;
global details_para;

% read and check the sigmoid apod para from the slider
temp = get(hObject,'Value');
if(ind_allvox == 1)
    ws_para.temp_SCSA(cur_vox)=temp - ws_para.SCSA(cur_vox);
    ws_para.SCSA(cur_vox) = temp;
else
    ws_para.temp_SCSA(details_para.selected_voxels(1,:))=temp- ws_para.SCSA(details_para.selected_voxels(1,:));
    ws_para.SCSA(details_para.selected_voxels(1,:)) = temp;  
end
set(preprocess_handles.ws_SCSA_edit,'string',num2str(ws_para.SCSA(cur_vox)));
%--------------------------------------------------------------------------------------------------------------------------------
function slider_HLSVD_Callback(hObject, eventdata)
%Callback for 'Sigmoid' slider 
global preprocess_handles;
global apod_para;
global ws_para;
global cur_vox;
global ind_allvox;
global details_para;

% read and check the sigmoid apod para from the slider
temp = get(hObject,'Value');
if(ind_allvox == 1)
    ws_para.temp_HLSVD(cur_vox)=temp - ws_para.HLSVD(cur_vox);
    ws_para.HLSVD(cur_vox) = temp;
else
    ws_para.temp_HLSVD(details_para.selected_voxels(1,:))=temp- ws_para.HLSVD(details_para.selected_voxels(1,:));
    ws_para.HLSVD(details_para.selected_voxels(1,:)) = temp;  
end
set(preprocess_handles.ws_HLSVD_edit,'string',num2str(ws_para.HLSVD(cur_vox)));
%--------------------------------------------------------------------------------------------------------------------------------
function HLSVD_ws_Callback(hObject, eventdata)
%Callback for 'Sigmoid' slider 
global preprocess_handles;
global apod_para;
global ws_para;
global cur_vox;
global ind_allvox;
global details_para;
global hsvd_para;


% read and check the sigmoid apod para from the slider
temp = str2double(get(hObject,'string'));
if((temp < 0)|| (temp > 30) || (isnan(temp)))
    errordlg('Invalid value');
    set(hObject,'string',num2str(ws_para.HLSVD(cur_vox)))
    return;
end
if(ind_allvox == 1)
    ws_para.HLSVD(cur_vox) = temp;
else
    ws_para.HLSVD(details_para.selected_voxels(1,:)) = temp;
end
set(preprocess_handles.ws_HLSVD_slide,'value',ws_para.HLSVD(cur_vox));
%--------------------------------------------------------------------------------------------------------------------------------
function SCSA_ws_Callback(hObject, eventdata)
%Callback for 'Sigmoid' slider 
global preprocess_handles;
global apod_para;
global ws_para;
global cur_vox;
global ind_allvox;
global details_para;
global scsa_para;


% read and check the sigmoid apod para from the slider
temp = str2double(get(hObject,'string'));
if((temp < 0)|| (temp > 30) || (isnan(temp)))
    errordlg('Invalid value');
    set(hObject,'string',num2str(ws_para.SCSA(cur_vox)))
    return;
end
if(ind_allvox == 1)
    ws_para.SCSA(cur_vox) = temp;
else
    ws_para.SCSA(details_para.selected_voxels(1,:)) = temp;
end
set(preprocess_handles.ws_SCSA_slide,'value',ws_para.SCSA(cur_vox));
%--------------------------------------------------------------------------------------------------------------------------------
function water_suppress_Callback(hObject, eventdata)

global data_value;
global details_para;
global preprocess_handles;
global preprocess_sel_vox;
global FID_FD_copy;
global pha_cor_data;
global ws_para;
global cur_vox;
global ind_allvox;
global Loop_Water_reduction Nh_WP_reduction;

% load(data_value.temp_filename);
freq_vector=[-1000,-90];

% set(preprocess_handles.fig,'Visible','off');

details_para.ws_count=details_para.ws_count+1;

if details_para.ws_count==1
    data_value.FID_TD_ws=data_value.FID_TD;
    data_value.FID_FD_ws=data_value.FID_FD;
end


cnt=1;

try
    if(ind_allvox == 1)
        if (ws_para.type(cur_vox)==1)
                cnt=cur_vox;
                Skip_W_Estim=0;                     % set to 1 in invitro case
                gm=0.5;fs=1; 
                WR_Ratio = 1 ;      % Water Reduction Ratio 
                MaxIter= 100;      % The maximum iteration for water residual reduction loop
                Loop_Water_reduction=4;          % 1 :  Reduce water residual
                Nh_WP_reduction= 3;   
%                 freq_vector=[-details_para.Fs/2,-50];
                org_data_unsuppressed_FD_SCSA=real(data_value.FID_FD_ws(:,cur_vox));
                min_water=min(org_data_unsuppressed_FD_SCSA);
                Nh_list=[15 24];
                if min_water<-20
                    min_water=0;
                end
                data_value.FID_FD_org1(:,cur_vox)=data_value.FID_FD_ws(:,cur_vox);
                data_value.FID_TD_org1(:,cur_vox)=ifft(ifftshift( data_value.FID_FD_org1(:,cur_vox)));
                yf_SCSA=SCSA_MRS_Water_Suppression_GUI(details_para.ppm_referenced',details_para.fres',freq_vector, org_data_unsuppressed_FD_SCSA,gm,fs,ws_para.SCSA(cur_vox));
                yf_SCSA=yf_SCSA+(min_water-min(yf_SCSA));
                data_value.FID_TD(:,cur_vox)=ifft(ifftshift(yf_SCSA'));
                data_value.FID_FD(:,cur_vox)=yf_SCSA';
%                 waitbar(1);
                details_para.ds_flag(:,cur_vox)=2;
        else
                
                cnt=cur_vox;
                org_data_unsuppressed=data_value.FID_TD_ws(:,cur_vox);
                data_value.FID_FD_org1(:,cur_vox)=data_value.FID_FD_ws(:,cur_vox);
                min_water=min(real(data_value.FID_FD_org1(:,cur_vox)));
                if min_water<-20
                    min_water=0;
                end
                wait_bar = waitbar(0.12,'Water peak estimation','Name','MRS water suppression using HLSVD');
                data_value.FID_TD_org1(:,cur_vox)=ifft(ifftshift( data_value.FID_FD_org1(:,cur_vox)));
                re_bl = ssa(real(fftshift(fft(org_data_unsuppressed))),40); % estimate the baseline for real part
                im_bl = ssa(imag(fftshift(fft(org_data_unsuppressed))),40); % estimate the baseline for imaginary part
                data_value.FID_FD_base(:,1)=re_bl + 1i*im_bl;
                data_value.FID_FD_b(:,1)=fftshift(fft(org_data_unsuppressed))-data_value.FID_FD_base(:,1);
                data_value.FID_TD_b(:,1)=ifft(ifftshift(data_value.FID_FD_b(:,1)));
                baseline(:,1)=data_value.FID_FD_base(:,1);
                baseline_t(:,1)=ifft(ifftshift(data_value.FID_FD_base(:,1)));
                [ fit_sig,amp,Ind_sig,freq,Z_sig,FWHM] = HSVDPK_new(((org_data_unsuppressed)),details_para.fres,details_para.Fs,ws_para.HLSVD(cur_vox));
                signal_fitted(:,1)=fftshift(fft(fit_sig));
                count=2;
                waitbar(0.56,wait_bar, 'Water peak parameters estimation');   
                while(1)
                    re_bl = ssa(real(fftshift(fft(((org_data_unsuppressed))))-signal_fitted(:,count-1)),30); % estimate the baseline for real part
                    im_bl = ssa(imag(fftshift(fft(((org_data_unsuppressed))))-signal_fitted(:,count-1)),30); % estimate the baseline for imaginary part
                    baseline(:,count)=re_bl + 1i*im_bl;
                    baseline_t(:,count)=ifft(ifftshift(baseline(:,count)));
                    [ fit_sig,amp,Ind_sig,freq,Z_sig,FWHM] = HSVDPK_new(((org_data_unsuppressed))-baseline_t(:,count),details_para.fres,details_para.Fs,ws_para.HLSVD(cur_vox));
                    signal_fitted(:,count)=fftshift(fft(fit_sig));
                    if (std(abs(baseline(:,count)-baseline(:,count-1)))<=1e-5 || count==50)
                        break;
                    end
                    count=count+1;
                end
                ind=find((freq>=-80 & freq<=80));
                hsvd_sig=sum(Ind_sig(:,ind),2);
                hsvd_sig=org_data_unsuppressed - hsvd_sig;
                hsvd_sig_FD=fftshift(fft((hsvd_sig)));
                data_value.FID_TD(:,cur_vox)=hsvd_sig;
                data_value.FID_FD(:,cur_vox)=hsvd_sig_FD+ (min_water-min(real(hsvd_sig_FD)));
                data_value.FID_TD(:,cur_vox)=ifft(ifftshift(data_value.FID_FD(:,cur_vox)));
%                 waitbar(1);
                details_para.ds_flag(:,cur_vox)=2;
                waitbar(1,wait_bar)  
                close(wait_bar)
        end
    
    else

    if (ws_para.type(details_para.selected_voxels(1,:))==1)
        for i = 1:details_para.num_vox
            cnt=i;
            Skip_W_Estim=0;                     % set to 1 in invitro case
            gm=0.5;fs=1; 
            WR_Ratio = 1 ;      % Water Reduction Ratio 
            MaxIter= 100;      % The maximum iteration for water residual reduction loop
            Loop_Water_reduction=4;          % 1 :  Reduce water residual
            Nh_WP_reduction= 3; 
%             freq_vector=[-details_para.Fs/2,-50];
            org_data_unsuppressed_FD_SCSA=real(data_value.FID_FD_ws(:,i));
            min_water=min(org_data_unsuppressed_FD_SCSA);
            Nh_list=[15 24];
            if min_water<-20
                    min_water=0;
            end
            data_value.FID_FD_org1(:,i)=data_value.FID_FD_ws(:,i);
            data_value.FID_TD_org1(:,i)=ifft(ifftshift( data_value.FID_FD_org1(:,i)));
            yf_SCSA=SCSA_MRS_Water_Suppression_GUI(details_para.ppm_referenced',details_para.fres',freq_vector, org_data_unsuppressed_FD_SCSA,gm,fs,ws_para.SCSA(cur_vox));
            yf_SCSA=yf_SCSA+(min_water-min(yf_SCSA));
            data_value.FID_TD(:,i)=ifft(ifftshift(yf_SCSA'));
            data_value.FID_FD(:,i)=yf_SCSA';
%             waitbar(i/details_para.num_vox);
            details_para.ds_flag(:,i)=2;
        end
    else
        
        for i = 1:details_para.num_vox
            cnt=i;
            org_data_unsuppressed=data_value.FID_TD_ws(:,i);
            data_value.FID_FD_org1(:,i)=data_value.FID_FD_ws(:,i);
            min_water=min(real(data_value.FID_FD_org1(:,i)));
            if min_water<-20
                    min_water=0;
            end
            wait_bar = waitbar(0.12,'Water peak estimation','Name','MRS water suppression using HLSVD');
            data_value.FID_TD_org1(:,i)=ifft(ifftshift( data_value.FID_FD_org1(:,i)));
            re_bl = ssa(real(fftshift(fft(org_data_unsuppressed))),40); % estimate the baseline for real part
            im_bl = ssa(imag(fftshift(fft(org_data_unsuppressed))),40); % estimate the baseline for imaginary part
            data_value.FID_FD_base(:,1)=re_bl + 1i*im_bl;
            data_value.FID_FD_b(:,1)=fftshift(fft(org_data_unsuppressed))-data_value.FID_FD_base(:,1);
            data_value.FID_TD_b(:,1)=ifft(ifftshift(data_value.FID_FD_b(:,1)));
            baseline(:,1)=data_value.FID_FD_base(:,1);
            baseline_t(:,1)=ifft(ifftshift(data_value.FID_FD_base(:,1)));
            [ fit_sig,amp,Ind_sig,freq,Z_sig,FWHM] = HSVDPK_new(((org_data_unsuppressed)),details_para.fres,details_para.Fs,ws_para.HLSVD(cur_vox));
            signal_fitted(:,1)=fftshift(fft(fit_sig));
            count=2;
            waitbar(0.56,wait_bar, 'Water peak parameters estimation'); 
            while(1)
                re_bl = ssa(real(fftshift(fft(((org_data_unsuppressed))))-signal_fitted(:,count-1)),30); % estimate the baseline for real part
                im_bl = ssa(imag(fftshift(fft(((org_data_unsuppressed))))-signal_fitted(:,count-1)),30); % estimate the baseline for imaginary part
                baseline(:,count)=re_bl + 1i*im_bl;
                baseline_t(:,count)=ifft(ifftshift(baseline(:,count)));
                [ fit_sig,amp,Ind_sig,freq,Z_sig,FWHM] = HSVDPK_new(((org_data_unsuppressed))-baseline_t(:,count),details_para.fres,details_para.Fs,ws_para.HLSVD(cur_vox));
                signal_fitted(:,count)=fftshift(fft(fit_sig));
                if (std(abs(baseline(:,count)-baseline(:,count-1)))<=1e-5 || count==50)
                    break;
                end
                count=count+1;
            end
            ind=find((freq>=-80 & freq<=80));
            hsvd_sig=sum(Ind_sig(:,ind),2);
            hsvd_sig=org_data_unsuppressed - hsvd_sig;
            hsvd_sig_FD=fftshift(fft((hsvd_sig)));
            data_value.FID_TD(:,i)=hsvd_sig;
            data_value.FID_FD(:,i)=hsvd_sig_FD++ (min_water-min(real(hsvd_sig_FD)));
            data_value.FID_TD(:,i)=ifft(ifftshift(data_value.FID_FD(:,i)));
%             waitbar(i/details_para.num_vox);
            details_para.ds_flag(:,i)=2;
            waitbar(1,wait_bar)  
            close(wait_bar)
        end
    end
    end
catch
%         close(h);
        details_para.ds_flag(cnt)=0;
        set(preprocess_handles.fig,'Visible','on');
        return;
end
set(preprocess_handles.fig,'Visible','on');
% FID_FD_copy=data_value.FID_FD;
pha_cor_data(:,details_para.selected_voxels(1,:))=data_value.FID_TD;
% seg = details_para.disp_seg;
% ppm = details_para.ppm_referenced;
% cla(preprocess_handles.axis);
% plot_var = real(data_value.FID_FD(seg(1):seg(2)));
% plot(preprocess_handles.axis,ppm(seg(1):seg(2))',plot_var);
% max_p = max(plot_var);
% mim_p = min(plot_var);
% axis(preprocess_handles.axis,[ppm(seg(1)), ppm(seg(2)),mim_p - 0.1*max_p, max_p + 0.1*max_p])
% set(preprocess_handles.axis,'XDir','reverse')
display_data()
%-----------------------------------------------------------------------------------------------------------
function ds_type_Callback(hObject, eventdata)
%Callback for 'Apodisation type' 
global preprocess_handles;
global apod_para;
global ds_para;
global ind_allvox;
global details_para;
global cur_vox;
% global data_value;

% get the selected apodization type
temp = get(hObject,'Value');
if(ind_allvox == 1)
    ds_para.type(cur_vox) = temp;
else
    ds_para.type(details_para.selected_voxels(1,:)) = temp;
end

% enable/disable based on the selected type
switch(temp)
  case 1 %SCSA
        set(preprocess_handles.ds_SCSA_edit, 'Enable', 'on');
        set(preprocess_handles.ds_HLSVD_edit, 'Enable', 'off');
        set(preprocess_handles.ds_SCSA_slide, 'Enable', 'on');
        set(preprocess_handles.ds_HLSVD_slide, 'Enable', 'off');
  case 2 %HLSVD
        set(preprocess_handles.ds_SCSA_edit, 'Enable', 'off');
        set(preprocess_handles.ds_HLSVD_edit, 'Enable', 'on');
        set(preprocess_handles.ds_SCSA_slide, 'Enable', 'off');
        set(preprocess_handles.ds_HLSVD_slide, 'Enable', 'on');
end
%-----------------------------------------------------------------------------------------------------------------
function slider_SCSA_ds_Callback(hObject, eventdata)
%Callback for 'Sigmoid' slider 
global preprocess_handles;
global apod_para;
global ds_para;
global cur_vox;
global ind_allvox;
global details_para;

% read and check the sigmoid apod para from the slider
temp = get(hObject,'Value');
if(ind_allvox == 1)
    ds_para.temp_SCSA(cur_vox)=temp - ds_para.SCSA(cur_vox);
    ds_para.SCSA(cur_vox) = temp;
else
    ds_para.temp_SCSA(details_para.selected_voxels(1,:))=temp- ds_para.SCSA(details_para.selected_voxels(1,:));
    ds_para.SCSA(details_para.selected_voxels(1,:)) = temp;  
end
set(preprocess_handles.ds_SCSA_edit,'string',num2str(ds_para.SCSA(cur_vox)));
%--------------------------------------------------------------------------------------------------------------------------------
function slider_HLSVD_ds_Callback(hObject, eventdata)
%Callback for 'Sigmoid' slider 
global preprocess_handles;
global apod_para;
global ds_para;
global cur_vox;
global ind_allvox;
global details_para;

% read and check the sigmoid apod para from the slider
temp = get(hObject,'Value');
if(ind_allvox == 1)
    ds_para.temp_HLSVD(cur_vox)=temp - ds_para.HLSVD(cur_vox);
    ds_para.HLSVD(cur_vox) = temp;
else
    ds_para.temp_HLSVD(details_para.selected_voxels(1,:))=temp- ds_para.HLSVD(details_para.selected_voxels(1,:));
    ds_para.HLSVD(details_para.selected_voxels(1,:)) = temp;  
end
set(preprocess_handles.ds_HLSVD_edit,'string',num2str(ds_para.HLSVD(cur_vox)));
%--------------------------------------------------------------------------------------------------------------------------------
function HLSVD_ds_Callback(hObject, eventdata)
%Callback for 'Sigmoid' slider 
global preprocess_handles;
global apod_para;
global ds_para;
global cur_vox;
global ind_allvox;
global details_para;
global hsvd_para;


% read and check the sigmoid apod para from the slider
temp = str2double(get(hObject,'string'));
if((temp < 0)|| (temp > 10) || (isnan(temp)))
    errordlg('Invalid value');
    set(hObject,'string',num2str(ds_para.HLSVD(cur_vox)))
    return;
end
if(ind_allvox == 1)
    ds_para.HLSVD(cur_vox) = temp;
else
    ds_para.HLSVD(details_para.selected_voxels(1,:)) = temp;
end
set(preprocess_handles.ds_HLSVD_slide,'value',ds_para.HLSVD(cur_vox));
%--------------------------------------------------------------------------------------------------------------------------------
function SCSA_ds_Callback(hObject, eventdata)
%Callback for 'Sigmoid' slider 
global preprocess_handles;
global apod_para;
global ds_para;
global cur_vox;
global ind_allvox;
global details_para;
global scsa_para;


% read and check the sigmoid apod para from the slider
temp = str2double(get(hObject,'string'));
if((temp < 0)|| (temp > 10) || (isnan(temp)))
    errordlg('Invalid value');
    set(hObject,'string',num2str(ds_para.SCSA(cur_vox)))
    return;
end
if(ind_allvox == 1)
    ds_para.SCSA(cur_vox) = temp;
else
   ds_para.SCSA(details_para.selected_voxels(1,:)) = temp;
end
set(preprocess_handles.ds_SCSA_slide,'value',ds_para.SCSA(cur_vox));
%--------------------------------------------------------------------------------------------------------------------------------
function denoise_Callback(hObject, eventdata)

global data_value;
global details_para;
global preprocess_handles;
global preprocess_sel_vox;
global FID_FD_copy;
global pha_cor_data;
global ds_para;
global a_var;
global ind_allvox;
global cur_vox;

% load(data_value.temp_filename);

freq_vector=[-1000,-90];
% set(preprocess_handles.fig,'Visible','off');

details_para.ds_count=details_para.ds_count+1;

% if details_para.ds_count==1
    data_value.FID_TD_ds=data_value.FID_TD;
    data_value.FID_FD_ds=data_value.FID_FD;
% end

global shif
shif=0;
Denoising_coeff=100;
gm=0.5;fs=1;
width_peaks=5;          %  The ratio w.r.t metabolite amplitude to set the threshold from which the peak can be selected
Th_peaks_ratio=2;        %  The with of the peak from its center(max values)

% h = waitbar(0,'Searching for optimum h value....');
cnt=1;
try
    if(ind_allvox == 1)
        if (ds_para.type(cur_vox)==1)
                cnt=cur_vox;
                org_data_unsuppressed=data_value.FID_TD_ds(:,cur_vox);
                data_value.FID_FD_org2(:,cur_vox)=data_value.FID_FD_ds(:,cur_vox);
                data_value.FID_TD_org2(:,cur_vox)=ifft(ifftshift( data_value.FID_FD_org2(:,cur_vox)));
                if (details_para.ds_flag(:,cur_vox)==0 || details_para.ds_flag(:,cur_vox)==2)
                    [yf_SCSA_Denoised, h_op, Nh]=SCSA_MRS_Denoising(500,details_para.ppm_referenced', real(data_value.FID_FD_ds(:,cur_vox)) , gm , fs , Th_peaks_ratio, width_peaks);
                    details_para.h_op(cur_vox)=h_op;
                    details_para.yf_SCSA_Denoised(:,cur_vox)=yf_SCSA_Denoised;
                end
                [yscsa,Nh]= Apply_SCSA_for_denoising(details_para.yf_SCSA_Denoised(:,cur_vox),details_para.ppm_referenced',details_para.h_op(cur_vox).*ds_para.SCSA(cur_vox),gm,fs)
                data_value.FID_FD(:,cur_vox)=yscsa;
                data_value.FID_TD(:,cur_vox)=ifft(ifftshift( data_value.FID_FD(:,cur_vox)));
%                 waitbar(1);
                details_para.ds_flag(:,cur_vox)=1;
        else
                cnt=cur_vox;
                org_data_unsuppressed=data_value.FID_TD_ds(:,cur_vox);
                data_value.FID_FD_org2(:,cur_vox)=data_value.FID_FD_ds(:,cur_vox);
                data_value.FID_TD_org2(:,cur_vox)=ifft(ifftshift( data_value.FID_FD_org2(:,cur_vox)));
                re_bl = ssa(real(fftshift(fft(org_data_unsuppressed))),ds_para.HLSVD(cur_vox)); % estimate the baseline for real part
                im_bl = ssa(imag(fftshift(fft(org_data_unsuppressed))),ds_para.HLSVD(cur_vox)); % estimate the baseline for imaginary part
                data_value.FID_FD(:,cur_vox)=re_bl + 1i*im_bl;
                data_value.FID_TD(:,cur_vox)=ifft(ifftshift( data_value.FID_FD(:,cur_vox)));
%                 waitbar(1);
                details_para.ds_flag(:,cur_vox)=1;
        end
        
    else

    if (ds_para.type(details_para.selected_voxels(1,:))==1)
        for i = 1:details_para.num_vox
            cnt=i;
            org_data_unsuppressed=data_value.FID_TD_ds(:,i);
            data_value.FID_FD_org2(:,i)=data_value.FID_FD_ds(:,i);
            data_value.FID_TD_org2(:,i)=ifft(ifftshift( data_value.FID_FD_org2(:,i)));
            if(details_para.ds_flag(:,i)==0 || details_para.ds_flag(:,i)==2)
                [yf_SCSA_Denoised, h_op, Nh]=SCSA_MRS_Denoising(500,details_para.ppm_referenced', real(data_value.FID_FD_ds(:,i)) , gm , fs , Th_peaks_ratio, width_peaks);
                 details_para.h_op(i)=h_op;
                 details_para.yf_SCSA_Denoised(:,i)=yf_SCSA_Denoised;
            end
            [yscsa,Nh]= Apply_SCSA_for_denoising(details_para.yf_SCSA_Denoised(:,i),details_para.ppm_referenced',details_para.h_op(i).*ds_para.SCSA(i),gm,fs);
            data_value.FID_FD(:,i)=yscsa;
            data_value.FID_TD(:,i)=ifft(ifftshift( data_value.FID_FD(:,i)));
%             waitbar(i/details_para.num_vox);
            details_para.ds_flag(:,i)=1;
        end
    else
        for i = 1:details_para.num_vox
            cnt=i;
            org_data_unsuppressed=data_value.FID_TD_ds(:,i);
            data_value.FID_FD_org2(:,i)=data_value.FID_FD_ds(:,i);
            data_value.FID_TD_org2(:,i)=ifft(ifftshift( data_value.FID_FD_org2(:,i)));
            re_bl = ssa(real(fftshift(fft(org_data_unsuppressed))),ds_para.HLSVD(i)); % estimate the baseline for real part
            im_bl = ssa(imag(fftshift(fft(org_data_unsuppressed))),ds_para.HLSVD(i)); % estimate the baseline for imaginary part
            data_value.FID_FD(:,i)=re_bl + 1i*im_bl;
            data_value.FID_TD(:,i)=ifft(ifftshift( data_value.FID_FD(:,i)));
%             waitbar(i/details_para.num_vox);
            details_para.ds_flag(:,i)=1;
        end
    end
    end
catch
%         close(h);
        details_para.ds_flag(cnt)=0;
        set(preprocess_handles.fig,'Visible','on');
        return;
end
% close(h);
set(preprocess_handles.fig,'Visible','on');
% FID_FD_copy=data_value.FID_FD;
pha_cor_data(:,details_para.selected_voxels(1,:))=data_value.FID_TD;
% seg = details_para.disp_seg;
% ppm = details_para.ppm_referenced;
% cla(preprocess_handles.axis);
% plot_var = real(data_value.FID_FD(seg(1):seg(2)));
% plot(preprocess_handles.axis,ppm(seg(1):seg(2))',plot_var);
% max_p = max(plot_var);
% mim_p = min(plot_var);
% axis(preprocess_handles.axis,[ppm(seg(1)), ppm(seg(2)),mim_p - 0.1*max_p, max_p + 0.1*max_p])
% set(preprocess_handles.axis,'XDir','reverse')
display_data()

a_var=data_value.FID_FD;
%-----------------------------------------------------------------------------------------------------------
function reload_callback(hObject, eventdata)
    
global details_para;
global data_value;
global pha_cor_data;
global preprocess_handles;
global apod_para;
global pha_para;
global ws_para;
global ds_para;
global ind_allvox;
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


try
    data_load_suppress_init();
catch
    return;
end
FID_FD_copy = data_value.FID_FD;
default_parameters.pha_0=0;
default_parameters.pha_1=0;
default_parameters.apod_para_lor=0;
default_parameters.apod_para_gaus=0;
default_parameters.apod_para_sigm=0;
default_parameters.apod_para_type=1;
default_parameters.ws_para_SCSA=17;
default_parameters.ws_para_HLSVD=10;
default_parameters.ws_para_type=1;
default_parameters.ds_para_SCSA=1;
default_parameters.ds_para_HLSVD=3;
default_parameters.ds_para_type=1;
ind_allvox=1;
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
pha_para.temp_pha0=default_parameters.pha_0*ones(details_para.num_vox,1);
pha_para.temp_pha1=default_parameters.pha_1*ones(details_para.num_vox,1);
pha_para.pha_0=default_parameters.pha_0*ones(details_para.num_vox,1);
pha_para.pha_1=default_parameters.pha_1*ones(details_para.num_vox,1);
apod_para.type = default_parameters.apod_para_type*ones(details_para.num_vox,1);
apod_para.lor = default_parameters.apod_para_lor*ones(details_para.num_vox,1);
apod_para.gaus = default_parameters.apod_para_gaus*ones(details_para.num_vox,1);
apod_para.sigm = default_parameters.apod_para_sigm*ones(details_para.num_vox,1);
apod_para.temp_lor = default_parameters.apod_para_lor*ones(details_para.num_vox,1);
apod_para.temp_gaus = default_parameters.apod_para_gaus*ones(details_para.num_vox,1);
apod_para.temp_sigm = default_parameters.apod_para_sigm*ones(details_para.num_vox,1);
ws_para.type = default_parameters.ws_para_type*ones(details_para.num_vox,1);
ws_para.SCSA = default_parameters.ws_para_SCSA*ones(details_para.num_vox,1);
ws_para.HLSVD = default_parameters.ws_para_HLSVD*ones(details_para.num_vox,1);
ds_para.type = default_parameters.ds_para_type*ones(details_para.num_vox,1);
ds_para.SCSA = default_parameters.ds_para_SCSA*ones(details_para.num_vox,1);
ds_para.HLSVD = default_parameters.ds_para_HLSVD*ones(details_para.num_vox,1);
pha_cor_data =  data_value.FID_TD;
hsvd_para.order=10;
scsa_para.order=17;
details_para.ds_flag=zeros(1,size(data_value.FID_TD,2));
if(isempty(apod_para)) % if pre-processing is not done previously load default pre-processing parameters
    
    apod_para.type = default_parameters.apod_para_type*ones(details_para.num_vox,1);
    apod_para.lor = default_parameters.apod_para_lor*ones(details_para.num_vox,1);
    apod_para.gaus = default_parameters.apod_para_gaus*ones(details_para.num_vox,1);
    apod_para.sigm = default_parameters.apod_para_sigm*ones(details_para.num_vox,1);
    ws_para.type = default_parameters.ws_para_type*ones(details_para.num_vox,1);
    ws_para.SCSA = default_parameters.ws_para_SCSA*ones(details_para.num_vox,1);
    ws_para.HLSVD = default_parameters.ws_para_HLSVD*ones(details_para.num_vox,1);
    ds_para.type = default_parameters.ds_para_type*ones(details_para.num_vox,1);
    ds_para.SCSA = default_parameters.ds_para_SCSA*ones(details_para.num_vox,1);
    ds_para.HLSVD = default_parameters.ds_para_HLSVD*ones(details_para.num_vox,1);
    apod_para.function_weight = ones(size(data_value.FID_FD));
    pha_para.pha_0 = default_parameters.pha_0*ones(details_para.num_vox,1);
    pha_para.pha_1 = default_parameters.pha_1*ones(details_para.num_vox,1);
    pha_para.temp_pha0=default_parameters.pha_0*ones(details_para.num_vox,1);
    pha_para.temp_pha1=default_parameters.pha_1*ones(details_para.num_vox,1);
    pha_cor_data =  data_value.FID_TD;
else % if pre-processing is done use the current pre-processing parameters 
    % and perform pre-processing using those values 
    for i = 1:length(details_para.selected_voxels)
        vox = details_para.selected_voxels(1,i);
%         pha_cor_data(:,vox) = phase_correction(FID_FD_copy(:,vox),pha_para.pha_0(vox),pha_para.pha_1(vox),details_para.fres); %phase correction is done on frequency domain non-fftshifted data
%         [temp2, function_weight] = apodize(pha_cor_data(:,vox),apod_para.type(vox),[details_para.Fs,apod_para.lor(vox),apod_para.gaus(vox),apod_para.sigm(vox)]);
%         sht_sig = fftshift(fft(temp2),1);
%         apod_para.function_weight(:, vox) = function_weight;
%         data_value.FID_TD(:,vox) = temp2;
%         data_value.FID_FD(:,vox) = sht_sig;
        apod_para.type=1*ones(details_para.num_vox,1);
        ws_para.type=1*ones(details_para.num_vox,1);
        ds_para.type=1*ones(details_para.num_vox,1);
    end
end
Pix_SS = get(0,'screensize');
x_max = Pix_SS(:,3)./2;
y_max = Pix_SS(:,4).*.75;
selected_vox = details_para.selected_voxels;
vox_no = cell(1,length(selected_vox));
for i = 1:length(selected_vox)
vox_no{i} =  num2str(selected_vox(1,i));
end
details_para.preprocess_going=1;
preprocess_handles.TF_FD_display = uibuttongroup(preprocess_handles.display_panel,'Units','normalized','Position',[0.04 0.76 .7 0.16],...
    'FontSize',9,'BackgroundColor',[0.68,0.92,.8],'BorderType','none','SelectionChangeFcn',@TD_FD_display_callback);
preprocess_handles.FD_disp_real = uicontrol(preprocess_handles.TF_FD_display,'Style','radiobutton','Units','normalized','Position',[.01 .1 .2 .9],'String',...
    'FD (real)','FontWeight','Bold','BackgroundColor',[0.68,0.92,.8]);
preprocess_handles.FD_disp_abs = uicontrol(preprocess_handles.TF_FD_display,'Style','radiobutton','Units','normalized','Position',[.3 .1 .2 .9],'String',...
    'FD (abs)','FontWeight','Bold','BackgroundColor',[0.68,0.92,.8]);
preprocess_handles.TD_disp = uicontrol(preprocess_handles.TF_FD_display,'Style','radiobutton','Units','normalized','String','FD (imag)','Position',...
    [.53 .1 .2 .9],'FontWeight','Bold','BackgroundColor',[0.68,0.92,.8]);
preprocess_handles.FD_disp_imag = uicontrol(preprocess_handles.TF_FD_display,'Style','radiobutton','Units','normalized','String','TD','Position',...
    [.75 .1 .2 .9],'FontWeight','Bold','BackgroundColor',[0.68,0.92,.8]);
% apodization panel
preprocess_handles.h2 = uipanel(preprocess_handles.fig,'Title','Apodization','Units','normalized','Position',[380/x_max 40/y_max (x_max-360)/(3.2*x_max) 320/y_max],...
    'BackgroundColor','yellow','FontSize',9, 'FontWeight','Bold');
preprocess_handles.Apodize = uicontrol(preprocess_handles.h2,'Style','text','Units','normalized','String','Apodization Type','Units','normalized',...
    'Position',[0.1 0.65 0.8 0.3],'BackgroundColor','yellow','FontSize',8,'FontWeight','Bold');
preprocess_handles.Apodize_type = uicontrol(preprocess_handles.h2,'Style','popup','String','Exponential|Gauss|Gauss-Exponential|sigmoid','Units','normalized','Position',...
    [0.1 0.59 0.8 0.3],'Value',apod_para.type(cur_vox), 'Callback',{@Apodize_type_Callback});

preprocess_handles.Apodize_Lor_text = uicontrol(preprocess_handles.h2, 'Style','text','Units','normalized','Position',[0.1 0.7 0.8 0.07],'String',...
    'order','FontSize',8,'FontWeight','Bold','BackgroundColor','yellow');
preprocess_handles.Apodize_Lor_edit = uicontrol(preprocess_handles.h2,'Style','edit','Units','normalized','Position',[0.1 0.65 0.8 0.07],'String',...
    num2str(apod_para.lor(cur_vox)),'Callback',{@lor_Callback}); 
preprocess_handles.Apodize_Lor_slide = uicontrol(preprocess_handles.h2,'Style', 'slider','Min',0,'Max',100,'Value',apod_para.lor(cur_vox),'Units','normalized',...
    'Position', [0.1 0.58 0.8 0.07],'SliderStep',[0.001,0.01],'Callback', {@slider_lor_Callback}); 
%     
preprocess_handles.Apodize_Gaus_text = uicontrol(preprocess_handles.h2,'Style','text','Units','normalized','Position',[0.1 0.45 0.8 0.07],'String',....
    'order','FontSize',8,'FontWeight','Bold','BackgroundColor','yellow');
preprocess_handles.Apodize_Gaus_edit = uicontrol(preprocess_handles.h2,'Style','edit','Units','normalized','Position',[0.1 0.4 0.8 0.07],...
    'String',num2str(apod_para.gaus(cur_vox)),'Callback',{@gaus_Callback});  
preprocess_handles.Apodize_Gaus_slide  = uicontrol(preprocess_handles.h2,'Style', 'slider','Min',0,'Max',1000,'Value',apod_para.gaus(cur_vox),'Units','normalized',...
    'Position', [0.1 0.33 0.8 0.07],'SliderStep',[0.0001,0.001],'Callback', {@slider_gaus_Callback});  

preprocess_handles.Apodize_Sig_text = uicontrol(preprocess_handles.h2,'Style','text','Units','normalized','Position',[0.1 0.2 0.8 0.07],'String',...
    'order','FontSize',8,'FontWeight','Bold','BackgroundColor','yellow');
preprocess_handles.Apodize_Sigm_edit = uicontrol(preprocess_handles.h2,'Style','edit','Units','normalized','Position',[0.1 0.15 0.8 0.07],...
    'String',num2str(apod_para.sigm(cur_vox)),'Callback',{@sigm_Callback});  
preprocess_handles.Apodize_Sigm_slide = uicontrol(preprocess_handles.h2,'Style', 'slider','Min',0,'Max',500,'Value',apod_para.sigm(cur_vox),'Units','normalized',...
    'Position', [0.1 0.08 0.8 0.07],'SliderStep',[0.001,0.01],'Callback', {@slider_sigm_Callback});

% ws panel
preprocess_handles.ws = uipanel(preprocess_handles.fig,'Title','Water suppression','Units','normalized','Position',[(380/x_max)+(x_max-360)/(3.2*x_max) 40/y_max (x_max-360)/(3.2*x_max) 320/y_max],...
    'BackgroundColor','yellow','FontSize',9, 'FontWeight','Bold');
preprocess_handles.ws_orient = uicontrol(preprocess_handles.ws,'Style','text','Units','normalized','String','WS Type','Units','normalized',...
    'Position',[0.1 0.65 0.8 0.3],'BackgroundColor','yellow','FontSize',8,'FontWeight','Bold');
preprocess_handles.ws_type = uicontrol(preprocess_handles.ws,'Style','popup','String','SCSA|HLSVD-PRO','Units','normalized','Position',...
    [0.1 0.59 0.8 0.3],'Value',ws_para.type(cur_vox), 'Callback',{@ws_type_Callback});

preprocess_handles.ws_SCSA_text = uicontrol(preprocess_handles.ws, 'Style','text','Units','normalized','Position',[0.1 0.7 0.8 0.07],'String',...
    'order','FontSize',8,'FontWeight','Bold','BackgroundColor','yellow');
preprocess_handles.ws_SCSA_edit = uicontrol(preprocess_handles.ws,'Style','edit','Units','normalized','Position',[0.1 0.65 0.8 0.07],'String',...
    num2str(ws_para.SCSA(cur_vox)),'Callback',{@SCSA_ws_Callback}); 
preprocess_handles.ws_SCSA_slide = uicontrol(preprocess_handles.ws,'Style', 'slider','Min',0,'Max',30,'Value',apod_para.lor(cur_vox),'Units','normalized',...
    'Position', [0.1 0.58 0.8 0.07],'SliderStep',[0.1,1],'Callback', {@slider_SCSA_Callback}); 
%     
preprocess_handles.ws_HLSVD_text = uicontrol(preprocess_handles.ws,'Style','text','Units','normalized','Position',[0.1 0.45 0.8 0.07],'String',....
    'order','FontSize',8,'FontWeight','Bold','BackgroundColor','yellow');
preprocess_handles.ws_HLSVD_edit = uicontrol(preprocess_handles.ws,'Style','edit','Units','normalized','Position',[0.1 0.4 0.8 0.07],...
    'String',num2str(ws_para.HLSVD(cur_vox)),'Callback',{@HLSVD_ws_Callback});  
preprocess_handles.ws_HLSVD_slide  = uicontrol(preprocess_handles.ws,'Style', 'slider','Min',0,'Max',30,'Value',apod_para.gaus(cur_vox),'Units','normalized',...
    'Position', [0.1 0.33 0.8 0.07],'SliderStep',[0.1,1],'Callback', {@slider_HLSVD_Callback}); 

% denoise panel
preprocess_handles.ds = uipanel(preprocess_handles.fig,'Title','Denoise','Units','normalized','Position',[(380/x_max)+(2*(x_max-360))/(3.2*x_max) 40/y_max (x_max-360)/(3.2*x_max) 320/y_max],...
    'BackgroundColor','yellow','FontSize',9, 'FontWeight','Bold');
preprocess_handles.ds_orient = uicontrol(preprocess_handles.ds,'Style','text','Units','normalized','String','Type','Units','normalized',...
    'Position',[0.1 0.65 0.8 0.3],'BackgroundColor','yellow','FontSize',8,'FontWeight','Bold');
preprocess_handles.ds_type = uicontrol(preprocess_handles.ds,'Style','popup','String','SCSA|SVD','Units','normalized','Position',...
    [0.1 0.59 0.8 0.3],'Value',ds_para.type(cur_vox), 'Callback',{@ds_type_Callback});

preprocess_handles.ds_SCSA_text = uicontrol(preprocess_handles.ds, 'Style','text','Units','normalized','Position',[0.1 0.7 0.8 0.07],'String',...
    'level','FontSize',8,'FontWeight','Bold','BackgroundColor','yellow');
preprocess_handles.ds_SCSA_edit = uicontrol(preprocess_handles.ds,'Style','edit','Units','normalized','Position',[0.1 0.65 0.8 0.07],'String',...
    num2str(ds_para.SCSA(cur_vox)),'Callback',{@SCSA_ds_Callback}); 
preprocess_handles.ds_SCSA_slide = uicontrol(preprocess_handles.ds,'Style', 'slider','Min',0.2,'Max',3,'Value',apod_para.lor(cur_vox),'Units','normalized',...
    'Position', [0.1 0.58 0.8 0.07],'SliderStep',[0.05,1],'Callback', {@slider_SCSA_ds_Callback}); 
%     
preprocess_handles.ds_HLSVD_text = uicontrol(preprocess_handles.ds,'Style','text','Units','normalized','Position',[0.1 0.45 0.8 0.07],'String',....
    'level','FontSize',8,'FontWeight','Bold','BackgroundColor','yellow');
preprocess_handles.ds_HLSVD_edit = uicontrol(preprocess_handles.ds,'Style','edit','Units','normalized','Position',[0.1 0.4 0.8 0.07],...
    'String',num2str(ds_para.HLSVD(cur_vox)),'Callback',{@HLSVD_ds_Callback});  
preprocess_handles.ds_HLSVD_slide  = uicontrol(preprocess_handles.ds,'Style', 'slider','Min',0,'Max',10,'Value',apod_para.gaus(cur_vox),'Units','normalized',...
    'Position', [0.1 0.33 0.8 0.07],'SliderStep',[0.1,1],'Callback', {@slider_HLSVD_ds_Callback});  

% phase correction panel
preprocess_handles.h3 = uipanel(preprocess_handles.fig,'Title','Phase correction','Units','normalized','Position',[20/x_max 160/y_max 360/x_max 200/y_max],...
    'BackgroundColor','yellow','FontSize',9, 'FontWeight','Bold');
preprocess_handles.pha_corr0_text = uicontrol(preprocess_handles.h3,'Style','text','Units','normalized','Position',[0.1 0.6 0.2 .2],'String',...
    'Zero-order phase correction ','FontSize',7.5,'FontWeight','Bold','BackgroundColor','yellow');
preprocess_handles.pha_corr1_text = uicontrol(preprocess_handles.h3,'Style','text','Units','normalized','Position',[0.1 0.3 0.2 0.2],'String',...
    'First-order phase correction','FontSize',7.5,'FontWeight','Bold','BackgroundColor','yellow');
preprocess_handles.pha_corr_value_text = uicontrol(preprocess_handles.h3,'Style','text','Units','normalized','Position',[0.31 0.61 0.2 .2],'String',...
    'Value','FontSize',8,'FontWeight','Bold','BackgroundColor','yellow');
preprocess_handles.pha_corr_slider_text = uicontrol(preprocess_handles.h3,'Style','text','Units','normalized','Position',[0.66 0.61 0.2 0.2],'String',...
    'Slider','FontSize',8,'FontWeight','Bold','BackgroundColor','yellow');
preprocess_handles.pha_corr0_value = uicontrol(preprocess_handles.h3,'Style','edit','Units','normalized','Position',[0.34 0.5 .15 .2],...
    'String',num2str(pha_para.pha_0(cur_vox)),'Callback',{@pha0_val_Callback});
preprocess_handles.pha_corr1_value = uicontrol(preprocess_handles.h3,'Style','edit','Units','normalized','Position',[0.34 0.2 .15 .2],...
    'String',num2str(pha_para.pha_1(cur_vox)),'Callback',{@pha1_val_Callback});
preprocess_handles.pha_corr0_slide = uicontrol(preprocess_handles.h3,'Style','slider','Min',-pi,'Max',pi,'Value',pha_para.pha_0(cur_vox),'Units','normalized','Position',...
    [0.55 0.5 .42 .2],'SliderStep',[0.001,0.01],'Callback', {@pha0_slider_Callback});  
preprocess_handles.pha_corr1_slide = uicontrol(preprocess_handles.h3,'Style','slider','Min',-10,'Max',10,'Value',pha_para.pha_1(cur_vox),'Units','normalized','Position',...
    [0.55 0.2 .42 .2],'SliderStep',[0.001,0.01],'Callback', {@pha1_slider_Callback});

preprocess_handles.water_suppress = uicontrol(preprocess_handles.fig,'Style', 'pushbutton','Units','normalized','Position',[((380/x_max)+(1*(x_max-360))/(3.0*x_max)) 70/y_max 80/x_max 25/y_max],...
'String', 'Water suppress','FontSize',7.5,'FontWeight','Bold', 'Interruptible','off','BusyAction','cancel','Callback',@water_suppress_Callback);
preprocess_handles.denoise = uicontrol(preprocess_handles.fig,'Style', 'pushbutton','Units','normalized','Position',[((380/x_max)+(2*(x_max-360))/(3.0*x_max)) 70/y_max 80/x_max 25/y_max],...
'String', 'Denoise','FontSize',7.5,'FontWeight','Bold', 'Interruptible','off','BusyAction','cancel','Callback',@denoise_Callback);
% preprocess_handles.wavelet = uicontrol(preprocess_handles.fig,'Style', 'pushbutton','Units','normalized','Position',[240/x_max 140/y_max 100/x_max 25/y_max],...
% 'String', 'Denoise','FontSize',9,'FontWeight','Bold', 'Interruptible','off','BusyAction','cancel','Callback',@denoise_Callback);
set(preprocess_handles.buttongroup,'SelectionChangeFcn',@ind_allvox_value);


if ind_allvox == 1
    set(preprocess_handles.buttongroup,'SelectedObject',preprocess_handles.ind_vox);
else
    set(preprocess_handles.buttongroup,'SelectedObject',preprocess_handles.all_vox);
end

set(preprocess_handles.ws_SCSA_slide,'value',ws_para.SCSA(cur_vox));
set(preprocess_handles.ws_HLSVD_slide,'value',ws_para.HLSVD(cur_vox));

set(preprocess_handles.ds_SCSA_slide,'value',ds_para.SCSA(cur_vox));
set(preprocess_handles.ds_HLSVD_slide,'value',ds_para.HLSVD(cur_vox));

if (apod_para.type(selected_vox(1,1)) == 1)
    set(preprocess_handles.Apodize_Lor_edit, 'Enable', 'on');
    set(preprocess_handles.Apodize_Lor_slide, 'Enable', 'on');
    set(preprocess_handles.Apodize_Sigm_edit, 'Enable', 'off');
    set(preprocess_handles.Apodize_Sigm_slide, 'Enable', 'off');
    set(preprocess_handles.Apodize_Gaus_edit, 'Enable', 'off');
    set(preprocess_handles.Apodize_Gaus_slide, 'Enable', 'off');
elseif (apod_para.type(selected_vox(1,1)) == 2)
    set(preprocess_handles.Apodize_Lor_edit, 'Enable', 'off');
    set(preprocess_handles.Apodize_Sigm_edit, 'Enable', 'off');
    set(preprocess_handles.Apodize_Lor_slide, 'Enable', 'off');
    set(preprocess_handles.Apodize_Sigm_slide, 'Enable', 'off');
    set(preprocess_handles.Apodize_Gaus_edit, 'Enable', 'on');
    set(preprocess_handles.Apodize_Gaus_slide, 'Enable', 'on');
elseif (apod_para.type(selected_vox(1,1)) == 3)
    set(preprocess_handles.Apodize_Lor_edit, 'Enable', 'on');
    set(preprocess_handles.Apodize_Sigm_edit, 'Enable', 'off');
    set(preprocess_handles.Apodize_Lor_slide, 'Enable', 'on');
    set(preprocess_handles.Apodize_Sigm_slide, 'Enable', 'off');
    set(preprocess_handles.Apodize_Gaus_edit, 'Enable', 'on');
    set(preprocess_handles.Apodize_Gaus_slide, 'Enable', 'on');
else
    set(preprocess_handles.Apodize_Lor_edit, 'Enable', 'off');
    set(preprocess_handles.Apodize_Sigm_edit, 'Enable', 'on');
    set(preprocess_handles.Apodize_Lor_slide, 'Enable', 'off');
    set(preprocess_handles.Apodize_Sigm_slide, 'Enable', 'on');
    set(preprocess_handles.Apodize_Gaus_edit, 'Enable', 'off');
    set(preprocess_handles.Apodize_Gaus_slide, 'Enable', 'off');
end

if (ws_para.type(selected_vox(1,1)) == 1)
    set(preprocess_handles.ws_SCSA_edit, 'Enable', 'on');
    set(preprocess_handles.ws_SCSA_slide, 'Enable', 'on');
    set(preprocess_handles.ws_HLSVD_edit, 'Enable', 'off');
    set(preprocess_handles.ws_HLSVD_slide, 'Enable', 'off');
elseif(ws_para.type(selected_vox(1,1)) == 2)
    set(preprocess_handles.ws_SCSA_edit, 'Enable', 'off');
    set(preprocess_handles.ws_SCSA_slide, 'Enable', 'off');
    set(preprocess_handles.ws_HLSVD_edit, 'Enable', 'on');
    set(preprocess_handles.ws_HLSVD_slide, 'Enable', 'on');
end

if (ds_para.type(selected_vox(1,1)) == 1)
    set(preprocess_handles.ds_SCSA_edit, 'Enable', 'on');
    set(preprocess_handles.ds_SCSA_slide, 'Enable', 'on');
    set(preprocess_handles.ds_HLSVD_edit, 'Enable', 'off');
    set(preprocess_handles.ds_HLSVD_slide, 'Enable', 'off');
elseif(ws_para.type(selected_vox(1,1)) == 2)
    set(preprocess_handles.ds_SCSA_edit, 'Enable', 'off');
    set(preprocess_handles.ds_SCSA_slide, 'Enable', 'off');
    set(preprocess_handles.ds_HLSVD_edit, 'Enable', 'on');
    set(preprocess_handles.ds_HLSVD_slide, 'Enable', 'on');
end
    

set(preprocess_handles.fig,'Visible','on');
set(preprocess_handles.display_seg_mim, 'Enable', 'on','string',num2str(details_para.disp_ppm_seg(1)))
set(preprocess_handles.display_seg_max, 'Enable', 'on','string',num2str(details_para.disp_ppm_seg(2)))
% set(preprocess_handles.main_val, 'Enable', 'off')


 
preprocess_display = 1;
set(preprocess_handles.fig,'CloseRequestFcn',@close_Callback);

% Display the selected voxel 
display_data()
waitfor(preprocess_handles.fig);

%--------------------------------------------------------------------------------------------------------------