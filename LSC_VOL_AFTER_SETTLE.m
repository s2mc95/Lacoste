clear;

%% change lines 9 and 60
%% CONTROL
load('26UN_NDX.mat');
[xx,yy,new_borders,NX]=border_values(unl_NDX,28,8,35);
% load('Control_growth_automated_Stopping_criteria_with_settle_with_step_by_step_wt_27_NS.mat')
load('Control_growth_automated_Stopping_criteria_with_settle_with_step_by_step_wt_51_1S.mat')% 
% load('D_BIG_phase2_just_grown_0.05%_Neurons.mat')
%wt_20 <><> 0.55,0.30;
%%
L4Z_max = 0.55;
L4Z_min = 0.30;
token = find(Rpos_master_settled(:,3)<L4Z_max & Rpos_master_settled(:,3)>=L4Z_min);
B_pos=Rpos_master_settled(token,:);

%%
[slice_details,gone_out,total_slce] = details(B_pos,token,xx,yy,NX,cp_ind);

%%
out_g= gone_out(1:nnz(gone_out(:,1)),:);
rem  =  find(out_g(:,3)==0);
if rem >0
    out_g(rem,:)=[];
else
    out_g=gone_out;
end

%%
L_R=zeros(size(slice_details,1),2);
r0=0.085;
kr=1;
radius_ref0=define_radius(r0,kr,Level_ID);
radius0=radius_ref0(Level_ID);
if size(slice_details,1)-1 > 0
    for lo=1:size(slice_details,1)
        L_R(lo,1)=sqrt(sum((Rpos_master_settled(slice_details(lo,1),:)-Rpos_master_settled(slice_details(lo,3),:)).^2,2));
    end
    
%     L_R(:,2)=radius0(slice_details(:,3));
end

%%
L_R_split = zeros(size(out_g,1),2);
L_R_split = Volume_proportion(out_g,xx,yy,Rpos_master_settled,NX,L_R_split,L4Z_max,L4Z_min);
% L_R_split(:,2)= radius0(out_g(:,3));

%%
Vol_T = L_R(:,1);
Vol_S = L_R_split(:,1);
%%
SUM = sum(Vol_S)+ sum(Vol_T);
Control_vol = ((0.6-0.0214)*(0.4071-0.2143)*(L4Z_max- L4Z_min));

Control_Branch_Density = total_slce/Control_vol
Control_Volume_Density = SUM/Control_vol
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LESIONED
load('new_lei_border.mat');
[xx_l,yy_l,new_borders_l,NX_l]=border_values(new_lei_border,28,8,35);
% load('Lesioned_growth_1.mat');
load('Lesioned_growth_automated_Stopping_criteria_with_settle_with_step_by_step_wt_51_1S.mat')

%%
L4Z_max = 0.55;
L4Z_min = 0.30;
token_l = find(Rpos_master_settled(:,3)<L4Z_max & Rpos_master_settled(:,3)>=L4Z_min);
B_pos_l=Rpos_master_settled(token_l,:);

%%
[slice_details_l,gone_out_l,total_slce_l] = details(B_pos_l,token_l,xx_l,yy_l,NX_l,cp_ind);

%%
out_g_l= gone_out_l(1:nnz(gone_out_l(:,1)),:);
rem  =  find(out_g_l(:,3)==0);
if rem >0
    out_g_l(rem,:)=[];
else
    out_g_l=gone_out_l;
end

%%
L_R_l=zeros(size(slice_details_l,1),2);

for lo=1:size(slice_details_l,1)
    L_R_l(lo,1)=sqrt(sum((Rpos_master_settled(slice_details_l(lo,1),:)-Rpos_master_settled(slice_details_l(lo,3),:)).^2,2));
end
r0=0.085;
kr=1;
radius_ref0=define_radius(r0,kr,Level_ID);
radius0=radius_ref0(Level_ID);
% L_R_l(:,2)=radius0(slice_details_l(:,3));
%%
L_R_split_l = zeros(size(out_g_l,1),2);
L_R_split_l = Volume_proportion(out_g_l,xx_l,yy_l,Rpos_master_settled,NX_l,L_R_split_l,L4Z_max,L4Z_min);
% L_R_split_l(:,2)= radius0(out_g_l(:,3));

%%
Vol_T_l = L_R_l(:,1);
Vol_S_l= L_R_split_l(:,1);

%%
SUM_l = sum(Vol_S_l)+ sum(Vol_T_l);
Lesioned_Volume = ((0.0215*(0.4071-0.1714))+(0.0214*(0.3857-0.1929))+((0.5571-0.0643)*(0.3643-0.2143))+(0.0214*(0.3857-0.1929))+(0.0214*(0.4071-0.1714)))*(L_MAX-L_MIN);

Lesioned_Branch_Density = total_slce_l/Lesioned_Volume
Lesioned_Volume_Density = SUM_l/Lesioned_Volume



Branch_Ratio =Control_Branch_Density /Lesioned_Branch_Density
Volume_Ratio =Control_Volume_Density /Lesioned_Volume_Density 

% keyboard;