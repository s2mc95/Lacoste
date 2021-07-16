clear;
clc;
%% CONTROL
% load('Control_growth_automated_Stopping_criteria_v0.mat')
% load('D_BIG_phase2_just_grown_0.043%_Neurons_4S.mat')
load('Lesioned_growth_automated_Stopping_criteria_with_settle_with_step_by_step_wt_51_1S.mat')
%% LESIONED
% load('Control_growth_automated_Stopping_criteria_with_settle_with_step_by_step_wt_42_1S.mat')
actv_inds = find(final_leaf_indx==1)';
leaf_pos = Rpos_master_settled(actv_inds,:);
spcf = find(leaf_pos(:,3)<0.55 & leaf_pos(:,3) >= 0.30);
spcf_pos = leaf_pos(spcf,:);
%%
% leaf_pos = Rpos_master(1:size(actv_inds,1),:);
% spcf = find(leaf_pos(:,3)<0.42 & leaf_pos(:,3) >= 0.38);
% spcf_pos = Rpos_master(spcf,:);


%%
these = find(Rpos_master_settled(:,3)<0.55 & Rpos_master_settled(:,3) >= 0.30);
these_pos = Rpos_master_settled(these,:);

these_cp_idx = cp_ind(these,:);

nm=size(these_cp_idx,1);
portion_cords = zeros(branching_factor,3,nm);

for e=1:nm
    portion_cords(1,:,e) = Rpos_master_settled(these_cp_idx(e,2),:) ;
    portion_cords(2,:,e) = Rpos_master_settled(these_cp_idx(e,1),:) ;
end

for ee=1:size(portion_cords,3)
    figure(35); plot3(portion_cords(:,1,ee),portion_cords(:,2,ee),portion_cords(:,3,ee),'-o','Color','r','MarkerSize',10,'MarkerFaceColor','#D9FFFF');
    hold on;
end
scatter3(spcf_pos(:,1),spcf_pos(:,2),spcf_pos(:,3),'b','filled')
title('LESIONED AFTER SETTLE');