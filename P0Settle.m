clear
load('D_BIG_phase2_just_grown_0.043%_Neurons.mat')
count=numel(Level_ID);
Rpos_Final=Rpos_master(1:count,:);
cp_final=cp_ind(1:count,:);
A_final=A(1:count,1:count);
%% Pathlength minimization
eta_wire = 1.15e-4;
eta_path= 1.05e-4;
eta_neuron = 2.5e-4;

path_epochs=5000;

rooted_inds=rooted_inds(find(rooted_inds(:,1)<= count));

[Rpos_master,DN2V,RMSE_mean,NVpair,RMSE_4m_delVals,eta_path,eta_wire,eta_neuron]=P0_settle_path_length(Height,Gpos,A_final,Rpos_Final,Level_ID,Npos,actv_inds,path_epochs,cp_final,eta_path,eta_wire,eta_neuron,rooted_inds,num_of_roots);
%                                                                                                 (Height,Gpos,A,Rpos_master,Level_ID,Npos,actv_inds,epochs,cp_ind,eta_path_max,eta_wire_max,eta_neuron_max,rooted_inds)
 save('D_BIG_phase2_just_grown_0.043%_Neurons_9S.mat');