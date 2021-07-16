clear
clc
load('D_BIG_phase2_just_grown_0.04%_Neurons.mat');
eta = 5e-5;
[Rpos_master,Level_ID,actv_inds,A,cp_ind,rooted_inds]=grow_tree_with_neural_act_2(Height,Npos,eta,branching_factor,Gpos,Rpos_master,Level_ID,actv_inds,A,cp_ind)
