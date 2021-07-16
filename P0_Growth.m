clear all;
 close all;
clc;
% cd F:\BHADRA\GitHub\CNN_ANVN\Phase1\Tsai_work\full_data;
file=csvread('10_d_neurons_stp_1.csv');
compare_vessel=csvread('10_d_vasculature_stp_1.csv');
non_neurons2=csvread('10_d_non_neurons_v2.csv');
% cd F:\BHADRA\GitHub\CNN_ANVN\Phase1\Tsai_work\Bhadra_version_Tsai_final;
microvasc=load('10_d_microvasculature_step_size_1.csv');

Xaxis=file(end,1);
xaxtemp=find(abs(non_neurons2(:,1)-Xaxis)<1e-5);
non_neurons=non_neurons2(1:xaxtemp,:);
Length = 0.6;% mm
Width  = 0.6 ;%mm
Height = 0.7;%mm
tissue_Height=0.05 ;%mm
th_mu_m=tissue_Height*1e3;

%% Unit adjustment of loaded file
given_resolution=1e-3*(file(2,1)-file(1,1)); %% value in file is in micrometer
required_resolution=round(tissue_Height/given_resolution);
Glia_given_resolution=1e-3*(non_neurons(2,1)-non_neurons(1,1)); %% value in file is in micrometer
Glia_required_resolution=round(tissue_Height/Glia_given_resolution);
file(:,2)=file(:,2)*1e5; % unit adjustment
non_neurons(:,2)=non_neurons(:,2)*1e5; % unit adjustment
%% Obtain the neural distribution
[Neuronpos,model_distribution,endN]=distribute_neurons(Length,Width,Height,tissue_Height,file,required_resolution);
[Gliapos,Glia_model_distribution,endG]=distribute_NON_Neu(Length,Width,Height,tissue_Height,non_neurons,Glia_required_resolution);


all_p_c = Neuronpos(1:end,:);
glia_p_c= Gliapos(1:end,:);
all_p = cell2mat(Neuronpos(1:end,:));
all_glia=cell2mat(Gliapos(1:end,:));
% taking a smaller box for debug purposes
small_L=0.6;small_W=0.6;
unit_box_vol=small_L*small_W*tissue_Height;
% % p=find(all_p(:,1)<small_L & all_p(:,2)<small_W & all_p(:,3)>0);
% % pg=find(all_glia(:,1)<small_L & all_glia(:,2)<small_W & all_glia(:,3)>0);
% % Npos = all_p(p,:);
% % Gpos=all_glia(pg,:);
% p=find(((Length/2)-(small_L/2))<all_p(:,1) & all_p(:,1) <((Length/2)+(small_L/2))& ((Width/2)-(small_L/2))<all_p(:,2)& all_p(:,2)<((Width/2)+(small_L/2))& all_p(:,3)>0);
% pg=find(((Length/2)-(small_L/2))<all_glia(:,1)& all_glia(:,1)<((Length/2)+(small_L/2))& ((Width/2)-(small_L/2))<all_glia(:,2)&all_glia(:,2)<((Width/2)+(small_L/2))& all_glia(:,3)>0); 
Npos = all_p;%(p,:);
Gpos=all_glia;%(pg,:);
%%
% figure(1);scatter3(Npos(:,1),Npos(:,2),Npos(:,3),'b','filled');
% figure(2);scatter3(Gpos(:,1),Gpos(:,2),Gpos(:,3),'b','filled');
% figurenum=10;
% % plot_neurons_in_box(small_L,small_W,Height,Npos,figurenum)
% plot_neurons_in_box_middle(Length,Width,Height,small_L,small_W,Npos,figurenum)
% figure(10);scatter3(Npos(:,1),Npos(:,2),Npos(:,3),'b','filled');
% title({['Small section of nodes, ',num2str(small_L),'mm x ',num2str(small_W),'mm x ',num2str(Height),'mm'],['No: neurons= ',num2str(size(Npos,1))]});
% figure(6);
%  plot(file(:,1),file(:,2)*1e-5,'b');
%  ylabel('Neural density (1e5 mm^{-3})');
%  xlabel('Depth below cortical surface (\mum)')
%  hold on;
%  plot(file(1:required_resolution:endN,1),model_distribution(:,1)*1e-5);
%  legend("Experimental","Model")
%  figure(7);
% plot(non_neurons(:,1),non_neurons(:,end)*1e-5,'b');
%  ylabel('Non Neural cell density (1e5 mm^{-3})');
%  xlabel('Depth below cortical surface (\mum)')
%  hold on;
%  plot(non_neurons(1:Glia_required_resolution:endG,1),Glia_model_distribution(:,1)*1e-5);
%  legend("Experimental","Model")
%% Tree growth using wirelength minimization

% root_node=[small_L/2,small_W/2,Height];
num_of_roots=44;
 xroot=0.001+(Length-0.001)*rand(num_of_roots,1);
 yroot=0.001+(Width-0.001)*rand(num_of_roots,1);
 zroot=Height*ones(num_of_roots,1);
 root_node=[xroot,yroot,zroot];
% root_node=[Length/4,Width/2,Height;Length-Length/4,Width/2,Height];
eta_wire=5e-5;
branching_factor=2;
stop_AT=round(0.045*size(Npos,1));
[Rpos_master,Level_ID,actv_inds,A,cp_ind,rooted_inds]=grow_tree_wire(stop_AT,Npos,root_node,eta_wire,branching_factor,Gpos);
Rpos_master_just_wire=Rpos_master;
A=sparse(A);
  save('D_BIG_phase2_just_grown_0.1%_Neurons.mat');
 %% Shaping up results
  count=numel(Level_ID);
 Rpos_Final=Rpos_master(1:count,:);
cp_final=cp_ind(1:count,:);
A_final=A(1:count,1:count);
% % % Pathlength minimization
% % eta_wire =1.5e-5;
% % eta_path=1.2e-5;
% % eta_neuron = 7e-5;
% % 
% % path_epochs=20000;
% % 
% % rooted_inds=rooted_inds(find(rooted_inds(:,1)<= count));
% % tic
% % [Rpos_Final,DN2V,RMSE_mean,NVpair,RMSE_4m_delVals,eta_path,eta_wire,eta_neuron]=settle_path_length(Height,Gpos,A_final,Rpos_Final,Level_ID,Npos,actv_inds,path_epochs,cp_final,eta_path,eta_wire,eta_neuron,rooted_inds);
% % 
% % toc
% % figure(20);  yyaxis left;
% % plot(RMSE_mean);hold on;
% % xlabel('epochs','FontSize',12);
% % ylabel('Delta derived RMS Error','FontSize',12);
% % hold on;yyaxis right;
% % plot(eta_wire);hold on; plot(eta_path);hold on; plot(eta_neuron);
% % ylabel("Variation of eta");
% % legend("RMS error","eta wire","eta path","eta neuron");
% % title({'RMS error'});
% % 
% % figure(21);  yyaxis left;
% % plot(RMSE_4m_delVals);hold on;
% % xlabel('epochs','FontSize',12);
% % ylabel('Delta derived RMS Error','FontSize',12);
% % hold on;yyaxis right;
% % plot(eta_wire);hold on; plot(eta_path);hold on; plot(eta_neuron);
% % ylabel("Variation of eta");
% % legend("Delta derived RMS error","eta wire","eta path","eta neuron");
% % title({'Delta derived RMS error'});
% % 
% % 
% %  
% % 
% % 
% % % length and radius calculation
% % Length_vasc = calculate_length(Rpos_Final,cp_ind,Level_ID);
% % r0=0.085;
% % kr=1;
% % radius_ref0=define_radius(r0,kr,Level_ID);
% % radius0=radius_ref0(Level_ID);
% % 
% % % Volume per box calculation
% % vasc_zcords=Rpos_Final(:,3);
% % if numel(find(vasc_zcords<0))>0
% %     warning('negative z coordinate, settling error');
% % end
% % [box_num,z_box]=find_box_no(vasc_zcords,Height,tissue_Height);
% % [~,BOX_INFO]=volume_per_box_temp(Rpos_Final,radius0,Length_vasc,cp_ind,box_num,z_box,Level_ID);
% % 
% % till=size(BOX_INFO,1);
% %  radius_ref=iterative_approx2(radius_ref0,BOX_INFO,compare_vessel([1,th_mu_m:th_mu_m:(till-1)*th_mu_m],2)*unit_box_vol);
% % %   radius_ref=grad_desc2(radius_ref0,BOX_INFO,compare_vessel([1,th_mu_m:th_mu_m:(till-1)*th_mu_m],2)*unit_box_vol);
% % radius=radius_ref(Level_ID);
% % clear BOX_INFO
% % [vol_per_box,BOX_INFO]=volume_per_box_temp(Rpos_Final,radius,Length_vasc,cp_ind,box_num,z_box,Level_ID);
% % 
% % density=vol_per_box(end:-1:1)/unit_box_vol;
% % ZBOX=z_box*1000;%=(Height-z_box(end:-1:2))*1000;
% % Rms_after_settle=rms(compare_vessel([1,th_mu_m:th_mu_m:(till-1)*th_mu_m],2)-density);
% % microvol=micro_vasc_volume_per_box(Rpos_Final,radius,Length_vasc,cp_ind,box_num,z_box,Level_ID);
% % micro_density=microvol(end:-1:1)/unit_box_vol;
% %     figure(5);
% %     % yyaxis right;
% %      plot(ZBOX,density(1:size(ZBOX),1),'-.k');
% %      ylabel('vascular volume (V/V)');
% %      hold on;
% %      plot(compare_vessel(1:end,1),compare_vessel(1:end,2),'r');hold on;
% %      plot(ZBOX,micro_density(1:size(ZBOX),1),'-.c');
% %      hold on;
% %      plot(microvasc(:,1),microvasc(:,2),'b');
% %      
% %      ylim([0,0.03])
% %       xlabel('Depth below cortical surface (\mum)');
% %      title({'Vascular densities along the depth below cortical surface RMSE = ',num2str(Rms_after_settle)});
% %     legend('Model Vasc','Experimental Vasc','Model Micro','exp micro'); hold on;
% %     
% % 
% %     
% % % figure(6);
% % %  plot(file(:,1),file(:,2)*1e-5,'b');
% % %  ylabel('Neural density (1e5 mm^{-3})');
% % %  xlabel('Depth below cortical surface (\mum)')
% % %  hold on;
% % %  plot(file(1:required_resolution:end,1),model_distribution(:,1)*1e-5);
% % %  legend("Experimental","Model")
% %  
% %  box_details_view(BOX_INFO,radius_ref,Npos,Rpos_Final,Height,tissue_Height,box_num,actv_inds,th_mu_m/1000);
%  temp_run;
%  %% plot_before_settling
%  plot_before_settling(th_mu_m,Npos,radius_ref0,Rpos_master_just_wire,radius0,cp_ind,Height,tissue_Height,Level_ID,compare_vessel,unit_box_vol,actv_inds)
% 
% 
% %% Connect to neurons and plot
% 
% NVpair(:,1)=1:size(NVpair,1);
% A=A_final;
% Aold=A; 
% pos_all(:,1)=nonzeros(Rpos_Final(:,1));
% pos_all(:,2)=nonzeros(Rpos_Final(:,2));
% pos_all(:,3)=nonzeros(Rpos_Final(:,3));
% 
% if size(A,1)<count+size(Npos,1) || size(A,2)<count+size(Npos,1)
%     A(count+1:count+size(Npos,1),count+1:count+size(Npos,1))=0;
% end
% 
% for de=1:size(Npos,1)
%     A(count+NVpair(de,1),NVpair(de,2))=1;
% end
% 
% dA = A(1:count+size(Npos,1),1:count+size(Npos,1));
% noOfNeurons = size(Npos,1);
% t.leaf_idx = actv_inds;
% t.leaf_pos = pos_all(actv_inds,:);
% X = [pos_all(:,1);Npos(:,1) ];%;X_random]
% Y = [pos_all(:,2); Npos(:,2)];%Y_random ]
% Z = [pos_all(:,3);Npos(:,3)];%Z_random ];
% t.dA = dA;
% t.X = X;
% t.Y = Y;
% t.Z = Z;
% 
%  adj = Aold(1:count,1:count);
%  check_r = sum(adj,1);
% leafnum = numel(find(check_r==0));
% Ntot=size(adj,2);
% 
% [plot_cords]= find_plot_cords(adj,X,Y,Z,Ntot,branching_factor,3,leafnum);
% 
% 
% hold on;
% figure(33)
% scatter3(Npos(:,1),Npos(:,2),Npos(:,3),'b','filled');%xlim([S_X(1)-0.5,S_X(end)+0.5]); ylim([S_Y(1)-0.5,S_Y(end)+0.5]);zlim([S_Z(1)-0.5,S_Z(end)+0.5]);
% hold on;
% 
% for i=1:size(plot_cords,3)
%     figure(33); plot3(plot_cords(:,1,i),plot_cords(:,2,i),plot_cords(:,3,i),'-o','Color','r','MarkerSize',10,'MarkerFaceColor','#D9FFFF');
%     hold on;
% end
% figure(33);legend('Neurons','Vessels');%'Root Node')
% zlabel({'Depth (in mm)'});
% ylabel({'Length (in mm)'});
% xlabel({'Width (in mm)'});
% title({'Tree growth for 0.1x0.1x0.8 mm3'});
% figure(33); scatter3(Rpos_Final(1,1),Rpos_Final(1,2),Rpos_Final(1,3),300,'r','filled');hold on;
% 
% figure(4);plot(digraph(dA),'Layout','layered');
% p=plot(digraph(dA),'Layout','layered');
% WW=count+1:size(dA,1);
% highlight(p,WW,'NodeColor','red');title('dA plot : check');
% %% Save data
% % save('Tree_new_1','ZBOX','density','Npos','Rpos_master','plot_cords','small_L','small_W','Height');
% % save('Tree_29_new');
% save('BIG_whisker_barrel_10_D_with_micro_20%Npos.mat');
%% plot tree

pos_all(:,1)=nonzeros(Rpos_master(1:count,1));
pos_all(:,2)=nonzeros(Rpos_master(1:count,2));
pos_all(:,3)=nonzeros(Rpos_master(1:count,3));
dA = A(1:count,1:count);
noOfNeurons = size(Npos,1);
t.leaf_idx = actv_inds;
t.leaf_pos = pos_all(1:size(actv_inds,1),:);
X = [pos_all(:,1) ];%;X_random]
Y = [pos_all(:,2)];%Y_random ]
Z = [pos_all(:,3)];%Z_random ];
t.dA = dA;
t.X = X;
t.Y = Y;
t.Z = Z;

 adj = A(1:count,1:count);
 check_r = sum(adj,1);
leafnum = numel(find(check_r==0));
Ntot=size(adj,2);

[plot_cords]= find_plot_cords(adj,X,Y,Z,Ntot,branching_factor,3,leafnum);


hold on;
figure(33)
scatter3(Npos(:,1),Npos(:,2),Npos(:,3),'b','filled');%xlim([S_X(1)-0.5,S_X(end)+0.5]); ylim([S_Y(1)-0.5,S_Y(end)+0.5]);zlim([S_Z(1)-0.5,S_Z(end)+0.5]);
hold on;

for i=1:size(plot_cords,3)
    figure(34); plot3(plot_cords(:,1,i),plot_cords(:,2,i),plot_cords(:,3,i),'-o','Color','r','MarkerSize',10,'MarkerFaceColor','#D9FFFF');
    hold on;
end

hold on;
for rootnd=1:size(root_node,1)
 figure(34); scatter3(Rpos_Final(rootnd,1),Rpos_Final(rootnd,2),Rpos_Final(rootnd,3),300,'r','filled');hold on;
end

[Rpos_master,Level_ID,actv_inds,A,cp_ind,rooted_inds]=grow_tree_with_neural_act(Npos,eta_wire,branching_factor,Gpos,Rpos_master,Level_ID,actv_inds,A,cp_ind);