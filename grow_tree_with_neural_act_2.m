function [Rpos_master,Level_ID,actv_inds,A,cp_ind,rooted_inds]=grow_tree_with_neural_act_2(Height,Npos,eta,branching_factor,Gpos,Rpos_master,Level_ID,actv_inds,A,cp_ind)

%   load('Lesi_itr80_Big_whisker.mat');
%  load('Big_whisker_SOM_sig4_Nrad20.mat');
Input_weights=2;
available_weights=10;
L_MAX = 0.55;
L_MIN = 0.30;

cd C:\Users\sarat\Documents\GitHub\Tree_growth_algorithm\Phase2_lacoste\SOM\som_final;
load('whisker_input.mat');
load(['wt_control_itr',num2str(Input_weights),'.mat']); %% always change this and line 232 together
% Yout=zeros(size(q,1),size(q,2));
Yout=zeros(40,40);
chosen=[6,7,8,11,12,13,16,17,18];
cd C:\Users\sarat\Documents\GitHub\Tree_growth_algorithm\Phase2_lacoste\tree_growth\som_functions;
for i=1:numel(chosen)
    inp_inx=chosen(i);
    [~,y_6] = movierespnew(wt,inp(inp_inx,:),1);%figure;imagesc(y_6);
    Yout=Yout+y_6;
end
Yout=Yout/numel(chosen);
Yout=Yout/max(max(Yout));
cd C:\Users\sarat\Documents\GitHub\Tree_growth_algorithm\Phase2_lacoste\tree_growth;
% stop_AT=round(0.6*size(Npos,1));
Nposx=Npos(:,1);
Nposy = Npos(:,2);
Nposz = Npos(:,3);

Gposx=Gpos(:,1);
Gposy = Gpos(:,2);
Gposz = Gpos(:,3);
Rpos=zeros(size(Rpos_master));
Rpos(actv_inds,:)=Rpos_master(actv_inds,:);

final_leaf_indx = zeros(size(Rpos,1),1);

count=numel(Level_ID);
k=branching_factor;

[Curr_share_N,Curr_share_G]=reallocate3d_initial_phase2(Gpos,Npos,actv_inds,Rpos,count);


i=0;
condition=1;

level_num=Level_ID;
extraL=size(Npos,1)+1;
Npos(extraL,:)=nan;

extraG=size(Gpos,1)+1;
Gpos(extraG,:)=nan;

Nposx(extraL)=nan;
Nposy(extraL)=nan;
Nposz(extraL)=nan;

Gposx(extraG)=nan;
Gposy(extraG)=nan;
Gposz(extraG)=nan;

[child_id,parent_id]=find(A(1:count,1:count));
%% Tree growth using wire Length minimization
while condition>0
    i=i+1;
    
    Curr_share_N(Curr_share_N==0)=extraL;
    Curr_share_G(Curr_share_G==0)=extraG;
    
    neural_act=som_resp(Npos,Yout,L_MAX,L_MIN);
    etaneur=0.0001*neural_act(Curr_share_N);

    distx=(etaneur+0.1*eta).*(Rpos(1:count,1)-Nposx(Curr_share_N));
    disty=(etaneur+0.1*eta).*(Rpos(1:count,2)-Nposy(Curr_share_N));
    distz =(etaneur+0.1*eta).*(Rpos(1:count,3)-Nposz(Curr_share_N));
    
    distx(isnan(distx))=0;
    disty(isnan(disty))=0;
    distz(isnan(distz))=0;
    
    Gdistx=(Rpos(1:count,1)-Gposx(Curr_share_G));
    Gdisty=(Rpos(1:count,2)-Gposy(Curr_share_G));
    Gdistz =(Rpos(1:count,3)-Gposz(Curr_share_G));
    
    
    
    Gdistx(isnan(Gdistx))=0;
    Gdisty(isnan(Gdisty))=0;
    Gdistz(isnan(Gdistz))=0;
    %
    %     delx=-1*((eta/10)*sum(distx,2)+(eta)*sum(Gdistx,2));
    %     dely=-1*((eta/10)*sum(disty,2)+(eta)*sum(Gdisty,2));
    %     delz= -1*((eta/10)*sum(distz,2)+(eta)*sum(Gdistz,2));
    delx=-1*(sum(distx,2)+(eta)*sum(Gdistx,2));
    dely=-1*(sum(disty,2)+(eta)*sum(Gdisty,2));
    delz= -1*(sum(distz,2)+(eta)*sum(Gdistz,2));
    
    
    Rpos(1:count,:)=Rpos(1:count,:)+[delx,dely,delz];
    
    chk_diff=abs(abs(delx)+abs(dely)+abs(delz));
    %     if i>=2 && restart==0
    %     delta_error(:,i)=chk_diff_prev-chk_diff;
    %     end
    chk_diff_prev=chk_diff;
    restart=0;
    
    %     if numel(find(abs(Rpos(:,1))>0))<numel(actv_inds)
    %         keyboard;
    %     end
    actv_inds=find(abs(Rpos(:,1))>0);
    
    rooted_inds = find(abs(Rpos(:,1))==0);
%     if numel(actv_inds)==0                       unc
%         keyboard;
%     end
    Rpos_master(actv_inds,:)= Rpos(actv_inds,:);
    disp(['Epochs no: ',num2str(i)])
    disp(['The number of active vasc nodes= ',num2str(numel(actv_inds))]);
    
    
    split_pts=actv_inds(chk_diff(actv_inds)<1e-9);
    %     if size(actv_inds,1)>=stop_AT
    %
    %         condition=0;
    %     end
    if size(actv_inds,1)<1
        
        condition=0;
    end
    tempn=Curr_share_N(split_pts,:);
    tempn(tempn==extraL)=0;
    %     tempn(isnan(tempn))=0;
    tempn2=sum(tempn~=0,2);
    if numel(split_pts)>0 && numel(tempn2)>0
        for levj=1:numel(find(child_id>0))
            level_num(child_id(levj))=level_num(parent_id(levj))+1;
        end
        for j=1:numel(split_pts)
            repeat=1;
            avg_neural_act_subreg=mean(neural_act(Curr_share_N(split_pts(j),:)));
            neural_act_threshold=0.06;                                                      %% PARAM
            if avg_neural_act_subreg<neural_act_threshold
                final_leaf_indx(split_pts(j))=1;
                Rpos(split_pts(j),:)=0;
                Curr_share_N(split_pts(j),:)=0;
                Curr_share_G(split_pts(j),:)=0;
                
            elseif tempn2(j)>1 && avg_neural_act_subreg>=neural_act_threshold
                A(count+1:count+k,split_pts(j))=1;
                if size(A,2)<size(A,1)
                    A(:,size(A,2)+1:size(A,1))=0;
                elseif size(A,2)>size(A,1)
                    disp('Impossible situation');
                    keyboard;
                end
                while repeat ==1
                    
                    incrx=-1+2*rand(k,3);
                    Rpos(count+1:count+k,:)=Rpos(split_pts(j),:)+0.001*incrx;
                    Rpos_master(count+1:count+k,:)=Rpos(count+1:count+k,:);
                    
                    [Curr_share_N,Curr_share_G,repeat] = reallocate3d_phase2(Npos,Gpos,Curr_share_N,Curr_share_G,split_pts(j),Rpos,count,k);
                    
                end
                
                
                
                
                Rpos(split_pts(j),:)=0;
                Curr_share_N(split_pts(j),:)=0;
                Curr_share_G(split_pts(j),:)=0;
                maximum_nonzero_N=max(sum(Curr_share_N~=0,2));
                maximum_nonzero_G=max(sum(Curr_share_G~=0,2));
                Curr_share_N=Curr_share_N(:,1:maximum_nonzero_N);
                Curr_share_G=Curr_share_G(:,1:maximum_nonzero_G);
                restart=1;
                count=count+k;
                clear child_id parent_id;
                [child_id,parent_id]=find(A(1:count,1:count));
                
                cp_ind(child_id,1)=child_id;
                cp_ind(child_id,2)=parent_id;
            end
        end
        for levj=1:numel(find(child_id>0))
            level_num(child_id(levj))=level_num(parent_id(levj))+1;
        end
    end
    
    %% for plotting uncomment section below
    %     if rem(i,100)==0
    %         xnode=[Rpos_master(1:count,1);Nposx];
    %         ynode=[Rpos_master(1:count,2);Nposy];
    %         znode=[Rpos_master(1:count,3);Nposz];
    %         c=ones(numel(xnode),1);
    %         s=50*ones(numel(xnode),1);
    %         %         c(count+1:end)=15;
    %         for kk=1:numel(actv_inds)
    %             c(actv_inds(kk))=10+2*kk;
    %             s(actv_inds(kk))=150;
    %             loc=find(BigNposx(actv_inds(kk),:)>0);
    %             c(count+loc)=10+2*kk;
    %         end
    %         s(count+1:end)=150;
    %         figure(2);scatter3(xnode(1:count),ynode(1:count),znode(1:count),s(1:count),c(1:count),'filled');hold on;
    %         scatter3(xnode(count+1:end),ynode(count+1:end),znode(count+1:end),s(count+1:end),c(count+1:end),'*');
    %         hold off;
    %         xlim([min(xnode),max(xnode)]); ylim([min(ynode),max(ynode)]);zlim([min(znode),max(znode)]);pause(0.1);
    %         figure(3);plot(digraph(A(1:count,1:count)))
    %     end
    %%
%     if rem(i,10000)==0 && sum(final_leaf_indx)>=1
%         Rpos_Final=Rpos_master(1:count,:);
%         cp_final=cp_ind(1:count,:);
%         A_final=A(1:count,1:count);
%         eta_wire =1.5e-5;
%         eta_path=1.2e-5;
%         eta_neuron = 7.5e-5;                       %%
%         Level_ID=level_num(1:count);
%         path_epochs=1000;                        %%
%         rooted_inds=rooted_inds(find(rooted_inds(:,1)<= count));
%         Leaf_Nodes=find(final_leaf_indx==1);
%         No_roots=find(cp_final(:,2)==1)-1;
%         [Rpos_master,DN2V,RMSE_mean,NVpair,RMSE_4m_delVals,eta_path,eta_wire,eta_neuron]=settle_path_length_witn_neu_act(No_roots,Height,Gpos,A_final,Rpos_Final,Level_ID,Npos,Leaf_Nodes,path_epochs,cp_final,eta_path,eta_wire,eta_neuron,rooted_inds,neural_act);
%         
%         Rpos(actv_inds,:)=Rpos_master(actv_inds,:);
%     end
    %%
    if rem(i,25)==0 && Input_weights<available_weights
        Input_weights=Input_weights+2;
        cd C:\Users\sarat\Documents\GitHub\Tree_growth_algorithm\Phase2_lacoste\SOM\som_final;
        load(['wt_control_itr',num2str(Input_weights),'.mat']);
        
        cd C:\Users\sarat\Documents\GitHub\Tree_growth_algorithm\Phase2_lacoste\tree_growth\som_functions;
        for blnd=1:numel(chosen)
            inp_inx=chosen(blnd);
            [~,y_6] = movierespnew(wt,inp(inp_inx,:),1);%figure;imagesc(y_6);
            Yout=Yout+y_6;
        end
        Yout=Yout/numel(chosen);
        Yout=Yout/max(max(Yout));
        cd C:\Users\sarat\Documents\GitHub\Tree_growth_algorithm\Phase2_lacoste\tree_growth;
    end
    
end
Npos_recovered=Npos;
clear Npos;
Npos=Npos_recovered(1:extraL-1,:);
save('Control_growth_automated_Stopping_criteria_with_settle_with_step_by_step_wt_33_NS.mat');
Rpos_Final=Rpos_master(1:count,:);
cp_final=cp_ind(1:count,:);
A_final=A(1:count,1:count);
eta_wire =1.5e-4;
eta_path=1.2e-4;
eta_neuron = 5e-4;                       %%
Level_ID=level_num(1:count);
path_epochs=20000;                        %%
rooted_inds=rooted_inds(find(rooted_inds(:,1)<= count));
Leaf_Nodes=find(final_leaf_indx==1);
No_roots=find(cp_final(:,2)==1)-1;
[Rpos_master_settled,DN2V,RMSE_mean,NVpair,RMSE_4m_delVals,eta_path,eta_wire,eta_neuron]=settle_path_length_witn_neu_act(No_roots,Height,Gpos,A_final,Rpos_Final,Level_ID,Npos,Leaf_Nodes,path_epochs,cp_final,eta_path,eta_wire,eta_neuron,rooted_inds,neural_act);


Level_ID=level_num(1:count);
% save('Control_growth_automated_Stopping_criteria_v1_starting0.01%neurons.mat');
save('Control_growth_automated_Stopping_criteria_with_settle_with_step_by_step_wt_33_1S.mat');
% save('Control_growth_automated_Stopping_criteria_with_settle_with_step_by_step_wt_25_1S.mat');
% yu=~isnan(BigNposx);
% for ju = 1:size(Npos,1)
%     Pair_nv(ju,1)= ju;
%     Pair_nv(ju,2) = find(yu(:,ju)==1);
% end