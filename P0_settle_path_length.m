function [Rpos_master,DN2V,RMSE_mean,NVpair,RMSE_4m_delVals,eta_path,eta_wire,eta_neuron]=P0_settle_path_length(Height,Gpos,A,Rpos_master,Level_ID,Npos,actv_inds,epochs,cp_ind,eta_path_max,eta_wire_max,eta_neuron_max,rooted_inds,No_roots)
RMSE_prev=1e15;
stopper=1;
stopper2=1;
Rpos_initial = Rpos_master;
for path_epochs=1:epochs
    path_epochs
    eta_path(path_epochs)=get_exp_eta(path_epochs,eta_path_max,eta_path_max,epochs);
    eta_wire(path_epochs)=get_exp_eta(path_epochs,eta_wire_max,eta_wire_max,epochs);
    eta_neuron(path_epochs)=get_exp_eta(path_epochs,eta_neuron_max,eta_neuron_max,epochs);
    %% move each path along parents to root and leaf to neurons
    Rpos_previous = Rpos_master;
    Dpath=zeros(size(Rpos_master));
    Dpl=zeros(size(Rpos_master));
    Dpl_up=zeros(size(Rpos_master));
    Dpl_down=zeros(size(Rpos_master));
    Dwl=zeros(size(Rpos_master));
    Dwl_up=zeros(size(Rpos_master));
    Dwl_down=zeros(size(Rpos_master));
    Dn=zeros(size(Rpos_master));
    for pji=1:numel(actv_inds)
        path_nodes=zeros(Level_ID(actv_inds(pji)),size(Npos,2));
        path_parent=zeros(Level_ID(actv_inds(pji)),size(Npos,2));
        
        Pid=zeros(Level_ID(actv_inds(pji)),1);
        Pid(end)=actv_inds(pji);
        path_parent(end,:)=Rpos_master(Pid(end),:);
        for pq=1:Level_ID(actv_inds(pji))-1
            Pid(Level_ID(actv_inds(pji))-pq)=cp_ind(Pid(Level_ID(actv_inds(pji))-pq+1),2);
            path_parent(Level_ID(actv_inds(pji))-pq,:)=Rpos_master( Pid(Level_ID(actv_inds(pji))-pq),:);
        end
        
        path_nodes(1:end-1,:)=path_parent(2:end,:);
        del_P=path_parent(1:end-1,:)-path_nodes(1:end-1,:);
        Rpos_master(Pid(2:end),:)= Rpos_master(Pid(2:end),:)+eta_path(path_epochs)*del_P;
        Dpath(Pid(2:end),:)=Dpath(Pid(2:end),:)+del_P;
        Dpl_up(Pid(2:end),:)=Dpl_up(Pid(2:end),:)+del_P;
        Dpl(Pid(2:end),:)= Dpl(Pid(2:end),:)+del_P;
        %%
        path_children=zeros(Level_ID(actv_inds(pji))-2,size(Npos,2));
        path_children=Rpos_master(Pid(2:end-1),:);
        del_C = path_children - path_parent(3:end,:);
        Rpos_master(Pid(2:end-1),:)= Rpos_master(Pid(2:end-1),:)+eta_path(path_epochs)*del_C;
        Dpath(Pid(2:end-1),:)=Dpath(Pid(2:end-1),:)+del_C;
        Dpl_down(Pid(2:end-1),:)=Dpl_down(Pid(2:end-1),:)+del_C;
        Dpl(Pid(2:end-1),:)= Dpl(Pid(2:end-1),:)+del_C;
        %%
        %     path_check{pji}=Pid;
        checkpoint=Level_ID(Pid)-(1:Level_ID(actv_inds(pji)))';
        if numel(find(checkpoint~=0))>0
            error('Levels are not in systematic order. Perform a debug');
        end
        %      neu_pair_loc=find(Pair_nv(:,2)==Pid(end));
        %     neu_pair=Pair_nv(neu_pair_loc,1);
        %    del_2N=Rpos_master(Pid(end),:)-Npos(neu_pair,:);
        %     Rpos_master(Pid(end),:)=Rpos_master(Pid(end),:)+eta_path*del_2N;
        
        
        clear Pid;
        
    end
    
for cot = No_roots+1:numel(rooted_inds)
        parent_loc=Rpos_master(rooted_inds(cot),:);
        [ch_no]=find(A(:,rooted_inds(cot))==1);
        ch_loc=Rpos_master(ch_no,:);
        delta_dis=sum((ch_loc-parent_loc),1);
        Rpos_master(rooted_inds(cot),:)=Rpos_master(rooted_inds(cot),:)+ eta_wire(path_epochs)*delta_dis;
        Dpath(rooted_inds(cot),:)=Dpath(rooted_inds(cot),:)+delta_dis;
        Dwl_down(rooted_inds(cot),:)=Dwl_down(rooted_inds(cot),:)+delta_dis;
        Dwl(rooted_inds(cot),:)=Dwl(rooted_inds(cot),:)+delta_dis;
    end
    
    for pot = No_roots+1:numel(rooted_inds)
        child_loc = Rpos_master(rooted_inds(pot),:);
        p_loc = Rpos_master(cp_ind(rooted_inds(pot),2),:);
        d_dis = p_loc - child_loc;
        Rpos_master(rooted_inds(pot),:)=Rpos_master(rooted_inds(pot),:)+ eta_wire(path_epochs)* d_dis;
        Dpath(rooted_inds(pot),:)=Dpath(rooted_inds(pot),:)+d_dis;
        Dwl_up(rooted_inds(cot),:)=Dwl_up(rooted_inds(cot),:)+d_dis;
        Dwl(rooted_inds(cot),:)=Dwl(rooted_inds(cot),:)+d_dis;
    end
    
    
    if size(find(isnan(Rpos_master(:,1))),1) >= 1
        keyboard;
    end
   %% Each neuron pulls its closest vessel
    for nid=1:size(Npos,1)
        eucl_N2V=euclidean(Npos(nid,:),Rpos_master(actv_inds,:));
        [DN2V(nid,1),Nearest_vessel_loc]=min(eucl_N2V);
        Nearest_vessel=actv_inds(Nearest_vessel_loc);
        del_2N=-(Rpos_master(Nearest_vessel,:)-Npos(nid,:));
        NVpair(nid,1)=nid;
        NVpair(nid,2)=Nearest_vessel;
        Rpos_master(Nearest_vessel,:)=Rpos_master(Nearest_vessel,:)+ eta_neuron(path_epochs)*del_2N;
        Dpath(Nearest_vessel,:)=Dpath(Nearest_vessel,:)+ del_2N ;
        Dn(Nearest_vessel,:)=Dn(Nearest_vessel,:)+ del_2N ;
    end
    if ~isnan(Gpos)
        for gid=1:size(Gpos,1)
            eucl_G2V=euclidean(Gpos(gid,:),Rpos_master(actv_inds,:));
            [DG2V(nid,1),G_Nearest_vessel_loc]=min(eucl_G2V);
            G_Nearest_vessel=actv_inds(G_Nearest_vessel_loc);
            del_2G=-(Rpos_master(G_Nearest_vessel,:)-Gpos(gid,:));
            GVpair(gid,2)=G_Nearest_vessel;
            Rpos_master(G_Nearest_vessel,:)=Rpos_master(G_Nearest_vessel,:)+ eta_neuron(path_epochs)*del_2G;
            Dpath(G_Nearest_vessel,:)=Dpath(G_Nearest_vessel,:)+ del_2G ;
            Dn(G_Nearest_vessel,:)=Dn(G_Nearest_vessel,:)+ del_2G ;
        end
    end
    outliers=find(Rpos_master(:,3)>Height);
    if numel(outliers)>0
        Rpos_master(outliers,3)=Height;
    end
    RMSE_mean(path_epochs,1) = sqrt(mean(mean((Rpos_master-Rpos_previous).^2)));
    RMSE_4m_delVals(path_epochs,1) =sqrt(mean(mean(Dpath).^2));
    DN2V_store(path_epochs,1)=mean(DN2V);
    
    RMSE_pl(path_epochs,1)=sqrt(mean(mean(Dpl).^2));
    RMSE_wl(path_epochs,1)=sqrt(mean(mean(Dwl).^2));
    RMSE_n(path_epochs,1)=sqrt(mean(mean(Dn).^2));
    
    
    RMSE_pl_up(path_epochs,1)=sqrt(mean(mean(Dpl_up).^2));
    RMSE_wl_up(path_epochs,1)=sqrt(mean(mean(Dwl_up).^2));
    RMSE_Leaf_up(path_epochs,1)=sqrt(mean(mean(Dpl_up(actv_inds)+Dwl_up(actv_inds)).^2));
    
    RMSE_pl_down(path_epochs,1)=sqrt(mean(mean(Dpl_down).^2));
    RMSE_wl_down(path_epochs,1)=sqrt(mean(mean(Dwl_down).^2));
    RMSE_Leaf_down(path_epochs,1)=sqrt(mean(mean(Dpl_down(actv_inds)+Dwl_down(actv_inds)+Dn(actv_inds)).^2));
    if (mean(DN2V)<0.15 && mean(DN2V)>0.13) || (RMSE_prev-RMSE_mean(path_epochs,1))<6e-9
        if (RMSE_prev-RMSE_mean(path_epochs,1))<6e-9
            stopper=stopper+1;
        else
            stopper=0;
        end
        if RMSE_prev<=RMSE_mean(path_epochs,1)
            stopper2=stopper2+1;
        else
            stopper2=0;
        end
        if stopper==100 || stopper2==100
            if stopper==100
                disp('converged');
            else
                disp('diverged');
            end
        break
        end
    end

    RMSE_prev=RMSE_mean(path_epochs,1);
end
