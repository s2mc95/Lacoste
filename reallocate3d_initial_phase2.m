function [Curr_share_N,Curr_share_G] = reallocate3d_initial_phase2(Gpos,Npos,strt_pts,Rpos,count)



for i=1:numel(strt_pts)
    euc(:,i)=euclidean(Rpos(strt_pts(i),:),Npos);
    eucG(:,i)=euclidean(Rpos(strt_pts(i),:),Gpos);
end

[~,assign_ind]=min(euc,[],2);



[~,assign_ind_G]=min(eucG,[],2);

[~,sharesizeN,~]=mode(assign_ind);
[~,sharesizeG,~]=mode(assign_ind_G);
 Curr_share_N=zeros(count,sharesizeN);
 Curr_share_G=zeros(count,sharesizeG);
    for i=1:numel(strt_pts)
        Npoints=find(assign_ind==i);
        Gpoints=find(assign_ind_G==i);
        Curr_share_N(strt_pts(i),1:numel(Npoints))=Npoints;
        Curr_share_G(strt_pts(i),1:numel(Gpoints))=Gpoints;
    end
end



