function [Curr_share_N,Curr_share_G] = reallocate3d_initial(Gpos,Npos,strt_pts,Rpos,count)



for i=1:numel(strt_pts)
    euc(:,i)=euclidean(Rpos(strt_pts(i),:),Npos);
    eucG(:,i)=euclidean(Rpos(strt_pts(i),:),Gpos);
end

[~,assign_ind]=min(euc,[],2);
[~,assign_ind_G]=min(eucG,[],2);
[~,sharesizeN,~]=mode(assign_ind);
[~,sharesizeG,~]=mode(assign_ind_G);

Curr_share_N=zeros(count+numel(strt_pts),sharesizeN);
 Curr_share_G=zeros(count+numel(strt_pts),sharesizeG);
    for i=1:numel(strt_pts)
        Npoints=find(assign_ind==i);
        Gpoints=find(assign_ind_G==i);
        Curr_share_N(count+i,1:numel(Npoints))=Npoints;
        Curr_share_G(count+i,1:numel(Gpoints))=Gpoints;
%         Nposx(strt_pts(i),assign_ind==i)=indsPx(assign_ind==i);
%         Nposy(strt_pts(i),assign_ind==i)=indsPy(assign_ind==i);
%         Nposz(strt_pts(i),assign_ind==i)=indsPz(assign_ind==i);
%          Gposx(strt_pts(i),assign_ind_G==i)=indsGx(assign_ind_G==i);
%         Gposy(strt_pts(i),assign_ind_G==i)=indsGy(assign_ind_G==i);
%         Gposz(strt_pts(i),assign_ind_G==i)=indsGz(assign_ind_G==i);
    end
end



