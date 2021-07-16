function [Curr_share_N,Curr_share_G,repeat] = reallocate3d_phase2(Npos,Gpos,Curr_share_N,Curr_share_G,split_pts,Rpos,count,k)
indsPx=Npos(Curr_share_N(split_pts,:),1);
indsPy=Npos(Curr_share_N(split_pts,:),2);
indsPz=Npos(Curr_share_N(split_pts,:),3);

nanlocs=find(Curr_share_N(split_pts,:)==size(Npos,1));

indsGx=Gpos(Curr_share_G(split_pts,:),1);
indsGy=Gpos(Curr_share_G(split_pts,:),2);
indsGz=Gpos(Curr_share_G(split_pts,:),3);

nanlocsG=find(Curr_share_G(split_pts,:)==size(Gpos,1));

for i=1:k
    euc(:,i)=euclidean(Rpos(count+i,:),[indsPx,indsPy,indsPz]);
    eucG(:,i)=euclidean(Rpos(count+i,:),[indsGx,indsGy,indsGz]);
end

[~,assign_ind]=min(euc,[],2);
assign_ind(nanlocs)=k+1;
chk1=ismember(1:k,assign_ind);

[~,assign_ind_G]=min(eucG,[],2);
assign_ind_G(nanlocsG)=k+1;
chk1G=ismember(1:k,assign_ind_G);

if nnz(chk1)< k && nnz(chk1G)< k
    repeat=1;
    disp('repeating random split')
else
    repeat=0;% keyboard;
    totalpoints=0;
    for i=1:k
         Npoints=find(assign_ind==i);
        Gpoints=find(assign_ind_G==i);
        Curr_share_N(count+i,1:numel(Npoints))=Curr_share_N(split_pts,Npoints);
        Curr_share_G(count+i,1:numel(Gpoints))=Curr_share_G(split_pts,Gpoints);
        if ismember(size(Npos,1),Curr_share_N(count+i,1:numel(Npoints)))
            keyboard;
        end
%         Nposx(count+i,assign_ind==i)=indsPx(assign_ind==i);
%         Nposy(count+i,assign_ind==i)=indsPy(assign_ind==i);
%         Nposz(count+i,assign_ind==i)=indsPz(assign_ind==i);
%          Gposx(count+i,assign_ind_G==i)=indsGx(assign_ind_G==i);
%         Gposy(count+i,assign_ind_G==i)=indsGy(assign_ind_G==i);
%         Gposz(count+i,assign_ind_G==i)=indsGz(assign_ind_G==i);

totalpoints=totalpoints+numel(Npoints);
    end
     parent=Curr_share_N(split_pts,:);
 parent(parent==size(Npos,1))=0;
 if nnz(parent)~=totalpoints
     keyboard;
 end
end

end


