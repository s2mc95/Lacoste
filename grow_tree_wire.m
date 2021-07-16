function [Rpos_master,Level_ID,actv_inds,A,cp_ind,rooted_inds]=grow_tree_wire(stop_AT,Npos,root_node,eta,branching_factor,Gpos)
Nposx=Npos(:,1);
Nposy = Npos(:,2);
Nposz = Npos(:,3);


% BigNposx=sparse(zeros(numel(Nposx)+1500,numel(Nposx)));
% BigNposy=sparse(BigNposx);
% BigNposz =sparse(BigNposx);

% BigNposx(2,:)=Nposx;
% BigNposy(2,:)=Nposy;
% BigNposz(2,:) = Nposz ;

Gposx=Gpos(:,1);
Gposy = Gpos(:,2);
Gposz = Gpos(:,3);

% BigGposx=sparse(zeros(numel(Nposx)+1500,numel(Gposx)));
% BigGposy=sparse(BigGposx);
% BigGposz =sparse(BigGposx);
%
% BigGposx(2,:)=Gposx;
% BigGposy(2,:)=Gposy;
% BigGposz(2,:)= Gposz ;

Rpos=zeros(size(Npos,1)+1500,3);
A=zeros(size(Rpos,1));
% A=zeros(size(root_node,1)+10);
% A(2,1)=1;
for rootnd=1:size(root_node,1)
    % Rpos(rootnd,:)=root_node(rootnd,:);
    Rpos_master(rootnd,:)=root_node(rootnd,:);
    Rpos(rootnd+size(root_node,1),:)=root_node(rootnd,:);
    A(rootnd+size(root_node,1),rootnd)=1;
    strt_pts(rootnd)=rootnd+size(root_node,1);
end
[Curr_share_N,Curr_share_G]=reallocate3d_initial(Gpos,Npos,strt_pts,Rpos,size(root_node,1));

count=size(Curr_share_N,1);
k=branching_factor;

movable=ones(size(Rpos,1),1);
movable(1:size(root_node,1))=0;
% Rpos_master=Rpos;
% Rpos(1,:)=zeros(1,3);
i=0;
condition=1;

level_num=zeros(size(A,1),1);

cp_ind=zeros(size(level_num,1),2);
[child_id,parent_id]=find(A);
cp_ind(child_id,1)=child_id;
cp_ind(child_id,2)=parent_id;
for rootnd=1:size(root_node,1)
    cp_ind(rootnd,1)=rootnd;
    level_num(rootnd)=1;
end

% level_num(1)=1;% root node is 1;
% Rpos=sparse(Rpos);
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

restart=0;
%% Tree growth using wire Length minimization
while condition>0
    i=i+1;
    
    Curr_share_N(Curr_share_N==0)=extraL;
    Curr_share_G(Curr_share_G==0)=extraG;
    
    distx=(movable(1:count).*Rpos(1:count,1)-Nposx(Curr_share_N));
    disty=(movable(1:count).*Rpos(1:count,2)-Nposy(Curr_share_N));
    distz =(movable(1:count).*Rpos(1:count,3)-Nposz(Curr_share_N));
    
    distx(isnan(distx))=0;
    disty(isnan(disty))=0;
    distz(isnan(distz))=0;
    
    Gdistx=(movable(1:count).*Rpos(1:count,1)-Gposx(Curr_share_G));
    Gdisty=(movable(1:count).*Rpos(1:count,2)-Gposy(Curr_share_G));
    Gdistz =(movable(1:count).*Rpos(1:count,3)-Gposz(Curr_share_G));
    
    
    
    Gdistx(isnan(Gdistx))=0;
    Gdisty(isnan(Gdisty))=0;
    Gdistz(isnan(Gdistz))=0;
    %
    delx=-1*((eta/10)*sum(distx,2)+(eta)*sum(Gdistx,2));
    dely=-1*((eta/10)*sum(disty,2)+(eta)*sum(Gdisty,2));
    delz= -1*((eta/10)*sum(distz,2)+(eta)*sum(Gdistz,2));
    
    
    Rpos(1:count,:)=Rpos(1:count,:)+[delx,dely,delz];
    
    chk_diff=abs(abs(delx)+abs(dely)+abs(delz));
    %     if i>=2 && restart==0
    %     delta_error(:,i)=chk_diff_prev-chk_diff;
    %     end
    chk_diff_prev=chk_diff;
    restart=0;
    
    
    actv_inds=find(abs(Rpos(:,1))>0);
    rooted_inds = find(abs(Rpos(:,1))==0);
    if numel(actv_inds)==0
        keyboard;
    end
    Rpos_master(actv_inds,:)= Rpos(actv_inds,:);
    disp(['Epochs no: ',num2str(i)])
    disp(['The number of active vasc nodes= ',num2str(numel(actv_inds))]);
    
    
    split_pts=actv_inds(chk_diff(actv_inds)<1e-9);
    if size(actv_inds,1)>=stop_AT
        
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
            
            if tempn2(j)>1
                A(count+1:count+k,split_pts(j))=1;
                while repeat ==1
                    
                    incrx=-1+2*rand(k,3);
                    Rpos(count+1:count+k,:)=Rpos(split_pts(j),:)+0.001*incrx;
                    Rpos_master(count+1:count+k,:)=Rpos(count+1:count+k,:);
                    
                    [Curr_share_N,Curr_share_G,repeat] = reallocate3d_phase2(Npos,Gpos,Curr_share_N,Curr_share_G,split_pts(j),Rpos,count,k);
                    
                end
                
                [child_id,parent_id]=find(A(1:count,1:count));
                
                cp_ind(child_id,1)=child_id;
                cp_ind(child_id,2)=parent_id;
                Rpos(split_pts(j),:)=0;
                Curr_share_N(split_pts(j),:)=0;
                Curr_share_G(split_pts(j),:)=0;
                maximum_nonzero_N=max(sum(Curr_share_N~=0,2));
                maximum_nonzero_G=max(sum(Curr_share_G~=0,2));
                Curr_share_N=Curr_share_N(:,1:maximum_nonzero_N);
                Curr_share_G=Curr_share_G(:,1:maximum_nonzero_G);
                restart=1;
                count=count+k;
            end
        end
        
    end
    for levj=1:numel(find(child_id>0))
        level_num(child_id(levj))=level_num(parent_id(levj))+1;
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
    
end
Level_ID=level_num(1:count);
% yu=~isnan(BigNposx);
% for ju = 1:size(Npos,1)
%     Pair_nv(ju,1)= ju;
%     Pair_nv(ju,2) = find(yu(:,ju)==1);
% end