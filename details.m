function[slice_details,gone_out,total_slce] = details(B_pos,token,xx,yy,NX,cp_ind)
oo=1;
for tm=1:2:size(NX,1)
    if NX(tm,1)== size(xx,1)
            slce_idx{oo,1} = find(B_pos(:,2)>= yy(NX(tm,1),NX(tm,2)) & B_pos(:,2)<= yy(NX(tm+1,1),NX(tm+1,2))  & B_pos(:,1)==xx(NX(tm,1),NX(tm,2)));

    else
    slce_idx{oo,1} = find(B_pos(:,2)>= yy(NX(tm,1),NX(tm,2)) & B_pos(:,2)<= yy(NX(tm+1,1),NX(tm+1,2)) & B_pos(:,1)>=xx(NX(tm,1),NX(tm,2))& B_pos(:,1)<=xx(NX(tm+2,1),NX(tm+2,2)) );
    oo=oo+1;
    end
                      
end


%% 
total_slce=0;
oo=1;
pp=1;
for sze=1:size(slce_idx,1)
    total_slce=total_slce+numel(slce_idx{sze,1});
end
 
% slice_details = zeros(total_slce,3);
% gone_out= zeros(total_slce,3);

oo=1;
pp=1;
for outc = 1:size(slce_idx,1)
    for innrc = 1:size(slce_idx{outc,1},1)
        if size(find(cp_ind(:,1)== token(slce_idx{outc,1}(innrc,1))),1) > 0
            fth= find(cp_ind(:,1)== token(slce_idx{outc,1}(innrc,1)));
             token_fth= find(token(:,1)==cp_ind(fth,2));
             if size(token_fth,1) > 0
                 slice_details(oo,:) = [token(slce_idx{outc,1}(innrc,1)),outc,token(token_fth)];
                 oo=oo+1;
             else
                 gone_out(pp,:) = [token(slce_idx{outc,1}(innrc,1)),outc,cp_ind(fth,2)];
                 pp=pp+1;
             end
        else
            keyboard;
        end
    end
end
if size(gone_out,1) == total_slce
    slice_details=0;
end
end