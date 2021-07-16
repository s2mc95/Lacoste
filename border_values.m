function [final_X,final_Y,new_border,NX] = border_values(border_indices,range,border_leftend,border_rightend)

borders=zeros(40,40);

for ty=1:size(border_indices,1)
    borders(border_indices(ty,1),border_indices(ty,2))=1;
end
clear ty;
[S_Ynew,S_Xnew] = meshgrid(0.6/range:0.6/range:0.6,0.6/range:0.6/range:0.6);

new_border = borders(border_leftend:border_rightend,border_leftend:border_rightend);

final_Y = S_Ynew.* new_border ;
final_X = S_Xnew.* new_border ;

oo=1;
for t=1:size(new_border,1)
    ixx = find(new_border(t,:)==1);
    if size(ixx,2)> 0
        for tt=1:size(ixx,2)
            NX(oo,:) = [t,ixx(1,tt)];
            oo=oo+1;
        end
    end
end
clear oo t tt ixx
end