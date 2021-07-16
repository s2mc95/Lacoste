function [L_R_split]= Volume_proportion(out,XX,YY,Rp_M,NX,L_R_split,L4Z_max,L4Z_min)
%BC - BORDER - CHILD DISTANCE
%BP = BORDER - PARENT DISTANCE
% NX = INDICES OF BORDERS WRT MESHGRID  

% out(:,1) = child 
%out(:,3)  = parent

%out(:,2) = the box, child is in 

for tt=1:size(out,1)
    if out(tt,3)~= 0
        act_len = sqrt(sum((Rp_M(out(tt,3),:)-Rp_M(out(tt,1),:)).^2,2));
        
        %%
        if Rp_M(out(tt,3),3) < L4Z_min                            
                                                    %% parent below box %%    
            BC = Rp_M(out(tt,1),3) - L4Z_min;                     
            BP = L4Z_min - Rp_M(out(tt,3),3);                       
            ratio = BC/(BC+BP);
            ln_contrib = ratio * act_len;
            if ln_contrib < 0
                keyboard;
            end
            L_R_split(tt,1) =  ln_contrib ;
            
           
        elseif Rp_M(out(tt,3),3) > L4Z_max  
                                                      %% parent above box             
            BP = Rp_M(out(tt,3),3) - L4Z_max;
            BC = L4Z_max - Rp_M(out(tt,1),3);
            ratio = BC/(BC+BP);
            ln_contrib = ratio * act_len;
            if ln_contrib < 0
                keyboard;
            end
            L_R_split(tt,1) =  ln_contrib ;
            
            %%
        else
            if Rp_M(out(tt,3),1) < XX(NX(((2*out(tt,2)-1)),1),NX(((2*out(tt,2)-1)),2))
                BC = Rp_M(out(tt,1),1)- XX(NX(((2*out(tt,2)-1)),1),NX(((2*out(tt,2)-1)),2));
                BP = XX(NX(((2*out(tt,2)-1)),1),NX(((2*out(tt,2)-1)),2)) - Rp_M(out(tt,3),1) ;
                ratio = BC/(BC+BP);
                ln_contrib = ratio * act_len;
                if ln_contrib < 0
                    keyboard;
                end
                L_R_split(tt,1) =  ln_contrib ;
            
            elseif Rp_M(out(tt,3),1) > XX(NX(((2*out(tt,2)+1)),1),NX(((2*out(tt,2)+1)),2))
                BC = XX(NX(((2*out(tt,2)+1)),1),NX(((2*out(tt,2)+1)),2)) - Rp_M(out(tt,1),1);
                BP = Rp_M(out(tt,3),1) - XX(NX(((2*out(tt,2)+1)),1),NX(((2*out(tt,2)+1)),2));
                ratio = BC/(BC+BP);
                ln_contrib = ratio * act_len;
                if ln_contrib < 0
                    keyboard;
                end

                L_R_split(tt,1) =  ln_contrib ; 
            
            else
                if Rp_M(out(tt,3),2)< YY(NX(((2*out(tt,2)-1)),1),NX(((2*out(tt,2)-1)),2))
                    BC = Rp_M(out(tt,1),2)- YY(NX(((2*out(tt,2)-1)),1),NX(((2*out(tt,2)-1)),2));
                    BP = YY(NX(((2*out(tt,2)-1)),1),NX(((2*out(tt,2)-1)),2)) - Rp_M(out(tt,3),2);
                    ratio = BC/(BC+BP);
                    ln_contrib = ratio * act_len;
                    if ln_contrib < 0
                        keyboard;
                    end
                    L_R_split(tt,1) =  ln_contrib ;
                else
                    BC = YY(NX(((2*out(tt,2))),1),NX(((2*out(tt,2))),2)) - Rp_M(out(tt,1),2);
                    BP = Rp_M(out(tt,3),2) - YY(NX(((2*out(tt,2))),1),NX(((2*out(tt,2))),2));
                    ratio = BC/(BC+BP);
                    ln_contrib = ratio * act_len;
                    if ln_contrib < 0
                        keyboard;
                    end
                    L_R_split(tt,1) =  ln_contrib ;
                end
            end
        end
    end
end
end