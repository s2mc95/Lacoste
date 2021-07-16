function [Neuronpos,model_distribution,endpt]=distribute_neurons(Length,Width,Height,tissue_Height,file,int)
volume_of_tissue = Length*Width*tissue_Height; %mm3
%%
X = [0;Length;Length;0;0];
Y = [0;0;Width;Width;0];
Z = [0;0;0;0;0];
Z_end =Z+Height;
findloc=abs(file(:,1)-Height*1000);
endpt=find(findloc==min(findloc));
% fifty=file((1:int:size(file,1))',:);
fifty=file((1:int:endpt)',:);
tsai_window_density = fifty(:,2);

tsai_window_number  = round(tsai_window_density * volume_of_tissue);
model_distribution = tsai_window_number / volume_of_tissue;
tot_neu=sum(tsai_window_number);

hold on;
figure(1);
hold on;
plot3(X,Y,Z,'Color', 'b', 'LineWidth', 1);
plot3(X,Y,Z_end,'Color', 'b', 'LineWidth', 1);
set(gca,'View',[-35,35]);
hold on;
for ll = 1:size(X,1)-1
    line([X(ll),X(ll)],[Y(ll),Y(ll)],[Z(ll),Z_end(ll)],'Color', 'b', 'LineWidth', 1);
    hold on;
end
xlabel('Tissue Length (mm)','fontweight','bold','fontsize',8)
ylabel('Tissue thickness(mm)','fontweight','bold','fontsize',8)
zlabel('Tissue depth(mm)')
title({[num2str(Length),'mm X ',num2str(Width),'mm X ',num2str(Height),'mm'];['Total neurons=',num2str(tot_neu)]});



%%

cook=1;
hold on;
figure(1);
for den =size(tsai_window_number,1)-1:-1:1
    CstNp(:,1) = Length*rand(tsai_window_number(cook,1),1);
    CstNp(:,2) = Width*rand(tsai_window_number(cook,1),1);
    CstNp(:,3) = ((tissue_Height*den)-(tissue_Height*(den-1)))*rand(tsai_window_number(cook,1),1)+(tissue_Height*(den-1));
    scatter3(CstNp(:,1),CstNp(:,2),CstNp(:,3));
    Neuronpos{cook,1} = [CstNp(:,1),CstNp(:,2),CstNp(:,3)];
    
    cook=cook+1;
    clear CstNp
end
end