function [neural_output]=som_resp(Npos,wt,inp,L_MAX,L_MIN)
neural_output=zeros(size(Npos,1),1);

% cd F:\Github_team\Tree_growth_algorithm\Phase2_lacoste\tree_growth\som_functions;
chosen=[6,7,8,11,12,13,16,17,18];

ci=randi(numel(chosen));
inp_inx=chosen(ci);

[~,y_6] = movierespnew(wt,inp(inp_inx,:),1);
y_trun=y_6(11:23,11:23);
y_oneD=reshape(y_trun,numel(y_trun),1);
% cd F:\Github_team\Tree_growth_algorithm\Phase2_lacoste\tree_growth;
Neu_L4_inds=Npos(:,3)<L_MAX & Npos(:,3)>=L_MIN;
Neu_L4=Npos(Neu_L4_inds,:);

[y_new,x_new]=meshgrid(0.6/size(y_trun,1):0.6/size(y_trun,1):0.6,0.6/size(y_trun,2):0.6/size(y_trun,2):0.6);
NNx=reshape(x_new,numel(x_new),1);
NNy=reshape(y_new,numel(y_new),1);
Xbig=repmat(NNx',size(Neu_L4,1),1);
Ybig=repmat(NNy',size(Neu_L4,1),1);
dist=sqrt((Neu_L4(:,1)-Xbig).^2+(Neu_L4(:,2)-Ybig).^2);
[minval,minloc]=min(dist,[],2);
act_val4m_SOM=exp(-minval).*y_oneD(minloc);
neural_output(Neu_L4_inds,1)=act_val4m_SOM;
