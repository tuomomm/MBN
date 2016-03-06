%drawmotres: Draw the single motif promotion figure.
%Tuomo M?ki-Marttunen, 2013-2016
function drawmotres()

load Mots_mean.mat
ps = 0.025:0.025:0.5
cols16 = colorsredtolila(17)*0.7;
styles = {'-','--','-.',':','-'}
uniquepatterns = [   0     0     0     0     0     0;...
                     0     0     0     0     0     1;...
                     0     0     0     0     1     1;...
                     0     0     0     1     0     1;...
                     0     0     0     1     1     0;...
                     0     0     1     0     1     0;...
                     0     0     0     1     1     1;...
                     0     0     1     0     1     1;...
                     0     1     0     1     0     1;...
                     0     1     1     0     0     1;...
                     0     0     1     1     1     1;...
                     0     1     0     1     1     1;...
                     0     1     1     0     1     1;...
                     0     1     1     1     1     0;...
                     0     1     1     1     1     1;...
                     1     1     1     1     1     1]
nedges = 1:6;
figure;

for myind = 3:3
 inds_compare = find(sum(uniquepatterns,2)==myind);
 cols = [.5 .5 .5;cols16(3+(myind-2):3:14,:); .5 .5 .5; 0 0 0];
 for j=1:length(inds_compare)
  SP(j) = subplot('position',[0.09 1-0.21*j 0.24 0.18]);
  MotsToCompare = zeros(length(ps),length(inds_compare),length(inds_compare));
  MotsToCompared = zeros(length(ps),length(inds_compare),length(inds_compare));
  for ip=1:length(ps)
    compare_inds = setdiff(1:length(inds_compare),j);
    RelMots = reshape(Motsm(1,inds_compare(j),ip,inds_compare,1),1,length(inds_compare));
    RelMotsd = reshape(Motsd(1,inds_compare(j),ip,inds_compare,1),1,length(inds_compare));
    RelRandMots = mean(RandMotsm(1,inds_compare(j),ip,:));
    MotsToCompare(ip,:,j) = RelMots%;/mean(RelRandMots);
    MotsToCompared(ip,:,j) = RelMotsd%;/mean(RelRandMots);
  end
  h=semilogy(ps,reshape(mean(RandMotsm(1,inds_compare(j),:,:),4),1,length(ps)),'color',[0.5 0.5 0.5],'linestyle',styles{1},'linewidth',2);
  hold on;
  reshape(mean(RandMotsm(1,inds_compare(j),:,:),4),1,length(ps))
  for k=[1:j-1 j+1:length(inds_compare) j]
    h=semilogy(ps,MotsToCompare(:,k,j),'color',cols(1+k,:),'linewidth',2);
    %set(h,'linestyle',styles{k+1}(1),'marker',styles{k+1}(2))
    set(h,'linestyle',styles{k+1})
    hold on;
    if k==j, set(h,'linewidth',4); MotsToCompare(:,k,j)', end
    if k==4, set(h,'marker','d','markersize',2); end
  end
  axis([0 0.5 1 1e5]);
  set(gca,'ytick',[1e0 1e1 1e2 1e3 1e4 1e5]);
  if j < length(inds_compare), set(gca,'xticklabel',[]); else, xlabel('degree of connectivity p'); end
  ylabel(['# Motifs ' num2str(j+6)]);
  box off;  
  
  SPb(j) = subplot('position',[0.385 1-0.21*j 0.16 0.18]);
  for k=1:length(inds_compare)
      %h=barwithq(k,mean(Motsmp(:,inds_compare(j),1,inds_compare(k))),mean(Motsmin(1,inds_compare(j),:,inds_compare(k))),mean(Motsmax(1,inds_compare(j),:,inds_compare(k))),0.75,0,'facecolor',[1 1 1]);
      h=barwithsd2(k,mean(Motsmp(:,inds_compare(j),1,inds_compare(k))),std(Motsmp(:,inds_compare(j),1,inds_compare(k))),0.75,0,'facecolor',[1 1 1]);
      if k==j, set(h,'linewidth',2); end
      hold on;
  end
  %barwithq(5,mean(RandMotsmp(:,inds_compare(j),1,inds_compare(k))),mean(RandMotsmin(1,inds_compare(j),:,inds_compare(k))),mean(RandMotsmax(1,inds_compare(j),:,inds_compare(k))),0.75,0,'facecolor',[0.5 0.5 0.5]);
  barwithsd2(5,mean(RandMotsmp(:,inds_compare(j),1,inds_compare(k))),std(RandMotsmp(:,inds_compare(j),1,inds_compare(k))),0.75,0,'facecolor',[0.5 0.5 0.5]);
  %barwithsd(5,mean(mean(RandMotsm(1,inds_compare(j),:,:),4)),mean(std(RandMotsm(1,inds_compare(j),:,:),[],4)),0.75,'facecolor',[0.5 0.5 0.5]);
  axis([0.45 5.55 0 35000]);
  set(gca,'xtick',[]);
  ylabel(['Average # Motifs ' num2str(j+6)]);
  %set(gca,'ytick',[1e4 2e4 3e4],'yticklabel',{'10k','20k','30k'});
  box off;
  drawonemotif(6+j,[0.43 1.132-0.21*j 0.075 0.075]);
  drawonemotif(6+j,[0.09 1.132-0.21*j 0.075 0.075]);
  drawonemotif(6+j,[0.387+(j-1)*0.031 0.122 0.031 0.032],0);
  if j==length(inds_compare), h=annotation('textbox',[0.382+j*0.031 0.122 0.031 0.033],'string','RN'); set(h,'linestyle','none'); end
end
end

l = axes('position',[0.09 0.01 0.24 0.1]);
ks = [2:5 1]
for ik=1:length(ks)
  k = ks(ik);
  hold on;
  xs = 2*ceil(ik/2);
  ys = 2-mod(ik,2);
  if ik==5, xs = 3; ys = 0; end
  h=plot(xs+[0.2 0.6 1],ys+[0 0 0],'color',cols(k,:),'linestyle',styles{k},'linewidth',2);
  if ik==4, set(h,'marker','d','markersize',2); end
  axis([1 6 -0.5 2.5]);
  if ik<5, 
    [vain,h2]=drawmymotif(5+k,[xs+1.5 ys],0.25); 
    set(h2,'color',[0 0 0]);
  else
    text(xs+1.2,ys,'RN');
  end
end

poss = [0.67 0.81 0.07 0.09;0.67 0.61 0.07 0.09;...
        0.62 0.43 0.16 0.09;0.62 0.31 0.16 0.09;0.62 0.19 0.16 0.09;0.62 0.07 0.16 0.09;...
        nan nan nan nan;nan nan nan nan;nan nan nan nan;nan nan nan nan;...
        0.83 0.87 0.16 0.09;0.83 0.75 0.16 0.09;0.83 0.63 0.16 0.09;0.83 0.51 0.16 0.09;...
        0.88 0.3 0.07 0.09;0.88 0.07 0.07 0.09];
%ylims = [1e5 1e5 1e5 1e5 1e5 1e5 nan nan nan nan 4e4 4e4 2e4 2e4 3e4 3e4];
ylims = [8e4 8e4 12e4 12e4 12e4 12e4 nan nan nan nan 6e4 6e4 6e4 6e4 2e4 2e4];
barcols = {{[1 1 1]},{[1 1 1]},{[1 1 1;0 0 0],[1 1 1;0 0 0;1 1 1;0 0 0;0 0 0;0 0 0],[0 0 0;1 1 1],[0 0 0; 0.5 0.5 0.5]},{},{[1 1 1;0 0 0],[1 1 1;0.75 0.75 0.75],[0.5 0.5 0.5;1 1 1],[0 0 0; 0.5 0.5 0.5]},{[1 1 1]},{[1 1 1]}};
barcldivs = {10,10,[10,4,20,10],[],[10,4,20,10],10,10};
pls = [0 0.01;0 0.01;0.04 -0.05;nan nan;0.04 -0.05;0 0.01;0 0.01];
for myind = [0:2 4 5 6]
if myind==0
  %cols = [.5 .5 .5;0.8 0 0; 0.4 0.5 0; 0.4 0 0.5; 0 0 0.8; .5 .5 .5; 0 0 0];
  cols = [.5 .5 .5;cols16(1,:); zeros(3,3); .5 .5 .5; 0 0 0];
elseif myind==1
  cols = [.5 .5 .5;cols16(2,:); zeros(3,3); .5 .5 .5; 0 0 0];
elseif myind<=4
  cols = [.5 .5 .5;cols16(3+(myind-2):3:14,:); .5 .5 .5; 0 0 0];
elseif myind==5
  cols = [.5 .5 .5;cols16(15,:); zeros(3,3); .5 .5 .5; 0 0 0];
elseif myind==6
  cols = [.5 .5 .5;cols16(16,:); zeros(3,3); .5 .5 .5; 0 0 0];
end
inds_compare = find(sum(uniquepatterns,2)==myind);
for j=1:length(inds_compare)
  SPc{myind+1}(j) = subplot('position',poss(inds_compare(j),:));
  for k=1:length(inds_compare)
      %h=barwithq(k,mean(Motsmp(:,inds_compare(j),1,inds_compare(k))),mean(Motsmin(1,inds_compare(j),:,inds_compare(k))),mean(Motsmax(1,inds_compare(j),:,inds_compare(k))),0.75,0,'facecolor',[1 1 1]);
      h=barwithsd2(k,mean(Motsmp(:,inds_compare(j),1,inds_compare(k))),std(Motsmp(:,inds_compare(j),1,inds_compare(k))),0.75,0,'facecolor',[1 1 1]);
      if k==j, set(h,'linewidth',2); end
      hold on;
  end
  %barwithq(length(inds_compare)+1,mean(RandMotsmp(:,inds_compare(j),1,inds_compare(k))),mean(RandMotsmin(1,inds_compare(j),:,inds_compare(k))),mean(RandMotsmax(1,inds_compare(j),:,inds_compare(k))),0.75,0,'facecolor',[0.5 0.5 0.5]);
  barwithsd2(length(inds_compare)+1,mean(RandMotsmp(:,inds_compare(j),1,inds_compare(k))),std(RandMotsmp(:,inds_compare(j),1,inds_compare(k))),0.75,0,'facecolor',[0.5 0.5 0.5]);
  axis([0.45 length(inds_compare)+1.55 0 ylims(inds_compare(j))]);
  set(gca,'xtick',[],'fontsize',8);
  ylabel(['Av. # Motifs ' num2str(inds_compare(j))]);
  box off;
  drawonemotif(inds_compare(j),[poss(inds_compare(j),1)+pls(myind+1,1) poss(inds_compare(j),2)+poss(inds_compare(j),4)+pls(myind+1,2) 0.075 0.075]);
  drawonemotif(inds_compare(j),[poss(inds_compare(j),1)+(j-1)*0.031 poss(inds_compare(end),2)-0.032 0.031 0.032],0); %don't take closer (upwards), would be overwritten
  %drawonemotif(6+j,[0.387+(j-1)*0.031 0.122 0.031 0.032],0);
  %set(gca,'ytick', ylims(inds_compare(j))*[0.5 1],'yticklabel',num2str(ylims(inds_compare(j))*[0.5 1]'/1000));
end
h=annotation('textbox',[poss(inds_compare(j),1)+j*0.031 poss(inds_compare(end),2)-0.032 0.031 0.033],'string','RN'); set(h,'linestyle','none');
end

%annotation(gcf,'line',[0.335 0.335],[0.98 0.02]);
annotation(gcf,'line',[0.557 0.557],[0.98 0.02]);
annotation(gcf,'line',[0.785 0.785],[0.98 0.02]);
annotation(gcf,'line',[0.557 0.785],[0.78 0.78]);
annotation(gcf,'line',[0.557 0.785],[0.58 0.58]);
annotation(gcf,'line',[0.785 0.98],[0.48 0.48]);
annotation(gcf,'line',[0.785 0.98],[0.25 0.25]);
set(gcf,'position',[680   127   734   851])
h=annotation('textbox',[0.06 0.01 0.1 0.1],'string','A'); set(h,'linestyle','none','fontsize',40);
h=annotation('textbox',[0.42 0.01 0.1 0.1],'string','B'); set(h,'linestyle','none','fontsize',40);
h=annotation('textbox',[0.55 0.89 0.1 0.1],'string','C'); set(h,'linestyle','none','fontsize',40);
h=annotation('textbox',[0.55 0.69 0.1 0.1],'string','D'); set(h,'linestyle','none','fontsize',40);
h=annotation('textbox',[0.55 0.49 0.1 0.1],'string','E'); set(h,'linestyle','none','fontsize',40);
h=annotation('textbox',[0.93 0.89 0.1 0.1],'string','F'); set(h,'linestyle','none','fontsize',40);
h=annotation('textbox',[0.81 0.38 0.1 0.1],'string','G'); set(h,'linestyle','none','fontsize',40);
h=annotation('textbox',[0.81 0.15 0.1 0.1],'string','H'); set(h,'linestyle','none','fontsize',40);

print('-depsc', 'motres.eps');
end

function a=drawonemotif(n,pos,on)
if nargin<3, on = 1; end
a=axes('position',pos);
[vain,h2]=drawmymotif(n);
set(h2,'color',[0 0 0]);
axis([-1 1.35 -1.2 1.2]);
if on
  plot([-1 1.35 1.35 -1 -1],1.2*[-1 -1 1 1 -1],'k-');
end
axis equal;
axis off;

end


function [h,nod]=drawmymotif(n,xy,c,ms,lw)
if nargin<2, xy=[0 0]; end
if nargin<3, c=1; end
if nargin<4, ms = 6; end
if nargin<5, lw = 1; end
uniquepatterns = [0     0     0     0     0     0;...
     0     0     0     0     0     1;...
     0     0     0     0     1     1;...
     0     0     0     1     0     1;...
     0     0     0     1     1     0;...
     0     0     1     0     1     0;...
     0     0     0     1     1     1;...
     0     0     1     0     1     1;...
     0     1     0     1     0     1;...
     0     1     1     0     0     1;...
     0     0     1     1     1     1;...
     0     1     0     1     1     1;...
     0     1     1     0     1     1;...
     0     1     1     1     1     0;...
     0     1     1     1     1     1;...
     1     1     1     1     1     1];
 
 for j=n
     M=zeros(3,3);
     M(~eye(3)) = uniquepatterns(j,:);
     [h,nod]=drawmygraph(M,[xy(1)+c*cos(2*pi/3:2*pi/3:2*pi)' xy(2)+c*sin(2*pi/3:2*pi/3:2*pi)'],[],0.75,lw,0.1);
 end
     
end


%function [h,nod]=drawgraph2(M,pos,acoeff,prc,lw,space)
%see also givegraphpos.m for optimal ordering of nodes
function [h,nod]=drawmygraph(M,pos,acoeff,prc,lw,space,plotdegrees)

N=size(M,1);

if nargin < 2 || isempty(pos)
    pos = giveposring(N,1); %N x 2 matrix of coordinates
end

if nargin < 3 %arrow shape
    acoeff = [];
end
if nargin < 4 %arrow shape
    prc = [];
end
if nargin < 5 %line width
    lw = [];
end
if nargin < 6 || isempty(space) %space between line and end points
    space = 0;
end
if nargin < 7 || isempty(plotdegrees)
    plotdegrees = 0;
end
Mdiag = diag(M);
M = M-diag(Mdiag);
Muni = ~~M;

nod=plot(pos(:,1),pos(:,2),'b.');
hold on;
if ~any(any(M))
    if any(Mdiag)
        drawselfarrow(pos(find(Mdiag),1),pos(find(Mdiag),2),min(max(pos(:,1))-min(pos(:,1)),max(pos(:,2))-min(pos(:,2)))/10,Mdiag(find(Mdiag)),lw);
        hold on;
    end
    h = nan;
    return
end

edges_x = zeros(sum(sum(Muni)),2);
edges_y = zeros(sum(sum(Muni)),2);
widths = zeros(sum(sum(Muni)),2);
nplaced = 0;
for i=1:N
    ind = find(M(i,:));
    n = length(ind);
    
    for j=1:n
        edges_x(nplaced+1,:) = [pos(i,1) pos(ind(j),1)];
        edges_y(nplaced+1,:) = [pos(i,2) pos(ind(j),2)];
        widths(nplaced+1,:) = M(i,ind(j));
        nplaced = nplaced + 1;
    end
end
if length(unique(M))==2
    [h,lens]=drawarrowpartial(edges_x,edges_y,acoeff,prc,lw,space);
else
    [h,lens]=drawarrowpartial(edges_x,edges_y,acoeff,prc,widths,space);
end
if any(Mdiag)
%    drawselfarrow(pos(find(Mdiag),1),pos(find(Mdiag),2),min(max(pos(:,1))-min(pos(:,1)),max(pos(:,2))-min(pos(:,2)))/10);
    drawselfarrow(pos(find(Mdiag),1),pos(find(Mdiag),2),min(lens)/5,Mdiag(find(Mdiag)),lw);
    hold on;
end
if plotdegrees
    for i=1:N
        text(pos(i,1),pos(i,2),[num2str(sum(M(:,i))) ',' num2str(sum(M(i,:)))],'color',[0.8 0 0]);
    end
end

axis equal;
axis off;
end    

function h=drawselfarrow(x,y,r,widths,lw)
N = length(x);
xall = [x(:)*ones(1,36)-r + r*repmat([cos((5:6:197)*2*pi/200) 0.7 cos(197*2*pi/200) 1.13],N,1) nan(N,1)];
yall = [y(:)*ones(1,36)   + r*repmat([sin((5:6:197)*2*pi/200) -0.28 sin(197*2*pi/200) -0.45],N,1) nan(N,1)];
if length(unique(widths))==1 && isempty(lw)
  plot(reshape(xall',N*37,1),reshape(yall',N*37,1),'k-');
elseif length(unique(widths))==1
  plot(reshape(xall',N*37,1),reshape(yall',N*37,1),'k-','linewidth',lw);
else
  for i=1:N
    plot(xall(i,:),yall(i,:),'k-','linewidth',widths(i));
    hold on;
  end
end 

end


%function h=drawarrow(x,y)
function [h,lens]=drawarrowpartial(x,y,acoeff,prc,lw,space,colors,linestyle,usefixedmink)

if nargin <3 || isempty(acoeff)
    acoeff = 1;
end
if nargin <4 || isempty(prc)
    prc = 0.9;
end
if nargin <5 || isempty(lw)
    lw = 1;
end
if nargin <6 ||isempty(space)
    space = 0; %maximum (0.5) would make the arrow length zero!
end

if size(x,2) ~= 2
    x = x';
end
if size(y,2) ~= 2
    y = y';
end

if nargin <7 ||isempty(colors)
    colors = zeros(size(x,1),3);
end
if nargin <8 ||isempty(linestyle)
    linestyle = '-'; %maximum (0.5) would make the arrow length zero!
end
if nargin <9 ||isempty(usefixedmink)
    usefixedmink = 0; %maximum (0.5) would make the arrow length zero!
end
d = [x(:,2)-x(:,1), y(:,2)-y(:,1)];
x = x+[space*d(:,1) -space*d(:,1)];
y = y+[space*d(:,2) -space*d(:,2)];
d = [x(:,2)-x(:,1), y(:,2)-y(:,1)];

k = sqrt(sum(d.^2,2));

d = d./(k*ones(1,2));

dperp = acoeff*[-d(:,2), d(:,1)];
if usefixedmink > 0
  lens = k-usefixedmink*(1-prc);
  perplen = 0.5*usefixedmink*(1-prc);
else
  lens = k-min(k)*(1-prc);
  perplen = 0.5*min(k)*(1-prc);
end

px = [x(:,1),x(:,2),x(:,1)+lens.*d(:,1)+perplen*dperp(:,1),x(:,1)+lens.*d(:,1)-perplen*dperp(:,1)];
py = [y(:,1),y(:,2),y(:,1)+lens.*d(:,2)+perplen*dperp(:,2),y(:,1)+lens.*d(:,2)-perplen*dperp(:,2)];

px = reshape([px(:,[1 2 3 2 4]) nan(size(x,1),1)]',size(x,1)*6,1);
py = reshape([py(:,[1 2 3 2 4]) nan(size(x,1),1)]',size(x,1)*6,1);
if length(lw)==1
  h=plot(px,py,'k-','linewidth',lw); set(h,'linestyle',linestyle);
else
  h = zeros(size(x,1),1);
  hold on;
  for i=1:size(x,1)
    h(i)=plot(px((1:5)+6*(i-1)),py((1:5)+6*(i-1)),'k-','linewidth',lw(i),'color',colors(i,:)); set(h(i),'linestyle',linestyle);
  end
end

end



function [h,h2]=barwithsd2(x,y,sd,width,tickwidth,varargin)

if nargin == 1
    sd = zeros(size(x));
    y = x;
    x = 1:length(u);
end
if nargin == 2
    sd = y;
    y = x;
    x = 1:length(y);
end
if isempty(x)
    x = 1:length(y);
end
if isempty(sd)
    sd = zeros(size(y));
end
if nargin < 5 || isempty(tickwidth)
    tickwidth=1/4;
end

h = bar(x,y,varargin{:});
dists = tril(abs(x(:)*ones(1,length(x))-ones(length(x),1)*x(:)'),-1);
mbw = width*tickwidth;
if nargin > 3 && ~isempty(width)
    set(h,'barwidth',width);
end

hold on;
n = length(x);
h2=plot(kron(x(:),ones(9,1))+repmat([0;0;-mbw;mbw;0;0;-mbw;mbw;nan],n,1),kron(y(:),ones(9,1))+reshape([zeros(1,n); repmat(sd(:)',4,1); repmat(-sd(:)',3,1);nan(1,n)],9*n,1),'k-');


end
