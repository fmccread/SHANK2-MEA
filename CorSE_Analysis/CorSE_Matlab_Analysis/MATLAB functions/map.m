
source_1=source;
target_1=target;
for i=1:length(source)
    temp=[source(i) target_1(i)];
    in1=find(target==temp(1,1));
    in2=find(source==temp(1,2));
    if isempty(intersect(in1,in2))~= 1
        source(i)=0;
        target(i)=0;
        weight(i)=0;
    end
end
source=(nonzeros(source))';
target=(nonzeros(target))';
weight=(nonzeros(weight))';
%%
% clear all
% clc
% close all
% cm = [0 1 1 0 0;1 0 0 1 1;1 0 0 0 0;0 0 0 0 1;1 0 1 0 0];
%         ids = {'M30931','L07625','K03454','M27323','M15390'};
%         bg2 = biograph(cm,ids)
%         get(bg2.nodes,'ID')
%  view(bg2)
clear channels
clear ids
j=1;
channels(j)=source(1);
for i=2:length(source)      %finds all source channels and add them to channels
    if source(i)~=channels
        j=j+1;
        channels(j)=source(i);
    end
end

for i=1:length(target)
    if target(i)~=channels
        j=j+1;
        channels(j)=target(i); %finds also all target channels and add them
    end                       %to channels if they are not added as source before
end

N=length(channels);
cm=zeros(N,N);

% for i=1:N
% formatSpec = 'channels: %2.0f and %2.0f \n';
% ids(i)= ('channel ',channels(i));
% end

for i=1:length(channels)
    [row,col]=find(source==channels(i)) %make connections from source to target
    if size(row) ~= 0
        x=target(row,col);
        x=x(1,:);
        for j=1: length(row)
            [row1,col1]=find(channels==x(j))
            cm(i,col1)=1;
        end
    end
end




for k=1:N
    ChannelName{k} = DataCell{channels(k),1}(end-2:end);
    ids{k}=['''CH' ChannelName{k} ''''];
end
% ids={'CH3','CH4','CH5','CH6','CH7','CH8','CH9','CH10','CH11','CH12','CH13','CH14','CH15','CH16'};

bg2=biograph(cm,ids,'NodeAutoSize','off','ShowArrows','off')

    get(bg2.nodes,'ID')
    set(bg2, 'EdgeType', 'curved')
    hbg= view(bg2);
% end

% Position of nodes!!!
% bg2.nodes(1).Position = [424 96];
% dolayout(bg2, 'Pathsonly', true)
% view(bg2)

%% Page size and node positioning EXAMPLE!!
%//////////////////////////////////////////////////////////////////////////
% page_Size = max(cell2mat(arrayfun(@(x) get(x,'Position'),...
%     get(hbg,'Nodes'),'Uniform',false)))
% for i = 1:28
%     set(hbg.Nodes(i),'Position',...
%         [center(1)+radius.*sin((i*2*pi/28)),...
%         center(2)+radius.*cos((i*2*pi/28))])
% end
%//////////////////////////////////////////////////////////////////////////
%% DO LAYOUT ACCORDING TO MEA SETUP
page_Size = [400 400];
% page_Size = max(cell2mat(arrayfun(@(x) get(x,'Position'),...
%     get(hbg,'Nodes'),'Uniform',false)));
XpageSize = page_Size(1);
YpageSize = page_Size(2);

for i = 1:length(ids)
nodeID = str2num(ids{i}(end-2:end-1));
tens = floor(nodeID / 10);
ones = mod(nodeID , 10);
xposition = floor(XpageSize / 8) * tens ;
yposition = floor(YpageSize / 8) * (9-ones); 
set(hbg.Nodes(i),'Position', [xposition, yposition]);
set(hbg.Nodes(i),'Color', [51/255 204/255 204/255]);
set(hbg.Nodes(i),'LineColor', [51/255 204/255 204/255]);
set(hbg.Nodes(i),'FontSize', 12);
set(hbg.Nodes(i),'Shape', 'circle');
set(hbg.Nodes(i),'Size', [30 30]);
end
%% SET PATH WIDTHS ACCORDING TO CORRELATION COEFFICIENTS i.e WEIGHTS

for i = 1:length(source)
set(hbg.Edges(i),'LineWidth', (10*(abs(weight(i))-treshold)+0.1 )); 
end
%%
set(hbg.Edges,'LineColor', [0/255 0/255 0/255]);
dolayout(hbg, 'Pathsonly', true)

