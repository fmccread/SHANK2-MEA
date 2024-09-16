clear row
clear col
clear source
clear target
%SELECT TRESHOLD VALUE with fixed CorSE value

% treshold=0.5; %for fixed threshold
% [row,col]=find(cross>treshold); %for fixed threshold

%% Threshold for maximum number of connections decided by user

val=sort(cross(:),'descend');
sortedval=sort(val,'descend');
sizeOfCross = size(cross,1)*size(cross,2);
sortedval=sortedval(1:2:sizeOfCross,1); 
maximum_N = 40; %write double the number of maximum connections, e.g., for maximum 20 channel write 40 
treshold= sortedval(maximum_N);
unq_sorted = -1*unique(-sortedval);
row=[];
col=[];
for r = 1 : maximum_N;
    [row_temp,col_temp]=find( abs(cross-unq_sorted(r)) < 1e-5);%if cross==unq_sorted(r)
    if length(row) + length(row_temp) > maximum_N
        row_temp = row_temp(1:maximum_N-length(row));
        col_temp = col_temp(1:maximum_N-length(col));
        row=[row; row_temp];
        col=[col; col_temp];
    else
        row=[row; row_temp];
        col=[col; col_temp];
    end
end

%% Calculating sources, targets and connection strengths(CorSEs)
j=0;

for i=1:length(row)
    
    

    j=j+1;

    formatSpec = 'channels: %2.0f and %2.0f \n';
    fprintf(formatSpec,col(i),row(i))
    source(j)=col(i);
    target(j)=row(i);
    weight(j)=cross(row(i),col(i));

    
end
N=length(source);

map;