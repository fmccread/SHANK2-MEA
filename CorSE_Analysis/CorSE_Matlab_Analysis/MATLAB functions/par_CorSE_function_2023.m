function cross = par_CorSE_function_2023(DataCell, fs, WS, overlap)

TotalDataSample = length(DataCell{1,7});%length data in samples
SElength = round(2*(TotalDataSample / WS) - 1);
chNO = height(DataCell);%number of channels to process
    
% 60 Hz notch filter to remove power line noise
%d = designfilt('bandstopiir','FilterOrder',2, ...
%                   'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
%                   'DesignMethod','butter','SampleRate',fs);

%60 Hz filter
[b,a]=rico(0.01,6,60,1/fs); 
               
%set up matrix
A=zeros(chNO, SElength-1);

%calculate SE in each window
parfor c = 1:chNO  
    x = DataCell{c,7};
    s1 = CorSE_body(x,fs,WS,overlap,SElength,b,a);
    A(c,:)=s1;    
end

%%CROSS-CORRELATION
cross=zeros(chNO,chNO);
delay=zeros(chNO,chNO);
CorrLength = (2*size(A,2))-1;
meancorr=zeros(chNO,CorrLength);


% ALL THE CHANNELS
for channel=1:chNO
    Correlations=zeros(59,CorrLength);
    
    y1=A(channel,:);
    j=0;
    valmax1=zeros(1,chNO-channel);
    lag1=zeros(1,chNO-channel);
    valmax2=zeros(1,chNO-channel);
    lag2=zeros(1,chNO-channel);
   
    for i=(channel+1):chNO
        j=j+1;
        y2=A(i,:);

        
        [xcf,lags,bounds]=crosscorr(y1,y2,(size(A,2)-1));
        
        Correlations(i,1:CorrLength)=xcf; %calculating cross-correlations' mean for each channel

        %
        maxDelayforCorr = 0;% max delay for max correlation in seconds
        zeroLagindex = find(lags==0);
        winSizeforDelay = (maxDelayforCorr/(WS/fs));
        
        valmax1(1,j)=max(xcf(zeroLagindex-winSizeforDelay:...
            zeroLagindex+winSizeforDelay));
        index=find(xcf==valmax1(j));
        lag1(1,j)=lags(index);
        
        [xcf,lags,bounds]=crosscorr(y2,y1,(size(A,2)-1));


        maxDelayforCorr = 0;% max delay for max correlation in seconds
        zeroLagindex = find(lags==0);
        winSizeforDelay = (maxDelayforCorr/(WS/fs));
        
        valmax2(1,j)=max(xcf(zeroLagindex-winSizeforDelay:...
            zeroLagindex+winSizeforDelay));
        index=find(xcf==valmax2(j));
        lag2(1,j)=lags(index);

    end
    cross((channel+1):chNO,channel)=valmax1;
    delay((channel+1):chNO,channel)=lag1;
    cross(channel,(channel+1):chNO)=valmax2;
    delay(channel,(channel+1):chNO)=lag2;
    meancorr(channel,1:CorrLength)=mean(Correlations);
end
    
