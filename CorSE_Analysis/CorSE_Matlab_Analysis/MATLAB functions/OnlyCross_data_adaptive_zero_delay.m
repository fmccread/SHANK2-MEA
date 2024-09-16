%%CROSS-CORRELATION
cross=zeros(chNO,chNO);
delay=zeros(chNO,chNO);
CorrLength = (2*size(A,2))-1;
meancorr=zeros(chNO,CorrLength);


% ALL THE CHANNELS
for channel=1:chNO
    Correlations=zeros(59,CorrLength);
    
    y1=A(channel,:);
    %         count=0;
    j=0;
    valmax1=zeros(1,chNO-channel);
    lag1=zeros(1,chNO-channel);
    valmax2=zeros(1,chNO-channel);
    lag2=zeros(1,chNO-channel);
    for i=(channel+1):chNO
        j=j+1;
        %             count=count+1;
        y2=A(i,:);
        %             A(count,:)=y2;
        %             [xcf,lags,bounds] = crosscorr(y1,y2);
        
        
        [xcf,lags,bounds]=crosscorr(y1,y2,(size(A,2)-1));
        
        Correlations(i,1:CorrLength)=xcf;%calculating cross-correlations' mean for each channel
        
        %             pause
        %             figure(1)
        %
        %             plot(lags,xcf);
        %             xlabel('time delay in number of samples');
        %             ylabel('cross correlation');
        %             title(['cross-correlation channels ',int2str(channel),' --> ',int2str(i)]);
        %
         maxDelayforCorr = 0;% max delay for max correlation in seconds
        zeroLagindex = find(lags==0);
        winSizeforDelay = (maxDelayforCorr/(WS/fs));
        
        valmax1(1,j)=max(xcf(zeroLagindex-winSizeforDelay:...
            zeroLagindex+winSizeforDelay));
        index=find(xcf==valmax1(j));
        lag1(1,j)=lags(index);
        
        [xcf,lags,bounds]=crosscorr(y2,y1,(size(A,2)-1));
        %             pause
        %             figure(2)
        %             plot(lags,xcf);
        %             xlabel('time delay in number of samples');
        %             ylabel('cross correlation');
        %             title(['cross-correlation channel ',int2str(i),' --> ',int2str(channel)]);
        %
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
    %         subplot(8,8,channel)
    %         plot(meancorr(channel,1:CorrLength))
    %         set(gca,'XTick',0:(CorrLength/5):CorrLength)
    %         set(gca,'XTickLabel',-(length(A)-1):(CorrLength/5):(size(A,2)-1))
end