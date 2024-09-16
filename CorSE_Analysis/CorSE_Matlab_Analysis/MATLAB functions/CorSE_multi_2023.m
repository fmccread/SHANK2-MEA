%% ENTROPY ALGORITHM 
% This algorithm calculates the correlated spectral entropy (CorSE) based functional connectivity between
% two time series. For more information, <http://journal.frontiersin.org/article/10.3389/fncom.2016.00112/full Kapucu et al.(2016)>. Please cite the paper if you have used or inspired by the method.

%% Introduction
% For this version, algorithm first calculates time variant spectral entropies (SE)
% from all the channels. Multichannel data, e.g., MEA data, has formatted
% as DataCell.

%% Parameters
fs = 12500; %sampling freq.
WS = fs*1; %window size to calculate Spectral Entropy in samples fs for 0.5s
overlap=0.5;%other than 0.5 makes it more complex and code needs to be changed
TotalDataSample = length(DataCell{1,7});%length data in samples
SElength = round(2*(TotalDataSample / WS) - 1);
chNO = length(DataCell);%number of channels to process

A=zeros(chNO, SElength-1);

for channel=1:chNO
    
    x=DataCell{channel,7}(1:end);

    %% Window Function
    
    for i=1:SElength-1 
        
        y=x((i-1)*WS*overlap +1:(i-1)*WS*overlap+WS);
        
        y=y-mean(y);%subtract mean from each window "preferable"



        %% Filtering

        [b,a]=rico(0.01,6,60,1/fs);  % rico(z,B,fc,T), con: z,B= dati; fc= 60 Hz; T=1/fs.
         
        z=filtfilt(b,a,y);         
        
        % High pass filter ~7Hz
        n=3;
        wn=7/(fs/2); %7HzHPF
        [b,a]=butter(n,wn,'high');
    
        
        z1=filtfilt(b,a,z);


        
        %% POWER SPECTRUM OF FILTERED SIGNAL
        
        nfft=2^15;
        w=hann(WS);
        o=0;
        [Pxx,f]=pwelch(z1,w,o,nfft,fs); %%Pxx lungo 16385
        
        pxx_norm1=Pxx/sum(Pxx);%normalizzo sulla somma
        
        
        
        %% SPECTRAL ENTROPY 
        
        k=0;
        s1(i)=0;
        n=length(pxx_norm1);
        for j=1:n
            if isfinite(log10(1/pxx_norm1(j)))==1 %%isfinite e' 1 se il numero e' finito
                s1(i)=s1(i)+ pxx_norm1(j)*log10(1/pxx_norm1(j));
            else
                k=1;
                s1(i)=s1(i);
            end
        end
        s1(i)=s1(i)/log10(n);
        
    end
    
    A(channel,:)=s1;
    
    
end

%% Plotting time variant SEs
% Plotting is done according to 8x8 MEA lay out. Comment out if plot is not needed
 asse_x=1:1:SElength-1;%temporary change added "-1"
 figure(1)
 for channel=1:chNO
     %CHname = str2num(DataCell{channel,1}(end-2:end));
     %ones = mod(CHname,10);
     %tens = floor(CHname/10);
     %CHlocation = ((ones - 1) * 8) + tens ;
     subplot(8,8,channel)
     plot(A(channel,:), 'b');
     %ylim([0.7 1]);
     %xlim([1 599]);
     axis([1 599 0.85 1]);
     %title(['Entropy channel ',int2str(CHname)]);
     h=get(gcf,'CurrentAxes');
     set(h,'FontSize',7)
 end
