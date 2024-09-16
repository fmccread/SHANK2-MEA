function s1 = CorSE_body(x, fs, WS, overlap, SElength, b,a)

for i=1:SElength-1
    
    y=x((i-1)*WS*overlap +1:(i-1)*WS*overlap+WS);
        
    y=y-mean(y);
    
    
    %% Filtering
         
    %z=filtfilt(d,y);         
    z=filtfilt(b,a,y);
        
    % 7 Hz high-pass filter  
    n=3;
    wn=7/(fs/2); 
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



