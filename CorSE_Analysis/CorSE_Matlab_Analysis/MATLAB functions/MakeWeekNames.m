function wknames = MakeWeekNames(Wstart,Wint,Wend)

wkseq = Wstart:Wint:Wend;
wknames = cell(1,length(wkseq));

for s = 1:length(wkseq)
   
   wknum = wkseq(s);
   tmp = ['Week ', num2str(wknum)];
   wknames{s} = tmp;
    
end