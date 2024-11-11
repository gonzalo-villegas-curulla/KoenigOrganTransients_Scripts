
signal_init = pipecurr(1:fix(fs*0.3));

Nfft   = 2^10; %2^15;
WinLen = 2*ceil(fs/f1estim);
win    = hann(WinLen, 'periodic');
HA     = WinLen;
OL     = HA-1;

[sp,fp,tp] = stft(signal_init,...
    fs,...
    'Window',win,...    
    'FFTLength',Nfft, ...
    'FrequencyRange',"onesided",...
    OverlapLength=OL);


Fcut = 5e3;
idxf = fp<=Fcut;

energy = sum(abs(sp(idxf,:)).^2,1);

figure(32); clf;

ah(1)=subplot(3,1,1);
plot(signal(find(tp)));
yyaxis right;
plot(energy);

ah(2)=subplot(3,1,[2,3]);
imagesc([], fp, db(abs(sp)));
ax=gca; ax.YDir = 'normal'; %colorbar();
linkaxes(ah, 'x'); axis tight;