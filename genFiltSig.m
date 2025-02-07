%
%
% genFiltSig.m
%
%

fname = '/usr/local/spklib/fek/fekxx012.tim';
speech =sread8(fname);
sound(speech)

[Pspeech,spFreq] = psd(speech,1024,8000,hanning(1024));
Pspeech = 10*log10(Pspeech);

%low Pass filter
[lpb,lpa] = cheby1(10,.25,.2);

%freqz(lpb,lpa,1024,8e3)
lpspeech = filter(lpb,lpa,speech);
lfP = 20*log10(abs(freqz(lpb,lpa,spFreq,8e3)));
lpP = 10*log10( psd(lpspeech,1024,8000,hanning(1024)));

%Band Pass filter
[bpb,bpa] = cheby1(10,.25,[.25 .56]);

%freqz(bpb,bpa,1024,8e3)
bpspeech = filter(bpb,bpa,speech);
bfP = 20*log10(abs(freqz(bpb,bpa,spFreq,8e3)));

bpP = 10*log10( psd(bpspeech,1024,8000,hanning(1024)));



%High Pass filter
[hpb,hpa] = cheby1(10,.25,.6,'high');

%freqz(hpb,hpa,1024,8e3)
hpspeech = filter(hpb,hpa,speech);
hfP = 20*log10(abs(freqz(hpb,hpa,spFreq,8e3)));

hpP = 10*log10( psd(hpspeech,1024,8000,hanning(1024)));

%Notch Filter

[npb,npa] = cheby1(10,.25,[.18  .26],'stop');

%freqz(npb,npa,1024,8e3)
ntspeech = speech+max(abs(speech))*sin(2*pi/8000*900*[1:length(speech)]');
npspeech = filter(npb,npa,ntspeech);
nfP = 20*log10(abs(freqz(npb,npa,spFreq,8e3)));

npP = 10*log10( psd(npspeech,1024,8000,hanning(1024)));

ntPspeech = 10*log10(psd(ntspeech,1024,8000,hanning(1024)));


save FiltSig speech Pspeech spFreq lpP lfP lpspeech hpspeech hpP ...
     bpb bpa bpspeech bpP hpa hpb hpspeech hpP npa npb npspeech npP ...
     nfP hfP bfP ntspeech ntPspeech
