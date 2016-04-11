%function[Trrc]=rrc(Fsampling,Fsymbol)
Fsymbol=1e6;
T=1/Fsymbol;
Fsampling=8e6;
beta=0.5;
RRCTaps=101;
fmax=Fsampling/2;
f=linspace(-fmax,fmax,RRCTaps);

pass=linspace(T,T,RRCTaps);
trans=T/2*(1+cos(pi*T/beta*(abs(f)-(1-beta)/(2*T))));
zero=linspace(0,0,RRCTaps);

part1=abs(f)<((1-beta)/(2*T));
part2=(abs(f)>=(1-beta)/(2*T))&(abs(f)<(1+beta)/(2*T));
part3=abs(f)>=(1+beta)/(2*T);

Hrc=pass.*part1+trans.*part2+zero.*part3;

Hrrc=sqrt(Hrc);
Trrc=ifft(ifftshift(Hrrc));

figure
plot(Hrc);


     

    