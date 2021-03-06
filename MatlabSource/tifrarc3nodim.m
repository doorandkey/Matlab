%==========================================================================
%
%   << (무차원) TICM을 이용한 사각(일정/변)단면 아크의 면내 자유진동해석 >>
%               회전관성과 전단변형을 고려 가능함
%              ( 매트랩 파일 : tifrarc3nodim.m )    2016. 03. 25
%
%==========================================================================
function tifrarc3nodim
%
clear all
%fid1=fopen('c:\tifrarc3nodim.out','w'); % 파일 출력
%
% <입력> 계산 모델
%        mt=0; 회전관성과 전단변형을 무시, mt=1; 회전관성과 전단변형을 고려
%        vs=0; 일정 단면, vs=1; 변단면(b: 일정, h: 선형 감소)
%
%it=1; mt=0; vs=0; % 일정단면 아크 (고정-고정) Ingenieur-Archiv (1979) 337-346 Fig.2
%it=1; mt=1; vs=0; % 일정단면 티모센코 아크 (고정-고정) Ingenieur-Archiv (1979) 337-346 Fig.2
%it=2; mt=0; vs=0; % 일정단면 아크 (고정-자유) Ingenieur-Archiv (1979) 337-346 Fig.3
%it=2; mt=1; vs=0; % 일정단면 티모센코 아크 (고정-자유) Ingenieur-Archiv (1979) 337-346 Fig.3
%it=2; mt=0; vs=1; % 변단면 아크 (고정-자유) Ingenieur-Archiv (1979) 337-346 Fig.3
%it=2; mt=1; vs=1; % 변단면 티모센코 아크 (고정-자유) Ingenieur-Archiv (1979) 337-346 Fig.3
%it=3; mt=1; vs=0; % 일정단면 티모센코 아크 (고정-고정) Ingenieur-Archiv (1979) 337-346 Fig.8a
%it=4; mt=1; vs=0; % 일정단면 티모센코 아크 (고정-고정) Ingenieur-Archiv (1979) 337-346 Fig.8b
it=5; mt=1; vs=0; % 일정단면 티모센코 아크 (고정-고정)
%
[FRQ,fno,modeYN]=analysis_data(it); % 자유진동해석 파라메타
[itn,ndof,nn,b0,h0,b1,h1_s,h1_d,h1_e,R,emod,roh,nu,kappa,sy0,ALPHA]=arc_data(it); % 아크 모델링
[KK]=support(ndof,nn,it); % 경계조건 및 지지조건
%
timestart=clock;
ii=0;
for alpha_degree=ALPHA(1):ALPHA(2):ALPHA(3)
   alpha_degree % 아크의 전체 각도 [degree]
   alpha=alpha_degree*pi/180; % [rad]
   Theta=linspace(0,alpha,nn); % 아크 각도의 분할
   %ii=ii+1; Alphas(ii)=alpha_degree; % [rad]
   for h1=h1_s:h1_d:h1_e
      %h1
      ii=ii+1; H1_h0s(ii)=h1/h0;
      [Lamda,XXs]=bsc(FRQ,fno,modeYN,itn,nn,ndof,Theta,sy0,KK,mt,emod,kappa,nu,b0,h0,b1,h1,alpha,vs);
      Lamdas(ii,1:fno)=Lamda
      Hz=(sqrt(emod/roh)*Lamda.^2/(R*sy0))/(2*pi); % 고유진동수 [Hz]
   end
end
timeend=clock;
computation_time=etime(timeend,timestart)
%
%plot(Alphas,Lamdas(:,1),'bo',Alphas,Lamdas(:,2),'bo',Alphas,Lamdas(:,3),'bo',Alphas,Lamdas(:,4),'bo')
%plot(H1_h0s,Lamdas(:,1),'bo',H1_h0s,Lamdas(:,2),'bo',H1_h0s,Lamdas(:,3),'bo',H1_h0s,Lamdas(:,4),'bo')
%axis([0.25 1.5 0 10])
%
%fprintf(fid1,'  Alpha    Lamda(1)       Lamda(2)       Lamda(3)       Lamda(4)\n');
%for jj=1:ii   
%   fprintf(fid1,' %3i :  %12.6e  %12.6e  %12.6e  %12.6e\n',Alphas(jj),Lamdas(jj,1),Lamdas(jj,2),Lamdas(jj,3),Lamdas(jj,4));
%end
%fclose(fid1);
%
%-----------------------------------------------------------------------
function [FRQ,fno,modeYN]=analysis_data(it)
%-----------------------------------------------------------------------
%
% 자유진동해석 파라메타
%
switch it
   case {1}
      fmin=1; % 검색 시작 고유치 [lamda]
      fdel=0.1; % 검색 초기 스텝 고유치 [lamda]
      fmax=30; % 검색 종료 고유치 [lamda]
      fno=4; % 검색할 고유치의 최고 차수
      feps=1e-9; % 고유치 검색 허용 오차 [lamda]
      modeYN=0; % 고유모드 계산 여부 ( Yes=1 / No=0 )
   case {2}
      fmin=1; % 검색 시작 고유치 [lamda]
      fdel=0.1; % 검색 초기 스텝 고유치 [lamda]
      fmax=10; % 검색 종료 고유치 [lamda]
      fno=4; % 검색할 고유진동수의 최고 차수
      feps=1e-9; % 고유치 검색 허용 오차 [lamda]
      modeYN=0; % 고유모드 계산 여부 ( Yes=1 / No=0 )
   case {3}
      fmin=1; % 검색 시작 고유치 [lamda]
      fdel=0.1; % 검색 초기 스텝 고유치 [lamda]
      fmax=10; % 검색 종료 고유치 [lamda]
      fno=4; % 검색할 고유진동수의 최고 차수
      feps=1e-7; % 고유치 검색 허용 오차 [lamda]
      modeYN=1; % 고유모드 계산 여부 ( Yes=1 / No=0 )
   case {4}
      fmin=1; % 검색 시작 고유치 [lamda]
      fdel=0.1; % 검색 초기 스텝 고유치 [lamda]
      fmax=10; % 검색 종료 고유치 [lamda]
      fno=4; % 검색할 고유진동수의 최고 차수
      feps=1e-7; % 고유치 검색 허용 오차 [lamda]
      modeYN=1; % 고유모드 계산 여부 ( Yes=1 / No=0 )
   case {5}
      fmin=1; % 검색 시작 고유치 [lamda]
      fdel=0.1; % 검색 초기 스텝 고유치 [lamda]
      fmax=10; % 검색 종료 고유치 [lamda]
      fno=4; % 검색할 고유진동수의 최고 차수
      feps=1e-7; % 고유치 검색 허용 오차 [lamda]
      modeYN=1; % 고유모드 계산 여부 ( Yes=1 / No=0 )      
   otherwise
end
%
FRQ=[fmin fdel fmax feps]; % [lamda]
%
%-----------------------------------------------------------------------
function [itn,ndof,nn,b0,h0,b1,h1_s,h1_d,h1_e,R,emod,roh,nu,kappa,sy0,ALPHA]=arc_data(it)
%-----------------------------------------------------------------------
%
% 아크 모델링
%
ndof=3; % 절점당 자유도수
%
switch it
   case {1}
      itn=50; % 필드행렬 적분 회수 (Runge-Kutta-Gill Method)
      nn=5; % 절점 수 (요소 수 +1)
      ALPHA(1)=180; % 아크의 각도 (계산 시작 각도) [degree]
      ALPHA(2)=5;  % 아크의 각도 (계산 증분 각도) [degree]
      ALPHA(3)=180; % 아크의 각도 (계산 종료 각도) [degree]
      R=1; % 아크의 곡률 반경 [m]
      kappa=0.823;
      sy0=50;
      h0=sqrt(12)*R/sy0; % 단면의 높이 (depth)
      b0=h0; % 단면의 폭 (width)
      %
      emod=2.0e11; % 종탄성계수 [N/m2]
      nu=0.33; % 프와송 비
      roh=7800; % 밀도 [kg/m3]
      % 무의미한 데이타(일정단면이므로)
      b1=b0;
      h1_s=h0;
      h1_d=h0;
      h1_e=h0;
   case {2}
      itn=100; % 필드행렬 적분 회수 (Runge-Kutta-Gill Method)
      nn=4; % 절점 수 (요소 수 +1)
      ALPHA(1)=90; % 아크의 각도 (계산 시작 각도) [degree]
      ALPHA(2)=1;  % 아크의 각도 (계산 증분 각도) [degree]
      ALPHA(3)=90; % 아크의 각도 (계산 종료 각도) [degree]
      R=1; % 아크의 곡률 반경 [m]
      kappa=0.85;
      sy0=50;
      h0=sqrt(12)*R/sy0; % 단면의 높이 (depth)
      b0=h0; % 단면의 폭 (width)
      b1=b0;
      h1_s=0.25*h0;
      h1_d=0.25*h0;
      h1_e=1.5*h0;
      %      
      emod=2.0e11; % 종탄성계수 [N/m2]
      nu=0.3; % 프와송 비
      roh=7800; % 밀도 [kg/m3]
   case {3}
      itn=1; % 필드행렬 적분 회수 (Runge-Kutta-Gill Method)
      nn=41; % 절점 수 (요소 수 +1)
      ALPHA(1)=90; % 아크의 각도 (계산 시작 각도) [degree]
      ALPHA(2)=1;  % 아크의 각도 (계산 증분 각도) [degree]
      ALPHA(3)=90; % 아크의 각도 (계산 종료 각도) [degree]
      R=1; % 아크의 곡률 반경 [m]
      kappa=0.85;
      sy0=100;
      h0=sqrt(12)*R/sy0; % 단면의 높이 (depth)
      b0=h0; % 단면의 폭 (width)
      %      
      emod=2.0e11; % 종탄성계수 [N/m2]
      nu=0.3; % 프와송 비
      roh=7800; % 밀도 [kg/m3]
      % 무의미한 데이타(일정단면이므로)
      b1=b0;
      h1_s=h0;
      h1_d=h0;
      h1_e=h0;
   case {4}
      itn=1; % 필드행렬 적분 회수 (Runge-Kutta-Gill Method)
      nn=41; % 절점 수 (요소 수 +1)
      ALPHA(1)=90; % 아크의 각도 (계산 시작 각도) [degree]
      ALPHA(2)=1;  % 아크의 각도 (계산 증분 각도) [degree]
      ALPHA(3)=90; % 아크의 각도 (계산 종료 각도) [degree]
      R=1; % 아크의 곡률 반경 [m]
      kappa=0.85;
      sy0=30;
      h0=sqrt(12)*R/sy0; % 단면의 높이 (depth)
      b0=h0; % 단면의 폭 (width)
      %      
      emod=2.0e11; % 종탄성계수 [N/m2]
      nu=0.3; % 프와송 비
      roh=7800; % 밀도 [kg/m3]
      % 무의미한 데이타(일정단면이므로)
      b1=b0;
      h1_s=h0;
      h1_d=h0;
      h1_e=h0;
   case {5}
      itn=1; % 필드행렬 적분 회수 (Runge-Kutta-Gill Method)
      nn=81; % 절점 수 (요소 수 +1)
      ALPHA(1)=180; % 아크의 각도 (계산 시작 각도) [degree]
      ALPHA(2)=1;  % 아크의 각도 (계산 증분 각도) [degree]
      ALPHA(3)=180; % 아크의 각도 (계산 종료 각도) [degree]
      R=1; % 아크의 곡률 반경 [m]
      kappa=0.85;
      sy0=50;
      h0=sqrt(12)*R/sy0; % 단면의 높이 (depth)
      b0=h0; % 단면의 폭 (width)
      %      
      emod=2.0e11; % 종탄성계수 [N/m2]
      nu=0.3; % 프와송 비
      roh=7800; % 밀도 [kg/m3]
      % 무의미한 데이타(일정단면이므로)
      b1=b0;
      h1_s=h0;
      h1_d=h0;
      h1_e=h0;
   otherwise
end
%
%-----------------------------------------------------------------------
function [KK]=support(ndof,nn,it)
%-----------------------------------------------------------------------
%
% 무차원 기초지지 스프링 입력 (Local Coordinate로 입력함)
%
KK=zeros(ndof,nn); % 절점 1에서 절점 n(마지막 절점)의 스프링상수
k_inf=1.0e10; % 무한 강성계수
%
switch it
   case {1,3,4} % 고정-고정
      %KK(i,j)=5.0; % 절점 j에서 자유도 i의 (무차원)스프링상수 (i=1: 접선 방향, i=2: 반경 방향, i=3: 회전 방향)
      KK(1,1)=k_inf; KK(2,1)=k_inf; KK(3,1)=k_inf; 
      KK(1,nn)=k_inf; KK(2,nn)=k_inf; KK(3,nn)=k_inf; 
   case {2} % 고정-자유
      %KK(i,j)=5.0; % 절점 j에서 자유도 i의 (무차원)스프링상수 (i=1: 접선 방향, i=2: 반경 방향, i=3: 회전 방향)
      KK(1,1)=k_inf; KK(2,1)=k_inf; KK(3,1)=k_inf; 
   case {5} % 고정-고정 (중앙에 스프링)
      node=(nn-1)/2+1;
      KK(2,node)=1e10; % 절점 j에서 자유도 i의 (무차원)스프링상수 (i=1: 접선 방향, i=2: 반경 방향, i=3: 회전 방향)
      KK(1,1)=k_inf; KK(2,1)=k_inf; KK(3,1)=k_inf; 
      KK(1,nn)=k_inf; KK(2,nn)=k_inf; KK(3,nn)=k_inf; 
   otherwise
end
%
%-----------------------------------------------------------------------
function [Lamda,XXs]=bsc(FRQ,fno,modeYN,itn,nn,ndof,Theta,sy0,KK,mt,emod,kappa,nu,b0,h0,b1,h1,alpha,vs)
%-----------------------------------------------------------------------
%
% 이분법
%
fmin=FRQ(1); fdel=FRQ(2); fmax=FRQ(3); feps=FRQ(4); % [lamda]
bmin=fmin;   lamda=fmin;
[isp]=ticmfree(lamda,itn,nn,ndof,Theta,sy0,KK,mt,emod,kappa,nu,b0,h0,b1,h1,alpha,vs);
%
fno_cal=0;
while ((lamda <= fmax) & (fno_cal < fno))
   lamda=lamda+fdel;
   [is]=ticmfree(lamda,itn,nn,ndof,Theta,sy0,KK,mt,emod,kappa,nu,b0,h0,b1,h1,alpha,vs);
   if is ~= isp
      bmax=lamda;
      while (bmax-bmin)/lamda > feps
         lamda=(bmin+bmax)/2;
         [is]=ticmfree(lamda,itn,nn,ndof,Theta,sy0,KK,mt,emod,kappa,nu,b0,h0,b1,h1,alpha,vs);
         if is == isp
            bmin=lamda;
         else
            bmax=lamda;
         end
      end
      fno_cal=fno_cal+1;
      lamda=(bmin+bmax)/2;
      Lamda(fno_cal)=lamda;
      %
      if modeYN == 0
          XXs=0; % No data
      else
          [XX]=ticmmode(lamda,itn,nn,ndof,Theta,sy0,KK,mt,emod,kappa,nu,b0,h0,b1,h1,alpha,vs); % 고유모드
          XXs(:,:,fno_cal)=XX;
      end
      %
      bmin=bmax;
      isp=-isp;
   else
      bmin=lamda;
   end
end
%
%-----------------------------------------------------------------------
function [is]=ticmfree(lamda,itn,nn,ndof,Theta,sy0,KK,mt,emod,kappa,nu,b0,h0,b1,h1,alpha,vs)
%-----------------------------------------------------------------------
%
% 전달영향계수법을 이용한 자유진동해석
%
is=1;
[KS]=make_pm(ndof,KK(:,1));
[PP,QQ,RR,SS]=make_pqrs(2,ndof,itn,lamda,Theta,sy0,mt,emod,kappa,nu,b0,h0,b1,h1,alpha,vs);
GG=RR+SS*KS; HH=PP+QQ*KS;
[KS]=make_pm(ndof,KK(:,2));
WW=GG+KS*HH;
if det(WW)<0; is=-is; end
TT=HH*inv(WW);
%
for ii=3:nn
   [PP,QQ,RR,SS]=make_pqrs(ii,ndof,itn,lamda,Theta,sy0,mt,emod,kappa,nu,b0,h0,b1,h1,alpha,vs);
   GG=RR*TT+SS; HH=PP*TT+QQ;
   [KS]=make_pm(ndof,KK(:,ii));
   WW=GG+KS*HH;
   if det(WW)<0; is=-is; end
   TT=HH*inv(WW);
end
%
%-----------------------------------------------------------------------
function [PK]=make_pm(ndof,Spring)
%-----------------------------------------------------------------------
%
% 포인트행렬 만들기 (기초지지 스프링)
%
PK=zeros(ndof);
PK(1,2)=Spring(1,1); % 접선 방향
PK(2,1)=Spring(2,1); % 반경 방향
PK(3,3)=Spring(3,1); % 회전 방향
%
%-----------------------------------------------------------------------
function [PP,QQ,RR,SS]=make_pqrs(ii,ndof,itn,lamda,Theta,sy0,mt,emod,kappa,nu,b0,h0,b1,h1,alpha,vs)
%-----------------------------------------------------------------------
%
% 필드행렬(일정단면/변단면 아크)로부터 부분동강성행렬 만들기
%
Y0=eye(ndof*2);
dTheta=Theta(ii)-Theta(ii-1);
th0=Theta(ii-1);
dh=dTheta/itn;
[FM]=multi_runge_kutta_gill(ndof,Y0,dh,itn,lamda,sy0,mt,emod,kappa,nu,b0,h0,b1,h1,alpha,th0,vs);
%
ndof1=ndof+1; ndof2=ndof*2;
PP=FM(1:ndof,1:ndof);      QQ=FM(1:ndof,ndof1:ndof2);
RR=FM(ndof1:ndof2,1:ndof); SS=FM(ndof1:ndof2,ndof1:ndof2);
%
%-----------------------------------------------------------------------
function [Y0]=multi_runge_kutta_gill(ndof,Y0,dh,itn,lamda,sy0,mt,emod,kappa,nu,b0,h0,b1,h1,alpha,th0,vs)
%-----------------------------------------------------------------------
%
% Runge-Kutta-Gill Method (일정단면/변단면 아크)
%
aa=1.0-sqrt(0.5);   bb=1.0+sqrt(0.5);
Q0=zeros(ndof*2);   Y(:,:,1)=Y0;
%
for ii=1:itn
   th=th0+(ii-1)*dh;
   [UM]=make_um(th,lamda,sy0,mt,emod,kappa,nu,b0,h0,b1,h1,alpha,vs);
   K1=dh*UM*Y0;    R1=(K1-2*Q0)/2;
   Y1=Y0+R1;       Q1=Q0+3*R1-K1/2;
   %
   [UM]=make_um(th+dh/2,lamda,sy0,mt,emod,kappa,nu,b0,h0,b1,h1,alpha,vs);
   K2=dh*UM*Y1;   R2=aa*(K2-Q1);
   Y2=Y1+R2;      Q2=Q1+3*R2-aa*K2;
   %   
   [UM]=make_um(th+dh/2,lamda,sy0,mt,emod,kappa,nu,b0,h0,b1,h1,alpha,vs);
   K3=dh*UM*Y2;   R3=bb*(K3-Q2);
   Y3=Y2+R3;      Q3=Q2+3*R3-bb*K3;
   %   
   [UM]=make_um(th+dh,lamda,sy0,mt,emod,kappa,nu,b0,h0,b1,h1,alpha,vs);
   K4=dh*UM*Y3;   R4=(K4-2*Q3)/6;
   Y0=Y3+R4;      Q0=Q3+3*R4-K4/2;
end
%
%-----------------------------------------------------------------------
function [U]=make_um(th,lamda,sy0,mt,emod,kappa,nu,b0,h0,b1,h1,alpha,vs)
%-----------------------------------------------------------------------
%
% 무차원 U 행렬 만들기 (일정단면/변단면 아크)
%
switch vs
   case {0} % 일정단면(b: 일정, h: 일정)
      a=1; % A(theta)/A0
      iy=1; % I2y(theta)/I2y0
   case {1} % 변단면(b: 일정, h: 선형)
      a=(1-(1-b1/b0)*(th/alpha)^0)*(1-(1-h1/h0)*(th/alpha)^1);
      iy=(1-(1-b1/b0)*(th/alpha)^0)*(1-(1-h1/h0)*(th/alpha)^1)^3;
   otherwise
end
%
if mt == 0 % 회전관성과 전단변형을 고려 안함
   k2=iy/(sy0^2*a);   l4s2=lamda^4/sy0^2;
   %
   U=[0        -1        1     0       0     0
      1         0        0     1/a     0     0
      0         0        0     0       0    -1/(a*k2)
      0        -l4s2*a   0     0       1     0
     -l4s2*a    0        0    -1       0     0
      0         0        0     0       1     0]; 
else % 회전관성과 전단변형을 고려함
   k2=iy/(sy0^2*a);   l4s2=lamda^4/sy0^2;
   k_12=2*k2;    k_22=k2;
   gmod=emod/(2*(1+nu)); 
   ekg=emod/(kappa*gmod);
   l4s2k2=l4s2*(1+k2);
   l4s2k_1=l4s2*k_12;
   l4s2k_2=l4s2*k_22;
   %
   U=[0          -1              1              0      ekg/a    0
      1           0              0              1/a    0        0
      0           0              0              0      0       -1/(a*k2)
      0          -l4s2k2*a      -l4s2k_1*a      0      1        0
     -l4s2k2*a    0              0             -1      0        0
      0           l4s2k_1*a      l4s2k_2*a      0      1        0]; 
end
%
%-----------------------------------------------------------------------
function [XX]=ticmmode(lamda,itn,nn,ndof,Theta,sy0,KK,mt,emod,kappa,nu,b0,h0,b1,h1,alpha,vs)
%-----------------------------------------------------------------------
%
% 전달영향계수법을 이용한 고유모드 계산
%
ndofm1=ndof-1;
II=eye(ndof);
TT=zeros(ndof,ndof,nn);
HH=zeros(ndof,ndof,nn-1);
XX=zeros(ndof,nn);
%
% 전달식
%
[KS]=make_pm(ndof,KK(:,1));
[PP,QQ,RR,SS]=make_pqrs(2,ndof,itn,lamda,Theta,sy0,mt,emod,kappa,nu,b0,h0,b1,h1,alpha,vs);
GG=RR+SS*KS; HH(:,:,1)=PP+QQ*KS;
[KS]=make_pm(ndof,KK(:,2));
WW=GG+KS*HH(:,:,1);
TT(:,:,2)=HH(:,:,1)*inv(WW);
%
for ii=3:nn-1
   iim1=ii-1;
   [PP,QQ,RR,SS]=make_pqrs(ii,ndof,itn,lamda,Theta,sy0,mt,emod,kappa,nu,b0,h0,b1,h1,alpha,vs);
   GG=RR*TT(:,:,iim1)+SS; HH(:,:,iim1)=PP*TT(:,:,iim1)+QQ;
   [KS]=make_pm(ndof,KK(:,ii));
   WW=GG+KS*HH(:,:,iim1);
   TT(:,:,ii)=HH(:,:,iim1)*inv(WW);
end
%
%ii=nn
nnm1=nn-1;
[PP,QQ,RR,SS]=make_pqrs(nn,ndof,itn,lamda,Theta,sy0,mt,emod,kappa,nu,b0,h0,b1,h1,alpha,vs);
GG=RR*TT(:,:,nnm1)+SS; HH(:,:,nnm1)=PP*TT(:,:,nnm1)+QQ;
TTb=HH(:,:,nnm1)*inv(GG);
[KS]=make_pm(ndof,KK(:,nn));
UU=KS*TTb+II;
TT(:,:,nn)=TTb*inv(UU);
%
% 모드 계산
%
%switch bc_no
%   case {1,4,7,0} % 우단: 자유
%      % 보완할 것
%   case {2,5,8} % 우단: 단순지지
%      % 보완할 것
%   case {3,6,9} % 우단: 고정
      FFb(3,1)=-1;
      FFb(1:ndofm1,1)=inv(TTb(1:ndofm1,1:ndofm1))*TTb(1:ndofm1,3);
      XX(:,nn)=zeros(3,1);
      XX(:,nn-1)=inv(RR-SS*inv(QQ)*PP)*FFb;
%   otherwise
%end
%
for ii=nn-2:-1:2
   XX(:,ii)=TT(:,:,ii)*inv(HH(:,:,ii))*XX(:,ii+1);
end
XX(:,1)=inv(HH(:,:,1))*XX(:,2);
%
[max_value,max_index]=max(abs(XX(1,:)));
XX=XX/XX(1,max_index);
U=XX(1,:);  W=XX(2,:);  Th=XX(3,:);  X=[1:nn];  O=zeros(1,nn);
%
%figure(2)
plot(X,U,'r-',X,W,'g--',X,Th/10,'b:',X,O,'k-');
%plot(X,U,'r',X,W,'g',X,-Th/10,'b',X,O,'k');
%title('TICM')
xlabel('Nodes')
axis([1 nn -1.2 1.2]);
pause
%