%====================================================================
%
%   << (무차원) 전달행렬법을 이용한 사각 일정단면 아크의 면내 자유진동해석 >>
%               회전관성과 전단변형을 고려 가능함
%              ( 매트랩 파일 : tmfrarc2nodim.m )    2016. 03. 19
%
%====================================================================
function tmfrarc2nodim
%
clear all
%fid1=fopen('C:\tmfrarc2nodim.out','w'); % 파일 출력
%
% <입력> 계산 모델
%it=1; mt=0; % 일정단면 아크 (회전관성과 전단변형을 무시, 고정-고정) Ingenieur-Archiv (1979) 337-346 Fig.2
%it=1; mt=1; % 일정단면 티모센코 아크 (회전관성과 전단변형을 고려, 고정-고정) Ingenieur-Archiv (1979) 337-346 Fig.2
%it=2; mt=0; % 일정단면 아크 (회전관성과 전단변형을 무시, 고정-자유) Ingenieur-Archiv (1979) 337-346 Fig.3
it=2; mt=1; % 일정단면 티모센코 아크 (회전관성과 전단변형을 고려, 고정-자유) Ingenieur-Archiv (1979) 337-346 Fig.3
%
[FRQ,fno,modeYN]=analysis_data(it); % 자유진동해석 파라메타
[itn,ndof,nn,b0,h0,R,emod,roh,nu,kappa,sy0,alpha_s,alpha_d,alpha_e]=arc_data(it); % 아크 모델링
[bc_no,KK]=support(ndof,nn,it); % 경계조건 및 지지조건 
%
%timestart=clock;
ii=0;
for alpha_degree=alpha_s:alpha_d:alpha_e
   alpha_degree % 아크의 전체 각도 [degree]
   alpha=alpha_degree*pi/180; % [rad]
   Theta=linspace(0,alpha,nn); % 아크 각도의 분할
   ii=ii+1; Alphas(ii)=alpha_degree; % [rad]
   [Lamda]=bsc(FRQ,fno,modeYN,itn,nn,ndof,Theta,sy0,bc_no,KK,mt,emod,kappa,nu)
   Lamdas(ii,1:fno)=Lamda;
   Hz=(sqrt(emod/roh)*Lamda.^2/(R*sy0))/(2*pi) % 고유진동수 [Hz]
end
%timeend=clock;
%computation_time=etime(timeend,timestart)
%
%plot(Alphas,Lamdas(:,1),'bo',Alphas,Lamdas(:,2),'bo',Alphas,Lamdas(:,3),'bo',Alphas,Lamdas(:,4),'bo')
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
   otherwise
end
%
FRQ=[fmin fdel fmax feps]; % [lamda]
%
%-----------------------------------------------------------------------
function [itn,ndof,nn,b0,h0,R,emod,roh,nu,kappa,sy0,alpha_s,alpha_d,alpha_e]=arc_data(it)
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
      alpha_s=180; % 아크의 각도 (계산 시작 각도) [degree]
      alpha_d=5;  % 아크의 각도 (계산 증분 각도) [degree]
      alpha_e=180; % 아크의 각도 (계산 종료 각도) [degree]
      R=1; % 아크의 곡률 반경 [m]
      kappa=0.823;
      sy0=50;
      h0=sqrt(12)*R/sy0; % 단면의 높이 (depth)
      b0=h0; % 단면의 폭 (width)
      %
      emod=2.0e11; % 종탄성계수 [N/m2]
      nu=0.33; % 프와송 비
      roh=7800; % 밀도 [kg/m3]
   case {2}
      itn=100; % 필드행렬 적분 회수 (Runge-Kutta-Gill Method)
      nn=4; % 절점 수 (요소 수 +1)
      alpha_s=90; % 아크의 각도 (계산 시작 각도) [degree]
      alpha_d=1;  % 아크의 각도 (계산 증분 각도) [degree]
      alpha_e=90; % 아크의 각도 (계산 종료 각도) [degree]
      R=1; % 아크의 곡률 반경 [m]
      kappa=0.85;
      sy0=50;
      h0=sqrt(12)*R/sy0; % 단면의 높이 (depth)
      b0=h0; % 단면의 폭 (width)
      %      
      emod=2.0e11; % 종탄성계수 [N/m2]
      nu=0.3; % 프와송 비
      roh=7800; % 밀도 [kg/m3]
   otherwise
end
%
%-----------------------------------------------------------------------
function [bc_no,KK]=support(ndof,nn,it)
%-----------------------------------------------------------------------
%
% 무차원 기초지지 스프링 및 경계조건 입력
%
KK=zeros(ndof,nn); % 절점 1에서 절점 n(마지막 절점)의 스프링상수
%
switch it
% <1> (좌단: 자유 & 우단: 자유),     <2> (좌단: 자유 & 우단: 단순지지),    <3> (좌단: 자유 & 우단: 고정)
% <4> (좌단: 단순지지 & 우단: 자유), <5> (좌단: 단순지지 & 우단: 단순지지), <6> (좌단: 단순지지 & 우단: 고정)
% <7> (좌단: 고정 & 우단: 자유),     <8> (좌단: 고정 & 우단: 단순지지),    <9> (좌단: 고정 & 우단: 고정)
% <0> (좌단: 탄성지지 & 우단: 탄성지지)
   case {1}
      bc_no=9;
      %KK(i,j)=5.0; % 절점 j에서 자유도 i의 (무차원)스프링상수 (i=1: 접선 방향, i=2: 반경 방향, i=3: 회전 방향)
   case {2}
      bc_no=7;
      %KK(i,j)=5.0; % 절점 j에서 자유도 i의 (무차원)스프링상수 (i=1: 접선 방향, i=2: 반경 방향, i=3: 회전 방향)
   otherwise
end
%
%-----------------------------------------------------------------------
function [Lamda]=bsc(FRQ,fno,modeYN,itn,nn,ndof,Theta,sy0,bc_no,KK,mt,emod,kappa,nu);
%-----------------------------------------------------------------------
%
% 이분법
%
fmin=FRQ(1); fdel=FRQ(2); fmax=FRQ(3); feps=FRQ(4); % [lamda]
bmin=fmin;   lamda=fmin;
[isp]=tmmfree(lamda,itn,nn,ndof,Theta,sy0,bc_no,KK,mt,emod,kappa,nu);
%
fno_cal=0;
while ((lamda <= fmax) & (fno_cal < fno))
   lamda=lamda+fdel;
   [is]=tmmfree(lamda,itn,nn,ndof,Theta,sy0,bc_no,KK,mt,emod,kappa,nu);
   if is ~= isp
      bmax=lamda;
      while (bmax-bmin)/lamda > feps
         lamda=(bmin+bmax)/2;
         [is]=tmmfree(lamda,itn,nn,ndof,Theta,sy0,bc_no,KK,mt,emod,kappa,nu);
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
          %XXs=0; % No data
      else
          %[XX]=tmmmode(); % 고유모드
          %XXs(:,:,fno_cal)=XX;
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
function [is]=tmmfree(lamda,itn,nn,ndof,Theta,sy0,bc_no,KK,mt,emod,kappa,nu)
%-----------------------------------------------------------------------
%
% 전달행렬법을 이용한 자유진동해석
%
ndof1=ndof+1;   ndof2=ndof*2;
%
% 전달식
%
if bc_no==0
   [TM]=make_pm(ndof,KK(:,1));
else
   TM=eye(ndof2);
end
%
for ii=2:nn
   [FM]=make_fm(ii,ndof,itn,lamda,Theta,sy0,mt,emod,kappa,nu);
   TM=FM*TM;
   [PM]=make_pm(ndof,KK(:,ii));
   TM=PM*TM;
end
%
% 진동수방정식
%
switch bc_no
   case {1,4,7,0}
      BC_row=[ndof1:ndof2];
   case {2,5,8}
      BC_row=[1 2 ndof2];
   case {3,6,9}
      BC_row=[1:ndof];
end
%
switch bc_no
   case {1,2,3,0}
      BC_col=[1:ndof];
   case {4,5,6}
      BC_col=[3 4 5];
   case {7,8,9}
      BC_col=[ndof1:ndof2];
end
%
TT=TM(BC_row,BC_col);
%
if det(TT) <0
   is=-1;
else
   is=1;
end
%
%-----------------------------------------------------------------------
function [PM]=make_pm(ndof,Spring)
%-----------------------------------------------------------------------
%
% 포인트행렬 만들기 (기초지지 스프링)
%
PM=eye(ndof*2);
PM(ndof+1,2)=Spring(1,1); % 접선 방향
PM(ndof+2,1)=Spring(2,1); % 반경 방향
PM(ndof+3,3)=Spring(3,1); % 회전 방향
%
%-----------------------------------------------------------------------
function [FM]=make_fm(ii,ndof,itn,lamda,Theta,sy0,mt,emod,kappa,nu)
%-----------------------------------------------------------------------
%
% 필드행렬 만들기 (일정단면 아크)
%
Y0=eye(ndof*2);
dTheta=Theta(ii)-Theta(ii-1);
dh=dTheta/itn;
[FM]=multi_runge_kutta_gill(ndof,Y0,dh,itn,lamda,sy0,mt,emod,kappa,nu);
%
%-----------------------------------------------------------------------
function [Y0]=multi_runge_kutta_gill(ndof,Y0,dh,itn,lamda,sy0,mt,emod,kappa,nu)
%-----------------------------------------------------------------------
%
% Runge-Kutta-Gill Method (일정단면 아크)
%
aa=1.0-sqrt(0.5);   bb=1.0+sqrt(0.5);
Q0=zeros(ndof*2);   Y(:,:,1)=Y0;
%
for ii=1:itn
   [UM]=make_um(lamda,sy0,mt,emod,kappa,nu);
   K1=dh*UM*Y0;    R1=(K1-2*Q0)/2;
   Y1=Y0+R1;       Q1=Q0+3*R1-K1/2;
   %
   [UM]=make_um(lamda,sy0,mt,emod,kappa,nu);
   K2=dh*UM*Y1;   R2=aa*(K2-Q1);
   Y2=Y1+R2;      Q2=Q1+3*R2-aa*K2;
   %   
   [UM]=make_um(lamda,sy0,mt,emod,kappa,nu);
   K3=dh*UM*Y2;   R3=bb*(K3-Q2);
   Y3=Y2+R3;      Q3=Q2+3*R3-bb*K3;
   %   
   [UM]=make_um(lamda,sy0,mt,emod,kappa,nu);
   K4=dh*UM*Y3;   R4=(K4-2*Q3)/6;
   Y0=Y3+R4;      Q0=Q3+3*R4-K4/2;
end
%
%-----------------------------------------------------------------------
function [U]=make_um(lamda,sy0,mt,emod,kappa,nu)
%-----------------------------------------------------------------------
%
% 무차원 U 행렬 만들기 (일정단면 아크)
%
if mt == 0 % 회전관성과 전단변형을 고려 안함
   k2=1/sy0^2;   l4s2=lamda^4/sy0^2;
   %
   U=[0      -1      1     0     0     0
      1       0      0     1     0     0
      0       0      0     0     0    -1/k2
      0      -l4s2   0     0     1     0
     -l4s2    0      0    -1     0     0
      0       0      0     0     1     0]; 
else % 회전관성과 전단변형을 고려함
   k2=1/sy0^2;   l4s2=lamda^4/sy0^2;
   k_12=2*k2;    k_22=k2;
   gmod=emod/(2*(1+nu)); 
   ekg=emod/(kappa*gmod);
   l4s2k2=l4s2*(1+k2);
   l4s2k_1=l4s2*k_12;
   l4s2k_2=l4s2*k_22;
   %
   U=[0        -1          1          0    ekg    0
      1         0          0          1    0      0
      0         0          0          0    0     -1/k2
      0        -l4s2k2    -l4s2k_1    0    1      0
     -l4s2k2    0          0         -1    0      0
      0         l4s2k_1    l4s2k_2    0    1      0]; 
end
%