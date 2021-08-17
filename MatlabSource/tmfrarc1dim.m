%====================================================================
%
%   << (������) ������Ĺ��� �̿��� �簢(����)�ܸ� ��ũ�� �鳻 ���������ؼ� >>
%               ȸ�������� ���ܺ����� ������
%              ( ��Ʈ�� ���� : tmfrarc1dim.m )    2016. 03. 19
%              
%====================================================================
function tmfrarc1dim
%
clear all
%fid1=fopen('C:\tmfrarc1dim.out','w'); % ���� ���
%
% <�Է�> ��� ��
%it=1; % �����ܸ� ��ũ (ȸ�������� ���ܺ����� ����, ����-����) Ingenieur-Archiv (1979) 337-346 Fig.2
it=2; % �����ܸ� ��ũ (ȸ�������� ���ܺ����� ����, ����-����) Ingenieur-Archiv (1979) 337-346 Fig.3
%
[FRQ,fno,modeYN]=analysis_data(it); % ���������ؼ� �Ķ��Ÿ
[itn,ndof,nn,b0,h0,R,emod,roh,nu,kappa,sy0,ALPHA]=arc_data(it); % ��ũ �𵨸�
[bc_no,KK]=support(ndof,nn,it); % ������� �� �������� 
%
%timestart=clock;
ii=0;
for alpha_degree=ALPHA(1):ALPHA(2):ALPHA(3)
   alpha_degree % ��ũ�� ��ü ���� [degree]
   alpha=alpha_degree*pi/180; % [rad]
   Theta=linspace(0,alpha,nn) % ��ũ ������ ����
   ii=ii+1; Alphas(ii)=alpha_degree; % [rad]
   [Hz]=bsc(FRQ,fno,modeYN,itn,nn,ndof,Theta,b0,h0,R,sy0,emod,roh,bc_no,KK)
   Hzs(ii,1:fno)=Hz;
end
%timeend=clock;
%computation_time=etime(timeend,timestart)
%
%plot(Alphas,Lamdas(:,1),'ro',Alphas,Lamdas(:,2),'ro',Alphas,Lamdas(:,3),'ro',Alphas,Lamdas(:,4),'ro')
%axis([0.25 1.5 0 8])
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
% ���������ؼ� �Ķ��Ÿ
%
switch it
   case {1}
      fmin=1; % �˻� ���� ����ġ [Hz]
      fdel=1; % �˻� �ʱ� ���� ������ [Hz]
      fmax=1000; % �˻� ���� ������ [Hz]
      fno=4; % �˻��� ������������ �ְ� ����
      feps=1e-9; % ���������� �˻� ��� ����
      modeYN=0; % ������� ��� ���� ( Yes=1 / No=0 )
   case {2}
      fmin=1; % �˻� ���� ����ġ [Hz]
      fdel=1; % �˻� �ʱ� ���� ������ [Hz]
      fmax=1000; % �˻� ���� ������ [Hz]
      fno=4; % �˻��� ������������ �ְ� ����
      feps=1e-9; % ���������� �˻� ��� ����
      modeYN=0; % ������� ��� ���� ( Yes=1 / No=0 )
   otherwise
end
%
FRQ=[fmin fdel fmax feps]*2*pi; % [rad/s]
%
%-----------------------------------------------------------------------
function [itn,ndof,nn,b0,h0,R,emod,roh,nu,kappa,sy0,ALPHA]=arc_data(it)
%-----------------------------------------------------------------------
%
% ��ũ �𵨸�
%
ndof=3; % ������ ��������
%
switch it
   case {1}
      itn=50; % �ʵ���� ���� ȸ�� (Runge-Kutta-Gill Method)
      nn=5; % ���� �� (��� �� +1)
      ALPHA(1)=180; % ��ũ�� ���� (��� ���� ����) [degree]
      ALPHA(2)=5;  % ��ũ�� ���� (��� ���� ����) [degree]
      ALPHA(3)=180; % ��ũ�� ���� (��� ���� ����) [degree]
      R=1; % ��ũ�� ��� �ݰ� [m]
      kappa=0.823;
      sy0=50;
      h0=sqrt(12)*R/sy0; % �ܸ��� ���� (depth)
      b0=h0; % �ܸ��� �� (width)
      %
      emod=2.0e11; % ��ź����� [N/m2]
      nu=0.33; % ���ͼ� ��
      roh=7800; % �е� [kg/m3]
   case {2}
      itn=100; % �ʵ���� ���� ȸ�� (Runge-Kutta-Gill Method)
      nn=4; % ���� �� (��� �� +1)
      ALPHA(1)=90; % ��ũ�� ���� (��� ���� ����) [degree]
      ALPHA(2)=1;  % ��ũ�� ���� (��� ���� ����) [degree]
      ALPHA(3)=90; % ��ũ�� ���� (��� ���� ����) [degree]
      R=1; % ��ũ�� ��� �ݰ� [m]
      kappa=0.85;
      sy0=50;
      h0=sqrt(12)*R/sy0; % �ܸ��� ���� (depth)
      b0=h0; % �ܸ��� �� (width)
      %      
      emod=2.0e11; % ��ź����� [N/m2]
      nu=0.3; % ���ͼ� ��
      roh=7800; % �е� [kg/m3]
   otherwise
end
%
%-----------------------------------------------------------------------
function [bc_no,KK]=support(ndof,nn,it)
%-----------------------------------------------------------------------
%
% ������ �������� ������ �� ������� �Է�
%
KK=zeros(ndof,nn); % ���� 1���� ���� n(������ ����)�� ���������
%
switch it
% <1> (�´�: ���� & ���: ����),     <2> (�´�: ���� & ���: �ܼ�����),    <3> (�´�: ���� & ���: ����)
% <4> (�´�: �ܼ����� & ���: ����), <5> (�´�: �ܼ����� & ���: �ܼ�����), <6> (�´�: �ܼ����� & ���: ����)
% <7> (�´�: ���� & ���: ����),     <8> (�´�: ���� & ���: �ܼ�����),    <9> (�´�: ���� & ���: ����)
% <0> (�´�: ź������ & ���: ź������)
   case {1}
      bc_no=9;
      %KK(i,j)=1.0e6; % ���� j���� ������ i�� (������)��������� (i=1: ���� ����, i=2: �ݰ� ����, i=3: ȸ�� ����)
   case {2}
      bc_no=7;
      %KK(i,j)=1.0e6; % ���� j���� ������ i�� (������)��������� (i=1: ���� ����, i=2: �ݰ� ����, i=3: ȸ�� ����)
   otherwise
end
%
%---------------------------------------------------------------------------------------------------------------------------
function [Hz]=bsc(FRQ,fno,modeYN,itn,nn,ndof,Theta,b0,h0,R,sy0,emod,roh,bc_no,KK)
%---------------------------------------------------------------------------------------------------------------------------
%
% �̺й�
%
fmin=FRQ(1);  fdel=FRQ(2);  fmax=FRQ(3);  feps=FRQ(4); % [rad/s]
bmin=fmin;  w=fmin;
[isp]=tmmfree(w,itn,nn,ndof,Theta,b0,h0,R,sy0,emod,roh,bc_no,KK);
%
fno_cal=0;
while ((w <= fmax) & (fno_cal < fno))
   w=w+fdel;
   [is]=tmmfree(w,itn,nn,ndof,Theta,b0,h0,R,sy0,emod,roh,bc_no,KK);
   if is ~= isp
      bmax=w;
      while (bmax-bmin)/w> feps
         w=(bmin+bmax)/2;
         [is]=tmmfree(w,itn,nn,ndof,Theta,b0,h0,R,sy0,emod,roh,bc_no,KK);
         if is == isp
            bmin=w;
         else
            bmax=w;
         end
      end
      fno_cal=fno_cal+1;
      w=(bmin+bmax)/2;
      Hz(fno_cal)=w/(2*pi);
      %
      if modeYN == 0
          %XXs=0; % No data
      else
          %[XX]=tmmmode(); % �������
          %XXs(:,:,fno_cal)=XX;
      end
      %
      bmin=bmax;
      isp=-isp;
   else
      bmin=w;
   end
end
%
%-----------------------------------------------------------------------
function [is]=tmmfree(w,itn,nn,ndof,Theta,b0,h0,R,sy0,emod,roh,bc_no,KK)
%-----------------------------------------------------------------------
%
% ������Ĺ��� �̿��� ���������ؼ�
%
ndof1=ndof+1;   ndof2=ndof*2;
%
% ���޽�
%
if bc_no==0
   [TM]=make_pm(ndof,KK(:,1));
else
   TM=eye(ndof2);
end
%
for ii=2:nn
   [FM]=make_fm(ii,ndof,itn,Theta,w,b0,h0,R,sy0,emod,roh);
   TM=FM*TM;
   [PM]=make_pm(ndof,KK(:,ii));
   TM=PM*TM;
end
%
% ������������
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
% ����Ʈ��� ����� (�������� ������)
%
PM=eye(ndof*2);
PM(ndof+1,2)=Spring(1,1); % ���� ����
PM(ndof+2,1)=Spring(2,1); % �ݰ� ����
PM(ndof+3,3)=Spring(3,1); % ȸ�� ����
%
%-----------------------------------------------------------------------
function [FM]=make_fm(ii,ndof,itn,Theta,w,b0,h0,R,sy0,emod,roh)
%-----------------------------------------------------------------------
%
% �ʵ���� ����� (�����ܸ� ��ũ)
%
Y0=eye(ndof*2);
dTheta=Theta(ii)-Theta(ii-1);
dh=dTheta/itn;
[FM]=multi_runge_kutta_gill(ndof,Y0,dh,itn,w,b0,h0,R,sy0,emod,roh);
%
%-----------------------------------------------------------------------
function [Y0]=multi_runge_kutta_gill(ndof,Y0,dh,itn,w,b0,h0,R,sy0,emod,roh)
%-----------------------------------------------------------------------
%
% Runge-Kutta-Gill Method (�����ܸ� ��ũ)
%
aa=1.0-sqrt(0.5);   bb=1.0+sqrt(0.5);
Q0=zeros(ndof*2);   Y(:,:,1)=Y0;
%
for ii=1:itn
   [UM]=make_um(w,b0,h0,R,sy0,emod,roh);
   K1=dh*UM*Y0;    R1=(K1-2*Q0)/2;
   Y1=Y0+R1;       Q1=Q0+3*R1-K1/2;
   %
   [UM]=make_um(w,b0,h0,R,sy0,emod,roh);
   K2=dh*UM*Y1;   R2=aa*(K2-Q1);
   Y2=Y1+R2;      Q2=Q1+3*R2-aa*K2;
   %   
   [UM]=make_um(w,b0,h0,R,sy0,emod,roh);
   K3=dh*UM*Y2;   R3=bb*(K3-Q2);
   Y3=Y2+R3;      Q3=Q2+3*R3-bb*K3;
   %   
   [UM]=make_um(w,b0,h0,R,sy0,emod,roh);
   K4=dh*UM*Y3;   R4=(K4-2*Q3)/6;
   Y0=Y3+R4;      Q0=Q3+3*R4-K4/2;
end
%
%-----------------------------------------------------------------------
function [U]=make_um(w,b0,h0,R,sy0,emod,roh)
%-----------------------------------------------------------------------
%
% ������ U ��� ����� (�����ܸ� ��ũ)
%
k2=1/sy0^2;    area=b0*h0;
EA=emod*area;  RARW2=roh*area*R*w^2;
%
U=[0      -1       R     0       0      0
   1       0       0     R/EA    0      0
   0       0       0     0       0     -1/(EA*k2*R)
   0      -RARW2   0     0       1      0
  -RARW2   0       0    -1       0      0
   0       0       0     0       R      0]; 
%