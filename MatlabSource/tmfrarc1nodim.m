%====================================================================
%
%   << (������) ������Ĺ��� �̿��� �簢 �����ܸ� ��ũ�� �鳻 ���������ؼ� >>
%               ȸ�������� ���ܺ����� ������
%              ( ��Ʈ�� ���� : tmfrarc1nodim.m )    2016. 03. 16
%
%====================================================================
function tmfrarc1nodim
%
clear all
%fid1=fopen('C:\tmfrarc1nodim.out','w'); % ���� ���
%
% <�Է�> ��� ��
%it=1; % �����ܸ� ��ũ (ȸ�������� ���ܺ����� ����, ����-����) Ingenieur-Archiv (1979) 337-346 Fig.2
it=2; % �����ܸ� ��ũ (ȸ�������� ���ܺ����� ����, ����-����) Ingenieur-Archiv (1979) 337-346 Fig.3
%
[FRQ,fno,modeYN]=analysis_data(it); % ���������ؼ� �Ķ��Ÿ
[itn,ndof,nn,b0,h0,R,emod,roh,nu,kappa,sy0,alpha_s,alpha_d,alpha_e]=arc_data(it); % ��ũ �𵨸�
[bc_no,KK]=support(ndof,nn,it); % ������� �� �������� 
%
%timestart=clock;
ii=0;
for alpha_degree=alpha_s:alpha_d:alpha_e
   alpha_degree % ��ũ�� ��ü ���� [degree]
   alpha=alpha_degree*pi/180; % [rad]
   Theta=linspace(0,alpha,nn); % ��ũ ������ ����
   ii=ii+1; Alphas(ii)=alpha_degree; % [rad]
   [Lamda]=bsc(FRQ,fno,modeYN,itn,nn,ndof,Theta,sy0,bc_no,KK)
   Lamdas(ii,1:fno)=Lamda;
   Hz=(sqrt(emod/roh)*Lamda.^2/(R*sy0))/(2*pi) % ���������� [Hz]
end
%timeend=clock;
%computation_time=etime(timeend,timestart)
%
%plot(Alphas,Lamdas(:,1),'ro',Alphas,Lamdas(:,2),'ro',Alphas,Lamdas(:,3),'ro',Alphas,Lamdas(:,4),'ro')
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
      fmin=1; % �˻� ���� ����ġ [lamda]
      fdel=0.1; % �˻� �ʱ� ���� ����ġ
      fmax=30; % �˻� ���� ����ġ
      fno=4; % �˻��� ����ġ�� �ְ� ����
      feps=1e-9; % ����ġ �˻� ��� ���� [lamda]
      modeYN=0; % ������� ��� ���� ( Yes=1 / No=0 )
   case {2}
      fmin=1; % �˻� ���� ����ġ [lamda]
      fdel=0.1; % �˻� �ʱ� ���� ����ġ
      fmax=10; % �˻� ���� ����ġ
      fno=4; % �˻��� ����ġ�� �ְ� ����
      feps=1e-9; % ����ġ �˻� ��� ���� [lamda]
      modeYN=0; % ������� ��� ���� ( Yes=1 / No=0 )
   otherwise
end
%
FRQ=[fmin fdel fmax feps]; % [lamda]
%
%-----------------------------------------------------------------------
function [itn,ndof,nn,b0,h0,R,emod,roh,nu,kappa,sy0,alpha_s,alpha_d,alpha_e]=arc_data(it)
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
      alpha_s=180; % ��ũ�� ���� (��� ���� ����) [degree]
      alpha_d=5;  % ��ũ�� ���� (��� ���� ����) [degree]
      alpha_e=180; % ��ũ�� ���� (��� ���� ����) [degree]
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
      alpha_s=90; % ��ũ�� ���� (��� ���� ����) [degree]
      alpha_d=1;  % ��ũ�� ���� (��� ���� ����) [degree]
      alpha_e=90; % ��ũ�� ���� (��� ���� ����) [degree]
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
      %KK(i,j)=5.0; % ���� j���� ������ i�� (������)��������� (i=1: ���� ����, i=2: �ݰ� ����, i=3: ȸ�� ����)
   case {2}
      bc_no=7;
      %KK(i,j)=5.0; % ���� j���� ������ i�� (������)��������� (i=1: ���� ����, i=2: �ݰ� ����, i=3: ȸ�� ����)
   otherwise
end
%
%-----------------------------------------------------------------------
function [Lamda]=bsc(FRQ,fno,modeYN,itn,nn,ndof,Theta,sy0,bc_no,KK);
%-----------------------------------------------------------------------
%
% �̺й�
%
fmin=FRQ(1); fdel=FRQ(2); fmax=FRQ(3); feps=FRQ(4); % [lamda]
bmin=fmin;   lamda=fmin;
[isp]=tmmfree(lamda,itn,nn,ndof,Theta,sy0,bc_no,KK);
%
fno_cal=0;
while ((lamda <= fmax) & (fno_cal < fno))
   lamda=lamda+fdel;
   [is]=tmmfree(lamda,itn,nn,ndof,Theta,sy0,bc_no,KK);
   if is ~= isp
      bmax=lamda;
      while (bmax-bmin)/lamda > feps
         lamda=(bmin+bmax)/2;
         [is]=tmmfree(lamda,itn,nn,ndof,Theta,sy0,bc_no,KK);
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
          %[XX]=tmmmode(); % �������
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
function [is]=tmmfree(lamda,itn,nn,ndof,Theta,sy0,bc_no,KK)
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
   [FM]=make_fm(ii,ndof,itn,lamda,Theta,sy0);
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
function [FM]=make_fm(ii,ndof,itn,lamda,Theta,sy0)
%-----------------------------------------------------------------------
%
% �ʵ���� ����� (�����ܸ� ��ũ)
%
Y0=eye(ndof*2);
dTheta=Theta(ii)-Theta(ii-1);
dh=dTheta/itn;
[FM]=multi_runge_kutta_gill(ndof,Y0,dh,itn,lamda,sy0);
%
%-----------------------------------------------------------------------
function [Y0]=multi_runge_kutta_gill(ndof,Y0,dh,itn,lamda,sy0)
%-----------------------------------------------------------------------
%
% Runge-Kutta-Gill Method (�����ܸ� ��ũ)
%
aa=1.0-sqrt(0.5);   bb=1.0+sqrt(0.5);
Q0=zeros(ndof*2);   Y(:,:,1)=Y0;
%
for ii=1:itn
   [UM]=make_um(lamda,sy0);
   K1=dh*UM*Y0;    R1=(K1-2*Q0)/2;
   Y1=Y0+R1;       Q1=Q0+3*R1-K1/2;
   %
   [UM]=make_um(lamda,sy0);
   K2=dh*UM*Y1;   R2=aa*(K2-Q1);
   Y2=Y1+R2;      Q2=Q1+3*R2-aa*K2;
   %   
   [UM]=make_um(lamda,sy0);
   K3=dh*UM*Y2;   R3=bb*(K3-Q2);
   Y3=Y2+R3;      Q3=Q2+3*R3-bb*K3;
   %   
   [UM]=make_um(lamda,sy0);
   K4=dh*UM*Y3;   R4=(K4-2*Q3)/6;
   Y0=Y3+R4;      Q0=Q3+3*R4-K4/2;
end
%
%-----------------------------------------------------------------------
function [U]=make_um(lamda,sy0)
%-----------------------------------------------------------------------
%
% ������ U ��� ����� (�����ܸ� ��ũ) ȸ�������� ���ܺ����� ��� ����
%
k2=1/sy0^2;   l4s2=lamda^4/sy0^2;
%
   U=[0      -1      1     0     0     0
      1       0      0     1     0     0
      0       0      0     0     0    -1/k2
      0      -l4s2   0     0     1     0
     -l4s2    0      0    -1     0     0
      0       0      0     0     1     0]; 
%