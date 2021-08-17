%====================================================================
%
%   << (������) ������Ĺ��� �̿��� �簢(����/��)�ܸ� ��ũ�� �鳻 ���������ؼ� >>
%               ȸ�������� ���ܺ����� ��� ������
%              ( ��Ʈ�� ���� : tmfrarc3dim.m )    2016. 03. 22
%              
%====================================================================
function tmfrarc3dim
%
clear all
%fid1=fopen('C:\tmfrarc3dim.out','w'); % ���� ���
%
% <�Է�> ��� ��
%        mt=0; ȸ�������� ���ܺ����� ����, mt=1; ȸ�������� ���ܺ����� ���
%        vs=0; ���� �ܸ�, vs=1; ���ܸ�(b: ����, h: ���� ����)
%
%it=1; mt=0; vs=0; % �����ܸ� ��ũ (����-����) Ingenieur-Archiv (1979) 337-346 Fig.2
%it=1; mt=1; vs=0; % �����ܸ� Ƽ���� ��ũ (����-����) Ingenieur-Archiv (1979) 337-346 Fig.2
%it=2; mt=0; vs=0; % �����ܸ� ��ũ (����-����) Ingenieur-Archiv (1979) 337-346 Fig.3
%it=2; mt=1; vs=0; % �����ܸ� Ƽ���� ��ũ (����-����) Ingenieur-Archiv (1979) 337-346 Fig.3
%it=2; mt=0; vs=1; % ���ܸ� ��ũ (����-����) Ingenieur-Archiv (1979) 337-346 Fig.3
%it=2; mt=1; vs=1; % ���ܸ� Ƽ���� ��ũ (����-����) Ingenieur-Archiv (1979) 337-346 Fig.3
it=3; mt=1; vs=0; % �����ܸ� Ƽ���� ��ũ (����-����) Ingenieur-Archiv (1979) 337-346 Fig.8a
%it=4; mt=1; vs=0; % �����ܸ� Ƽ���� ��ũ (����-����) Ingenieur-Archiv (1979) 337-346 Fig.8b
%
[FRQ,fno,modeYN]=analysis_data(it); % ���������ؼ� �Ķ��Ÿ
[itn,ndof,nn,b0,h0,b1,h1_s,h1_d,h1_e,R,emod,roh,nu,kappa,sy0,ALPHA]=arc_data(it); % ��ũ �𵨸�
[bc_no,KK]=support(ndof,nn,it); % ������� �� �������� 
%
%timestart=clock;
ii=0;
for alpha_degree=ALPHA(1):ALPHA(2):ALPHA(3)
   alpha_degree % ��ũ�� ��ü ���� [degree]
   alpha=alpha_degree*pi/180; % [rad]
   Theta=linspace(0,alpha,nn); % ��ũ ������ ����
   %ii=ii+1; Alphas(ii)=alpha_degree; % [rad]
   for h1=h1_s:h1_d:h1_e
      h1
      ii=ii+1; H1_h0s(ii)=h1/h0;
      [Hz,XXs]=bsc(FRQ,fno,modeYN,itn,nn,ndof,Theta,sy0,bc_no,KK,mt,emod,kappa,nu,roh,R,b0,h0,b1,h1,alpha,vs);
      Hzs(ii,1:fno)=Hz;
   end   
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
      fmin=1; % �˻� ���� ������ [Hz]
      fdel=1; % �˻� �ʱ� ���� ������ [Hz]
      fmax=1000; % �˻� ���� ������ [Hz]
      fno=4; % �˻��� ������������ �ְ� ����
      feps=1e-9; % ���������� �˻� ��� ����
      modeYN=0; % ������� ��� ���� ( Yes=1 / No=0 )
   case {2}
      fmin=1; % �˻� ���� ������ [Hz]
      fdel=1; % �˻� �ʱ� ���� ������ [Hz]
      fmax=1000; % �˻� ���� ������ [Hz]
      fno=4; % �˻��� ������������ �ְ� ����
      feps=1e-9; % ���������� �˻� ��� ����
      modeYN=0; % ������� ��� ���� ( Yes=1 / No=0 )
   case {3}
      fmin=10; % �˻� ���� ������ [Hz]
      fdel=10; % �˻� �ʱ� ���� ������ [Hz]
      fmax=1000; % �˻� ���� ������ [Hz]
      fno=4; % �˻��� ������������ �ְ� ����
      feps=1e-9; % ���������� �˻� ��� ����
      modeYN=1; % ������� ��� ���� ( Yes=1 / No=0 )
   case {4}
      fmin=10; % �˻� ���� ������ [Hz]
      fdel=10; % �˻� �ʱ� ���� ������ [Hz]
      fmax=2000; % �˻� ���� ������ [Hz]
      fno=4; % �˻��� ������������ �ְ� ����
      feps=1e-9; % ���������� �˻� ��� ����
      modeYN=1; % ������� ��� ���� ( Yes=1 / No=0 )
   otherwise
end
%
FRQ=[fmin fdel fmax feps]*2*pi; % [rad/s]
%
%-----------------------------------------------------------------------
function [itn,ndof,nn,b0,h0,b1,h1_s,h1_d,h1_e,R,emod,roh,nu,kappa,sy0,ALPHA]=arc_data(it)
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
      % ���ǹ��� ����Ÿ(�����ܸ��̹Ƿ�)
      b1=b0;
      h1_s=h0;
      h1_d=h0;
      h1_e=h0;      
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
      b1=b0;
      h1_s=0.25*h0;
      h1_d=0.25*h0;
      h1_e=1.5*h0;
      %      
      emod=2.0e11; % ��ź����� [N/m2]
      nu=0.3; % ���ͼ� ��
      roh=7800; % �е� [kg/m3]
   case {3}
      itn=50; % �ʵ���� ���� ȸ�� (Runge-Kutta-Gill Method)
      nn=51; % ���� �� (��� �� +1)
      ALPHA(1)=90; % ��ũ�� ���� (��� ���� ����) [degree]
      ALPHA(2)=1;  % ��ũ�� ���� (��� ���� ����) [degree]
      ALPHA(3)=90; % ��ũ�� ���� (��� ���� ����) [degree]
      R=1; % ��ũ�� ��� �ݰ� [m]
      kappa=0.85;
      sy0=100;
      h0=sqrt(12)*R/sy0; % �ܸ��� ���� (depth)
      b0=h0; % �ܸ��� �� (width)
      %      
      emod=2.0e11; % ��ź����� [N/m2]
      nu=0.3; % ���ͼ� ��
      roh=7800; % �е� [kg/m3]
      % ���ǹ��� ����Ÿ(�����ܸ��̹Ƿ�)
      b1=b0;
      h1_s=h0;
      h1_d=h0;
      h1_e=h0;
   case {4}
      itn=50; % �ʵ���� ���� ȸ�� (Runge-Kutta-Gill Method)
      nn=51; % ���� �� (��� �� +1)
      ALPHA(1)=90; % ��ũ�� ���� (��� ���� ����) [degree]
      ALPHA(2)=1;  % ��ũ�� ���� (��� ���� ����) [degree]
      ALPHA(3)=90; % ��ũ�� ���� (��� ���� ����) [degree]
      R=1; % ��ũ�� ��� �ݰ� [m]
      kappa=0.85;
      sy0=30;
      h0=sqrt(12)*R/sy0; % �ܸ��� ���� (depth)
      b0=h0; % �ܸ��� �� (width)
      %      
      emod=2.0e11; % ��ź����� [N/m2]
      nu=0.3; % ���ͼ� ��
      roh=7800; % �е� [kg/m3]
      % ���ǹ��� ����Ÿ(�����ܸ��̹Ƿ�)
      b1=b0;
      h1_s=h0;
      h1_d=h0;
      h1_e=h0;    
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
   case {1,3,4}
      bc_no=9;
      %KK(i,j)=1.0e6; % ���� j���� ������ i�� (������)��������� (i=1: ���� ����, i=2: �ݰ� ����, i=3: ȸ�� ����)
   case {2}
      bc_no=7;
      %KK(i,j)=1.0e6; % ���� j���� ������ i�� (������)��������� (i=1: ���� ����, i=2: �ݰ� ����, i=3: ȸ�� ����)
   otherwise
end
%
%---------------------------------------------------------------------------------------------------------------------------
function [Hz,XXs]=bsc(FRQ,fno,modeYN,itn,nn,ndof,Theta,sy0,bc_no,KK,mt,emod,kappa,nu,roh,R,b0,h0,b1,h1,alpha,vs)
%---------------------------------------------------------------------------------------------------------------------------
%
% �̺й�
%
fmin=FRQ(1);  fdel=FRQ(2);  fmax=FRQ(3);  feps=FRQ(4); % [rad/s]
bmin=fmin;  w=fmin;
[isp]=tmmfree(w,itn,nn,ndof,Theta,b0,h0,R,sy0,emod,roh,bc_no,KK,mt,kappa,nu,h1,alpha,vs);
%
fno_cal=0;
while ((w <= fmax) & (fno_cal < fno))
   w=w+fdel;
   [is]=tmmfree(w,itn,nn,ndof,Theta,b0,h0,R,sy0,emod,roh,bc_no,KK,mt,kappa,nu,h1,alpha,vs);
   if is ~= isp
      bmax=w;
      while (bmax-bmin)/w> feps
         w=(bmin+bmax)/2;
          [is]=tmmfree(w,itn,nn,ndof,Theta,b0,h0,R,sy0,emod,roh,bc_no,KK,mt,kappa,nu,h1,alpha,vs);
         if is == isp
            bmin=w;
         else
            bmax=w;
         end
      end
      fno_cal=fno_cal+1;
      w=(bmin+bmax)/2;
      Hz(fno_cal)=w/(2*pi)
      %
      if modeYN == 0
          XXs=0; % No data
      else
          [XX]=tmmmode(w,itn,nn,ndof,Theta,b0,h0,R,sy0,emod,roh,bc_no,KK,mt,kappa,nu,h1,alpha,vs); % �������
          XXs(:,:,fno_cal)=XX;
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
function [is]=tmmfree(w,itn,nn,ndof,Theta,b0,h0,R,sy0,emod,roh,bc_no,KK,mt,kappa,nu,h1,alpha,vs)
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
   [FM]=make_fm(ii,ndof,itn,Theta,w,b0,h0,R,sy0,emod,roh,mt,kappa,nu,h1,alpha,vs);
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
function [FM]=make_fm(ii,ndof,itn,Theta,w,b0,h0,R,sy0,emod,roh,mt,kappa,nu,h1,alpha,vs)
%-----------------------------------------------------------------------
%
% �ʵ���� ����� (�����ܸ�/���ܸ� ��ũ)
%
Y0=eye(ndof*2);
dTheta=Theta(ii)-Theta(ii-1);
th0=Theta(ii-1);
dh=dTheta/itn;
[FM]=multi_runge_kutta_gill(ndof,Y0,dh,itn,w,b0,h0,R,sy0,emod,roh,mt,kappa,nu,h1,alpha,th0,vs);
%
%-----------------------------------------------------------------------
function [Y0]=multi_runge_kutta_gill(ndof,Y0,dh,itn,w,b0,h0,R,sy0,emod,roh,mt,kappa,nu,h1,alpha,th0,vs)
%-----------------------------------------------------------------------
%
% Runge-Kutta-Gill Method (�����ܸ�/���ܸ� ��ũ)
%
aa=1.0-sqrt(0.5);   bb=1.0+sqrt(0.5);
Q0=zeros(ndof*2);   Y(:,:,1)=Y0;
%
for ii=1:itn
   th=th0+(ii-1)*dh;
   [UM]=make_um(th,w,b0,h0,R,sy0,emod,roh,mt,kappa,nu,h1,alpha,vs);
   K1=dh*UM*Y0;    R1=(K1-2*Q0)/2;
   Y1=Y0+R1;       Q1=Q0+3*R1-K1/2;
   %
   [UM]=make_um(th+dh/2,w,b0,h0,R,sy0,emod,roh,mt,kappa,nu,h1,alpha,vs);
   K2=dh*UM*Y1;   R2=aa*(K2-Q1);
   Y2=Y1+R2;      Q2=Q1+3*R2-aa*K2;
   %   
   [UM]=make_um(th+dh/2,w,b0,h0,R,sy0,emod,roh,mt,kappa,nu,h1,alpha,vs);
   K3=dh*UM*Y2;   R3=bb*(K3-Q2);
   Y3=Y2+R3;      Q3=Q2+3*R3-bb*K3;
   %   
   [UM]=make_um(th+dh,w,b0,h0,R,sy0,emod,roh,mt,kappa,nu,h1,alpha,vs);
   K4=dh*UM*Y3;   R4=(K4-2*Q3)/6;
   Y0=Y3+R4;      Q0=Q3+3*R4-K4/2;
end
%
%-----------------------------------------------------------------------
function [U]=make_um(th,w,b0,h0,R,sy0,emod,roh,mt,kappa,nu,h1,alpha,vs)
%-----------------------------------------------------------------------
%
% ������ U ��� ����� (�����ܸ�/���ܸ� ��ũ)
%
switch vs
   case {0} % �����ܸ�(b: ����, h: ����)
      b=b0;
      h=h0;
   case {1} % ���ܸ�(b: ����, h: ����)
      b=b0;
      h=h0-(h0-h1)*(th/alpha)^1;
   otherwise
end
area=b*h; % ���� 
i2y=b*h^3/12; % ����2�����Ʈ
k2=i2y/(R^2*area);
%
if mt == 0 % ȸ�������� ���ܺ����� ��� ����
   %k2=1/sy0^2;    area=b0*h0;
   EA=emod*area;  RARW2=roh*area*R*w^2;
   %
   U=[0      -1       R     0       0      0
      1       0       0     R/EA    0      0
      0       0       0     0       0     -1/(EA*k2*R)
      0      -RARW2   0     0       1      0
     -RARW2   0       0    -1       0      0
      0       0       0     0       R      0]; 
else % ȸ�������� ���ܺ����� �����
   %k2=1/sy0^2;    area=b0*h0;
   EA=emod*area;  RARW2=roh*area*R*w^2;
   k_12=2*k2;    k_22=k2;
   gmod=emod/(2*(1+nu));
   kga=kappa*gmod*area;   
   RARW2K2=RARW2*(1+k2);
   RAR2W2K12=RARW2*R*k_12;
   RAR3W2K22=RARW2*R^2*k_22;
   %
   U=[0        -1            R              0       R/kga     0
      1         0            0              R/EA    0         0
      0         0            0              0       0        -1/(EA*k2*R)
      0        -RARW2K2     -RAR2W2K12      0       1         0
     -RARW2K2   0            0             -1       0         0
      0         RAR2W2K12    RAR3W2K22      0       R         0]; 
end
%
%-----------------------------------------------------------------------
function [Z]=tmmmode(w,itn,nn,ndof,Theta,b0,h0,R,sy0,emod,roh,bc_no,KK,mt,kappa,nu,h1,alpha,vs)
%-----------------------------------------------------------------------
%
% ������Ĺ��� �̿��� ���������ؼ�
%
ndof1=ndof+1;   ndof2=ndof*2;
TM=zeros(ndof2,ndof2,nn);
%
% ���޽�
%
if bc_no==0
   [TM(:,:,1)]=make_pm(ndof,KK(:,1));
else
   TM(:,:,1)=eye(ndof2);
end
%
for ii=2:nn
   [FM]=make_fm(ii,ndof,itn,Theta,w,b0,h0,R,sy0,emod,roh,mt,kappa,nu,h1,alpha,vs);
   TM(:,:,ii)=FM*TM(:,:,ii-1);
   [PM]=make_pm(ndof,KK(:,ii));
   TM(:,:,ii)=PM*TM(:,:,ii);
end
%
% ��� ���
%
switch bc_no
   case {1,4,7,0} % ���: ����
      BC_row=[ndof1:ndof2];
   case {2,5,8} % ���: �ܼ�����
      BC_row=[1 2 ndof2];
   case {3,6,9} % ���: ����
      BC_row=[1:ndof];
end
%
switch bc_no
   case {1,2,3,0} % �´�: ����
      BC_col=[1:ndof];
   case {4,5,6} % �´�: �ܼ�����
      BC_col=[3 4 5];
   case {7,8,9} % �´�: ����
      BC_col=[ndof1:ndof2];
end
%
TT=TM(BC_row,BC_col,nn);
%
Z=zeros(ndof2,nn);
ndofm1=ndof-1;
switch bc_no
   case {1,2,3,0} % �´�: ����
      Z(1:ndofm1,1)=inv(TT(1:ndofm1,1:ndofm1))*TT(1:ndofm1,ndof);
      Z(ndof,1)=-1;
   case {4,5,6} % �´�: �ܼ�����
      Z(3:4,1)=inv(TT(1:ndofm1,1:ndofm1))*TT(1:ndofm1,ndof);
      Z(5,1)=-1;
   case {7,8,9} % �´�: ����
      Z(ndof1:ndof2-1,1)=inv(TT(1:ndofm1,1:ndofm1))*TT(1:ndofm1,ndof);
      Z(ndof2,1)=-1;
   otherwise
end
%
for ii=2:nn
    Z(:,ii)=TM(:,:,ii)*Z(:,1);
end
U=Z(1,:);
[max_value,max_index]=max(abs(U));
Z=Z/U(max_index);
U=Z(1,:);
W=Z(2,:);
Th=Z(3,:);
X=[1:nn];
O=zeros(1,nn);
%plot(X,-U,'r',X,-W,'g',X,Th/10,'b',X,O,'k');
plot(X,U,'r',X,W,'g',X,-Th/10,'b',X,O,'k');
axis([1 nn -1.2 1.2]);
pause
%