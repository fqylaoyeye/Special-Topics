clear all
close all

video = false; % if want to save a movie, then let video = true
if(video)
    writerObj = VideoWriter('wheel_on_ground.mp4','MPEG-4');
    writerObj.FrameRate = 30;
    open(writerObj);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Variables for the mashgrid
n1=30;
n2=30;
h1=0.5;
h2=0.5;
x=1:h1:n1;
y=1:h2:n2;
r_m=0.5;

lmax_m = (n1-1)*n2+(n2-1)*n1;       % linkmunber
kmax_m = n1*n2;                     % piontmunber
R0_m = r_m*ones(lmax_m,1);                % rest length
S_m = 1000000*ones(lmax_m,1);           % Stiffness
D_m = 0*ones(lmax_m,1);              % Damping constant
m = 0.5;                               % Mass of the points
M_m = m*ones(kmax_m,1);               % Mass array

g = 9.8;                            % gravity constant


%link_h = 1:(n1-1)*n2;
%link_v = (n1-1)*n2+1;
nskip = 5;
X_m = zeros(kmax_m,3);              % setup position of the points on the meshgrid
for i = 1:n1
    for j = 1:n2
index=(i-1)*n1+j;
X_m(index,:)=[x(i),y(j),0];
    end
end 

% Indexing all points on the meshgrid with indeces of links
jj_m=zeros(lmax_m,1);kk_m=zeros(lmax_m,1);
for k=1:n1-1
    for l=1:n2
        jj_m(k+(n1-1)*(l-1))=k+n1*(l-1);  kk_m(k+(n1-1)*(l-1))=k+n1*(l-1)+1;
    end
end
for l=1:n2-1
    for k=1:n1
        jj_m((n1-1)*n2+l+(n2-1)*(k-1))=k+n1*(l-1);
        kk_m((n1-1)*n2+l+(n2-1)*(k-1))=k+n1*l;
    end
end

% Variables for Ball and links bewteen ball and the meshgrid
r_b = 1;
xpb=10;
ypb=10;
zpb=3;
R0_b = zeros(kmax_m,1);
X_b = [xpb,ypb,zpb];
S_b = zeros(kmax_m,1);
D_b = zeros(kmax_m,1);
M_b = 5;

% Merge all the variables

kmax = kmax_m + 1;                  % number of points on the meshgrid + the ball
lmax = lmax_m + kmax_m;             
% number of the links on the meshgrid + all the links that link ball and points on the meshgrid
jj = zeros(lmax,1);                   % index of points of one end of any link
kk = zeros(lmax,1);                   % index of points of the other end of any link
for p = 1:lmax_m
    jj(p)=jj_m(p);
    kk(p)=kk_m(p);
end
for p = (lmax_m+1):lmax
    jj(p) = p-lmax_m;
    kk(p) = kmax;
end
X = zeros(kmax,3);                   % set up position for all the points on the meshgrid
for k=1:kmax-1
    X(k,:)=X_m(k,:);                 % set up position for the ball specifically
end
X(kmax,:)=X_b;

R0 = [R0_m;R0_b];
S = [S_m;S_b];
D = [D_m;D_b];
M = [M_m;M_b];

%G=m*g;
%U_b = zeros(zpb,3);

%T=sqrt(2*zpb/g); % 设置下落时间
%H0=zpb;
%dt=0.001;%高度更新时间周期

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Set up Animation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lm = 1:lmax_m;                       % indeces of links on the meshgrid
lb = (lmax_m+1):lmax;                % indeces of links between ball and points on the meshgrid
pm = 1:kmax_m;                       % indeces of points on the meshgrid
figure(1);
fig = figure(1);
    phb=plot3(X_b(1,1),X_b(1,2),X_b(1,3),".","Markersize",40,"Color","r");  % plot the ball
    hold on
    phm=plot3(X(pm,1),X(pm,2),X(pm,3),"O","Markersize",6,"Color","b");         % plot points of the meshgrid
%    lhm=plot3([X(jj(lm),1),X(kk(lm),1)]',[X(jj(lm),2),X(kk(lm),2)]',...
%        [X(jj(lm),3),X(kk(lm),3)]',"Color","b","LineWidth",1);            % plot links on the meshgrid
%    lhb=plot3([X(jj(lb),1),X(kk(lb),1)]',[X(jj(lb),2),X(kk(lb),2)]',...
%        [X(jj(lb),3),X(kk(lb),3)]',"Color","g","LineWidth",1); 
    
 % plot links between ball and points on the meshgrid
 hold on
    axis equal
    axis manual
    axis([0,n1/2+1,0,n2/2+1,-5,5]);
 drawnow


tmax = 5;   % duration of simulation (s)
clockmax = 150000; %number of time steps
dt = tmax/clockmax; %(s)
X_save=zeros(clockmax,3);
t_save=zeros(clockmax,1);
E_save = zeros(clockmax,1);
Xm_z = zeros(kmax_m,1);                      % array to store the height of points on the meshgrid
Hmax_save = zeros(clockmax,1);               % store maximum height of the mesh
Hmin_save = zeros(clockmax,1);               % store minimum height of the mesh
Hb_save = zeros(clockmax,1);                 % store height of ball

fixed = [1,n1,n1*n2-n1+1,n1*n2];   % Set up fixed points

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% original codes for directly apply a force to the meshgrid %%%%%%%%
%forced = 15; %Set up where we apply the force

t_ext_start = 5; %time that external force starts (s)
%t_ext_stop =  1; %time that external force stops (s)
%F_ext = -10*m*g; %external force (Kg.m/s^2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%initialize velocity
U = zeros(kmax,3);
U(kmax,:) = [0,0,-20];

%initialize difference between the ball and any point on the meshgrid
X_mb = zeros(kmax_m,3);
d_mb = zeros(kmax_m,1);

for clock=1:clockmax
  t = clock*dt;

%%%%%%% compare distance between the ball and points on the meshgrid %%%%%%
  for k=1:kmax_m
    X_mb(k,:) = X(kmax,:) - X(k,:);
    d_mb(k,1) = sqrt(sum(X_mb(k,:).^2,2)); %compute distance between each point on 
                                           %the meshgrid and ball
  end
  d_min = min(d_mb);

% reset the links between ball and points on the meshgrid (disconnect them)
for l = (lmax_m+1):lmax
    R0(l) = 0;
    S(l) = 0;
    D(l) = 0;
end

% throw away all points between ball and meshgrid if distance is big
  if (d_min > r_b+r_m)
    l_mb = [];
  end

% select points with distance smaller than r_b+r_m
  if (d_min <= r_b+r_m)
    %d_mb2 = d_mb - d_min*ones(kmax_m,1);
    %l_mb = lmax_m + find(~d_mb2);         % find all the indices of the 
                                          % links that we want to connect 
    l_mb = lmax_m + find(d_mb <= r_b+r_m);
    
  end
  
% connect ball and selected points 
  for l = l_mb
    R0(l) = r_m+r_b;
    S(l) = 5000;
    D(l) = 0;
  end

%%%%%%%%%%%%%%%%%%%%% start to compute the motion %%%%%%%%%%%%%%%%%%%%%%%%%
  DX = X(jj,:) - X(kk,:); %link vectors
  DU = U(jj,:) - U(kk,:); %link velocity difference vectors
  R = sqrt(sum(DX.^2,2)); %link lengths
  T = S.*(R-R0) + (D./R).*sum(DX.*DU,2); %link tensions
  TR=T./R; %link tensions divided by link lengths
  FF=[TR,TR,TR].*DX; %link force vectors

  F=zeros(kmax,3); %initialize force array for mass points 
  F(1:kmax_m,3) = - m*g*ones(kmax_m,1);    %apply force of gravity to each meshgrid point
  F(kmax,3) = - M_b*g;       %apply force of gravity to the ball

  % For each link, add its force with correct sign
  % to the point at each end of the link:
  for link=1:lmax
    F(kk(link),:)=F(kk(link),:)+FF(link,:);
    F(jj(link),:)=F(jj(link),:)-FF(link,:);
  end

  disp(norm(U(kmax,:)));
  disp(norm(U(500,:)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% original codes for directly apply a force to the meshgrid %%%%%%%%
  %apply external force during specified time interval:
  %if((t_ext_start < t) && (t < t_ext_stop))
  %  F(forced,:) = F(forced,:) + F_ext;
  %end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  U = U + dt*F./[m,m,m]; %update velocities of all points,
  %but if the index of the point is on the list "fixed",
  %set its velocity equal to zero:
  U(fixed,:)=0; 

  X = X + dt*U; %update positions of all points

  %store some results for future plotting
  %X_save(clock,:)=X(n+2,:);
  t_save(clock) = t;


%   E_save(clock) = ...                                         % Save current total energy
%         1 / 2 * sum(M .* sqrt(sum(U .^ 2, 2)) .^ 2) + ...       % Total kinetic energy
%         sum(M .* X(:, 3) * g) + ...                             % Total gravitational potential energy
%         1 / 2 * sum(S .* (R-R0) .^ 2);                         % Total elastic potential energy

  Hb_save(clock) = X(kmax,3);
  
  Xm_z(1:kmax_m) = X(1:kmax_m,3);
  Hmin_save(clock) = min(Xm_z);
  Hmax_save(clock) = max(Xm_z);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%   ANIMATION   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf(fig)

    phb=plot3(X(kmax,1),X(kmax,2),X(kmax,3),".","Markersize",40,"Color","r");  % plot the ball
    hold on
    phm=plot3(X(pm,1),X(pm,2),X(pm,3),"O","Markersize",6,"Color","b");         % plot points of the meshgrid
    
    lhmb = plot3([X(jj(l_mb),1),X(kk(l_mb),1)]',[X(jj(l_mb),2),X(kk(l_mb),2)]',...
        [X(jj(l_mb),3),X(kk(l_mb),3)]',"Color","b","LineWidth",1);
%    lhm=plot3([X(jj(lm),1),X(kk(lm),1)]',[X(jj(lm),2),X(kk(lm),2)]',...
%        [X(jj(lm),3),X(kk(lm),3)]',"Color","b","LineWidth",1);
%        % plot links on the meshgrid
%    lhb=plot3([X(jj(lb),1),X(kk(lb),1)]',[X(jj(lb),2),X(kk(lb),2)]',...
%        [X(jj(lb),3),X(kk(lb),3)]',"Color","g","LineWidth",1);
    
% plot links between ball and points on the meshgrid
 hold on
    axis equal
    axis manual
    axis([0,n1/2+1,0,n2/2+1,-5,5]);
 drawnow

    if video
            frame = getframe(gcf);
            writeVideo(writerObj,frame);
    end
end


% figure(2)
% plot(t_save',E_save');

figure(3)
plot(t_save',Hmax_save',t_save',Hmin_save');

figure(4)
plot(t_save',Hb_save');

if  video 
    close(writerObj);
end


%for t=0:dt:T%循环开始
%   if(norm(X(jj,:) - X_b) > R)
%       H0=H0-1/2*g*t^2;
%   end
%h=plot3(xpb,ypb,H0,'O', 'markersize', 20);%绘制实时小球高度图
%   drawnow;%实时更新坐标图窗口
%   hold on
%end

