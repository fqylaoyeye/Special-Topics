
clear all
close all

video = false; % if want to save a movie, then let video = true
if(video)
    writerObj = VideoWriter('wheel_on_ground.mp4','MPEG-4');
    writerObj.FrameRate = 30;
    open(writerObj);
end

% Variables for the mashgrid
n1=30;
n2=30;
h1=1;
h2=1;
x=1:h1:n1;
y=1:h2:n2;

lmax = (n1-1)*n2+(n2-1)*n1;       % linkmunber
kmax = n1*n2;                     % piontmunber
R0 = ones(lmax,1);                % rest length
S = 1000000*ones(lmax,1);           % Stiffness
D = 0*ones(lmax,1);              % Damping constant
m = 0.5;                            % Mass of the points
M = m*ones(kmax,1);               % Mass array

g = 9.8;


link_h = 1:(n1-1)*n2;
link_v = (n1-1)*n2+1;
nskip = 5;
X = zeros(kmax,3);
for i = 1:n1
    for j = 1:n2
index=(i-1)*n1+j;
X(index,:)=[x(i),y(j),0];
    end
end

jj=zeros(lmax,1);kk=zeros(lmax,1);
for k=1:n1-1
    for l=1:n2
        jj(k+(n1-1)*(l-1))=k+n1*(l-1);  kk(k+(n1-1)*(l-1))=k+n1*(l-1)+1;
    end
end
for l=1:n2-1
    for k=1:n1
        jj((n1-1)*n2+l+(n2-1)*(k-1))=k+n1*(l-1);
        kk((n1-1)*n2+l+(n2-1)*(k-1))=k+n1*l;
    end
end

figure(1);
fig = figure(1);
    h1 = plot3(5,5,5,"O","Markersize",15,"Color","r");
    hold on
    h=plot3(X(:,1),X(:,2),X(:,3),"O","Markersize",3,"Color","b");
    hold off
    axis equal
    axis manual
    axis([0,n1,0,n2,-3,5]);
drawnow

tmax = 10;   % duration of simulation (s)
clockmax = 100000; %number of time steps
dt = tmax/clockmax; %(s)
Ek=zeros(clockmax,1);
Epg=zeros(clockmax,1);
Eps=zeros(clockmax,1);
E_save=zeros(clockmax,1);
t_save=zeros(clockmax,1);

fixed = [1,n1,n1*n2-n1+1,n1*n2];   % Set up fixed points
forced = 400; %Set up where we apply the force

t_ext_start = 1; %time that external force starts (s)
t_ext_stop =  2; %time that external force stops (s)
F_ext = -10*m*g; %external force (Kg.m/s^2)

%initial velocity
U = zeros(kmax,3);

for clock=1:clockmax
  t = clock*dt;
  DX = X(jj,:) - X(kk,:); %link vectors
  DU = U(jj,:) - U(kk,:); %link velocity difference vectors
  R = sqrt(sum(DX.^2,2)); %link lengths
  T = S.*(R-R0) + (D./R).*sum(DX.*DU,2); %link tensions
  TR=T./R; %link tensions divided by link lengths
  FF=[TR,TR,TR].*DX; %link force vectors

  F=zeros(kmax,3); %initialize force array for mass points 
  F(:,3) = - m*g;    %apply force of gravity to each link point

  % For each link, add its force with correct sign
  % to the point at each end of the link:
  for link=1:lmax
    F(kk(link),:)=F(kk(link),:)+FF(link,:);
    F(jj(link),:)=F(jj(link),:)-FF(link,:);
  end

  %apply external force during specified time interval:
  if((t_ext_start < t) && (t < t_ext_stop))
    F(forced,:) = F(forced,:) + F_ext;
  end

  U = U + dt*F./[m,m,m]; %update velocities of all points,
  %but if the index of the point is on the list "fixed",
  %set its velocity equal to zero:
  U(fixed,:)=0; 

  X = X + dt*U; %update positions of all points

  %store some results for future plotting
  %X_save(clock,:)=X(n+2,:);
  t_save(clock) = t;

  for k = 1:kmax
      Ek(clock) = Ek(clock) + 1/2*M(k)*(U(k,1)^2+U(k,2)^2+U(k,3)^2);
      Epg(clock) = Epg(clock) + M(k)*g*X(k,3);
  end
  for l = 1:lmax
      Eps(clock) = Eps(clock) + 1/2* S(l) *(R(l)-R0(l)) ^2;
  end
  E_save(clock) = Ek(clock)+Epg(clock)+Eps(clock);

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   ANIMATION   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if mod(clock,nskip)==0
    clf(fig)
        h1 = plot3(5,5,5,"O","Markersize",15,"Color","r");
        hold on
        h=plot3(X(:,1),X(:,2),X(:,3),"O","Markersize",3,"Color","b");
        hold off
        axis equal
        axis manual
        axis([0,n1,0,n2,-3,5]);
    drawnow
end
%}

    if video
            frame = getframe(gcf);
            writeVideo(writerObj,frame);
    end
end

figure(2)
  plot(t_save',E_save);

if  video 
    close(writerObj);
end