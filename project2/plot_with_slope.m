%write animition to mp4
video = true;
if(video)
    if exist('savings','var') == 1
        writerObj = VideoWriter('sav_lnr.mp4','MPEG-4');
    else
        writerObj = VideoWriter('ori_lnr.mp4','MPEG-4');
    end
    writerObj.FrameRate = 30;
    open(writerObj);
end


%histogram of ln(r) with time
% figure('Name','histo_lnr');
% for i = 1:t_max
%     edges = linspace(-7, 7, 56);
%     graph_r = log(r_total(:,i));
%     if exist('savings','var') == 1
%         graph_r(graph_r==-Inf)=-6;
%     end
%     histogram(graph_r,'BinEdges',edges);
%     drawnow;
%     frame = getframe(gcf);
%     writeVideo(writerObj,frame);
%     pause(0.02)
% end

%linear garph
figure('Name','r');
plot(r_total');

figure('Name','log_r');
r_tmp = log(r_total');
if exist('savings','var') == 1
    r_tmp(r_tmp==-inf)=-6;
end
plot(r_tmp);
xlabel('$r$','interpreter','latex','fontsize',16)
ylabel('$\ln{r}$','interpreter','latex','fontsize',16)

figure('Name','logsum_r');
plot(log(sum(r_total)));
xlabel('$sum(r_total)$','interpreter','latex','fontsize',16)
ylabel('$\log(sum(r_total))$','interpreter','latex','fontsize',16)


figure('Name','p');
plot(p_total');
xlabel('$t$','interpreter','latex','fontsize',16)
ylabel('$p_total$','interpreter','latex','fontsize',16)

figure('Name','E');
plot(excess_total');
xlabel('$t$','interpreter','latex','fontsize',16)
ylabel('$excess_total$','interpreter','latex','fontsize',16)

slopes = diff(log(sum(r_total)), 1);
figure('Name','slope');
plot(slopes');
xlabel('$t$','interpreter','latex','fontsize',16)
ylabel('$log(sum(r_{total}))$','interpreter','latex','fontsize',16)

acc= diff(slopes, 1)
figure('Name','acc');
plot(acc');
xlabel('$t$','interpreter','latex','fontsize',16)
ylabel('$acc))$','interpreter','latex','fontsize',16)
% if exist('savings','var') == 1
%     figure('Name','savings');
%     plot(savings');
% end

crp=zeros(2,t_max+1);
% these computations does not take savings into consideration, so we have
% to do it manually to make sure it is conserved.
for t=1:t_max+1
    crp(1,t)=(r_total(:,t+1)')*A*p_total(:,t) ;
    if exist('savings','var') == 1
        if t==1
            crp(1,t)=(r_total(:,t+1)')*A*p_total(:,t) + sum(savings(:,t));
            crp(2,t)=(r_total(:,t)')*B*p_total(:,t);
        else
            crp(2,t)=(r_total(:,t)')*B*p_total(:,t) - sum(savings(:,t)) + sum(savings(:,t-1));
        end
    else
        crp(2,t)=(r_total(:,t)')*B*p_total(:,t);
    end
end

if  video 
    close(writerObj);
end