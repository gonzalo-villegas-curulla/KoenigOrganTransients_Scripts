buff = median(1./MX(:,:,21),1,'omitnan');
buff = 1./MX(:,:,21);
figure(1);clf;
plot(A, buff,'d');ax=gca;ax.YLim(1)=0;
%
figure(2);clf;
plot(B, buff,'d');ax=gca;ax.YLim(1)=0;
%
figure(3);clf;
plot(C, buff,'d');
xlim([0 2.2e-3]);ax=gca;ax.YLim(1)=0;
%
figure(4);clf;
plot(D, buff,'d');
xlim([0 1.7e-3]);ax=gca;ax.YLim(1)=0;
%
figure(5);clf;
plot(sigMa, buff,'d');ax=gca;ax.YLim(1)=0;

%%
figure();
plot(1./MX(:,:,26), MX(:,:,21), 'o');
ax=gca;ax.YLim(1) = 0;

%%
figure(2);clf;
bof = 1.*MX(:,:,45);
plot(fax, bof, 'o');
grid on; ax=gca; ax.YLim(1) = 0;