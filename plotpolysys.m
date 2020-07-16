function plotpolysys(polysys, d)
%
%
M=getM(polysys,d);
imagesc(abs(M))
caxis([0,1]);
myColorMap = jet(2); %jet(256);
myColorMap(1,:) = 1;
myColorMap(end,:)=0;
colormap(myColorMap);
xlabel('Column');
ylabel('Row');
