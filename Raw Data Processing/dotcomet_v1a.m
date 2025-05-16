figure
dotsize = 30;
dt = 0.00000000000000000001;
xlim([-250 250])
ylim([-200 200])
axis manual;
hold on
plot(dat.Y1(1),dat.Z1(1),'.','Color','#990000','MarkerSize',dotsize)
plot(dat.Y2(1),dat.Z2(1),'.','Color','#009900','MarkerSize',dotsize)
plot(dat.Y3(1),dat.Z3(1),'.','Color','#000099','MarkerSize',dotsize)
for i = 2:length(dat.Y1)
    plot(dat.Y1(i-1),dat.Z1(i-1),'.','Color','#FFCCCC','MarkerSize',dotsize)
    plot(dat.Y2(i-1),dat.Z2(i-1),'.','Color','#CCFFCC','MarkerSize',dotsize)
    plot(dat.Y3(i-1),dat.Z3(i-1),'.','Color','#CCCCFF','MarkerSize',dotsize)
    plot(dat.Y1(i),dat.Z1(i),'.','Color','#990000','MarkerSize',dotsize)
    plot(dat.Y2(i),dat.Z2(i),'.','Color','#009900','MarkerSize',dotsize)
    plot(dat.Y3(i),dat.Z3(i),'.','Color','#000099','MarkerSize',dotsize)
    pause(dt)
end