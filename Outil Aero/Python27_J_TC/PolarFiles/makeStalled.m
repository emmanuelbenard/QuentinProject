lift = load('Stalled_lift.txt');
lift(1:14,2) = interp1(lift([1,14],1),lift([1,14],2),lift(1:14,1))
drag = load('Stalled_Drag.txt');
mom = load('Stalled_Mom.txt');
nb = max([size(lift,1),size(drag,1),size(mom,1)]);
data = zeros(nb,5);
data(:,1) = linspace(0,360,nb)';
data(:,2) = interp1(lift(1:4:end,1),lift(1:4:end,2),data(:,1),'spline');
data(:,3) = interp1(drag(1:5:end,1),drag(1:5:end,2),data(:,1),'spline');
data(:,4) = interp1(mom(1:3:end,1),mom(1:3:end,2),data(:,1),'spline');
figure, hold on, 
plot(data(:,1),data(:,2));
plot(data(:,1),data(:,3));
plot(data(:,1),data(:,4));
data(51:end,1) = data(51:end,1)-360;
[~,index] = sort(data(:,1),'ascend');
data = data(index,:);
figure, hold on, 
plot(data(:,1),data(:,2));
plot(data(:,1),data(:,3));
plot(data(:,1),data(:,4));

polar = data;
fileID = fopen('Stalled.txt','w');
for i = 1:size(polar,1)
    fprintf(fileID,'%9.1f %9.4f %9.4f %9.4f %9.4f \n',polar(i,1),polar(i,2),polar(i,3),polar(i,5),polar(i,4));
end
fclose(fileID);
