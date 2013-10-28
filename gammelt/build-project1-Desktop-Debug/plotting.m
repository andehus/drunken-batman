fileID=fopen('/home/anders/Project1/build-project1-Desktop-Debug/pdata100.dat');
X = textscan(fileID, '%f %f %f');
fclose(fileID);
plot(X{1},X{2});