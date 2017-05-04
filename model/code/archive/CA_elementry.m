clear variables;
clc;

% 7     6   5   4   3   2   1   0
% 128   64  32  16  8   4   2   1

len=501; 
GRID=zeros(len,len);
GRID(1,round(len/2))=1;

up=[2:len 1]; 
down=[len 1:len-1];
colormap(gray(2));
for i=1:len-1
    
    neighbours= 4*GRID(down,down) + 2*GRID(down,:) + GRID(down,up);
    GRID = GRID==1 | neighbours==6 | neighbours==5 | neighbours==4 | neighbours==3 | neighbours==2 | neighbours==1;

    image(GRID*2); 
    pause(0.01);
end