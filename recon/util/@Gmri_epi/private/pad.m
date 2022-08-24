function y = pad(x,left,right,above,below,value)

y = value.*ones(size(x,1)+above+below,size(x,2)+left+right);
y(above+1:above+size(x,1),left+1:left+size(x,2))=x;
