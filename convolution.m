function [h]=convolution(f,g)
    [n, m, c] = size(f);
    [x,y]=size(g);
    h=zeros(n,m,c);
    padx=floor(x/2);
    pady=floor(y/2);
    paddedf=padarray(f,[padx,pady],'replicate');

    for channel=1:c
        imgChannel = paddedf(:,:,channel);
        imshow(imgChannel)
        for i=2:n-1
            for j=2:m-1
                tmp=imgChannel(i:i+x-1,j:j+y-1);
                h(i,j)=sum(sum(tmp.*g));
            end
        end
    end
    
end

% img = imread("img\camera.bmp");
% g=[-1 -1 -1; -1 8 -1; -1 -1 -1];
% h=uint8(convolution(double(img),double(g)));
% 
% 
% imshow(h)
% hm=uint8(conv2(double(img), double(g)));
% figure,imshow(hm)