function [h]=convolution(f,g)
    [n, m, c] = size(f);
    h=zeros(n,m,c);
    for channel=1:c
        imgChannel = f(:,:,channel);
        imshow(imgChannel)
        for i=2:n-1
            for j=2:m-1
                h(i,j)= imgChannel(i-1,j-1)*g(1,1)+imgChannel(i-1,j)*g(1,2)+imgChannel(i-1,j+1)*g(1,3)+imgChannel(i,j-1)*g(2,1)+imgChannel(i,j)*g(2,2)+imgChannel(i,j+1)*g(2,3)+imgChannel(i+1,j-1)*g(3,1)+imgChannel(i+1,j)*g(3,2)+imgChannel(i+1,j+1)*g(3,3);
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