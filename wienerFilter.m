function [blurred, result]=wienerFilter(img, length, theta, nsr)
    if length == 0
        blurred = img;
        result = img;
        return;
    end

    [blurred, kernel] = motionBlur(img, length, theta, nsr);

    [x, y, z] = size(img);
    result = zeros(size(img));

    for i = 1:z
        blurred_channel = blurred(:, :, i);
        blurred_fft = fft2(blurred_channel);
        kernel_padded = padarray(kernel, [x - length, y - length], 'post');
        kernel_fft = fft2(kernel_padded);

        H = conj(kernel_fft) ./ (abs(kernel_fft).^2 + nsr / 10);

        result_channel_fft = blurred_fft .* H;
        result_channel = ifft2(result_channel_fft);
        result(:, :, i) = real(result_channel);
    end
end

function [result, kernel]=motionBlur(img, length, theta, nsr)
    kernel = zeros(length, length);
    center = length / 2;
    for i = 1:length
        x = round(center + (i - center) * cosd(theta));
        y = round(center + (i - center) * sind(theta));
        kernel(y, x) = 1;
    end
    kernel = kernel + nsr * randn(size(kernel));
    kernel = kernel / sum(kernel(:));

    [x, y, z] = size(img);
    result = zeros(size(img));

    for i = 1:z
        channel = img(:, :, i);
        channel_fft = fft2(channel);
        kernel_padded = padarray(kernel, [x - length, y - length], 'post');
        kernel_fft = fft2(kernel_padded);
        blur = channel_fft .* kernel_fft;
        blur = ifft2(blur);
        result(:, :, i) = real(blur);
    end
end