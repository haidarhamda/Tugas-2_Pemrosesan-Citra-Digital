function [noise, result]=periodicNoiseRemoval(img, rotation, frequency, process2dims)  
    noise = insertNoise(img, rotation, frequency, process2dims);
    
    [x, y, z] = size(img);
    notch_filter = ones(x, y);
    radius = 5;
    center1 = [100, 150];
    center2 = [200, 250];

    [a, b] = meshgrid(1:x, 1:y);
    notch1 = sqrt((a - center1(2)).^2 + (b - center1(1)).^2) <= radius;
    notch2 = sqrt((a - center2(2)).^2 + (b - center2(1)).^2) <= radius;
    notch_filter = notch_filter .* ~notch1 .* ~notch2;

    result = zeros(size(img));

    for c = 1:z
        channel = img(:, :, c);
        channel_fft = fft2(channel);
        channel_fftshift = fftshift(channel_fft);

        filtered_fftshift = channel_fftshift .* notch_filter;

        filtered_fft = ifftshift(filtered_fftshift);
        filtered = ifft2(filtered_fft);
        result(:, :, c) = real(filtered);
    end
end

function noise=insertNoise(img, rotation, frequency, process2dims)
    [x, y, z] = size(img);

    if process2dims
        [a, b] =  meshgrid(1:x, 1:y);
        n = 0.2 * sin(deg2rad(rotation) * frequency * a / y) + 0.2 * sin(deg2rad(rotation) * frequency * b / x);
    else
        a = meshgrid(1:x, 1:y);
        n = 0.2 * sin(deg2rad(rotation) * frequency * a / y);
    end

    noise = img;
    for c = 1:z
        noise(:, :, c) = img(:, :, c) + n;
    end
end