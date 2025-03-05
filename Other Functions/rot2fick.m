function fick = rot2fick(rot)%Purpose Converts rotation vectors into Fick angles%Algorithm from a C-program written by Th. Haslwanter in 1990.% Essentially does rot2mat() then mat2fick(), but we only% compute the 3 elements of the mat that we need for Fick.x = rot(:,1);y = rot(:,2);z = rot(:,3);clear rot;scalar = sqrt(x.*x + y.*y + z.*z);alpha = 2*atan(scalar);% Avoid divide-by-zero errors.zr = scalar == 0;scalar(zr) = 1;nx = x ./ scalar; %normalizeny = y ./ scalar;nz = z ./ scalar;clear x y zsn = sin(alpha);cs = cos(alpha);clear alpha%nx2 = nx.*nx;%ny2 = ny.*ny;%nz2 = nz.*nz;%f11 = nx2 + cs .* (ny2 + nz2);%f22 = ny2 + cs .* (nz2 + nx2);%f33 = nz2 + cs .* (nx2 + ny2);%f12 = nx .* ny .* (1-cs) - nz .* sn;%f23 = ny .* nz .* (1-cs) - nx .* sn;f31 = nz .* nx .* (1-cs) - ny .* sn; %f13 = nx .* nz .* (1-cs) + ny .* sn;f21 = ny .* nx .* (1-cs) + nz .* sn;f32 = nz .* ny .* (1-cs) + nx .* sn;clear sn csclear nx ny nzvert = asin(f31);horz = asin(f21 ./ cos(vert));tors = asin(f32 ./ cos(vert));fick = [horz -vert tors] .* 180 / pi;