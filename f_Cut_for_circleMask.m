function E_out = f_Cut_for_circleMask(E_in,Center_x,Center_y,Cir_Radiu,A,Flag)
% 函数功能：将输入矩阵按照像素值为单位进行圆形区域截取
% E_in：输入矩阵
% Center_x：截取圆形区域的中心坐标x
% Center_y：截取圆形区域的中心坐标x
% Cir_Radiu：截取圆形区域的截取半径像素值
% 注意：这三个参数均是像素值为单位
% A：截取后的其余部分赋值
% Flag:判断截取内部还是外部
[xnums,ynums] = size(E_in);
E_out = E_in;
if Flag
    TAm=1;
    TAn=A;
else
    TAm=A;
    TAn=1;
end
for nx = 1:xnums
    for ny = 1:ynums
        if abs((nx - Center_x) + 1i*(ny - Center_y)) < Cir_Radiu
            E_out(nx, ny) = TAm*E_in(nx, ny);
        else
            E_out(nx, ny) = TAn*E_in(nx, ny);
        end
    end
end
end

