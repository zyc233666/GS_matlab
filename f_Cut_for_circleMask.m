function E_out = f_Cut_for_circleMask(E_in,Center_x,Center_y,Cir_Radiu,A,Flag)
% �������ܣ����������������ֵΪ��λ����Բ�������ȡ
% E_in���������
% Center_x����ȡԲ���������������x
% Center_y����ȡԲ���������������x
% Cir_Radiu����ȡԲ������Ľ�ȡ�뾶����ֵ
% ע�⣺������������������ֵΪ��λ
% A����ȡ������ಿ�ָ�ֵ
% Flag:�жϽ�ȡ�ڲ������ⲿ
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

