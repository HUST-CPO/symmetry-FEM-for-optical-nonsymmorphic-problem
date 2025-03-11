function out=BF_Triangle(type,order,num,u,v)
%三角网格基函数
out=zeros(3,1);

if type==0 %Ez
    if order==1
        if num==1
            out(3)=1. - u - v;
        elseif num==2
            out(3)=u;
        elseif num==3
            out(3)=v;
        else
            disp('超出基函数最大数目');
        end
    elseif order==2
        if num==1 %200
            out(3)=1.0 - 3.0 * u + 2.0 * u * u - 3.0 * v + 4.0 * u * v + 2.0 * v * v;
        elseif num==2 %020
            out(3)=-u + 2.0 * u * u;
        elseif num==3 %002
            out(3)=-v + 2.0 * v * v;
        elseif num==4 %110
            out(3)=4.0 * u - 4.0 * u * u - 4.0 * u * v;
        elseif num==5 %101
            out(3)=4.0 * v - 4.0 * u * v - 4.0 * v * v;
        elseif num==6 %011
            out(3)=4.0 * u * v;
        else
            disp('超出基函数最大数目');
        end
    else
        disp('基函数阶数错误');
    end
elseif type==1 %gradEz
    if order==1
        if num==1
            out(1)=-1;
            out(2)=-1;
        elseif num==2
            out(1)=1;
            out(2)=0;
        elseif num==3
            out(1)=0;
            out(2)=1;
        else
            disp('超出基函数最大数目');
        end
    elseif order==2
        if num==1 %200
            out(1)=-3. + 4. * u + 4. * v;
            out(2)=-3. + 4. * u + 4. * v;
        elseif num==2 %020
            out(1)=-1. + 4. * u;
            out(2)=0;
        elseif num==3 %002
            out(1)=0;
            out(2)=-1. + 4. * v;
        elseif num==4 %110
            out(1)=4. - 8. * u - 4. * v;
            out(2)=-4. * u;
        elseif num==5 %101
            out(1)=-4. * v;
            out(2)=4. - 4. * u - 8. * v;
        elseif num==6 %011
            out(1)=4. * v;
            out(2)=4. * u;
        else
            disp('超出基函数最大数目');
        end
    else
        disp('基函数阶数错误');
    end
elseif type==2 %curlEz
    if order==1
        if num==1
            out(1)=-1;
            out(2)=1;
        elseif num==2
            out(1)=0;
            out(2)=-1;
        elseif num==3
            out(1)=1;
            out(2)=0;
        else
            disp('超出基函数最大数目');
        end
    elseif order==2
        if num==1 %200
            out(1)=-3. + 4. * u + 4. * v;
            out(2)=3. - 4. * u - 4. * v;
        elseif num==2 %020
            out(1)=0;
            out(2)=1. - 4. * u;
        elseif num==3 %002
            out(1)=-1. + 4. * v;
            out(2)=0;
        elseif num==4 %110
            out(1)=-4. * u;
            out(2)= -4. + 8. * u + 4. * v;
        elseif num==5 %101
            out(1)=4. - 4. * u - 8. * v;
            out(2)=4. * v;
        elseif num==6 %011
            out(1)=4. * u;
            out(2)=-4. * v;
        else
            disp('超出基函数最大数目');
        end
    else
        disp('基函数阶数错误');
    end
elseif type==3 %Et
    if order==1
        if num==1
            out(1)=1. - v;
            out(2)=u;
        elseif num==2
            out(1)=v;
            out(2)=-u + 1.;
        elseif num==3
            out(1)=-v;
            out(2)=u;
        else
            disp('超出基函数最大数目');
        end
    elseif order==2
        if num==1 %3-120 1-2
            out(1)=-1. + 3. * u + v - 3. * u * v;
            out(2)=-u + 3. * u * u;
        elseif num==2 %3-210  1-2
            out(1)=2. - 3. * u - 5. * v + 3. * u * v + 3. * v * v;
            out(2)=2. * u - 3. * u * u - 3. * u * v;
        elseif num==3 %2-102  1-3
            out(1)= -v + 3. * v * v;
            out(2)=-1. + u + 3. * v - 3. * u * v;
        elseif num==4 %2-201 1-3
            out(1)=2. * v - 3. * u * v - 3. * v * v;
            out(2)=2. - 5. * u + 3. * u * u - 3. * v + 3. * u * v;
        elseif num==5 %1-012  2-3
            out(1)=v - 3. * v * v;
            out(2)=-u + 3. * u * v;
        elseif num==6 %1-021  2-3
            out(1)=v - 3. * u * v;
            out(2)=-u + 3. * u * u;
        elseif num==7 %3-111 1-2
            out(1)=1.5 * (3. * v - 3. * v * v);
            out(2)=1.5 * (3. * u * v);
        elseif num==8 %2-111 1-3
            out(1)=1.5 * (3. * u * v);
            out(2)=1.5 * (3. * u - 3. * u * u);
        else
            disp('超出基函数最大数目');
        end
    else
        disp('基函数阶数错误');
    end
elseif type==4 %curlEt
    if order==1
        if num==1
            out(3)=2;
        elseif num==2
            out(3)=-2;
        elseif num==3
            out(3)=2;
        else
            disp('超出基函数最大数目');
        end
    elseif order==2
        if num==1 %3-120 1-2
            out(3)=-2. + 9. * u;
        elseif num==2 %3-210  1-2
            out(3)= 7. - 9. * u - 9. * v;
        elseif num==3 %2-102  1-3
            out(3)= 2 - 9. * v;
        elseif num==4 %2-201 1-3
            out(3)= -7. + 9. * u + 9. * v;
        elseif num==5 %1-012  2-3
            out(3)=-2. + 9. * v;
        elseif num==6 %1-021  2-3
            out(3)=-2. + 9. * u;
        elseif num==7 %3-111 1-2
            out(3)=1.5 * (-3. + 9. * v);
        elseif num==8 %2-111 1-3
            out(3)=1.5 * (3. - 9. * u);
        else
            disp('超出基函数最大数目');
        end
    else
        disp('基函数阶数错误');
    end
else
    disp('基函数类型错误');
end