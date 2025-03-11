function TriIndex=findTris(flag,mesh)

%不存在flag，直接返回
if isempty(flag)
    TriIndex=[];
    return
end

%找到对应边界面
index=find(mesh.Trisflag==flag(1));
if length(flag)>1
    for i=2:length(flag)
        index=[index;find(mesh.Trisflag==flag(i))];
    end
end
index=sort(index);
boundary=mesh.Tris(:,index)';
boundary=sort(boundary,2);
el2no=boundary';
n1=el2no([1 1 2],:);
n2=el2no([2 3 3],:);
el_ed2no_array=[n1(:) n2(:)];
[boundary,~]=unique(el_ed2no_array,'rows');
%找到边界边的全局编码
[~,index]=ismember(boundary,mesh.Edge,'rows');


%数据拆分
n=length(index);
TriIndex=zeros(n,1);
TriIndex(:,1)=index;


end