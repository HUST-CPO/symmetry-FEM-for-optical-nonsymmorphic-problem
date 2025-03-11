function edgeIndex=findEdge(flag,mesh)
%n*2 1-网格 2-网格边

%不存在flag，直接返回
if isempty(flag)
    edgeIndex=[];
    return
end

%找到对应边界边
index=find(mesh.Edgesflag==flag(1));
if length(flag)>1
    for i=2:length(flag)
        index=[index;find(mesh.Edgesflag==flag(i))];
    end
end
index=sort(index);
boundary=mesh.Edges(index,:);
boundary=sort(boundary,2);

%找到边界边的全局编码
[~,index]=ismember(boundary,mesh.Edge,'rows');

%找到边界边对应网格编码
[~,index]=ismember(index,mesh.EdgeOfElements,'rows');
%数据拆分
n=length(index);
edgeIndex=zeros(n,2);
edgeIndex(:,1)=fix((index-1)/3)+1;
edgeIndex(:,2)=index-(edgeIndex(:,1)-1)*3;


end