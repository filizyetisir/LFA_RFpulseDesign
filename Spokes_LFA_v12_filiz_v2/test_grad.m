




function test_grad(fhandle,x)
% function test_gradient(fhandle,x)

[f df]=fhandle(x);

% STA
%eps=1e-6;
% LFA
eps = 1e-3;

nunk=size(df,1);
nconst=size(df,2);

for j=1:nconst
    fprintf('f%d  ',j);
end

for i=1:nunk
    x2=x;
    x2(i)=x2(i)+eps;
    [f2 df2]=fhandle(x2);
    fprintf('\n');
    for j=1:nconst
        fprintf('#%d: %e (%e) (%e)  ',i,(f2(j)-f(j))/eps,df(i,j),df2(i,j) );
        diff(i,j) = abs((f2(j)-f(j))/eps-df(i,j))./df(i,j);
    end    
end

ind = abs(diff)>1e9;
diff(ind) = 0;
max(diff(:))*100

fprintf('\n');



