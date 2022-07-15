




function test_gradient_fn(fhandle,x)
% function test_gradient(fhandle,x)

[f feq df dfeq]=fhandle(x);

eps=1e-6;

nunk=size(df,1);
nconst=size(df,2);

% for j=1:nconst
%     fprintf('f%3.0d\t\t\t\t\t\t\t\t\t',j);
% end

for i=1:nunk
    x2=x;
    x2(i)=x2(i)+eps;
    [f2 f2eq df2 df2eq]=fhandle(x2);
    fprintf('\n');
    for j=1:nconst
            fprintf('#%2.0d: %3.3e (%3.4e) \t\t\t ',i,(f2(j)-f(j))/eps,df(i,j) );
        diff(i,j) =( ((f2(j)-f(j))/eps-df(i,j))./df(i,j) )*100;
    end
end

fprintf('\n Max diff: %f percent\n', max(diff(:)));



