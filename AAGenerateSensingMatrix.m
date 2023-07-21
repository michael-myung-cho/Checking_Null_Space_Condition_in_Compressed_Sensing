function [rss,pnt]=AAGenerateSensingMatrix(imenu,n,m);
% This code was extracted from that released by Juditsky & Nemirovski
% Generates random CS matrices knwon to satisfy the RIP property
% 1    strF='Fourier';
% 3    str1='+1/-1';
% 4    strG='Gauss';
% 9    strV='DeVore';
rss.A=[];
rss.m=0;
rss.n=0;
Scaling=0;
if imenu==9,
    clear rss;
    disp('Primes: 3 5 7 11 13 17 19 23 29 31 17');
    p=input('m=p^2 for a prime p. p > ');
    rss.prime=p;
    m=p^2;
    rss.m=m;
    r=input('p^2 < n <= p^(r+1) for 1 < r < p. r > ');
    if r<=1,
        r=2;
    end;
    if r>=p,
        r=p-1;
    end;
    n=input(sprintf('%1d <= n <= %d. n > ',p^2+1,p^(r+1)));
    if n<p^2+1,
        n=p^2+1;
    end;
    if n>p^(r+1),
        n=p^(r+1);
    end;
    rss.r=r;
    rss.n=n;
    %disp(sprintf('%1dx%1d matrix',m,n));
    [rss.A,rss.RIP]=DeVore(p,r,n);
    disp(sprintf('RIP sparsity: %1d',rss.RIP.s));
    if Scaling,
        for i=1:n,
            rss.A(:,i)=rss.A(:,i)/norm(rss.A(:,i));
        end;
    end;
    pnt=9;
    rss.pnt=pnt;
end;
if imenu==1,
    clear rss;
    hm=floor(m/2);
    hn=floor(n/2);
    hm=min(hm,hn-2);
    [tmp,ind]=sort(randn(hn-2,1));
    freq=sort(ind(1:hm));
    rss.freq=freq;
    rss.A=zeros(m,n);
    rss.m=m;
    rss.n=n;
    for i=1:hm,
        if freq(i)~=0,
            rss.A(2*i-1,:)=cos(2*pi*freq(i)*[0:1:n-1]/n)*sqrt(0.5);
            rss.A(2*i,:)=sin(2*pi*freq(i)*[0:1:n-1]/n)*sqrt(0.5);
        else
            rss.A(2*i-1,:)=1;
        end;
    end;
    [U,D,V]=svd(rss.A);
    dm=max(diag(D));
    mm=m;
    for i=1:m,
        if D(i,i)<=1.e-8*dm,
            mm=i-1;
            break;
        end;
    end;
    if mm<m,
        D=D(1:mm,:);
        rss.A=D*V';
        rss.m=mm;
        m=mm;
        %disp(sprintf('%1dx%1d matrix',m,n));
    end;
    if Scaling,
        for i=1:n,
            rss.A(:,i)=rss.A(:,i)/norm(rss.A(:,i));
        end;
    end;
    pnt=1;
    rss.pnt=pnt;
end;
if imenu==3,
    clear rss;
    rss.m=m;
    rss.n=n;
    rss.A=sign(randn(m,n));
    [U,D,V]=svd(rss.A);
    dm=max(diag(D));
    mm=m;
    for i=1:m,
        if D(i,i)<=1.e-8*dm,
            mm=i-1;
            break;
        end;
    end;
    if mm<m,
        D=D(1:mm,:);
        rss.A=D*V';
        rss.m=mm;
        m=mm;
        %disp(sprintf('%1dx%1d matrix',m,n));
    end;
    pnt=4;
    rss.pnt=pnt;
end;
if imenu==4,
    clear rss;
    rss.m=m;
    rss.n=n;
    rss.A=randn(m,n);
    [U,D,V]=svd(rss.A);
    dm=max(diag(D));
    mm=m;
    for i=1:m,
        if D(i,i)<=1.e-8*dm,
            mm=i-1;
            break;
        end;
    end;
    if mm<m,
        D=D(1:mm,:);
        rss.A=D*V';
        rss.m=mm;
        m=mm;
        %disp(sprintf('%1dx%1d matrix',m,n));
    end;
    if Scaling,
        for i=1:n,
            rss.A(:,i)=rss.A(:,i)/norm(rss.A(:,i));
        end;
    end;
    pnt=5;
    rss.pnt=pnt;
end;
