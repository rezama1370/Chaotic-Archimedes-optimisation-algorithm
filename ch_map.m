function chmap=ch_map(chnum,chstart,Maxit)


chmap(1)=chstart;
%Chebyshev map
if chnum == 1
for i=1:Maxit
    chmap(i+1)=cos(i*acos(chmap(i)));
%     % G(i)=((chmap(i)+1)*Value)/2;
end
elseif chnum == 2
%Circle map
a=0.5;
b=0.2;
for i=1:Maxit
    chmap(i+1)=mod(chmap(i)+b-(a/(2*pi))*sin(2*pi*chmap(i)),1);
%     % G(i)=chmap(i)*Value;
end


elseif chnum == 3
%Iterative map
a=0.7;
for i=1:Maxit
    chmap(i+1)=sin((a*pi)/chmap(i));
    % G(i)=((chmap(i)+1)*Value)/2;
end



elseif chnum == 4
%Gauss/mouse map
for i=1:Maxit
    if chmap(i)==0
        chmap(i+1)=0;
    else
        chmap(i+1)=mod(1/chmap(i),1);
    end
%     % G(i)=chmap(i)*Value;
end



elseif chnum == 5
%Logistic map
a=4;
for i=1:Maxit
    chmap(i+1)=a*chmap(i)*(1-chmap(i));
    % G(i)=chmap(i)*Value;
end



elseif chnum == 6
%Piecewise map
P=0.4;
for i=1:Maxit
    if chmap(i)>=0 && chmap(i)<P
        chmap(i+1)=chmap(i)/P;
    end
    if chmap(i)>=P && chmap(i)<0.5
        chmap(i+1)=(chmap(i)-P)/(0.5-P);
    end
    if chmap(i)>=0.5 && chmap(i)<1-P
        chmap(i+1)=(1-P-chmap(i))/(0.5-P);
    end
    if chmap(i)>=1-P && chmap(i)<1
        chmap(i+1)=(1-chmap(i))/P;
    end    
    % G(i)=chmap(i)*Value;
end



elseif chnum == 7
 %Sinusoidal map
 for i=1:Maxit
     chmap(i+1) = 2.3*chmap(i)^2*sin(pi*chmap(i));
 end

 

elseif chnum == 8
 %Singer map 
 u=1.07;
 for i=1:Maxit
     chmap(i+1) = u*(7.86*chmap(i)-23.31*(chmap(i)^2)+28.75*(chmap(i)^3)-13.302875*(chmap(i)^4));
     % G(i)=(chmap(i))*Value;
 end
 
 

 
elseif chnum == 9
 %Tent map
%  chmap(1)=0.6;
 for i=1:Maxit
     if chmap(i)<0.7
         chmap(i+1)=chmap(i)/0.7;
     end
     if chmap(i)>=0.7
         chmap(i+1)=(10/3)*(1-chmap(i));
     end
 end

 
elseif chnum == 10
%Sine map
for i=1:Maxit
     chmap(i+1) = sin(pi*chmap(i));
     % G(i)=(chmap(i))*Value;
end 
 
 
 
end

end