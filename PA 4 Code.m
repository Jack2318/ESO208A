clear;
close all;
clc
%Aviral Agarwal
%180167 - Tutorial Section - J6
%Computer Assignment 3
prompt = "Which method do you want to run first?\n1) Linear spline\n2)  Quadratic spline\n3)  Natural cubic spline\n4) Not-a-knot cubic spline\n5)  Periodic cubic spline\n6)Clamped cubic spine \nEnter list as [] consisting of space separated values for all the method that you want to run?: ";
methods = input(prompt);
file_name = "input_file"+".txt";
prompt = "Enter the absolute path of the folder in which you want the output files ?\n(See the pdf for examples):-";
out = input(prompt);
while isempty(out)
    prompt = "Programme will not proceed until you enter a path:- ";
    out = input(prompt);
end
output_file_name = out + "\output_file" + ".txt";
fileoID = fopen(output_file_name,'wt');
fclose(fileoID);
fid = fopen(file_name);
tline = fgetl(fid);
k = 0;
x_train = [];
y_train = [];
x_pred = [];
while ischar(tline)
    k = k+1;
    if(k==1)
        n = str2num(tline);
        x_train = zeros(1,n);
        y_train = zeros(1,n);
    elseif (k==2)
        m = str2num(tline);
    elseif (k<=n+2)
        curr = str2num(tline);
        x_train(k-2) = curr(1);
        y_train(k-2) = curr(2);
    elseif(k<=(m+n+2))
        x_pred(k-n-2) = str2num(tline);
    else
        if(any(methods(:)==6) && k==m+n+3)
            curr = str2num(tline);
            s_0 = curr(1);
            s_n = curr(2);
        end
    end
    tline = fgetl(fid);
end
fclose(fid);
for z = 1:length(methods)
    method_no = methods(z);
    if(method_no == 1)
        %Linear Spline
        y_pred = zeros(1,m);
        for j = 1:m
            for i = 2:n
                if(x_pred(j) <= x_train(i) && x_pred(j) >= x_train(i-1))
                    y_pred(j) = (((x_pred(j)-x_train(i-1))*y_train(i))/(x_train(i)-x_train(i-1))) + (((x_pred(j)-x_train(i))* y_train(i-1))/(x_train(i-1)-x_train(i)));
                end
            end
        end
        %Printing Output In file
        fileoID = fopen(output_file_name,'a');
        fprintf(fileoID,'Linear Spline\n\n');
        for i = 1:m
            fprintf(fileoID,'%f %f\n',x_pred(i),y_pred(i));
        end
        fprintf(fileoID,'\n\n');
        fclose(fileoID);
        hold on;
        grid on;
        x = [];
        y = [];
        for i = 2:n
            x1 = x_train(i-1):0.001:x_train(i);
            y1 = ((x1 - x_train(i-1))./(x_train(i) - x_train(i-1))).*y_train(i) + ((x1 - x_train(i))./(x_train(i-1) - x_train(i))).*y_train(i-1);
            x = [x, x1];
            y = [y, y1];
        end
        plot(x,y,'-r*','MarkerSize',5,'MarkerIndices',1:20:length(y),'DisplayName','Linear spline');
        legend
    elseif(method_no == 2)
        %Quadratic Spline
        A = zeros(3*(n-1),3*(n-1));
        i = 1;
        for l = 1:(n-1)
            j = 3*(l-1)+1;
            A(i,j) = x_train(l)^2;
            A(i,j+1) = x_train(l);
            A(i,j+2) = 1;
            A(i+1,j) = x_train(l+1)^2;
            A(i+1,j+1) = x_train(l+1);
            A(i+1,j+2) = 1;
            i = i+2;
        end
        j = 1;
        l = 2;
        for i = 2*(n-1)+1:3*(n-1)-1
            A(i,j) = 2*x_train(l);
            A(i,j+1) = 1;
            A(i,j+2) = 0 ;
            A(i,j+3) = -2*x_train(l);
            A(i,j+4) = -1;
            A(i,j+5) = 0;
            l = l+1;
            j = j+3;
        end
        A(3*(n-1),1) = 1;
        B = zeros(3*(n-1),1);
        B(1) = y_train(1);
        j = 2;
        for i = 2:n-1
            B(j) = y_train(i);
            B(j+1) = y_train(i);
            j = j+2;
        end
        B(2*(n-1)) = y_train(n);
        B(2*(n-1)+1:3*(n-1)) = 0;
        sigma_curr = linsolve(A,B);
        l = 1;
        for j = 1:m
            for i = 2:n
                if(x_pred(j) <= x_train(i) && x_pred(j) >= x_train(i-1))
                    y_pred(j) = sigma_curr(l)*x_pred(j)^2 + sigma_curr(l+1)*x_pred(j) + sigma_curr(l+2);
                    l = l+3;
                end
            end
        end
        fileoID = fopen(output_file_name,'a');
        fprintf(fileoID,'Quadratic Spline\n\n');
        for i = 1:m
            fprintf(fileoID,'%f %f\n',x_pred(i),y_pred(i));
        end
        fprintf(fileoID,'\n\n');
        fclose(fileoID);
        hold on;
        grid on;
        x =[];
        y = [];
        l = 1;
        for i = 2:n
            x1 = x_train(i-1):0.001:x_train(i);
            y1 = sigma_curr(l).*x1.^2 + sigma_curr(l+1).*x1 + sigma_curr(l+2);
            l = l+3;
            x = [x,x1];
            y = [y,y1];
        end 
        pb = plot(x,y,'-k','DisplayName','Quadratic spline');
        legend
    elseif(method_no == 3)
        %Natural Cubic Spline
        h = zeros(n,1);
        g = zeros(n,1);
        for i = 2:n
            h(i) = x_train(i) - x_train(i-1);
            g(i) = (y_train(i) - y_train(i-1))/h(i);
        end
        sigma_curr = zeros(n,1);
        B = zeros(n,1);
        for i = 2:n-1
            B(i) = 6*(g(i+1)-g(i));
        end
        A = zeros(n,n);
        for i=2:n-1
            for j = 2:n-1
                if(i==j)
                    A(i,j) = 2*(h(i)+h(i+1));
                elseif(i==j-1)
                    A(i,j) = h(i+1);
                elseif(i==j+1)
                    A(i,j) = h(i);
                end
            end
        end
        sigma_curr(2:n-1) = linsolve(A(2:n-1,2:n-1),B(2:n-1));
        y_pred = zeros(1,m);
        for j = 1:m
            for i = 2:n
                if(x_pred(j) <= x_train(i) && x_pred(j) >= x_train(i-1))
                    A_new = sigma_curr(i)/(6*h(i));
                    B_new = sigma_curr(i-1)/(6*h(i));
                    C = (y_train(i)/h(i)) - ((sigma_curr(i)/6)*h(i));
                    D = (y_train(i-1)/h(i)) - ((sigma_curr(i-1)/6)*h(i));
                    y_pred(j) = A_new*power(x_pred(j)-x_train(i-1),3) - B_new*power(x_pred(j)-x_train(i),3) + C*(x_pred(j)-x_train(i-1)) - D*(x_pred(j)-x_train(i)); 
                end
            end
        end
        %Printing Output In file
        fileoID = fopen(output_file_name,'a');
        fprintf(fileoID,'Natural Cubic Spline\n\n');
        for i = 1:m
            fprintf(fileoID,'%f %f\n',x_pred(i),y_pred(i));
        end
        fprintf(fileoID,'\n\n');
        fclose(fileoID);
        hold on;
        grid on;
        x = [];
        y = [];
        for i = 1:n-1
            x1 = x_train(i):0.001:x_train(i+1);
            A_new = sigma_curr(i+1)/(6*h(i+1));
            B_new = sigma_curr(i)/(6*h(i+1));
            C = (y_train(i+1)/h(i+1)) - ((sigma_curr(i+1)/6)*h(i+1));
            D = (y_train(i)/h(i+1)) - ((sigma_curr(i)/6)*h(i+1));
            y1 = A_new.*((x1-x_train(i)).^3)-B_new.*((x1-x_train(i+1)).^3)+C.*(x1-x_train(i))-D.*(x1-x_train(i+1));
            x = [x,x1];
            y = [y,y1];
        end
        pf = plot(x,y,'-g','MarkerSize',5,'MarkerIndices',1:20:length(y),'DisplayName','Natural Cubic Spline');
        legend
    elseif(method_no == 4)
        %Not a Knot Cubic Spline
        h = zeros(n,1);
        g = zeros(n,1);
        for i = 2:n
            h(i) = x_train(i) - x_train(i-1);
            g(i) = (y_train(i) - y_train(i-1))/h(i);
        end
        sigma_curr = zeros(n,1);
        B = zeros(n,1);
        for i = 2:n-1
            B(i) = 6*(g(i+1)-g(i));
        end
        A = zeros(n,n);
        A(1,1) = h(3);
        A(1,2) = -(h(2)+h(3));
        A(1,3) = h(2);
        for i=2:n-1
            for j = 1:n
                if(i==j)
                    A(i,j) = 2*(h(i)+h(i+1));
                elseif(i==j-1)
                    A(i,j) = h(i+1);
                elseif(i==j+1)
                    A(i,j) = h(i);
                end
            end
        end
        A(n,n-2) = h(n);
        A(n,n-1) = -(h(n-1)+h(n));
        A(n,n) = h(n-1);
        sigma_curr = linsolve(A,B);
        y_pred = zeros(1,m);
        for j = 1:m
            for i = 2:n
                if(x_pred(j) <= x_train(i) && x_pred(j) >= x_train(i-1))
                    A_new = sigma_curr(i)/(6*h(i));
                    B_new = sigma_curr(i-1)/(6*h(i));
                    C = (y_train(i)/h(i)) - ((sigma_curr(i)/6)*h(i));
                    D = (y_train(i-1)/h(i)) - ((sigma_curr(i-1)/6)*h(i));
                    y_pred(j) = A_new*power(x_pred(j)-x_train(i-1),3) - B_new*power(x_pred(j)-x_train(i),3) + C*(x_pred(j)-x_train(i-1)) - D*(x_pred(j)-x_train(i)); 
                end
            end
        end
        
        %Printing Output In file
        fileoID = fopen(output_file_name,'a');
        fprintf(fileoID,'Not-a-Knot Cubic Spline\n\n');
        for i = 1:m
            fprintf(fileoID,'%f %f\n',x_pred(i),y_pred(i));
        end
        fprintf(fileoID,'\n\n');
        fclose(fileoID);
        hold on;
        grid on;
        x = [];
        y = [];
        for i = 1:n-1
            x1 = x_train(i):0.001:x_train(i+1);
            A_new = sigma_curr(i+1)/(6*h(i+1));
            B_new = sigma_curr(i)/(6*h(i+1));
            C = (y_train(i+1)/h(i+1)) - ((sigma_curr(i+1)/6)*h(i+1));
            D = (y_train(i)/h(i+1)) - ((sigma_curr(i)/6)*h(i+1));
            y1 = A_new.*((x1-x_train(i)).^3)-B_new.*((x1-x_train(i+1)).^3)+C.*(x1-x_train(i))-D.*(x1-x_train(i+1));
            x = [x,x1];
            y = [y,y1];
        end
        pf = plot(x,y,'-b','MarkerSize',5,'MarkerIndices',1:20:length(y),'DisplayName','Not-a-Knot Cubic Spline');
        legend
    elseif(method_no == 5)
        %Periodic Cubic Spline
        h = zeros(n,1);
        g = zeros(n,1);
        for i = 2:n
            h(i) = x_train(i) - x_train(i-1);
            g(i) = (y_train(i) - y_train(i-1))/h(i);
        end
        sigma_curr = zeros(n,1);
        B = zeros(n,1);
        B(1) = -6*((y_train(n)-y_train(n-1))/h(n)) + (6*(y_train(2) - (y_train(1)))/h(2));
        for i = 2:n-1
            B(i) = 6*(g(i+1)-g(i));
        end
        A = zeros(n,n);
        A(1,1) = 2*h(2);
        A(1,2) = h(2);
        A(1,n-1) = h(n);
        A(1,n) = 2*h(n,1);
        for i=2:n-1
            for j = 1:n
                if(i==j)
                    A(i,j) = 2*(h(i)+h(i+1));
                elseif(i==j-1)
                    A(i,j) = h(i+1);
                elseif(i==j+1)
                    A(i,j) = h(i);
                end
            end
        end
        A(n,1) = 1;
        A(n,n) = -1;   
        sigma_curr = linsolve(A,B);
        y_pred = zeros(1,m);
        for j = 1:m
            for i = 2:n
                if(x_pred(j) <= x_train(i) && x_pred(j) >= x_train(i-1))
                    A_new = sigma_curr(i)/(6*h(i));
                    B_new = sigma_curr(i-1)/(6*h(i));
                    C = (y_train(i)/h(i)) - ((sigma_curr(i)/6)*h(i));
                    D = (y_train(i-1)/h(i)) - ((sigma_curr(i-1)/6)*h(i));
                    y_pred(j) = A_new*power(x_pred(j)-x_train(i-1),3) - B_new*power(x_pred(j)-x_train(i),3) + C*(x_pred(j)-x_train(i-1)) - D*(x_pred(j)-x_train(i)); 
                end
            end
        end
        %Printing Output in File
        fileoID = fopen(output_file_name,'a');
        fprintf(fileoID,'Periodic Cubic Spline\n\n');
        for i = 1:m
            fprintf(fileoID,'%f %f\n',x_pred(i),y_pred(i));
        end
        fprintf(fileoID,'\n\n');
        fclose(fileoID);
        hold on;
        grid on;
        x = [];
        y = [];
        for i = 1:n-1
            x1 = x_train(i):0.001:x_train(i+1);
            A_new = sigma_curr(i+1)/(6*h(i+1));
            B_new = sigma_curr(i)/(6*h(i+1));
            C = (y_train(i+1)/h(i+1)) - ((sigma_curr(i+1)/6)*h(i+1));
            D = (y_train(i)/h(i+1)) - ((sigma_curr(i)/6)*h(i+1));
            y1 = A_new.*((x1-x_train(i)).^3)-B_new.*((x1-x_train(i+1)).^3)+C.*(x1-x_train(i))-D.*(x1-x_train(i+1));
            x = [x,x1];
            y = [y,y1];
        end
        pf = plot(x,y,'-c','MarkerSize',5,'MarkerIndices',1:20:length(y),'DisplayName','Periodic Cubic Spline');
        legend
    elseif(method_no == 6)
        %Clamped Cubic Spline
        h = zeros(n,1);
        g = zeros(n,1);
        for i = 2:n
            h(i) = x_train(i) - x_train(i-1);
            g(i) = (y_train(i) - y_train(i-1))/h(i);
        end
        sigma_curr = zeros(n,1);
        B = zeros(n,1);
        B(1) = 6*(((y_train(2)-y_train(1))/h(2)) - s_0);
        B(n) = 6*(((y_train(n-1)-y_train(n))/h(n)) + s_n);
        for i = 2:n-1
            B(i) = 6*(g(i+1)-g(i));
        end
        A = zeros(n,n);
        A(1,1) = 2*h(2);
        A(1,2) = h(2);
        for i=2:n-1
            for j = 1:n
                if(i==j)
                    A(i,j) = 2*(h(i)+h(i+1));
                elseif(i==j-1)
                    A(i,j) = h(i+1);
                elseif(i==j+1)
                    A(i,j) = h(i);
                end
            end
        end
        A(n,n) = 2*h(n);
        A(n,n-1) = h(n);
        sigma_curr = linsolve(A,B);
        y_pred = zeros(1,m);
        for j = 1:m
            for i = 2:n
                if(x_pred(j) <= x_train(i) && x_pred(j) >= x_train(i-1))
                    A_new = sigma_curr(i)/(6*h(i));
                    B_new = sigma_curr(i-1)/(6*h(i));
                    C = (y_train(i)/h(i)) - ((sigma_curr(i)/6)*h(i));
                    D = (y_train(i-1)/h(i)) - ((sigma_curr(i-1)/6)*h(i));
                    y_pred(j) = A_new*power(x_pred(j)-x_train(i-1),3) - B_new*power(x_pred(j)-x_train(i),3) + C*(x_pred(j)-x_train(i-1)) - D*(x_pred(j)-x_train(i)); 
                end
            end
        end
        %Printing Output in File
        fileoID = fopen(output_file_name,'a');
        fprintf(fileoID,'Clamped Cubic Spline\n\n');
        for i = 1:m
            fprintf(fileoID,'%f %f\n',x_pred(i),y_pred(i));
        end
        fprintf(fileoID,'\n\n');
        fclose(fileoID);
        hold on;
        grid on;
        x = [];
        y = [];
        for i = 1:n-1
            x1 = x_train(i):0.001:x_train(i+1);
            A_new = sigma_curr(i+1)/(6*h(i+1));
            B_new = sigma_curr(i)/(6*h(i+1));
            C = (y_train(i+1)/h(i+1)) - ((sigma_curr(i+1)/6)*h(i+1));
            D = (y_train(i)/h(i+1)) - ((sigma_curr(i)/6)*h(i+1));
            y1 = A_new.*((x1-x_train(i)).^3)-B_new.*((x1-x_train(i+1)).^3)+C.*(x1-x_train(i))-D.*(x1-x_train(i+1));
            x = [x,x1];
            y = [y,y1];
        end
        pf = plot(x,y,'-mo','MarkerSize',5,'MarkerIndices',1:20:length(y),'DisplayName','Clamped Cubic Spline');
        legend
    else
        dips("Method number out of range");
    end
end
clear