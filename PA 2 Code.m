clear;
close all;
clc
%Aviral Agarwal
%180167 - Tutorial Section - J6
%Computer Assignment 1
prompt = "Which method do you want to run first?\n1) Gauss elimination (GE; without pivoting)\n2) GE (with pivoting)\n3) GE (with scaling and pivoting)\n4) LU decomposition by using GE (without pivoting)\n5) LU decomposition by using GE (with pivoting)\n6) LU decomposition by using Crout method (without pivoting)\n7) Cholesky decomposition (for symmetric positive definite matrix)\nEnter value from 1 to 7: ";
method_no = input(prompt);
file_name = "input_file_method_"+string(method_no)+".txt";
prompt = "Enter the absolute path of the folder in which you want the output files ?\n(See the pdf for examples):-";
out = input(prompt);
while isempty(out)
    prompt = "Programme will not proceed until you enter a path:- ";
    out = input(prompt);
end
output_file_name = out + "\output_file_method_" + string(method_no) + ".txt";
fid = fopen(file_name);
tline = fgetl(fid);
n = str2num(tline);
m = n+1;
count = 0;
fclose(fid);
AUG = readmatrix(file_name);
if(method_no == 1)
    %Gaussian Elimination without Pivoting or Scaling
    for i = 1:n
        for j = i+1:n
            multiply_factor = AUG(j,i)/AUG(i,i);
            for k = 1:m
                AUG(j,k) = AUG(j,k) - (multiply_factor * AUG(i,k));
            end
        end
    end
    %Back substitution
    X = [];
    for i = n:-1:1
        sum = 0;
        for j = m-1:-1:m-(n-i)
            sum = sum + X(j)*AUG(i,j);
        end
        X(i) = (AUG(i,m) - sum)/AUG(i,m-n+i-1);
    end
    U = AUG(:,1:n);
    %Printing Output In file
    fileoID = fopen(output_file_name,'wt');
    fprintf(fileoID,'X\n\n');
    for i = 1:n
        fprintf(fileoID,'%f\n',X(i));
    end
    fprintf(fileoID,'\n');
    fprintf(fileoID,'U for GE with no pivot and no scaling\n\n');
    for i = 1:n
        for j = 1:n
            fprintf(fileoID,'%f ',U(i,j));
        end
        fprintf(fileoID,'\n');
    end
    fclose(fileoID);
elseif(method_no == 2)
    %Gaussian Elimination with Pivoting (Total) and no Scaling
    row_permutation = eye(n);
    column_permutation = eye(n);
    for i = 1:n
        %Pivoting
        max_element = 0;
        row = 0;
        column = 0;
        for k = i:n
            for j = i:m-1
                if abs(AUG(k,j)) > max_element
                    max_element = abs(AUG(k,j));
                    row = k;
                    column = j;
                end 
            end
        end
        %row pivoting
        AUG([i,row],:) = AUG([row,i],:);
        row_permutation([i,row],:) = row_permutation([row,i],:);
        %column pivoting
        AUG(:,[i,column]) = AUG(:,[column,i]);
        column_permutation(:,[i,column]) = column_permutation(:,[column,i]);

        %Applying Gaussian Elimination after pivoting
        for j = i+1:n
            multiply_factor = AUG(j,i)/AUG(i,i);
            for k = 1:m
                AUG(j,k) = AUG(j,k) - (multiply_factor * AUG(i,k));
            end
        end
    end

    %Back substitution
    X = [];
    for i = n:-1:1
        sum = 0;
        for j = m-1:-1:m-(n-i)
            sum = sum + X(j)*AUG(i,j);
        end
        X(i) = (AUG(i,m) - sum)/AUG(i,m-n+i-1);
    end
    U = AUG(:,1:n);
    %rearranging X because during column pivoting they got rotated
    X = X * (column_permutation^-1);
    %Printing Output In file
    fileoID = fopen(output_file_name,'wt');
    fprintf(fileoID,'X\n\n');
    for i = 1:n
        fprintf(fileoID,'%f\n',X(i));
    end
    fprintf(fileoID,'\n');
    fprintf(fileoID,'U for GE pivot and no scaling\n\n');
    for i = 1:n
        for j = 1:n
            fprintf(fileoID,'%f ',U(i,j));
        end
        fprintf(fileoID,'\n');
    end
    fprintf(fileoID,'\n');
    fprintf(fileoID,'Row Permutation Matrix for GE with pivot\n\n');
    for i = 1:n
        for j = 1:n
            fprintf(fileoID,'%f ',row_permutation(i,j));
        end
        fprintf(fileoID,'\n');
    end
    fprintf(fileoID,'\n');
    fprintf(fileoID,'Column Permutation Matrix for GE with pivot\n\n');
    for i = 1:n
        for j = 1:n
            fprintf(fileoID,'%f ',column_permutation(i,j));
        end
        fprintf(fileoID,'\n');
    end
    fclose(fileoID);
elseif(method_no == 3)
    %Gaussian Elimination with Pivoting (Total) and Scaling
    row_permutation = eye(n);
    column_permutation = eye(n);
    for i = 1:n
        %Scaling
        AUG_COPY = AUG;
        max_each_row = max(abs(AUG),[],2);
        %Scaling every row
        for k = i:n
            for j = i:m
                AUG_COPY(i,j) = AUG_COPY(i,j)/max_each_row(i);
            end
        end

        %Pivoting
        max_element = 0;
        row = 0;
        column = 0;
        for k = i:n
            for j = i:m-1
                if abs(AUG_COPY(k,j)) > max_element
                    max_element = abs(AUG_COPY(k,j));
                    row = k;
                    column = j;
                end 
            end
        end
        %row pivoting
        AUG([i,row],:) = AUG([row,i],:);
        row_permutation([i,row],:) = row_permutation([row,i],:);
        %column pivoting
        AUG(:,[i,column]) = AUG(:,[column,i]);
        column_permutation(:,[i,column]) = column_permutation(:,[column,i]);

        %Applying Gaussian Elimination after pivoting
        for j = i+1:n
            multiply_factor = AUG(j,i)/AUG(i,i);
            for k = 1:m
                AUG(j,k) = AUG(j,k) - (multiply_factor * AUG(i,k));
            end
        end
    end

    %Back substitution
    X = [];
    for i = n:-1:1
        sum = 0;
        for j = m-1:-1:m-(n-i)
            sum = sum + X(j)*AUG(i,j);
        end
        X(i) = (AUG(i,m) - sum)/AUG(i,m-n+i-1);
    end
    U = AUG(:,1:n);
    %rearranging X because during column pivoting they got rotated
    X = X * (column_permutation^-1);
    %Printing Output In file
    fileoID = fopen(output_file_name,'wt');
    fprintf(fileoID,'X\n\n');
    for i = 1:n
        fprintf(fileoID,'%f\n',X(i));
    end
    fprintf(fileoID,'\n');
    fprintf(fileoID,'U for GE pivot and scaling\n\n');
    for i = 1:n
        for j = 1:n
            fprintf(fileoID,'%f ',U(i,j));
        end
        fprintf(fileoID,'\n');
    end
    fprintf(fileoID,'\n');
    fprintf(fileoID,'Row Permutation Matrix for GE with pivot and scaling\n\n');
    for i = 1:n
        for j = 1:n
            fprintf(fileoID,'%f ',row_permutation(i,j));
        end
        fprintf(fileoID,'\n');
    end
    fprintf(fileoID,'\n');
    fprintf(fileoID,'Column Permutation Matrix for GE with pivot and scaling\n\n');
    for i = 1:n
        for j = 1:n
            fprintf(fileoID,'%f ',column_permutation(i,j));
        end
        fprintf(fileoID,'\n');
    end
    fclose(fileoID);
elseif(method_no == 4)
    %LU decomposition using Gaussian Elimination
    A = AUG(:,1:m-1);
    B = AUG(:,m);
    %A = LU where U is upper triangular matrix found using GE
    %Applying Gaussian Elimination for fidning U
    U = A;
    L = eye(n);
    for i = 1:n
        for j = i+1:n
            multiply_factor = U(j,i)/U(i,i);
            L(j,i) = multiply_factor;
            for k = 1:m-1
                U(j,k) = U(j,k) - (multiply_factor * U(i,k));
            end
        end
    end
    %Forward substitution
    Y = [];
    for i = 1:n
        sum = 0;
        for j = 1:i-1
            sum = sum + Y(j)*L(i,j);
        end
        Y(i) = (B(i) - sum)/L(i,i);
    end
    %Back substitution
    X = [];
    for i = n:-1:1
        sum = 0;
        for j = m-1:-1:m-(n-i)
            sum = sum + X(j)*U(i,j);
        end
        X(i) = (Y(i) - sum)/U(i,m-n+i-1);
    end
    %Printing Output In file
    fileoID = fopen(output_file_name,'wt');
    fprintf(fileoID,'X\n\n');
    for i = 1:n
        fprintf(fileoID,'%f\n',X(i));
    end
    fprintf(fileoID,'\n');
    fprintf(fileoID,'L for LU with GE\n\n');
    for i = 1:n
        for j = 1:n
            fprintf(fileoID,'%f ',L(i,j));
        end
        fprintf(fileoID,'\n');
    end
    fprintf(fileoID,'\n');
    fprintf(fileoID,'U for LU with GE\n\n');
    for i = 1:n
        for j = 1:n
            fprintf(fileoID,'%f ',U(i,j));
        end
        fprintf(fileoID,'\n');
    end
    fclose(fileoID);
elseif(method_no == 5)
    %LU decomposition using Gaussian Elimination with pivoting (Row)
    row_permutation = eye(n);
    column_permutation = eye(n);
    L = eye(n);
    B = AUG(:,m);
    for i = 1:n
        %Row Pivoting
        max_element = 0;
        row = 0;
        column = 0;
        for k = i:n
            if abs(AUG(k,i)) > max_element
                max_element = abs(AUG(k,i));
                row = k;
                column = j;
            end 
        end
        %row pivoting
        AUG([i,row],:) = AUG([row,i],:);
        row_permutation([i,row],:) = row_permutation([row,i],:);
        B([i,row],:) = B([row,i],:);
        %Applying Gaussian Elimination after pivoting
        for j = i+1:n
            multiply_factor = AUG(j,i)/AUG(i,i);
            L(j,i) = multiply_factor;
            for k = 1:m
                AUG(j,k) = AUG(j,k) - (multiply_factor * AUG(i,k));
            end
        end
    end
    A = AUG(:,1:m-1);
    U = A;
    
    %Forward substitution to find Y from LY = B
    Y = [];
    for i = 1:n
        sum = 0;
        for j = 1:i-1
            sum = sum + Y(j)*L(i,j);
        end
        Y(i) = (B(i) - sum)/L(i,i);
    end

    %Back substitution to find X from UX = Y
    X = [];
    for i = n:-1:1
        sum = 0;
        for j = m-1:-1:m-(n-i)
            sum = sum + X(j)*U(i,j);
        end
        X(i) = (Y(i) - sum)/U(i,m-n+i-1);
    end
    %rearranging X because during column pivoting they got rotated
    X = X * (column_permutation^-1);
    %Printing Output In file
    fileoID = fopen(output_file_name,'wt');
    fprintf(fileoID,'X\n\n');
    for i = 1:n
        fprintf(fileoID,'%f\n',X(i));
    end
    fprintf(fileoID,'\n');
    fprintf(fileoID,'L for LU with GE pivot\n\n');
    for i = 1:n
        for j = 1:n
            fprintf(fileoID,'%f ',L(i,j));
        end
        fprintf(fileoID,'\n');
    end
    fprintf(fileoID,'\n');
    fprintf(fileoID,'U for LU with GE pivot\n\n');
    for i = 1:n
        for j = 1:n
            fprintf(fileoID,'%f ',U(i,j));
        end
        fprintf(fileoID,'\n');
    end
    fprintf(fileoID,'\n');
    fprintf(fileoID,'Permutation Matrix for LU with GE pivot\n\n');
    for i = 1:n
        for j = 1:n
            fprintf(fileoID,'%f ',row_permutation(i,j));
        end
        fprintf(fileoID,'\n');
    end
    fclose(fileoID);
elseif(method_no == 6)
    %LU decomposition using Crout Method
    clear sum
    A = AUG(:,1:m-1);
    B = AUG(:,m);
    %A = LU where U is upper triangular matrix with U_ii = 1 in crout method
    U = eye(n);
    L = zeros(n);
    for i= 1:n
        for j= i:n
            num = sum(L(j,1:i-1)*U(1:i-1,i));
            L(j,i)= A(j,i) - num;
        end
        for j= i+1:n
            num = sum(L(i,1:i-1)*U(1:i-1,j));
            U(i,j) = (A(i,j) - num)/L(i,i);
        end
    end
    %Forward substitution to find Y from LY = B
    Y = [];
    for i = 1:n
        sum = 0;
        for j = 1:i-1
            sum = sum + Y(j)*L(i,j);
        end
        Y(i) = (B(i) - sum)/L(i,i);
    end
    %Back substitution to find X from UX = Y
    X = [];
    for i = n:-1:1
        sum = 0;
        for j = m-1:-1:m-(n-i)
            sum = sum + X(j)*U(i,j);
        end
        X(i) = (Y(i) - sum)/U(i,m-n+i-1);
    end
    %Printing Output In file
    fileoID = fopen(output_file_name,'wt');
    fprintf(fileoID,'X\n\n');
    for i = 1:n
        fprintf(fileoID,'%f\n',X(i));
    end
    fprintf(fileoID,'\n');
    fprintf(fileoID,'L for Crouts Method\n');
    for i = 1:n
        for j = 1:n
            fprintf(fileoID,'%f ',L(i,j));
        end
        fprintf(fileoID,'\n');
    end
    fprintf(fileoID,'\n');
    fprintf(fileoID,'U for Crouts Method\n\n');
    for i = 1:n
        for j = 1:n
            fprintf(fileoID,'%f ',U(i,j));
        end
        fprintf(fileoID,'\n');
    end
    fclose(fileoID);
else
    %Cholesky decomposition
    L = zeros(n);
    L(1,1) = sqrt(AUG(1,1));
    for i = 2:n
        L(i,1) = AUG(i,1)/L(1,1);
    end
    for j=2:n
        sum = 0;
        for k = 1:j-1
            sum = sum + (L(j,k))^2;
        end
        L(j,j) = sqrt(AUG(j,j)-sum);
        for i = j+1:n
            sum = 0;
            for k = 1:j-1
                 sum = sum + L(i,k)*L(j,k);
            end
        L(i,j) = (AUG(i,j)-sum)/L(j,j);
        end
    end  
    %U = L transpose
    U = L';
    B = AUG(:,m);
    %Forward substitution
    Y = [];
    for i = 1:n
        sum = 0;
        for j = 1:i-1
            sum = sum + Y(j)*L(i,j);
        end
        Y(i) = (B(i) - sum)/L(i,i);
    end
    %Back substitution
    X = [];
    for i = n:-1:1
        sum = 0;
        for j = m-1:-1:m-(n-i)
            sum = sum + X(j)*U(i,j);
        end
        X(i) = (Y(i) - sum)/U(i,m-n+i-1);
    end
    fileoID = fopen(output_file_name,'wt');
    fprintf(fileoID,'X\n\n');
    for i = 1:n
        fprintf(fileoID,'%f\n',X(i));
    end
    fprintf(fileoID,'\n');
    fprintf(fileoID,'Cholesky factor L_c for Cholesky Decomposition\n\n');
    for i = 1:n
        for j = 1:n
            fprintf(fileoID,'%f ',L(i,j));
        end
        fprintf(fileoID,'\n');
    end
    fclose(fileoID);
end