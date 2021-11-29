clear;
close all;
clc
%Aviral Agarwal
%180167 - Tutorial Section - J6
%Computer Assignment 3
prompt = "Which method do you want to run first?\n1) Direct power method\n2) Inverse power method\n3)  Shifted-power method\n4) QR method\nEnter value from 1 to 4: ";
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
A = readmatrix(file_name);
A = A([1:n],:);
k = 0;
while ischar(tline)
    k = k+1;
    if(k==m+1)
        max_iter = str2num(tline);
    elseif (k==m+2)
        max_e = str2num(tline);
    else
        if(k==m+3 && method_no == 3)
            shift = str2num(tline);
        end
    end
    tline = fgetl(fid);
end
fclose(fid);

if(method_no == 1)
    %Direct Power Method
    X = ones(n,1);
    iter_no = 0;
    eig_vec = ones(n,1);
    eig_val = 1;
    tot_iter_no = 0;
    while(iter_no < max_iter)
        Y = A*X;
        s = max(abs(Y));
        if iter_no > 0
            curr_e = abs((s-prev_s)/s)*100;
            if curr_e < max_e
                eig_val = s;
                eig_vec = (X/sqrt(sum(X.^2)));
                tot_iter_no = iter_no + 1;
                break
            end
        end
        prev_s = s;
        X = Y/s;
        iter_no = iter_no + 1;
    end
    if iter_no >= max_iter
        eig_val = s;
        eig_vec = (X/sqrt(sum(X.^2)));
        tot_iter_no = iter_no;
    end
    %Printing Output In file
    fileoID = fopen(output_file_name,'wt');
    fprintf(fileoID,'Direct Power Method\n\n');
    fprintf(fileoID,'Eigenvalue\n\n');
    fprintf(fileoID,'%f\n\n',eig_val);
    fprintf(fileoID,'Eigenvector\n\n');
    for i = 1:n
        fprintf(fileoID,'%f\n',eig_vec(i));
    end
    fprintf(fileoID,'\n');
    fprintf(fileoID,'Iterations\n\n');
    fprintf(fileoID,'%d\n\n',tot_iter_no);
    fclose(fileoID);
elseif(method_no == 2)
    %Inverse Power Method
    A = inv(A);
    X = ones(n,1);
    iter_no = 0;
    eig_vec = ones(n,1);
    eig_val = 1;
    tot_iter_no = 0;
    while(iter_no < max_iter)
        Y = A*X;
        s = max(abs(Y));
        if iter_no > 0
            curr_e = abs((s-prev_s)/s)*100;
            if curr_e < max_e
                eig_val = 1/s;
                eig_vec = (X/sqrt(sum(X.^2)));
                tot_iter_no = iter_no + 1;
                break
            end
        end
        prev_s = s;
        X = Y/s;
        iter_no = iter_no + 1;
    end
    if iter_no >= max_iter
        eig_val = 1/s;
        eig_vec = (X/sqrt(sum(X.^2)));
        tot_iter_no = iter_no;
    end
    %Printing Output In file
    fileoID = fopen(output_file_name,'wt');
    fprintf(fileoID,'Inverse Power Method\n\n');
    fprintf(fileoID,'Eigenvalue\n\n');
    fprintf(fileoID,'%f\n\n',eig_val);
    fprintf(fileoID,'Eigenvector\n\n');
    for i = 1:n
        fprintf(fileoID,'%f\n',eig_vec(i));
    end
    fprintf(fileoID,'\n');
    fprintf(fileoID,'Iterations\n\n');
    fprintf(fileoID,'%d\n\n',tot_iter_no);
    fclose(fileoID);
elseif(method_no == 3)
    %Shifted Power method based on Gershgorin disc
    A = A - shift*eye(n);
    A = inv(A);
    X = ones(n,1);
    iter_no = 0;
    eig_vec = ones(n,1);
    eig_val = 1;
    tot_iter_no = 0;
    while(iter_no < max_iter)
        Y = A*X;
        s = max(abs(Y));
        if iter_no > 0
            curr_e = abs((s-prev_s)/s)*100;
            if curr_e < max_e
                eig_val = (1/s) + shift;
                eig_vec = (X/sqrt(sum(X.^2)));
                tot_iter_no = iter_no + 1;
                break
            end
        end
        prev_s = s;
        X = Y/s;
        iter_no = iter_no + 1;
    end
    if iter_no >= max_iter
        eig_val = (1/s) + shift;
        eig_vec = (X/sqrt(sum(X.^2)));
        tot_iter_no = iter_no;
    end
    %Printing Output In file
    fileoID = fopen(output_file_name,'wt');
    fprintf(fileoID,'Shifted Power Method\n\n');
    fprintf(fileoID,'Eigenvalue\n\n');
    fprintf(fileoID,'%f\n\n',eig_val);
    fprintf(fileoID,'Eigenvector\n\n');
    for i = 1:n
        fprintf(fileoID,'%f\n',eig_vec(i));
    end
    fprintf(fileoID,'\n');
    fprintf(fileoID,'Iterations\n\n');
    fprintf(fileoID,'%d\n\n',tot_iter_no);
    fclose(fileoID);
elseif(method_no == 4)
    %QR method to find all the eigen values
    iter_no = 0;
    eig_val = [];
    tot_iter_no = 0;
    while(iter_no < max_iter)
        Q = zeros(n,n);
        A_new = zeros(n,n);
        % Finding q1,q2,q3...
        for i = 1:n
            if i == 1
                Q([1:n],i) = A([1:n],i)/sqrt(sum(A([1:n],i).^2));
            else
                A_new([1:n],i) = A([1:n],i);
                for j = 1:i-1
                    A_new([1:n],i) = A_new([1:n],i) - ((Q([1:n],j)'*A([1:n],i))*Q([1:n],j));
                end
                Q([1:n],i) = A_new([1:n],i)/sqrt(sum(A_new([1:n],i).^2));
            end
        end
        % Finding R
        R = zeros(n,n);
        for i = 1:n
            q_trans = Q([1:n],i)';
            for j = i:n
                R(i,j) = q_trans*A([1:n],j);
            end
        end
        %New value of A
        disp(Q)
        A_new = R*Q;
        s = max(abs(A_new(:)));
        if iter_no > 0
            curr_e = abs((s-prev_s)/s)*100;
            if curr_e < max_e
                eig_val = diag(A_new);
                tot_iter_no = iter_no + 1;
                break
            end
        end
        prev_s = s;
        %Again decomposing new A into Q and R and repeating the steps
        A = A_new;
        iter_no = iter_no + 1;
    end
    if iter_no >= max_iter
        eig_val = diag(A);
        tot_iter_no = iter_no;
    end
    %Printing Output In file
    fileoID = fopen(output_file_name,'wt');
    fprintf(fileoID,'QR Method\n\n');
    fprintf(fileoID,'Eigenvalues\n\n');
    for i = 1:n
        fprintf(fileoID,'%f\n',eig_val(i));
    end
    fprintf(fileoID,'\n');
    fprintf(fileoID,'Iterations\n\n');
    fprintf(fileoID,'%d\n\n',tot_iter_no);
    fclose(fileoID);
else
    print("Method number out of range please rerun the program and enter the correct value");
end
clear