% Practice on sampling & reconstruction

clc;
clear; % clears all variables
close all % close all windows
% Creating "analog" signal
t = 0:.1:20;
F1 = input('Please insert first frequency component:');
F2 = input('Please insert second frequency component:');
F3 = input('Please insert third frequency component:');
B = max([F1,F2,F3]);
x = sin(2*pi*F1*t)+sin(2*pi*F2*t)+sin(2*pi*F3*t);

% Sampling
Fs = input('Please insert sampling frequency:');
Ts = 1/Fs;        % Sampling time
x_samples=x(1:10*Ts:201);

% Plotting
figure(1);
subplot(2,1,1);
plot(t,x,'k-');
title(['Original signal: B=',num2str(B)])
xlabel('t');
ylabel('x(t)');
grid on

% % Creating dialog box with explanations
% l1=[blanks(10),'Sample by sample reconstruction.'];
% l2='Blue dots: Input samples.';
% l3='Blue curve: reconstructed signal.';
% l4='Red curve: contribution to output sample from current sample.';
% l5='Click or press any key to update with 1 iteration.';
% l6='(You can keep this window open while watching the reconstruction)';
% information ={l1,'',l2,l3,l4,'',l5,'',l6};
% messagebox=msgbox(information,'Information','help');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Starting reconstruction process
subplot(2,1,2);
x_recon=0;
for k=0:length(x_samples)-1
    stem(0:Ts:20,x_samples,'filled');
    if k==length(x_samples)-1
       title('Reconstruction finished');
       xlabel('t');
       ylabel('x_r_e_c_o_n(t)');
    elseif k==0
        title(['Sampled signal: F_s =',num2str(Fs)]); 
        xlabel('n');
        ylabel('x_s(n)');
    else
       title('Sample by sample reconstruction'); 
       xlabel('t');
       ylabel('x_r_e_c_o_n(t)');
    end
    grid on;

    l=k:-.1/Ts:k-20/Ts;
    x_recon=x_recon+x_samples(k+1)*sinc(l);
    hold;
    plot(t,x_samples(k+1)*sinc(l),'r')
    plot(t,x_recon);
    hold off;
    waitforbuttonpress; 
    clc
end