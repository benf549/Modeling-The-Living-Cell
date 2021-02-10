%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Benjmain Fry (bfry2)
% 01/28/21, coded on MATLAB R2020b
% 
% Fixed fix_me_script() function. Takes in numeric frequency and dampening constant
% and a positive integer for figure number. Plots a damp oscillator for t between 0 to 10 for
% the given frequency and dampening constant.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function[]=fix_me_script(freq, damp, figure_num)
    t=0:.1:10;
    f1=exp(-damp*t);
    f2=sin(freq*t);
    
    %define a damped oscillator
    y=f1.*f2; %performs elementwise multiplication rather than matrix multiplication.
    
    figure(figure_num);
    plot(t,y,'r-','LineWidth',3)
    title("Plot of Damped Oscillator as a Function of Time")
    xlabel("t")
    ylabel("y")
end