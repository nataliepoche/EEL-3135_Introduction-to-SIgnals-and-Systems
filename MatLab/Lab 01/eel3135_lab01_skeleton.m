%% PREAMBLE
% DO NOT REMOVE THE LINE BELOW
clear;




%% QUESTION 1: COMMENTING
% =======================

% Copy and comment every line of the following MATLAB script. Say what
% each line is doing in your comment. Explain each MATLAB line by using
% no more than one comment line, as done in the first line below. Run and
% publish the script:
a=zeros(1,10) % Generate and print a 1x10 row vector of zeros
b=ones(4,2)  
c=size(b);      
abs([-2.2 , 3])  
floor(1.6)       
d=[1:-2.5:-9];   
f=d(2); g=sin(pi/2);     
K=[1.4, 2.3; 5.1, 7.8];     
m=K(1,2);                   
n=K(:,2);                  
p=K(1,2);                   
comp = 10+40i;              
real(comp)                 
imag(comp)                 
abs(comp)                   
angle(comp)                
disp('haha, MATLAB is fun');    
3^2                             
4==4                           
[2==8 3~=5]                     
x=[1:3:10];                    
y=[5 9 6 8];                   
tic; pause(0.2); toc          

q = zeros(10,1);               
for ii = 1:10                  
    q(ii) = ii^2;              
end                             
figure(129);                  
stem(x,y)                      
hold on;                       
plot(x,y, 'k', 'linewidth', 2)     
plot(x,y,'+r', 'markersize', 20); 
hold off;                   
xlabel('Horizontal Axis')     
ylabel('Vertical Axis')        




%% QUESTION 2: PLOTTING
% =======================

%% 2(a) PLOT RESULT



%% 2(b) PLOT RESULT



%% 2(c) PLOT RESULT



%% QUESTION 3: COMPLEX ROOTS
% =======================

%% 3(a) WRITE FUNCTION IN SEPARATE FILE (TEMPLATE PROVIDED)
type('myroots.m')


%% 3(b) ANSWER QUESTION



%% 3(c) OUTPUT RESULTS



