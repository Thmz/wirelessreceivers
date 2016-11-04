function output = matched_filter_np(input,os_factor,filterlength)
rolloff = 0.22;

pulse  = rrc(os_factor,rolloff,filterlength);
output = zeros(1,length(input)-2*filterlength);
data  = input(:);

for i=1:length(input)-2*filterlength
    
    segment   = data(i:i+2*filterlength);
    output(i) = segment.'*pulse;
    
end

return