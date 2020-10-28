
% EDIT THE FOLLOWING LINE IN CASE OF DIFFERENT INPUT FILE...

FILENAME='plotdata.dat';     %Remember to remove the commented header in the .dat file!

step=144;					 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plotdata=load(FILENAME);

ran=[1:step:(size(plotdata,2)-1)];
outdata=ones(size(plotdata,1), size(ran,2)+1); outdata(:,1)=plotdata(:,1);

for index1=1:size(ran,2) 
    buf=zeros(size(plotdata,1),1);
    for index2=1:(step-1)
        buf=buf+plotdata(:,(ran(index1)+index2));
    end
    %data=[plotdata(:,1) buf];
    %newname=['data_',num2str(ran(index1)),'_',num2str(ran(index1)+(step-1)),'.dat'];
    %save newname data -ascii
    %
    outdata(:,1+index1)=buf./step;
end

NEW_FILENAME=[FILENAME,'_merged_data.dat'];
save NEW_FILENAME outdata -ascii


% Coded in quick and dirty way... For Info, Matteo x8498