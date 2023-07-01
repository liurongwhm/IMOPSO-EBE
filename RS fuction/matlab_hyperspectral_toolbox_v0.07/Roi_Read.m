%读取ENVI的ROI
function [C Name rgb npts Train]=Roi_Read()
[File_Name, File_Path] = uigetfile('.txt', 'Select image ROI to open');
Fid_roi=fopen([File_Path File_Name],'r');
elements={'; Number of ROIs' '; ROI name' '; ROI rgb value' '; ROI npts'};
% Number of ROIs 类别数量
% ROI name 类别名称
% ROI rgb value 类别彩色分量
% ROI npts 类别训练样本数
Name='';
spa=' ';
flag0=0;
while 1
    fline=fgetl(Fid_roi);
    [first,second]=strtok(fline,':');
    flag=0;%用于判断
    switch first
        case elements(1)
              [f,num]=strtok(second);
               classes=str2num(num);
               C=classes;
               flag=1;
               flag0=1;
               c=0;
        case elements(2)
                 [f,Nam]=strtok(second);
                 if c<9
                 na=strcat(Nam,spa);
                 Name=strcat(Name,na);
                 else
                     Name=strcat(Name,Nam);
                 end
        case elements(3)
                 [f,str]=strtok(second);
                 num=length(str);
                 rgb_str=str(3:num-1);
                 [f1,f2]=strtok(rgb_str,',');
                 rgb(c,1)=str2num(f1);
                 [f1,f2]=strtok(f2);
                 [f1,f2]=strtok(f2);
                 rgb(c,2)=str2num(f1);
                 [f1,f2]=strtok(f2);
                 rgb(c,3)=str2num(f1);
        case elements(4)
                 [f,num]=strtok(second);
                 npts(c,1)=str2num(num);
                 flag=1;
    end
    if flag==1
        c=c+1;
    end
     if flag0==1 & c==classes+1
         break;
     end
end
fline=fgetl(Fid_roi);
    for i=1:sum(npts)
        if mod(i,1000)==1
         tic
        end
        fline=fgetl(Fid_roi);
        if isempty(fline)
            fline=fgetl(Fid_roi);
        end
        Train(i,:)=str2num(fline);
        if mod(i,1000)==0
            i
            toc
        end
      
    end
fclose(Fid_roi);
end

                 
    