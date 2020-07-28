%-----------------------------------------------------------------------------|
%This is a GUI created for the calculation of the force distribution on beam
%             written by Hanfeng Zhai and Shengjie Weng
% 
%            For the SHU summer internship program, 2018
%-----------------------------------------------------------------------------|

function varargout = finalwork(varargin)
%FINALWORK MATLAB code file for finalwork.fig
%      FINALWORK, by itself, creates a new FINALWORK or raises the existing
%      singleton*.
%
%      H = FINALWORK returns the handle to a new FINALWORK or the handle to
%      the existing singleton*.
%
%      FINALWORK('Property','Value',...) creates a new FINALWORK using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to finalwork_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      FINALWORK('CALLBACK') and FINALWORK('CALLBACK',hObject,...) call the
%      local function named CALLBACK in FINALWORK.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%a
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help finalwork

% Last Modified by GUIDE v2.5 04-Jul-2019 19:53:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @finalwork_OpeningFcn, ...
                   'gui_OutputFcn',  @finalwork_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before finalwork is made visible.
function finalwork_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for finalwork
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes finalwork wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = finalwork_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function img3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to img3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: place code in OpeningFcn to populate img3
   



function photo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to photo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% Hint: place code in OpeningFcn to populate photo


% --- Executes during object creation, after setting all properties.
function img1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to img1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: place code in OpeningFcn to populate img1


% --- Executes during object creation, after setting all properties.
function img4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to img4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate img4


% --- Executes during object creation, after setting all properties.
function img2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to img2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate img2


% --- Executes during object creation, after setting all properties.

% --- Executes during object creation, after setting all properties.
function uipanel1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function uipanel2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Value=get(handles.popupmenu1,'Value');
switch Value
    case 1%左端固定
        F=str2num(get(handles.force,'String'));
        a=str2num(get(handles.fplace,'String'));
        b=str2num(get(handles.mplace,'String'));
        L=str2num(get(handles.length,'String'));
        q=str2num(get(handles.load,'String'));
        M=str2num(get(handles.moment,'String'));
        E=str2num(get(handles.E,'String'));
        I=str2num(get(handles.I,'String'));
        c=str2num(get(handles.start,'String'));
        d=str2num(get(handles.stop,'String'));
        x=linspace(0,L,200*L+1);
        dx=1/200;
        %jianli
        F1=[(F)*ones(1,200*a),zeros(1,200*(L-a))];
        F2=zeros(1,200*L);
        Fb=[];
        for i=200*c+1:1:200*d
            Fa=q*(d-i/200);
            Fb=[Fb,Fa];
        end
        F3=[(q*(d-c))*ones(1,200*c),Fb,zeros(1,200*(L-d))];
        Fs=F1+F2+F3;
        %wanju
        Mb=[];
        for i=1:1:200*a
            Ma=-F*(a-i/200);
            Mb=[Mb,Ma];
        end
        M1=[Mb,zeros(1,200*(L-a))];
        M2=[(-M)*ones(1,200*b),zeros(1,200*(L-b))];

        Md=[];
        for i=1:1:200*c
            Mc=((d+c)*0.5-i/200)*q*(d-c);
            Md=[Md,-Mc];
        end
        Mf=[];
        for i=200*c+1:1:200*d
            Me=(d-i/200)^2*q;
            Mf=[Mf,-Me];
        end
        M3=[Md,Mf,zeros(1,200*(L-d))];
        Ms=M1+M2+M3;
        %zhuanjiao
        for i=200*L:-1:1
            t=0;
            for j=1:1:(200*L-i+1)
                t=t+Ms(j)*dx/(E*I);
            end
            A(i)=t;
        end
        %naodu
        for i=200*L:-1:1
            r=0;
            for j=1:1:(200*L-i+1)
                r=r+A(j)*dx;
            end
            Y(i)=r;
        end
        %draw
        y=1/200:1/200:L;
        %subplot(2,2,1),bar(y,Fs),grid
        bar(handles.img1,y,Fs);
        %subplot(2,2,2),bar(y,Ms),grid
        bar(handles.img2,y,Ms);
        %subplot(2,2,3),bar(y,A),grid
        bar(handles.img3,y,A);
        %subplot(2,2,4),bar(y,Y),grid
        bar(handles.img4,y,Y);
        FQmax=max(abs(Fs)); %求最大剪力
        Mmax=max(abs(Ms)); %求最大弯矩
        Amax=max(abs(A)); %求最大转角
        ymax=max(abs(y)); %求最大挠度 
        set(handles.maxmoment,'String',Mmax);
        set(handles.maxdeflection,'String',ymax);
        set(handles.maxforce,'String',FQmax);
        set(handles.maxangle,'String',Amax);
    case 2%右端固定
        F=str2num(get(handles.force,'String'));
        a=str2num(get(handles.fplace,'String'));
        b=str2num(get(handles.mplace,'String'));
        L=str2num(get(handles.length,'String'));
        q=str2num(get(handles.load,'String'));
        M=str2num(get(handles.moment,'String'));
        E=str2num(get(handles.E,'String'));
        I=str2num(get(handles.I,'String'));
        c=str2num(get(handles.start,'String'));
        d=str2num(get(handles.stop,'String'));
        x=linspace(0,L,200*L+1);
        dx=1/200;
        %jianli
        F1=[zeros(1,200*a),(-F)*ones(1,200*(L-a))];
        F2=zeros(1,200*L);
        Fb=[];
        for i=200*c+1:1:200*d
            Fa=q*(c-i/200);
            Fb=[Fb,Fa];
        end
        F3=[zeros(1,200*c),Fb,q*(c-d)*ones(1,200*(L-d))];
        Fs=F1+F2+F3;
        %wanju
        Mb=[];
        for i=200*a+1:1:200*L
            Ma=F*(a-i/200);
            Mb=[Mb,Ma];
        end
        M1=[zeros(1,200*a),Mb];
        M2=[zeros(1,200*b),(M)*ones(1,200*(L-b))];

        Md=[];
        for i=200*c+1:1:200*d
            Mc=-0.5*q*(i/200-c)^2;
            Md=[Md,Mc];
        end
        Mf=[];
        for i=200*d+1:1:200*L
            Me=-q*(d-c)*(i/200-0.5*(d+c));
            Mf=[Mf,Me];
        end
        M3=[zeros(1,200*c),Md,Mf];
        Ms=M1+M2+M3;
        %zhuanjiao
        for i=200*L:-1:1
            t=0;
            for j=1:1:(200*L-i+1)
                t=t+Ms(j)*dx/(E*I);
            end
            A(i)=t;
        end
        %naodu
        for i=200*L:-1:1
            r=0;
            for j=1:1:(200*L-i+1)
                r=r+A(j)*dx;
            end
            Y(i)=r;
        end
        %draw
        y=1/200:1/200:L;
        %subplot(2,2,1),bar(y,Fs),grid
        bar(handles.img1,y,Fs);
        %subplot(2,2,2),bar(y,Ms),grid
        bar(handles.img2,y,Ms);
        %subplot(2,2,3),bar(y,A),grid
        bar(handles.img3,y,A);
        %subplot(2,2,4),bar(y,Y),grid
        bar(handles.img4,y,Y);
        FQmax=max(abs(Fs)); %求最大剪力
        Mmax=max(abs(Ms)); %求最大弯矩
        Amax=max(abs(A)); %求最大转角
        ymax=max(abs(y)); %求最大挠度 
        set(handles.maxmoment,'String',Mmax);
        set(handles.maxdeflection,'String',ymax);
        set(handles.maxforce,'String',FQmax);
        set(handles.maxangle,'String',Amax);
    case 3%左端外伸
        x=0:0.1:100;
        y=0:0.1:100;
        F=str2num(get(handles.force,'String'));
        a=str2num(get(handles.start,'String'));
        b=str2num(get(handles.stop,'String'));
        L=str2num(get(handles.length,'String'));
        q=str2num(get(handles.load,'String'));
        Me=str2num(get(handles.moment,'String'));
        E=str2num(get(handles.E,'String'));
        I=str2num(get(handles.I,'String'));
        %function y=test1(F,a,b,L,q,Me,E,I)
        %左外伸梁
        %F=20000;a=1;b=4;L=5;q=10000;Me=40000;E=200e9;I=2e-5;%输入已知条件
        FRA =(F*L-Me)/b+0.5*q*b;
        FRB =F+q*b-FRA; %求约束反力
        x=linspace(0,L,1001); %在 0 与 L 之间产生 1001 个点
        FQ1=(-F)*ones([1,200]); %第一段剪力
        FQ2=-F+FRA-q*(x(201:1001)-a);%第二段剪力
        FQ=[FQ1,FQ2];%全梁剪力
        dx=L/1000; %步长为 L/1000
        M3=cumtrapz(FQ1)*dx;
        M4=cumtrapz(FQ2)*dx; %对剪力积分求弯矩
        C1=-M3(1);
        C2=-M4(1001-200); %确定弯矩积分常数
        M1=M3+C1;
        M2=M4+C2;
        M=[M1, M2]; %全梁弯矩方程
        A0=cumtrapz(M)*dx/(E*I);
        y0=cumtrapz(A0)*dx; %积分求转角、挠度
        B=[a,1;L,1]\[-y0(201);-y0(1001)];
        C3=B(1);
        C4=B(2);%确定转角、挠度积分常数
        A=A0+C3;
        y=y0+C3*x+C4; %全梁转角、挠度方程

        %subplot(3,2,1),plot(x,FQ),grid %画剪力图
        bar(handles.img1,x,FQ);
        %subplot(3,2,2),plot(x,M),grid %画弯矩图
        bar(handles.img2,x,M);
        %subplot(3,2,3),plot(x,A),grid %画转角图
        bar(handles.img3,x,A);
        %subplot(3,2,4),plot(x,y),grid %画挠度图
        bar(handles.img4,x,y);
        FQmax=max(abs(FQ)); %求最大剪力
        Mmax=max(abs(M)); %求最大弯矩
        Amax=max(abs(A)); %求最大转角
        ymax=max(abs(y)); %求最大挠度 
        set(handles.maxmoment,'String',Mmax);
        set(handles.maxdeflection,'String',ymax);
        set(handles.maxforce,'String',FQmax);
        set(handles.maxangle,'String',Amax);
        %D=FQmax;
        %D=set(handles.uitable1;'Data','FQmax');
    case 4%右端外伸
        x=0:0.1:100;
        y=0:0.1:100;
        p=str2num(get(handles.force,'String'));
        l3=str2num(get(handles.fplace,'String'));
        l1=str2num(get(handles.mplace,'String'));
        L=str2num(get(handles.length,'String'));
        q=str2num(get(handles.load,'String'));
        m=str2num(get(handles.moment,'String'));
        E=str2num(get(handles.E,'String'));
        I=str2num(get(handles.I,'String'));
        l2=L-l1-l3;
        Fra=(m+0.5*q*(l2+l3)*(l2-l3)-p*l3)/(l1+l2);
        %Frb = (q*(12+ 13)* (l1+0.5 *(l2+l3))+p*(l1+l2+l3)-m)/(l1+l2);% 求支反力 
        x=linspace(0,L,121);
        dx=L/120; 
        M1=Fra*x(1:20); 
        Fq1=Fra*ones([1,20]) 
        M2=(Fra*x(21:100)-m-0.5*q*(x(21:100)-2).^2); 
        Fq2=Fra-q*(x(21:100)-l1) 
        M3=(-p*(L-x(101:121))-0.5*q*(L-x(101:121)).^2); 
        Fq3=p+q*(L-x(101:121)) 
        Fq=[Fq1,Fq2,Fq3];
        M=[-M1,-M2,-M3];% 分段列剪力、 弯矩方程 
        A0=cumtrapz(-M)*dx/(E*I);% 积分一次求转角
        Y0=cumtrapz(A0)*dx;
        C=[0,1;10,1]\[-Y0(1);-Y0(101)];
        Ca=C(1),Cy=C(2);
        A=A0+Ca,Y=Y0+Ca*x+Cy;
        % subplot(3,2,1),plot(x,Fq),grid
        bar(handles.img1,x,Fq);
        % subplot(3,2,2),plot(x,M),grid
        bar(handles.img2,x,M);
        % subplot(3,2,3),plot(x,A),grid
        bar(handles.img3,x,A);
        % subplot(3,2,4),plot(x,Y),grid
        bar(handles.img4,x,Y);
        FQmax=max(abs(Fq)); %求最大剪力
        Mmax=max(abs(M)); %求最大弯矩
        Amax=max(abs(A)); %求最大转角
        ymax=max(abs(y)); %求最大挠度 
        set(handles.maxmoment,'String',Mmax);
        set(handles.maxdeflection,'String',ymax);
        set(handles.maxforce,'String',FQmax);
        set(handles.maxangle,'String',Amax);
    case 5%简支梁
        F=str2num(get(handles.force,'String'));
        b=str2num(get(handles.fplace,'String'));
        a=str2num(get(handles.mplace,'String'));
        L=str2num(get(handles.length,'String'));
        q=str2num(get(handles.load,'String'));
        M=str2num(get(handles.moment,'String'));
        E=str2num(get(handles.E,'String'));
        I=str2num(get(handles.I,'String'));
        c=str2num(get(handles.start,'String'));
        d=str2num(get(handles.stop,'String'));
        Fr1=-M/L;
        Fr2=-F*(L-b)/L;
        Fr3=q*(d-c)*(L-0.5*(c+d))/L;
        x=linspace(0,L,200*L+1);dx=1/200;
        %jianli
        F1=(Fr1)*ones(1,200*L);
        F2=[(Fr2)*ones(1,200*b),(Fr2-F)*ones(1,200*(L-b))];
        Fb=[];
        for i=200*c+1:1:200*d
            Fa=Fr3-q*(i/200-c);
            Fb=[Fb,Fa];
        end
        F3=[(Fr3)*ones(1,200*c),Fb,(Fr3-q*(d-c))*ones(1,200*(L-d))];
        Fs=F1+F2+F3;
        %wanju
        Mb=[];
        for i=1:1:200*a
            Ma=Fr1*i/200;
            Mb=[Mb,Ma];
        end
        Md=[];
        for i=200*a+1:1:200*L
            Mc=Fr1*i/200+M;
            Md=[Md,Mc];
        end
        M1=[Mb,Md];

        Mf=[];
        for i=1:1:200*b
            Me=Fr2*i/200;
            Mf=[Mf,Me];
        end
        Mh=[];
        for i=200*b+1:1:200*L
            Mg=Fr2*i/200-F*(i/200-b);
            Mh=[Mh,Mg];
        end
        M2=[Mf,Mh];

        Mj=[];
        for i=1:1:200*c
            Mi=Fr3*i/200;
            Mj=[Mj,Mi];
        end
        Ml=[];
        for i=200*c+1:1:200*d
            Mk=Fr3*i/200-0.5*q*(i/200-c)^2;
            Ml=[Ml,Mk];
        end
        Mn=[];
        for i=200*d+1:1:200*L
            Mm=Fr3*i/200-0.5*q*(i/200-c)^2+0.5*q*(i/200-d)^2;
            Mn=[Mn,Mm];
        end
        M3=[Mj,Ml,Mn];
        Ms=M1+M2+M3;
        A0=cumtrapz(Ms)*dx/(E*I);
        Y0=cumtrapz(A0)*dx;
        %draw
        y=1/200:1/200:L;
        C=[0 1;L 1]\[-Y0(1);-Y0(200*L)];
        Ca=C(1),Cy=C(2);
        A=A0+Ca;
        Y=Y0+Ca.*y+Cy;
        %subplot(2,2,1),bar(y,Fs),grid
        bar(handles.img1,y,Fs);
        %subplot(2,2,2),bar(y,Ms),grid
        bar(handles.img2,y,Ms);
        %subplot(2,2,3),bar(y,A),grid
        bar(handles.img3,y,A);
        %subplot(2,2,4),bar(y,Y),grid
        bar(handles.img4,y,Y);
        FQmax=max(abs(Fs)); %求最大剪力
        Mmax=max(abs(Ms)); %求最大弯矩
        Amax=max(abs(A)); %求最大转角
        ymax=max(abs(y)); %求最大挠度 
        set(handles.maxmoment,'String',Mmax);
        set(handles.maxdeflection,'String',ymax);
        set(handles.maxforce,'String',FQmax);
        set(handles.maxangle,'String',Amax);
end



% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.E,'String','');
set(handles.I,'String','');
set(handles.force,'String','');
set(handles.load,'String','');
set(handles.moment,'String','');
set(handles.start,'String','');
set(handles.stop,'String','');
set(handles.maxmoment,'String','');
set(handles.maxdeflection,'String','');
set(handles.maxforce,'String','');
set(handles.maxangle,'String','');
set(handles.fplace,'String','');
set(handles.mplace,'String','');
set(handles.length,'String','');
try
    delete(allchild(handles.img1));
    delete(allchild(handles.img2));
    delete(allchild(handles.img3));
    delete(allchild(handles.img4));
end



% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function load_Callback(hObject, eventdata, handles)
% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of load as text
%        str2double(get(hObject,'String')) returns contents of load as a double


% --- Executes during object creation, after setting all properties.
function load_CreateFcn(hObject, eventdata, handles)
% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function stop_Callback(hObject, eventdata, handles)
% hObject    handle to stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stop as text
%        str2double(get(hObject,'String')) returns contents of stop as a double


% --- Executes during object creation, after setting all properties.
function stop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function start_Callback(hObject, eventdata, handles)
% hObject    handle to start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of start as text
%        str2double(get(hObject,'String')) returns contents of start as a double


% --- Executes during object creation, after setting all properties.
function start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function length_Callback(hObject, eventdata, handles)
% hObject    handle to length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of length as text
%        str2double(get(hObject,'String')) returns contents of length as a double


% --- Executes during object creation, after setting all properties.
function length_CreateFcn(hObject, eventdata, handles)
% hObject    handle to length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function E_Callback(hObject, eventdata, handles)
% hObject    handle to E (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of E as text
%        str2double(get(hObject,'String')) returns contents of E as a double


% --- Executes during object creation, after setting all properties.
function E_CreateFcn(hObject, eventdata, handles)
% hObject    handle to E (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function I_Callback(hObject, eventdata, handles)
% hObject    handle to I (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of I as text
%        str2double(get(hObject,'String')) returns contents of I as a double


% --- Executes during object creation, after setting all properties.
function I_CreateFcn(hObject, eventdata, handles)
% hObject    handle to I (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxmoment_Callback(hObject, eventdata, handles)
% hObject    handle to maxmoment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxmoment as text
%        str2double(get(hObject,'String')) returns contents of maxmoment as a double


% --- Executes during object creation, after setting all properties.
function maxmoment_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxmoment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxangle_Callback(hObject, eventdata, handles)
% hObject    handle to maxangle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxangle as text
%        str2double(get(hObject,'String')) returns contents of maxangle as a double


% --- Executes during object creation, after setting all properties.
function maxangle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxangle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxforce_Callback(hObject, eventdata, handles)
% hObject    handle to maxforce (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxforce as text
%        str2double(get(hObject,'String')) returns contents of maxforce as a double


% --- Executes during object creation, after setting all properties.
function maxforce_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxforce (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxdeflection_Callback(hObject, eventdata, handles)
% hObject    handle to maxdeflection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxdeflection as text
%        str2double(get(hObject,'String')) returns contents of maxdeflection as a double


% --- Executes during object creation, after setting all properties.
function maxdeflection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxdeflection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function force_Callback(hObject, eventdata, handles)
% hObject    handle to force (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of force as text
%        str2double(get(hObject,'String')) returns contents of force as a double


% --- Executes during object creation, after setting all properties.
function force_CreateFcn(hObject, eventdata, handles)
% hObject    handle to force (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fplace_Callback(hObject, eventdata, handles)
% hObject    handle to fplace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fplace as text
%        str2double(get(hObject,'String')) returns contents of fplace as a double


% --- Executes during object creation, after setting all properties.
function fplace_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fplace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function moment_Callback(hObject, eventdata, handles)
% hObject    handle to moment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of moment as text
%        str2double(get(hObject,'String')) returns contents of moment as a double


% --- Executes during object creation, after setting all properties.
function moment_CreateFcn(hObject, eventdata, handles)
% hObject    handle to moment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mplace_Callback(hObject, eventdata, handles)
% hObject    handle to mplace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mplace as text
%        str2double(get(hObject,'String')) returns contents of mplace as a double


% --- Executes during object creation, after setting all properties.
function mplace_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mplace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in show.
function show_Callback(hObject, eventdata, handles)
% hObject    handle to show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Value=get(handles.popupmenu1,'Value');
switch Value
    case 1
       imshow(imread('pic1.PNG'));
    case 2
       imshow(imread('pic2.PNG'));
    case 3
       imshow(imread('photo.PNG'));
    case 4
       imshow(imread('yws.PNG'));
    case 5
       imshow(imread('pic3.PNG'));
end


% --- Executes on button press in disp.
function disp_Callback(hObject, eventdata, handles)
% hObject    handle to disp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Value=get(handles.popupmenu1,'Value');
switch Value
    case 1
       imshow(imread('pic1.PNG'));
    case 2
       imshow(imread('pic2.PNG'));
    case 3
       imshow(imread('photo.PNG'));
    case 4
       imshow(imread('yws.PNG'));
    case 5
       imshow(imread('pic3.PNG'));
end


% --- Executes during object creation, after setting all properties.
