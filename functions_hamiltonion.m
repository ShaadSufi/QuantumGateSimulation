function varargout = GuideSQUID(varargin)


gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GuideSQUID_OpeningFcn, ...
                   'gui_OutputFcn',  @GuideSQUID_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end
if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
function GuideSQUID_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;

current_dir=pwd;
newdir= uigetdir(current_dir,'Select the folder you want to see the IV curves of:');
set(handles.textpath,'String',newdir);
[I,V,Vctrl]=fileimp(newdir);
[points,N]= size(I);
Rctrl=5550;
Rleads=250;
nrange=10;
traceup=floor(points/4):-1:1; 
retraceup=floor(points/4)+1:floor(points/2); 
tracedown=floor(points/2)+1:floor(3*points/4);
retracedown=points:-1:floor(3*points/4)+1; 
trace=[traceup,tracedown];
retrace=[retraceup,retracedown];
Imod=zeros(size(I));
Iexc=zeros(N,2);
RN=zeros(N,2);
for k=1:N
    pU=polyfit(V(floor(points/4-points/20):floor(points/4+points/20),k),I(floor(points/4-points/20):floor(points/4+points/20),k),1);
    pD=polyfit(V(floor(3.*points/4-points/20):floor(3.*points/4+points/20),k),I(floor(3.*points/4-points/20):floor(3.*points/4+points/20),k),1);
    Imod(:,k)=[I(1:floor(points/2),k)-polyval(pU,V(1:floor(points/2),k));I(floor(points/2)+1:points,k)-polyval(pD,V(floor(points/2)+1:points,k))];
    Iexc(k,:)=[pU(2),pD(2)];
    RN(k,:)=[1./pU(1),1./pD(1)];
end
G=LSQ(I,V,nrange);

Ij=Vctrl./(Rctrl+Rleads);
Ijm=zeros(points,N);
for k=1:points
    Ijm(k,:)=Ij;
end
k=1;

handles.Rctrl=Rctrl;
handles.Rleads=Rleads;
handles.nrange=nrange;
handles.Ij=Ij;
handles.Ijm=Ijm;
handles.k=k;
handles.I0= I;
handles.V0= V;
handles.I= I;
handles.V= V;
handles.Vctrl= Vctrl;
handles.N=N;
handles.G=G;
handles.points=points;
handles.trace=trace;
handles.retrace=retrace;
handles.Imod=Imod;
handles.Iexc=Iexc;
handles.RN=RN;
handles.flag=0;
handles.offset=[0,0]; %Offset=[Voff,Ioff];

set(handles.kslider,'Max',N);
set(handles.kslider,'Min',1);
set(handles.kslider,'Value',1);
set(handles.kslider,'SliderStep',[1./N .05]);
set(handles.leadr,'String',num2str(Rleads));
set(handles.crtlr,'String',num2str(Rctrl));
set(handles.range,'String',num2str(nrange));
set(handles.ioffset,'String',num2str(0));
set(handles.voffset,'String',num2str(0));
guidata(hObject, handles);
visualize(hObject,handles);
textupdate(hObject,handles);
function varargout = GuideSQUID_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;
function kslider_Callback(hObject, eventdata, handles)

k = get(hObject,'Value');
k=ceil(k);
handles.k=k;

guidata(hObject, handles);
visualize(hObject,handles);
textupdate(hObject,handles);

function PlotSTYLE_Callback(hObject, eventdata, handles)
setslider(hObject, handles)
visualize(hObject,handles);
textupdate(hObject,handles);
function plotstyle_Callback(hObject, eventdata, handles)
plotstyle=get(handles.plotstyle,'Value');
if plotstyle == 3  | plotstyle == 4  
    set(handles.kslider,'Value',1)
    handles.k=1;
    % Update handles structure
    guidata(hObject, handles);
end
visualize(hObject,handles);
textupdate(hObject,handles);

function figureextract_Callback(hObject, eventdata, handles)
figure,
visualize(hObject,handles);

Vctrl=handles.Vctrl;
Rctrl=handles.Rctrl;
Rleads=handles.Rleads;
RN=handles.RN;
G=handles.G;
I=handles.I;
V=handles.V;
Ijm=handles.Ijm;
Ij=handles.Ij;
Iexc=handles.Iexc;
datapath=get(handles.textpath,'String');
uisave({'Vctrl','Rctrl','Rleads','G','I','V','Ijm','Ij','RN','Iexc','datapath'})
function leadr_Callback(hObject, eventdata, handles)

Buffer=str2double(get(handles.leadr,'String'));
Current=handles.Rleads;
if (Current-Buffer)~=0 && ~(isequalwithequalnans(Buffer,NaN))
    handles.Rleads=Buffer;
    handles.flag=max([1,handles.flag]);
else
    set(handles.leadr,'String', Current);
end
guidata(hObject, handles);

function crtlr_Callback(hObject, eventdata, handles)

Buffer=str2double(get(handles.crtlr,'String'));
Current=handles.Rctrl;
if (Current-Buffer)~=0 && ~(isequalwithequalnans(Buffer,NaN))
    handles.Rctrl=Buffer;
    handles.flag=max([1,handles.flag]);
else
    set(handles.crtlr,'String', Current);
end
% Update handles structure
guidata(hObject, handles);
% hObject    handle to voffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of voffset as text
%        str2double(get(hObject,'String')) returns contents of voffset as a double
Buffer=str2double(get(handles.voffset,'String'));
Current=handles.offset;
if (Current(1)-Buffer)~=0 && ~(isequalwithequalnans(Buffer,NaN))
    Current(1)=Buffer;
    handles.offset=Current;
    handles.flag=max([1,handles.flag]);
else
    set(handles.crtlr,'String', Current);
end
% Update handles structure
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ioffset_Callback(hObject, eventdata, handles)
% hObject    handle to ioffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of ioffset as text
%        str2double(get(hObject,'String')) returns contents of ioffset as a double
Buffer=str2double(get(handles.voffset,'String'));
Current=handles.offset;
if (Current(2)-Buffer)~=0 && ~(isequalwithequalnans(Buffer,NaN))
    Current(2)=Buffer;
    handles.offset=Current;
    handles.flag=max([1,handles.flag]);
else
    set(handles.crtlr,'String', Current);
end
% Update handles structure
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function range_Callback(hObject, eventdata, handles)
% hObject    handle to range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Buffer=ceil(str2double(get(handles.range,'String')));
Current=handles.nrange;
if (Current-Buffer)~=0 && ~(isequalwithequalnans(Buffer,NaN))
    handles.nrange=Buffer;
    handles.flag=max([2,handles.flag]);
else
    set(handles.nrange,'String', Current);
end

guidata(hObject, handles);

function getoffset_Callback(hObject, eventdata, handles)

set(handles.plotstyle,'Value',1)
visualize(hObject,handles);
offset=mean(ginput); offset=offset.*[1e-3,1e-6];
set(handles.voffset,'String',num2str(offset(1)));
set(handles.ioffset,'String',num2str(offset(2)));
handles.offset=offset;
handles.flag=1;

guidata(hObject, handles);


function update_Callback(hObject, eventdata, handles)

flag=handles.flag;
if flag~=0
    
    I=handles.I0;
    V=handles.V0;
    
   
    offset=handles.offset;
    Vctrl=handles.Vctrl;
    Rctrl=handles.Rctrl;
    Rleads=handles.Rleads;
    points=handles.points;
    N=handles.N;
    nrange=handles.nrange;
    Ij=Vctrl./(Rctrl+Rleads);
    Ijm=zeros(points,N);
    for k=1:points
        Ijm(k,:)=Ij;
    end

    handles.Ij=Ij;
    handles.Ijm=Ijm;
    if flag==2
        G=LSQ(I,V,nrange);
        handles.G=G;
    end
    
    I=I-offset(2);
    V=V-offset(1);
    
    Imod=zeros(size(I));
    Iexc=zeros(N,2);
    RN=zeros(N,2);
    for k=1:N
        pU=polyfit(V(floor(points/4-points/20):floor(points/4+points/20),k),I(floor(points/4-points/20):floor(points/4+points/20),k),1);
        pD=polyfit(V(floor(3.*points/4-points/20):floor(3.*points/4+points/20),k),I(floor(3.*points/4-points/20):floor(3.*points/4+points/20),k),1);
        Imod(:,k)=[I(1:floor(points/2),k)-polyval(pU,V(1:floor(points/2),k));I(floor(points/2)+1:points,k)-polyval(pD,V(floor(points/2)+1:points,k))];
        Iexc(k,:)=[pU(2),pD(2)];
        RN(k,:)=[1./pU(1),1./pD(1)];
    end
    
    
    handles.I=I;
    handles.V=V;
    handles.Imod=Imod;
    handles.Iexc=Iexc;
    handles.RN=RN;
    handles.flag=0;

    guidata(hObject, handles);
end
visualize(hObject,handles);
textupdate(hObject,handles);

function kslider_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function trange_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function range_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function crtlr_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function leadr_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ioffset_CreateFcn(hObject, eventdata, handles)
groundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes during object creation, after setting all properties.
function voffset_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function plotstyle_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes during object creation, after setting all properties.
function getoffset_CreateFcn(hObject, eventdata, handles)

function setslider(hObject, handles)
N=handles.N;
set(handles.kslider,'Max',N);
set(handles.kslider,'Min',1);
set(handles.kslider,'Value',1);
set(handles.kslider,'SliderStep',[1./N .05]);
% Update handles structure
guidata(hObject, handles);
% --- Text update function ---
function textupdate(hObject,handles)
ptype = get(handles.plotstyle,'Value');
switch ptype
    case {1,2}
        k=handles.k;
        Ij=handles.Ij;
        Iexc=handles.Iexc;
       handles.RN;
        
        set(handles.xcmaint,'String','Excess Current, [?A]');
        set(handles.xcrt,'String','Right Branch :');
        set(handles.xclt,'String','Left Branch :');
        set(handles.excl,'String',num2str(Iexc(k,1).*1e6,'%6.4g'));
        set(handles.excr,'String',num2str(Iexc(k,2).*1e6,'%6.4g'));
        
        set(handles.nrmaint,'String','Normal Resistence, [ohm] :');
        set(handles.rnrt,'String','Right Branch :');
        set(handles.rnlt,'String','Left Branch :');
        set(handles.rnl,'String',num2str(RN(k,1),'%6.4g'));
        set(handles.rnr,'String',num2str(RN(k,2),'%6.4g'));
        
        set(handles.injt,'String','Injection Current, [?A]:');
        set(handles.injcurr,'String',num2str(Ij(k).*1e6,'%6.4g'));
        
        set(handles.valuek,'String',num2str(k));
    case {3,4,5,6}
        k=handles.k;
         
        % changing indication values
        set(handles.xcmaint,'String',' ');
        set(handles.xcrt,'String',' ');
        set(handles.xclt,'String',' ');
        set(handles.excl,'String',' ');
        set(handles.excr,'String',' ');
        set(handles.nrmaint,'String',' ');
        set(handles.rnrt,'String',' ');
        set(handles.rnlt,'String',' ');
        set(handles.rnl,'String',' ');
        set(handles.rnr,'String',' ');
        set(handles.injt,'String',' ');
        set(handles.injcurr,'String',' ');
        set(handles.valuek,'String',num2str(k));
end        

guidata(hObject, handles);


function visualize(hObject,handles)
ptype = get(handles.plotstyle,'Value');
switch ptype
    case 1
        V=handles.V;
        G=handles.G;
        I=handles.I;
        trace=handles.trace;
        retrace=handles.retrace;
        
        k=handles.k;
        
        % Plotting IV
        subplot(1,2,1),
        plot(V(trace,k).*1e3,I(trace,k).*1e6,V(retrace,k).*1e3,I(retrace,k).*1e6);
        xlabel('Voltage, [mV]'); ylabel('Current, [\muA]'); grid on;
        axis([-max(max(1e3.*V)) max(max(1e3.*V)) -max(max(I.*1e6)) max(max(I.*1e6))]);
                
        % COnductance plot
        subplot(1,2,2),
        plot(V(trace,k).*1e3,G(trace,k),V(retrace,k).*1e3,G(retrace,k));
        xlabel('Voltage, [mV]'); ylabel('Conductance, [S]'); grid on;
        axis([-max(max(1e3.*V)) max(max(1e3.*V)) -max(max(G)) max(max(G))]);
        hold off   
    case 2
        V=handles.V;
        G=handles.G;
        I=handles.I;
        
        trace=handles.trace;
        retrace=handles.retrace;
        
        Iexc=handles.Iexc;
        RN=handles.RN;
        Imod=handles.Imod;
        
        points=handles.points;
        
        k=handles.k;
        
        % Linear interpolation IV plot
        subplot(1,2,1),
        plot(V(trace,k).*1e3,I(trace,k).*1e6,V(retrace,k).*1e3,I(retrace,k).*1e6), hold on;
        %plot((I(1:floor(points/4),k)-Iexc(k,1)).*RN(k,1).*1e3,I(1:floor(points/4),k).*1e6,'r-.',...
        %	(I(floor(3*points/4):points,k)-Iexc(k,2)).*RN(k,2).*1e3,I(floor(3*points/4):points,k).*1e6,'r-.'), hold off;
        plot(V(1:floor(points/4),k).*1e3,(V(1:floor(points/4),k)./RN(k,1)+Iexc(k,1)).*1e6,'r-.',...
            V(floor(3*points/4):points,k).*1e3,(V(floor(3*points/4):points,k)./RN(k,2)+Iexc(k,2)).*1e6,'r-.'), hold off;
        xlabel('Voltage, [mV]'); ylabel('Current, [\muA]'); grid on;
        axis([-max(max(1e3.*V)) max(max(1e3.*V)) -max(max(I.*1e6)) max(max(I.*1e6))]);
           
        % Excess current plot
        subplot(2,2,2),
        plot(V(trace,k).*1e3,Imod(trace,k).*1e6,V(retrace,k).*1e3,Imod(retrace,k).*1e6)
        axis([min(min(V.*1e3)),max(max(V.*1e3)),min(min(Imod.*1e6)),max(max(Imod.*1e6))])
        ylabel('Supercurrent [\muA]'); xlabel('Squid Voltage [mV]'); grid on;
                
        % Conductance plot
        subplot(2,2,4),
        plot(V(trace,k).*1e3,G(trace,k),V(retrace,k).*1e3,G(retrace,k));
        xlabel('Voltage, [mV]'); ylabel('Conductance, [S]'); grid on;
        axis([-max(max(1e3.*V)) max(max(1e3.*V)) -max(max(G)) max(max(G))]);
        
     case 3
        V=handles.V;
        I=handles.I;
        Ijm=handles.Ijm;
        trace=handles.trace;
        retrace=handles.retrace;
         
        k=handles.k;
         
        subplot(1,2,1),contour(Ijm(trace,:)*1e6,V(trace,:).*1e3,I(trace,:).*1e6,k,'LineWidth',2),colorbar,colormap jet;
        ylabel('Squid Voltage [mV]'); xlabel('Injection Current [\muA]'); grid on;
        title('Current Contour :: TRACE')
        subplot(1,2,2),contour(Ijm(retrace,:)*1e6,V(retrace,:).*1e3,I(retrace,:).*1e6,k,'LineWidth',2),colorbar,colormap jet;
        ylabel('Squid Voltage [mV]'); xlabel('Injection Current [\muA]'); grid on;
        title('Current Contour :: RETRACE')
    case 4
        V=handles.V;
        I=handles.I;
        Ijm=handles.Ijm;
        trace=handles.trace;
        retrace=handles.retrace;
        
        k=handles.k;
        subplot(1,2,1),contour(Ijm(trace,:)*1e6,I(trace,:).*1e6,V(trace,:).*1e3,k,'LineWidth',2),colorbar,colormap jet;
        ylabel('Squid Courrent [\muA]'); xlabel('Injection Current [\muA]'); grid on;
        title('Voltage Contour :: TRACE')
        subplot(1,2,2),contour(Ijm(retrace,:)*1e6,I(retrace,:).*1e6,V(retrace,:).*1e3,k,'LineWidth',2),colorbar,colormap jet;
        ylabel('Squid Courrent [\muA]'); xlabel('Injection Current [\muA]'); grid on;
        title('Voltage Contour :: RETRACE')
     case 5
        V=handles.V;
        G=handles.G; G=10.*log10(abs(G));
        I=handles.I;
        Ijm=handles.Ijm;
        trace=handles.trace;
        retrace=handles.retrace;
        subplot(1,2,1),
        pcolor(Ijm(trace,:)*1e6,I(trace,:).*1e6,G(trace,:)), shading interp, colorbar,colormap gray; %#ok<DUALC>
        title('Trace Logaritmic Conductance [dBS]')
        ylabel('Squid Current [\muA]'); xlabel('Injection Current [\muA]');
        subplot(1,2,2),
        pcolor(Ijm(retrace,:)*1e6,I(retrace,:).*1e6,G(trace,:)), shading interp, colorbar,colormap gray; %#ok<DUALC>
        ylabel('Squid Current [\muA]'); xlabel('Injection Current [\muA]');
        title('Retrace Logaritmic Conductance [dBS]')
        
     case 6
        V=handles.V;
        G=handles.G; G=10.*log10(abs(G));
        I=handles.I;
        Ijm=handles.Ijm;
        trace=handles.trace;
        retrace=handles.retrace;
        subplot(1,2,1),
        pcolor(Ijm(trace,:)*1e6,V(trace,:).*1e3,G(trace,:)), shading interp, colorbar,colormap gray; 
        ylabel('Squid Voltage [mV]'); xlabel('Injection Current [\muA]');
        title('Trace Logaritmic Conductance [dBS]')
        subplot(1,2,2),
        pcolor(Ijm(retrace,:)*1e6,V(retrace,:).*1e3,G(trace,:)), shading interp, colorbar,colormap gray;
        ylabel('Squid Voltage [mV]'); xlabel('Injection Current [\muA]');
        title('Retrace Logaritmic Conductance [dBS]')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- File import function ---
function [I, V, B]=fileimp(newdir)
current_dir=pwd;
cd(newdir);
s='*[*]';
filelist=dir(s);
m=length(filelist);
Bus=zeros(m,1);
for k=1:m;
    A=dlmread(filelist(k).name,'\t',1,0);
    gV=find(abs(A(:,2))>100);
    A(gV,2)=0; A(gV,1)=0;
    gI=find(abs(A(:,1))>100);
    A(gI,2)=0; A(gI,1)=0;
    Vus(:,k)=A(:,2); Ius(:,k)=A(:,1); %#ok<AGROW>
    istart=findstr(filelist(k).name,'[');
    istop=findstr(filelist(k).name,']');
    Bus(k)=str2double(filelist(k).name(istart+1:istop-1));
end;
[Y,In]=sort(Bus);
V=Vus(:,In); I=Ius(:,In);
B=Bus(In);
cd(current_dir);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Conductance extraction function ---
function G=LSQ(I,V,range)
range=range-mod(range,2)+1;
[points,N]=size(I);
G=zeros(points,N);
w=(range-1)./2;
uni=ones(range,1);
for h=1:N
    v=V(:,h);  i=I(:,h);
    for k=1:points
        inizio=mod(k-w,points)+points.*(k-w==0);
        fine=mod(k+w,points)+points.*(k+w==0);
        if(fine-inizio)> 0
                window=inizio:fine;
        end
        if(fine-inizio)< 0
                window=[1:fine,inizio:points];
        end
        if(fine-inizio)==0
                window=k;
        end
        Svi=v(window)'*i(window);       Svv=v(window)'*v(window);
        Sv=uni'*v(window);              Si=uni'*i(window);
        G(k,h)=(range.*Svi-Sv.*Si)./(range.*Svv-Sv.*Sv);
    end
end
