% Punto 1

Gx=7.8*10^-3; %Gradiente di campo [T/m]
Gz=12.8*10^-3; %Gradiente di campo [T/m]

B0=3; %Campo statico
BW=18*10^3; %Banda del segnale di selezione
n=11; %Numero avvolgimenti spira

ss=18; %Slice Selection

Lx=7.2*10^-2; %[m]
Ly=Lx;
Px=600; %Pixel 
Py=Px; 
h=1.055*10^-34; %Costante di Planck [J/s]
gamma=42.62*2*pi*10^6; %Rapporto giromagnetico di un protone di idrogeno [Hz/T]
f0=gamma*B0/(2*pi); %Frequenza centrale [Hz]
BWx=Gx*Lx*(gamma/(2*pi)); %Banda del segnale lungo x [Hz]

delta_x=Lx/Px; %Dimensione della cella lungo x [m]
delta_y=delta_x;%Dimensione della cella lungo y [m]
delta_z=(2*pi*BW)/(gamma*Gz);%Dimensione della cella lungo z [m]

%I risultati ottenuti sulla dimensione della cella sono:
%delta_x=delta_y=1.2e-04 m
%delta_z=0.033 m

% Il sistema di riferimento statico per tutte le immagini MRI prevede l’asse z allineato con
% la lunghezza del corpo umano, l’asse x corrispondente a destra–sinistra, e quello y antero–
% posteriore.
% 
% Il valore di ∆z corrisponde alla risoluzione della immagine lungo z, che
% solitamente è dell'ordine di grandezza di qualche mm. Di conseguenza, il
% risultato ottenuto può essere ritenuto accettabile.

% Punto 2

% Segmentare le immagini biomedicali fornite nel proprio teams con gli algoritmi:
% -sogliatura e metodo di Otsu;
% -clustering, usando l’algoritmo k-means;
% -l’algoritmo di Region Growing.
% Analizzare e spiegare i risultati ottenuti e le differenze con le diverse metodologie

% Caricamento e visualizzazione delle immagini originali con relativo
% istogramma

I132=dicomread('132.dcm');
I133=dicomread('133.dcm');
I134=dicomread('134.dcm');

figure
subplot(2,3,1)
imshow(I132,'DisplayRange',[])
title('Immagine 132 originale')

subplot(2,3,2)
imshow(I133,'DisplayRange',[])
title('Immagine 133 originale')

subplot(2,3,3)
imshow(I134,'DisplayRange',[])
title('Immagine 134 originale')

subplot(2,3,4)
imhist(I132)

subplot(2,3,5)
imhist(I133)

subplot(2,3,6)
imhist(I134)

% Equalizzazione per migliorare le immagini da segmentare e visualizzazione
% delle immagini equalizzate con relativo istogramma

J132=histeq(I132);
J133=histeq(I133);
J134=histeq(I134);

figure
subplot(2,3,1)
imshow(J132,'DisplayRange',[])
title('Immagine 132 equalizzata')

subplot(2,3,2)
imshow(J133,'DisplayRange',[])
title('Immagine 133 equalizzata')

subplot(2,3,3)
imshow(J134,'DisplayRange',[])
title('Immagine 134 equalizzata')

subplot(2,3,4)
imhist(J132)

subplot(2,3,5)
imhist(J133)

subplot(2,3,6)
imhist(J134)

%Sogliatura e metodo di Otsu con segmentazione manuale

figure('Name','Segmentazione manuale')

J1 = J132;
subplot(1,3,1)
BW = roipoly(J1);
J1(BW==0)=1;
imshowpair(J132,J1)
title('Immagine 132')

J2 = J133;
subplot(1,3,2)
BW = roipoly(J2);
J2(BW==0)=1;
imshowpair(J133,J2)
title('Immagine 133')

J3 = J134;
subplot(1,3,3)
BW = roipoly(J3);
J3(BW==0)=1;
imshowpair(J134,J3)
title('Immagine 134')

% Metodo di Otsu con funzione graythresh (alternativa alla segmentazione
% manuale)

figure('Name','Metodo Otsu con graythresh')

level = graythresh(J132);
BW = imbinarize(J132,level);
subplot(1,3,1)
imshowpair(J132,BW)
title('Immagine 132')

level = graythresh(J133);
BW = imbinarize(J133,level);
subplot(1,3,2)
imshowpair(J133,BW)
title('Immagine 133')

level = graythresh(J134);
BW = imbinarize(J134,level);
subplot(1,3,3)
imshowpair(J134,BW)
title('Immagine 134')

% Clustering, usando l'algoritmo k-means

Nclass = 10; % Numero di classi scelto in base all'aspetto dell'immagine

figure('Name','Clustering - k-Means')
subplot(1,3,1)
L = imsegkmeans(J132,Nclass);
B = labeloverlay(J132,L);
imshow(B)
title('Immagine 132')

Nclass2=11;  
subplot(1,3,2)
L = imsegkmeans(J133,Nclass2);
B = labeloverlay(J133,L);
imshow(B)
title('Immagine 133')

subplot(1,3,3)
L = imsegkmeans(J134,Nclass2);
B = labeloverlay(J134,L);
imshow(B)
title('Immagine 134')

% Algoritmo di Region Growing I132

% Definizione valori soglia
T1 = 100;  % Seed threshold
T2 = 200;  % Growing threshold

% Trovare tutti i seed points e creare una maschera binaria dove i seeds
% sono considerati come l’oggetto. Visualizzare l’immagine originale con i
% seeds rappresentati in bianco.
seeds = find(I132 >= T1);

% Vettore ToProcess e maschera di segmentazione K2
ToProcess = seeds;
% ToProcess_u=ToProcess;
K2 = zeros(size(I132));
% in alternativa K2(seeds)=1; 
% Ciclo finchè ci sono ancora elementi nel vettore ToProcess

[Nrows,Ncols]=size(I132);
while isempty(ToProcess)==0
    
    % Inserimento del pixel da processare in questo ciclo
    current=ToProcess(1);
    
    % Aggiornamento della maschera di segmentazione
    K2(current)=1;
    
    % Insieme 8-connesso del pixel corrente 
    neighbors8 = conneigh8(Nrows,Ncols,current);
    
    % Controllo se i pixel dell'insieme 8-connesso sono tra T1 e T2, in
    % caso positivo sono seed, quindi li metto in coda al vettore ToProcess
    A = neighbors8(I132(neighbors8)<T1);
    B = neighbors8(I132(neighbors8)>=T2);
    C = neighbors8(K2(neighbors8)~=1);
    ToProcess=[ToProcess; intersect(intersect(A,B),C)'];
    K2(intersect(A,B))=1;
    %ToProcess_u=[ToProcess_u; neighbors8];
    
    % Elimino il primo elemento del vettore ToProcess (cioè il pixel appena analizzato in questo ciclo)
    ToProcess(1)=[];
        
    % Rimozione degli elementi ripetuti all'interno del vettore ToProcess (i.e. unique())    
    ToProcess=unique(ToProcess);
    
     
end

% Immagine iniziale con vasi segmentati 
I_rg(K2==1)=1;

figure
% imshow(I132);
imshowpair(I132,I_rg)
title('Region Growing Segmentation: I132')


%Region Growing I133

seeds = find(I133 >= T1);

ToProcess = seeds;
% ToProcess_u=ToProcess;
K2 = zeros(size(I133));
% K2(seeds)=1;

[Nrows,Ncols]=size(I133);
while isempty(ToProcess)==0

    current=ToProcess(1);
   
    K2(current)=1;
    
    neighbors8 = conneigh8(Nrows,Ncols,current);
   
    A = neighbors8(I133(neighbors8)<T1);
    B = neighbors8(I133(neighbors8)>=T2);
    C = neighbors8(K2(neighbors8)~=1);
    ToProcess=[ToProcess; intersect(intersect(A,B),C)'];
    K2(intersect(A,B))=1;
    % ToProcess_u=[ToProcess_u; neighbors8];
 
    ToProcess(1)=[];

    ToProcess=unique(ToProcess);
end

% Immagine iniziale con vasi segmentati 
I_rg(K2==1)=1;

figure
% imshow(I33);
imshowpair(I133,I_rg)
title('Region Growing Segmentation: I133')
    
%Region Growing I134

seeds = find(I134 >= T1);

ToProcess = seeds;
% ToProcess_u=ToProcess;
K2 = zeros(size(I134));
% K2(seeds)=1;

[Nrows,Ncols]=size(I134);
while isempty(ToProcess)==0
   
    current=ToProcess(1); 
    
    K2(current)=1;
     
    neighbors8 = conneigh8(Nrows,Ncols,current);
   
    A = neighbors8(I134(neighbors8)<T1);
    B = neighbors8(I134(neighbors8)>=T2);
    C = neighbors8(K2(neighbors8)~=1);
    ToProcess=[ToProcess; intersect(intersect(A,B),C)'];
    K2(intersect(A,B))=1;
    % ToProcess_u=[ToProcess_u; neighbors8];
    
    ToProcess(1)=[];
       
    ToProcess=unique(ToProcess);
    
     
end

% Immagine iniziale con vasi segmentati 
I_rg(K2==1)=1;

figure
% imshow(I134);
imshowpair(I134,I_rg)
title('Region Growing Segmentation: I134')




