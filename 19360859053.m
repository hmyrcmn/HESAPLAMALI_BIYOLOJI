sekans=getgenbank('NC_001416','SequenceOnly',true);%48.502 bazdan oluşan canlıya ait dizilim genbank adlı siteden elde edildi
HPA2='CCGG'; % CCGG RESTRİKCİONS ENZİMİ DİZİLİMDE KESME İŞLEMİNİ YAPACAK (KESİLECEK YERLERİ TANIYACAK OLAN ENZİM )
H=strfind(sekans,HPA2(1,:));

fragmanlar=cell(1,length(H)); %boş hücreler tanımlandı
uzunluk=zeros(1,length(H));  %elemanlar sıfır olarak tanımlandı 

for i=1:length(H)-1 
    fragmanlar(i)={sekans(H(i)+3:H(i+1)-1)};
    uzunluk(i)=length(fragmanlar{i});
end
fragmanlar(end)={sekans(H(end)+3:end)};
uzunluk(end)=length(fragmanlar{end});



% [1, 100, 200, 300, 400, 500, 600, 48502]; 
%ifadesini elde etmek için kod blogu

fprintf( "intergal belirli aralıkları martis(1,8) olarak saklanıyor \n"); 

belirliintegral_araliklari=zeros(1,8);
i=1;
j=0;
while i<8
    belirliintegral_araliklari(i)=j;
    j=j+100; % j+=100 yaptım kabul etmedi ! 
    i=i+1;
  
end
belirliintegral_araliklari(8)=48502;% döngüdeki artık miktarına uymayan özel değer.
fprintf("*******integral aralık dizisi tanımlandı !! ************")
fprintf("\n \n********** Hesaplanan Degerler Sonuc Tablosu **********");



P=0.0039; 
integral_formul = @(x) P * exp(-P*x);   %integrali alınacak fonksiyon 

total_chi_square=0;
kac_kez=0;
a=1;
while a<8
    tahmin= integral(integral_formul,belirliintegral_araliklari(a), belirliintegral_araliklari(a+1)) * length(H);

    real=0;
    for k = uzunluk 
        if k >= belirliintegral_araliklari(a)
            if k < belirliintegral_araliklari(a+1)   
                kac_kez= kac_kez + 1;
                real = real + 1;
            end
        end
    end
    chi_square = (real - tahmin).^2 / tahmin ;      
    
     total_chi_square = total_chi_square + chi_square ;
     fprintf("\n")
     
    fprintf("%3.d ---%5.d |  %.3d  \t  |  %4.5f \t |    %4.6f  \t|\n", belirliintegral_araliklari(a), belirliintegral_araliklari(a+1), real, tahmin, chi_square);
   
    
   
    a=a+1;
end
       
    
     
    
