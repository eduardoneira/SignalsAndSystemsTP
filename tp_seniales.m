#TP de seniales

pkg load signal

a = 1; 										#para que no piense que es un archivo de funcion

function espectrograma(y,fs,ms)
	step = fix((fs*ms)/8);
	window = fix(fs*ms);	  				# 2 ms data window
	fftn = 2^nextpow2(window);

	printf("Usando ventana de %d ms\n",ms*1000);
	printf("params\n\tfftn: %d\n\tfs: %d\n\twindow: %d\n\tstep: %d\n",fftn,fs,window,step);

	graphics_toolkit gnuplot;
	specgram(y, fftn, fs, hanning(window),step);

endfunction

function ejercicio1(y,fs)

	# Espectograma de banda ancha
	printf("ESPECTOGRAMA DE BANDA ANCHA\n");
	ms = 2/1000;
	filename = sprintf("ej1_espectograma_banda_ancha_%dms.png",ms*1000);
	espectrograma(y,fs,ms);

	set (gca, "xlim", [0, 40]);
	print(filename,"-dpng");

	# Espectograma de banda angosta
	printf("ESPECTOGRAMA DE BANDA ANGOSTA\n");
	ms = 160/1000;
	filename = sprintf("ej1_espectograma_banda_angosta_%dms.png",ms*1000);
	espectrograma(y,fs,ms);

	set (gca, "xlim", [0, 40]);
	print(filename,"-dpng");

endfunction

function ejercicio2(y,fs)
	printf("Ejercicio 2: Espectograma de notas\n");
	ms = 100/1000;
	filename = sprintf("ej2_espectograma_banda_angosta_%dms.png",ms*1000);
	espectrograma(y,fs,ms);

	set (gca, "ylim", [500, 2000]);
	set (gca, 'ytick', 500:25:2000);
	set (gca, "xlim", [0, 40]);
	print(filename,"-dpng");

endfunction

function new_y = hold_orden_cero(y)
  	r = rows(y);

	size = (2*r)-1;

	new_y = zeros(size,1);
	new_y(1:2:size) = y(:);
	new_y(2:2:size) = y(1:r-1);
 
  
endfunction

function new_y = hold_orden_uno(y)
  	r = rows(y);

	size = (2*r)-1;

	new_y = zeros(size,1);
	new_y(1:2:size) = y(:);
	new_y(2:2:size) = (y(1:r-1) + y(2:r)) /2; 
  
endfunction

function ret = sinc_windowed(samples,fs,fc)
	t = (-samples/2:samples/2 -1);
	ws = (fc/fs);
	_sinc = 2 * ws * sinc(2* ws *t);

	han = blackman(length(_sinc));

	_sinc(:) = _sinc(:) .* han(:);

	#plot(t, _sinc);
  	#print("sinc(t).png","-dpng");
  	#_sinc_fft = fftshift(abs(fft(_sinc)));
  
	#_sinc_f = fs*(-samples/2:samples/2-1)/samples;
	#plot(_sinc_f, _sinc_fft);
	#print("sinc_fft.png","-dpng");

	ret = _sinc;
endfunction

function new_y = interpolar_zeros_y_filtrar(y,fs) 
	r = rows(y);

	size = (2*r)-1;

	y_0=zeros(size,1);
	y_0(1:2:size) = y(:);

	_sinc = sinc_windowed(length(y_0),fs,4000);

	y_filtrada = fftconv(y_0,_sinc);
  
  	new_y = y_filtrada(length(y_filtrada)/4:(3/4)*length(y_filtrada));

endfunction

function ejercicio3(y,fs)

	new_y = hold_orden_cero(y);

	ms = 100/1000;
	espectrograma(new_y,fs,ms);

	set(gca, "xlim", [0, 80]);
	print("espectrograma_hold0","-dpng");

	new_y = hold_orden_uno(y);

	espectrograma(new_y,fs,ms);

	set(gca, "xlim", [0, 80]);
	print("espectrograma_hold1","-dpng");

	new_y = interpolar_zeros_y_filtrar(y,fs);

	espectrograma(new_y,fs,ms);

	set(gca, "xlim", [0, 80]);
	print("espectrograma_interpolar_0_y_sinc","-dpng");

	audiowrite("audio_ej3.wav", new_y, fs);

endfunction

function ejercicio4(y,fs)
	
	filter = sinc_windowed(length(y),fs,4000);
	y_filtrada = fftconv(y,filter);
	y_filtrada = y_filtrada(length(y_filtrada)/4:(3/4)*length(y_filtrada));

	new_y = zeros(ceil(length(y)/2),1);
	new_y(:) = y_filtrada(1:2:length(y_filtrada));

	ms = 100/1000;
	espectrograma(new_y,fs,ms);

	set(gca, "xlim", [0, 20]);
	print("espectrograma_decimacion","-dpng");

	audiowrite("audio_ej4.wav", new_y, fs);
endfunction

#Comienza programa

[y, fs] = audioread("Audio.wav");
printf ("Duracion = %.2f s\n", rows (y)/fs);
printf ("Sampling rate = %.2f Hz\n", fs);

ejercicio4(y,fs);
