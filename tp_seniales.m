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

function ejercicio3(y,fs)
	r = rows(y);

	size = (2*r)-1;

	new_y = zeros(size,1);
	new_y(1:2:size) = y(:);
	new_y(2:2:size) = y(1:r-1);

	ms = 100/1000;
	espectrograma(new_y,fs,ms);

	set(gca, "xlim", [0, 80]);
	print("espectrograma_hold0","-dpng");

	espectrograma(y,fs,ms);

	set(gca, "xlim", [0, 40]);
	print("espectrograma_comun","-dpng");

	#falta filtrar


endfunction

#Comienza programa

[y, fs] = audioread("Audio.wav");
printf ("Duracion = %.2f s\n", rows (y)/fs);
printf ("Sampling rate = %.2f Hz\n", fs);

ejercicio3(y,fs);
