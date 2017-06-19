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

function x = istft(stft, wlen, step, nfft, fs)

	coln = size(stft, 2)
	xlen = wlen + (coln-1)*step;
	x = zeros(1, xlen);

	% synthesis window
	win = hanning(wlen);

	% initialize the signal time segment index
	indx = 0;

	% perform ISTFT (via IFFT and Weighted-OLA)
	if rem(nfft, 2)                     % odd nfft excludes Nyquist point
	    for col = 1:coln
	        % extract FFT points
	        X = stft(:, col);
	        X = [X; conj(X(end:-1:2))];
	        
	        % IFFT
	        xprim = real(ifft(X));
	        xprim = xprim(1:wlen);
	        
	        % weighted-OLA
	        x((indx+1):(indx+wlen)) = x((indx+1):(indx+wlen)) + (xprim.*win)';
	        
	        % update the index
	        indx = indx + step;
	    end
	else                                % even nfft includes Nyquist point
	    for col = 1:coln
	        % extract FFT points
	        X = stft(:, col);
	        X = [X; conj(X(end-1:-1:2))];
	        
	        % IFFT
	        xprim = real(ifft(X));
	        xprim = xprim(1:wlen);
	        
	        % weighted-OLA
	        x((indx+1):(indx+wlen)) = x((indx+1):(indx+wlen)) + (xprim.*win)';
	        
	        % update the index
	        indx = indx + step;
	    end
	end

	% scale the signal
	W0 = sum(win.^2);
	scaleFactor = step/W0                  
	x = x.*scaleFactor;                      

end


function ejercicio5(y,fs)

	ms = 100/1000;
	window = fix(fs*ms);	  				
	fftn = 2^nextpow2(window);
	hop = window/4;

	specgram(y,fftn,fs,hanning(window),window-hop);
	print("espectrograma_ej5_stft.png","-dpng");
	
	[s,f,t] = specgram(y,fftn,fs,hanning(window),window-hop);

	new_y = istft(s,window,hop,fftn,fs);

	error = 0;

	for i = 1:length(new_y)
		error+=(y(i) - new_y(i));
	endfor

	printf("Error cuadrático medio : %.2f\n",error);
	specgram(new_y,fftn,fs,hanning(window),window-hop);
	print("espectrograma_ej5_istft.png","-dpng");

	audiowrite("audio_ej5.wav", new_y, fs);
endfunction

function ejercicio6(y,fs)
	ms = 100/1000;
	window = fix(fs*ms);	  				
	fftn = 2^nextpow2(window);
	hop = window/4;

	[s,f,t] = specgram(y,fftn,fs,hanning(window),window-hop);
	#nueva matriz
	new_s = zeros(rows(s),2*columns(s)-1);
	new_s(:,1:2:columns(new_s)) = s(:,:);

	fc = 4000;	

	for i = 1:rows(new_s)
		[s1,f1,t1] = specgram(new_s(i,:),fftn,fs,hanning(window),window-hop);
		min = 10000;
		found = 0;

		for j=1:rows(s1)
			aux = sum(abs(s1(j,:)));
			if (aux < min)
				min = aux;
				found = j;
			endif
		endfor

		fc = 8000*found/rows(s1);
		_sinc = sinc_windowed(columns(new_s),fs,fc);

		row = fftconv(new_s(i,:),_sinc);

		row_no_padding = row(fix(length(row)/4):fix((3/4)*length(row)));
		new_s(i,:) = row_no_padding(:);

	endfor

	new_y = istft(new_s,window,hop,fftn,fs);

	specgram(new_y,fftn,fs,hanning(window),window-hop);
	print("espectrograma_ej6.png","-dpng");

	audiowrite("audio_ej6_mejor.wav",new_y,fs);	

endfunction

function ejercicio7(y,fs)
	ms = 100/1000;
	window = fix(fs*ms);	  				
	fftn = 2^nextpow2(window);
	hop = window/4;

	[s,f,t] = specgram(y,fftn,fs,hanning(window),window-hop);
	#nueva matriz
	new_s = zeros(rows(s),2*columns(s)-1);
	new_s(:,1:2:columns(new_s)) = s(:,:);	

	for i = 1:rows(new_s)
		filename = sprintf("ej6/plot_row_%d_fft.png",i);
		specgram(s(i,:),fftn,fs,hanning(window/2),window/2-hop);
		print(filename,"-dpng");

		filename = sprintf("ej6_interpolado/plot_row_%d_fft.png",i);
		specgram(new_s(i,:),fftn,fs,hanning(window/2),window/2-hop);
		print(filename,"-dpng");
		
		[s1,f1,t1] = specgram(new_s(i,:),fftn,fs,hanning(window),window-hop);
		min = 10000;
		found = 0;

		for j=1:rows(s1)
			aux = sum(abs(s1(j,:)));
			if (aux < min)
				min = aux;
				found = j;
			endif
		endfor

		fc = 8000*found/rows(s1);

		_sinc = sinc_windowed(columns(new_s),fs,fc);

		row = fftconv(new_s(i,:),_sinc);

		row_no_padding = row(fix(length(row)/4):fix((3/4)*length(row)));
		new_s(i,:) = row_no_padding(:);

	endfor

	new_y = istft(new_s,window,hop,fftn,fs);

	specgram(new_y,fftn,fs,hanning(window),window-hop);
	print("espectrograma_ej6_sin_filtrar.png","-dpng");

	audiowrite("audio_ej6.wav",new_y,fs);
endfunction


#Comienza programa

[y, fs] = audioread("Audio.wav");
printf ("Duracion = %.2f s\n", rows (y)/fs);
printf ("Sampling rate = %.2f Hz\n", fs);

ejercicio5(y,fs);