fftconvol: FFT-based\ Convolution/fftconvol.c
	cd FFT-based\ Convolution && \
	gcc -lfftw3 fftconvol.c -o fftconvol.exec && \
	./fftconvol.exec && \
	rm fftconvol.exec

pmesh: Particle\ Mesh\ Force\ Law/pmesh.c
	cd Particle\ Mesh\ Force\ Law && \
	gcc -lfftw3 pmesh.c -o pmesh.exec && \
	./pmesh.exec && \
	python pmesh.py && \
	rm pmesh.exec *.dat
