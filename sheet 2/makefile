stiffeq: Stiff\ Equation/stiffeq.c Stiff\ Equation/stiffeq.py
	cd Stiff\ Equation && \
	gcc stiffeq.c -o stiffeq.exec && \
	./stiffeq.exec && \
	python stiffeq.py && \
	rm stiffeq.exec *.dat

doublepend: Double\ Pendulum/doublepend.c Double\ Pendulum/error.py Double\ Pendulum/movie.py
	cd Double\ Pendulum && \
	gcc doublepend.c -o doublepend.exec && \
	./doublepend.exec && \
	python error.py && \
	mkdir -p Frame\ Buffer && \
	python movie.py && \
	cd Frame\ Buffer && \
		ffmpeg -r 24 -i frame_%03d.png -pix_fmt yuv420p ../double-pendulum.mp4 && \
	rm -rf Frame\ Buffer doublepend.exec *.dat

3: exercise_3.c
	gcc -o exercise_3 exercise_3.c
	./exercise_3
	python exercise_3cd.py
	mkdir -p Pendulum
	python exercise_3e.py
	cd Pendulum && \
		ffmpeg -r 24 -i frame_%03d.png -pix_fmt yuv420p ../double-pendulum.mp4
	rm -rf Pendulum, *.dat
