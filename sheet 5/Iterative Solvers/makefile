all: jacobi gaussseidel cb_gaussseidel multigrid resdecay density potential
.PHONY : all

.NOTPARALLEL :

jacobi: Jacobi/jacobi.c
	cd Jacobi && \
	gcc jacobi.c -o jacobi.exec && \
	./jacobi.exec && \
	rm -f jacobi.exec


gaussseidel: Gauss-Seidel/gaussseidel.c
	cd Gauss-Seidel && \
	gcc gaussseidel.c -o gaussseidel.exec && \
	./gaussseidel.exec && \
	rm -f gaussseidel.exec


cb_gaussseidel: Chessboard\ Gauss-Seidel/cb_gaussseidel.c
	cd Chessboard\ Gauss-Seidel && \
	gcc cb_gaussseidel.c -o cb_gaussseidel.exec && \
	./cb_gaussseidel.exec && \
	rm -f cb_gaussseidel.exec


multigrid: Multigrid/multigrid.c
	cd Multigrid && \
	gcc multigrid.c -o multigrid.exec && \
	./multigrid.exec && \
	rm -f multigrid.exec

resdecay: res_decay.py
	python res_decay.py

density: density.py
	python density.py

potential: potential.py
	python potential.py

clean:
	find . -type f -name '*.dat' -delete
	find . -type f -name '*.txt' -delete
