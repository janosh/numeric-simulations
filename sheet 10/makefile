rejection: Rejection/rejection.c Rejection/rejection.py
	cd Rejection && \
	gcc rejection.c -o rejection.exec && \
	./rejection.exec && \
	python rejection.py && \
	rm rejection.exec *.txt

markovchain: Markov\ Chain/markovchain.c Markov\ Chain/markovchain.py
	cd Markov\ Chain && \
	gcc markovchain.c -o markovchain.exec && \
	./markovchain.exec && \
	python markovchain.py && \
	rm markovchain.exec *.txt

ising: Ising/ising.c Ising/ising.py
	cd Ising && \
	gcc ising.c -o ising.exec && \
	./ising.exec && \
	python ising.py && \
	rm ising.exec lattice*.txt

hardcoded: hardcoded.c hardcoded.py
	gcc hardcoded.c -o hardcoded
	./hardcoded
	python hardcoded.py
	rm *.txt
