numgen: Random\ Number\ Generation/numgen.c Random\ Number\ Generation/plot_points.py
	cd Random\ Number\ Generation/ && \
	gcc numgen.c -o numgen.exec && \
	./numgen.exec && \
	python plot_points.py && \
	rm numgen.exec *.txt

montecarlo: Monte\ Carlo/montecarlo.c
	cd Monte\ Carlo/ && \
	gcc montecarlo.c -o montecarlo.exec && \
	./montecarlo.exec && \
	rm montecarlo.exec
