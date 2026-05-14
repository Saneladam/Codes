all: figures pdf visualizer

figures: tex/figures/poisson.pdf

tex/figures/poisson.pdf: code/poisson.py
	python code/poisson.py

pdf:
	cd tex && pdflatex main.tex
	cd tex && bibtex main
	cd tex && pdflatex main.tex
	cd tex && pdflatex main.tex

visualizer:
	zathura tex/main.pdf


