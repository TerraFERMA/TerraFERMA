execute_process(
     COMMAND pdflatex tftutorials.tex && pdflatex tftutorials.tex && bibtex tftutorials.aux && bibtex tftutorials.aux && pdflatex
tftutorials.tex && pdflatex tftutorials.tex
     )


