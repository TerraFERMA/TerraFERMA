execute_process(COMMAND rm tftutorials.aux tftutorials.toc tftutorials.out tftutorials.pdf tftutorials.log tftutorials.blg tftutorials.bbl)
execute_process(COMMAND pdflatex tftutorials.tex)
execute_process(COMMAND pdflatex tftutorials.tex)
execute_process(COMMAND bibtex tftutorials.aux)
execute_process(COMMAND bibtex tftutorials.aux)
execute_process(COMMAND pdflatex tftutorials.tex)
execute_process(COMMAND pdflatex tftutorials.tex)


