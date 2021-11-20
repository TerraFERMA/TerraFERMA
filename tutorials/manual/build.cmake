execute_process(COMMAND latexmk -C)
execute_process(COMMAND latexmk -pdf -interaction=nonstopmode tftutorials.tex)


