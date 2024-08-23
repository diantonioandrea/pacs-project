.PHONY: presentation.pdf all clean

# All.
all: presentation.pdf

%.tex: %.raw
	@./raw2tex $< > $@

%.tex: %.dat
	@./dat2tex $< > $@

# Presentation.
presentation.pdf: presentation.tex
	@latexmk -pdf -pdflatex="pdflatex -interaction=nonstopmode" -use-make presentation.tex

# Clean.
clean:
	@latexmk -CA
	@$(RM) -r *.bbl
	@$(RM) -r *.nav
	@$(RM) -r *.snm