.PHONY: report.pdf all clean

# All.
all: report.pdf

%.tex: %.raw
	@./raw2tex $< > $@

%.tex: %.dat
	@./dat2tex $< > $@

# Report.
report.pdf: report.tex
	@latexmk -pdf -pdflatex="pdflatex -interaction=nonstopmode" -use-make report.tex

# Clean.
clean:
	@latexmk -CA
	@$(RM) -r *.bbl
	@$(RM) -r *.xml