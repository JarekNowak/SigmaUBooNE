# Build the selection report. Figures are pulled from ../figs (regenerate with
# PlotSigma.C). Requires pdflatex (no siunitx dependency).
selection.pdf: selection.tex
	pdflatex -interaction=nonstopmode selection.tex
	pdflatex -interaction=nonstopmode selection.tex
clean:
	rm -f selection.aux selection.log selection.out build.log
