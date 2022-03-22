MAIN=latex
WINDOWS=1

SVGDIR=svg
SVGOUTDIR=figs_svg

SVGTEXDIR=svgtex

SVGTEXOUTDIR=figs_svgtex
SVGTEXS=$(notdir $(wildcard $(SVGTEXDIR)/*.svg))
SVGTEXS2PDFS=$(patsubst %,$(SVGTEXOUTDIR)/%,$(SVGTEXS:.svg=.pdf))
SVGS=$(notdir $(wildcard $(SVGDIR)/*.svg))
SVGS2PDFS=$(patsubst %,$(SVGOUTDIR)/%,$(SVGS:.svg=.pdf))



INKSCAPE=inkscape
RM=rm -f
MV=mv
SLASH="/"
# Windows specific
ifeq ($(WINDOWS),1)
   INKSCAPE="C:/Program Files/Inkscape/inkscape.exe"
   RM=del
   MV=move /Y 
   SLASH="\\"
else
   INKSCAPE=inkscape
   RM=rm
   MV=mv
   SLASH="/"
endif

all: figspdf pdf

pdf: figspdf
	pdflatex -synctex=1 --file-line-error-style --shell-escape --interaction=nonstopmode $(MAIN)
	bibtex $(MAIN).aux

bibtex:
	bibtex $(MAIN).aux

figspdf: $(SVGOUTDIR) $(SVGTEXOUTDIR) $(SVGS2PDFS) $(SVGTEXS2PDFS)

$(SVGTEXOUTDIR)/%.pdf:$(SVGTEXDIR)/%.svg 
	$(INKSCAPE) -z -D --file="$<" --export-pdf="$(SVGTEXOUTDIR)/$*.pdf" --export-latex
	$(MV) "$(SVGTEXOUTDIR)$(SLASH)$*.pdf_tex" "$(SVGTEXOUTDIR)$(SLASH)$*.tex"

$(SVGOUTDIR)/%.pdf: $(SVGDIR)/%.svg
	$(INKSCAPE) -D -A "$@" "$<"
 
clean:
	$(RM) *.log *.aux *.out *.bbl *.blg *.toc
