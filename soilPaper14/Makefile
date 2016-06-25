# modified from http://robjhyndman.com/researchtips/makefiles/

# Usually, only these lines need changing
TEXFILE= soil_main
PDFFILE = Wutzler14_seam
RDIR= ../inst/soilPaper14R
FIGDIR= ./fig
# Don't put any pdf in this directory what should not be cropped

# list R files
RFILES := $(wildcard $(RDIR)/*.R)

# Indicator files to show R file has run
OUT_FILES:= $(RFILES:.R=.Rout)

# pdf figures created by R
PDFFIGS := $(wildcard $(FIGDIR)/*f)

# Indicator files to show pdfcrop has run
CROP_FILES:= $(PDFFIGS:.pdf=.pdfcrop)

all: $(PDFFILE).pdf $(OUT_FILES) $(CROP_FILES)

# May need to add something here if some R files depend on others.

# RUN EVERY R FILE
#   R CMD BATCH '--args paperBGC' $<
#   need to run from projects root directory, then need to change the filename that includes the path (remove ../)
$(RDIR)/%.Rout: $(RDIR)/%.R $(RDIR)/functions.R
	$ cd ..; echo $< | sed 's/\.\.\///' | xargs -n 1 R CMD BATCH --vanilla '--args paperBGC'

    
# CROP EVERY PDF FIG FILE
$(FIGDIR)/%.pdfcrop: $(FIGDIR)/%.pdf
	#pdfcrop --margins 5 $< $< && touch $@
	pdfcrop $< $< && touch $@

# copy the current bib-database from z and run latex several times
# pdflatex main -job-name=Wutzler13_ModelError -output-directory=tmp
# bibtex main 
#	pdflatex -quiet $(TEXFILE) -job-name=$(PDFFILE) -output-directory=tmp
#	cd tmp; bibtex  $(PDFFILE) 
#	cp tmp/$(PDFFILE).bbl $(TEXFILE).bbl
#	pdflatex -quiet $(TEXFILE) -job-name=$(PDFFILE) -output-directory=tmp
#	pdflatex -quiet $(TEXFILE) -job-name=$(PDFFILE) -output-directory=tmp
#	mv tmp/$(PDFFILE).pdf .
bibtex:
	cp /cygdrive/m/people/twutz/journals/twutz_txt.bib tmp
	pdflatex -quiet $(TEXFILE) -output-directory=tmp
	cd tmp; bibtex  $(TEXFILE) 
	pdflatex -quiet $(TEXFILE) -output-directory=tmp
	pdflatex -quiet $(TEXFILE) -output-directory=tmp
	mv tmp/$(TEXFILE).pdf ./$(PDFFILE).pdf
	
# Compile main tex file
$(PDFFILE).pdf: $(TEXFILE).tex $(OUT_FILES) $(CROP_FILES) *.tex
	pdflatex -quiet $(TEXFILE) -output-directory=tmp
	mv tmp/$(TEXFILE).pdf $(PDFFILE).pdf

# Run R files
R: $(OUT_FILES)

# prepare graphics
# R scripts should generate figures in subdir /fig
# Figures that are not derived should stay in subdir figOrig
figs: $(PDFFIGS)
	cp figOrig/*.pdf fig    # copy figOrig to fig compile dir
	
# crop pdf figures	
cropfigs: figs $(CROP_FILES)

# View main tex file
view: $(PDFFILE).pdf
	sumatra $(PDFFILE).pdf &

# Clean up stray files
clean:
	rm -fv $(OUT_FILES)
	rm -fv $(CROP_FILES)
	rm -fv $(FIGDIR)/*
	rm -fv *.aux *.log *.toc *.blg *.bbl *.synctex.gz *.out *.bcf *blx.bib *.run.xml
	rm -fv *.fdb_latexmk *.fls
	rm -fv $(TEXFILE).pdf

.PHONY: all clean

