help:
	@echo "Available makes:"
	@echo "manuscript - compile the latex files to build the manuscript"
	@echo "figs - rename, convert to tiff and set colorpace, dpi and size of figures"
	@echo "clean - remove unused files"

all:
	make manuscript
	make figs
	make clean

# Compile the manuscript
manuscript:
	pdflatex manuscript.tex
	bibtex references
	pdflatex manuscript.tex
	pdflatex manuscript.tex

# Convert and rename figures
figs:
	convert Fig/Figure1.png -colorspace CMYK -units PixelsPerInch Fig/Figure1.tiff
	convert Fig/Figure2.png -colorspace CMYK -units PixelsPerInch -resample 600 -resize 2000 Fig/Figure2.tiff
	convert Fig/Figure3.png -colorspace CMYK -units PixelsPerInch -resample 600 -resize 2000 Fig/Figure3.tiff
	convert Fig/Figure4.png -colorspace CMYK -units PixelsPerInch -resample 600 -resize 2000 Fig/Figure4.tiff
	convert Fig/Figure5.png -colorspace CMYK -units PixelsPerInch -resample 600 -resize 2000 Fig/Figure5.tiff
	convert Fig/Figure6.png -colorspace CMYK -units PixelsPerInch -resample 600 -resize 2000 Fig/Figure6.tiff
	convert Fig/Figure7.png -colorspace CMYK -units PixelsPerInch -resample 600 -resize 2000 Fig/Figure7.tiff
	convert Fig/Figure8.png -colorspace CMYK -units PixelsPerInch -resample 600 -resize 2000 Fig/Figure8.tiff
	convert Fig/Figure9.png -colorspace CMYK -units PixelsPerInch -resample 600 -resize 2000 Fig/Figure9.tiff
	convert Fig/Figure10.png -colorspace CMYK -units PixelsPerInch -resample 600 -resize 2000 Fig/Figure10.tiff
	convert Fig/Figure11.png -colorspace CMYK -units PixelsPerInch -resample 600 -resize 2000 Fig/Figure11.tiff
	convert Fig/Figure12.png -colorspace CMYK -units PixelsPerInch -resample 600 -resize 2000 Fig/Figure12.tiff
	convert Fig/Figure13.png -colorspace CMYK -units PixelsPerInch -resample 600 -resize 2000 Fig/Figure13.tiff
	convert Fig/Figure14.png -colorspace CMYK -units PixelsPerInch -resample 600 -resize 2000 Fig/Figure14.tiff

	

# Remove unused files
clean:
	rm *.aux *.bbl *.blg *.log *.pdf *.lof *.lot *.synctex.gz *.out

