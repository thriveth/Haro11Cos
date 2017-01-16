#!/bin/sh
# Copy text files
# cp Paper/Paper.tex ToSubmit/RiveraThorsen.tex
cp Paper/main.bib ToSubmit/
cp Paper/table1.tex ToSubmit/
# Copy figures
cp Figs/Haroslit.pdf ToSubmit/
cp Figs/FullSpec.pdf ToSubmit
cp Figs/EffResol.pdf ToSubmit/
cp Figs/HISLISProfiles.pdf ToSubmit/
cp Figs/AOD-details-example.pdf ToSubmit/
cp Figs/Fc_haro11c_LIS.pdf ToSubmit/
cp Figs/Fc_haro11c-SiIV.pdf ToSubmit/
cp Figs/LyACoverfracs.pdf ToSubmit/
cp Figs/1260-fit-fluor.pdf ToSubmit/
cp Figs/CIIdoubletmodel.pdf ToSubmit/
# Correct figure paths
sed s,../Figs/,./, Paper/Paper.tex > ToSubmit/RiveraThorsen.tex
# Build document
cd ToSubmit
pdflatex RiveraThorsen.tex
pdflatex RiveraThorsen.tex
bibtex RiveraThorsen
pdflatex RiveraThorsen.tex
cd ..

