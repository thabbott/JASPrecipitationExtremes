cd manuscript
echo "Compiling figures"
pdflatex Figure1.tex
pdflatex Figure2.tex
pdflatex Figure3.tex
pdflatex Figure4.tex
pdflatex Figure5.tex
pdflatex Figure6.tex
pdflatex Figure7.tex
pdflatex Figure8.tex
pdflatex Figure9.tex
pdflatex Figure10.tex
pdflatex FigureA1.tex
pdflatex FigureA2.tex
pdflatex FigureA3.tex
pdflatex FigureA4.tex
pdflatex FigureA5.tex
echo "Finished compiling figures"
echo "Compiling manuscript"
pdflatex JASPrecipitationExtremes.tex
pdflatex JASPrecipitationExtremes.tex
bibtex JASPrecipitationExtremes.tex
bibtex JASPrecipitationExtremes.tex
pdflatex JASPrecipitationExtremes.tex
pdflatex JASPrecipitationExtremes.tex
echo "Compiled manuscript at manuscript/JASPrecipitationExtremes.pdf"
cd ..
