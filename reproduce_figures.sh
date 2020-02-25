echo "Reproducing figures"
cd figures
matlab -nodesktop -nosplash -wait -r "Figure1; exit"
matlab -nodesktop -nosplash -wait -r "Figure2; exit"
matlab -nodesktop -nosplash -wait -r "Figure3; exit"
matlab -nodesktop -nosplash -wait -r "Figure4; exit"
matlab -nodesktop -nosplash -wait -r "Figure5; exit"
matlab -nodesktop -nosplash -wait -r "Figure6; exit"
matlab -nodesktop -nosplash -wait -r "Figure7; exit"
matlab -nodesktop -nosplash -wait -r "Figure8and9; exit"
matlab -nodesktop -nosplash -wait -r "Figure10; exit"
matlab -nodesktop -nosplash -wait -r "FigureA1; exit"
matlab -nodesktop -nosplash -wait -r "FigureA5; exit"
cd ..
cd figures_rcemip
matlab -nodesktop -nosplash -wait -r "Figure5; exit"
matlab -nodesktop -nosplash -wait -r "Figure8and9; exit"
cd ..
cd figures_rcemip_m2005
matlab -nodesktop -nosplash -wait -r "Figure5; exit"
matlab -nodesktop -nosplash -wait -r "Figure8and9; exit"
cd ..
cd figures_rcemip_tke15
matlab -nodesktop -nosplash -wait -r "Figure5; exit"
matlab -nodesktop -nosplash -wait -r "Figure8and9; exit"
cd ..
cd figures_rcemip_m2005_tke15
matlab -nodesktop -nosplash -wait -r "Figure5; exit"
matlab -nodesktop -nosplash -wait -r "Figure8and9; exit"
cd ..
echo "Finished"
