echo "Reproducing figures"
cd figures
matlab -r Figure1
matlab -r Figure2
matlab -r Figure3
matlab -r Figure4
matlab -r Figure5
matlab -r Figure6
matlab -r Figure7
matlab -r Figure8and9
matlab -r Figure10
matlab -r FigureA1
matlab -r FigureA5
cd ..
cd figures/rcemip
matlab -r Figure5
matlab -r Figure8and9
cd ..
cd figures/rcemip_m2005
matlab -r Figure5
matlab -r Figure8and9
cd ..
cd figures/rcemip_tke15
matlab -r Figure5
matlab -r Figure8and9
cd ..
echo "Finished"
