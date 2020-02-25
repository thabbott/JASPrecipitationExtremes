echo "Recalculating values in Table 1"
cd figures
echo "Channel, default"
matlab -nodesktop -nosplash -wait -r "TableA1"
cd ..
cd figures_rcemip
echo "Small-domain, default"
matlab -nodesktop -nosplash -wait -r "TableA1"
cd ..
cd figures_rcemip_m2005
echo "Small-domain, M2005"
matlab -nodesktop -nosplash -wait -r "TableA1"
cd ..
cd figures_rcemip_tke15
echo "Small-domain, TKE15"
matlab -nodesktop -nosplash -wait -r "TableA1"
cd ..
cd figures_rcemip_m2005_tke15
echo "Small-domain, M2005 + TKE15"
matlab -nodesktop -nosplash -wait -r "TableA1"
cd ..
echo "Finished"
