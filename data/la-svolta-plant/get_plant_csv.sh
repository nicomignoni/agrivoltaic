cat <(echo 'pos_east,pos_north,width,depth') <(inkscape --query-all plant_layout_only.csv | tail -n +4 | cut -d, -f 2-) > plant_2.csv
