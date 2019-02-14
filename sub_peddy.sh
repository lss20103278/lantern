#########################################################################
# File Name: sub_peddy.sh
# Author: sslyu
# mail: 1573077420@qq.com
# Created Time:Fri 12 Oct 2018 09:53:08 AM CST
#########################################################################
#!/bin/bash

cd peddy
sed -i "s/^sample.*$/sample=peddy/;s/^#SBATCH -J .*$/#SBATCH -J peddy/" trio_cmd.sh
sbatch trio_cmd.sh
cd ..
