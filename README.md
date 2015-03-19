# atk-scripts
Various scripts to deal with Atomistic data

There's the collection of scripts I've used for basic analysis of the Bloch states. The current version of the collection is shown in file version in that folder. Whenever you need these scripts you can copy them to your working directory

$ cp /home/common/shared/misha/scripts_ksec/* /home/my/current/working/directory

Next follows somewhat longish text, its TL;DR version is in the WORKFLOW section below

stanalysis.py
        Usage:  atkpython stanalysis.py -c Name_of_File_with_DFT_results.nc -k [K | G | Kbar] --size_up Up_Num --size_down Down_Num

        It takes the results of the DFT calculations and for the chosen k-point (K, G or Kbar) goes through bands from Up_Num below the Fermi level to Down_Num above the Fermi level. For each band it calculates the center of mass of the Bloch state and average values of the momentum and angular momentum. More precisely, <r>, <\nabla>, <r \times \nabla>. It should be noted that G, K, Kbar actually mean vectors [0,0,0], [1/3, 1/3, 0], [-1/3, -1/3, 0]. So, they work for arbitrarily shaped supercells.

        The recommended way of using it is to create in the folder with the ATK data the file main.txt, which must contain a single line with the name of the file where the ATK wrote the results of the DFT calculations. This name is specified in VNL and if you don't remember it, you can always look it up in the .py-file produced by VNL in the line like nlsave("FileName.nc", bulk_configuration)

        When main.txt exists it's usually sufficient to provide the k-point only. In this case it will take 17 bands above and 17 bands below.

It may take a long time for stanalysis.py to finish, depending on the number of atoms in the elementary cell. It may easily be days if the number of atoms is more than couple hundreds.

Therefore three scripts are provided

subStan-G.sh
subStan-Kbar.sh
subStan-K.sh

for submission to Torque through qsub. For example, $ qsub subStan-G.sh


Important! These scripts assume that the name of the file with the results of DFT calculations is in the main.txt

When these scripts are done they will produce a bunch of .dat files with the results of processing, Bloch_Something.nc files with respective Bloch states,  and log files

outStG.log, outStK.log, outStKbar.log

These logs are worth looking into and (mostly for historical reasons, meaning at one point I've screwed up with dat files and had to deal with logs instead) they are of big importance.

These log files must be parsed into more computer friendly form and for this purpose there's the script

parselog.py

All three log files are processed by the shell script

./passlogs.sh

Invoking this shell script will produce (if no errors occurred, of course) three .json-files These files are (somewhat) human readable and are made for further automatic analysis.

At the moment such "further analysis" is simple plotting. For this there's the script

processjson.py

It's envoked by (the usual python is sufficient, no need for atkpython)

python processjson.py FileName.json

For example, python processjson.py outStG.dat.json

It will produce three .png figures showing <r>, <p>, <L> - <r> \times <p>

Finally, in order to see Bloch states one needs to use the script

visBloch-red.py

It requires the atkpython and the name of a file with the Bloch states, for instance

atkpython visBloch-red.py Bloch_something.nc

It will probably take about an hour to finish for a system of decent size. Unfortunately, at the moment it cannot be used with Torque.

When visBloch-red.py is done it will make the whole bunch of png files with graphical representations of Bloch states, which then should be sorted by hands. Files 'fig_0_' show the probability distribution in the plane of the sheet at z = 0.25 (fractional coordinate), files 'fig_1_' show the distribution integrated over z-axis.

WORKFLOW

After DFT calculations are finished a file File_with_DFT_results.nc is created.

1. We create the main.txt file

$ echo -n "File_with_DFT_results.nc" > main.txt

2. Submit to Torque

$ qsub subStan-G.sh
$ qsub subStan-K.sh
$ qsub subStan-Kbar.sh

3. Process log files

$ ./passlogs.sh

4. Plot results

$ python processjson.py outStK.log.json
$ python processjson.py outStKbar.log.json
$ python processjson.py outStKbar.log.json

5. Move the resultant png files into a special folder for easier access

$ mkdir analysis
$ mv *.png analysis

6. Plot Band states for the G point

$ atkpython visBloch-red.py Bloch_G_something.nc

7. Move figures to special folder

$ mkdir BSFigs_G
$ mv *.png BSFigs_G

8. Repeat 6 and 7 for K and Kbar

If there are any questions and especially if there are any errors during any of these steps, please, let me know
