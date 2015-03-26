# atk-scripts
Various scripts to deal with Atomistic data

There's the collection of scripts I've used for basic analysis of the Bloch states. The current version of the collection is shown in file version in that folder. Whenever you need these scripts you can copy them to your working directory.

    $ cp /home/common/shared/misha/scripts_ksec/* /home/my/current/working/directory

Next follows somewhat longish text, its TL;DR version is in the WORKFLOW section below

stanalysis.py

Usage:  stanalysis.py [-c Name_of_File_with_DFT_results.nc] [-k [K | G | Kbar]] [--size_up Up_Num] [--size_down Down_Num]

It takes the results of the DFT calculations and for the chosen k-point (K, G or Kbar) goes through bands from Up_Num below the Fermi level to Down_Num above the Fermi level. For each band it calculates the center of mass of the Bloch state and average values of the momentum and angular momentum. More precisely, <r>, <\nabla>, <r \times \nabla>. It should be noted that G, K, Kbar actually mean vectors [0,0,0], [1/3, 1/3, 0], [-1/3, -1/3, 0]. So, they work for arbitrarily shaped supercells.

The recommended way of using it is to create in the folder with the ATK data the file main.txt, which must contain a single line with the name of the file where the ATK wrote the results of the DFT calculations. This name is specified in VNL and if you don't remember it, you can always look it up in the .py-file produced by VNL in the line like nlsave("FileName.nc", bulk_configuration)

When main.txt exists it's usually sufficient to provide the k-point only. In this case it will take 17 bands above and 17 bands below. Thus, if you simply run

    $ stanalysis.py -k K

it will try to take the name of the file with the results of the DFT calculations from the file 'main.txt' and will perform calculations for the K point for 17 bands above and 17 bands below the Fermi level.

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

It will produce three .png figures showing < r >, < p >, < L > - < r > \times < p >

Finally, in order to see Bloch states one needs to use the script

visBloch-red.py

For instance

    visBloch-red.py Bloch_something.nc

A bit more complete information is obtained by expected

    $ visBloch-red.py -- help

    usage: visBloch-red.py [-h] [-s S] [-c C] [--section SECTION] [-t]
                           blochFileName

    positional arguments:
      blochFileName      The name of the .nc-file with Bloch functions

    optional arguments:
      -h, --help         show this help message and exit
      -s S               Projection of the spin
      -c C               The name of the .nc-file with the configuration
      --section SECTION  The list of sections to plot. Expected format is
                         --section [s1, s2, ...]. Parameters 0 < s_k < 1 specify
                         the relative z-coordinate of the (x, y)-plane. If s_k is
                         outside of this range the distribution integrated over z
                         is plotted.
      -t                 Provide to plot the distribution in the z-direction.

It will probably take about an hour to finish for a system of decent size and it's the good idea to submit it to Torque (see subBlochview-?.sh)

When visBloch-red.py is done it will make the whole bunch of png files with graphical representations of Bloch states, which then should be sorted by hands. Files 'fig_N_', where N is a number, show the probability distribution in the plane of the sheet at the respective z-coordinate, determined by the position (started from zero) of N in the list of sections specified by the option --section. If -t is provided, then the z-dependence will be in files with the largest N.

# WORKFLOW

After DFT calculations are finished a file File_with_DFT_results.nc is created.

1. We create the main.txt file

        $ echo -n "File_with_DFT_results.nc" > main.txt

2. If only Bloch functions and their pictures are need then do the following, otherwise go to 6
        $ cp $HOME/progs/atk-scripts/subBlochview-G.sh subBlochview-G.sh

3. Correct kPoint and bands in subBlochview-G.sh or if all you need is G point, 30 bands below and above the Fermi level, and their basic plots, leave it as it is

4. Submit to Torque

        $ qsub subBlochview-G.sh

5. Go to 12
        
6. If you need to perform the basic analysis of the Bloch states. Then (don't miss the point at the end)

        $ cp $HOME/progs/atk-scripts/subStan-*.sh .

7. Submit to Torque

        $ qsub subStan-G.sh
        $ qsub subStan-K.sh

8. Process log files

        $ ./passlogs.sh

9. Plot results

        $ processjson.py outStG.log.json
        $ processjson.py outStK.log.json

10. Move the resultant png files into a special folder for easier access

        $ mkdir analysis
        $ mv *.png analysis

11. Plot Band states (for example for the G point)

        $ visBloch-red.py Bloch_G_something.nc

12. Move figures to special folder

        $ mkdir BSFigs_G
        $ mv *.png BSFigs_G

This documentation is with necessity incomplete. For better understanding of functionality of the scripts it's worth looking into them.

If there are any questions and especially if there are any errors during any of these steps, please, let me know
