#!/bin/bash

# "Installs" or updates atk scripts from github
# Usage:
#	install_atk_scripts.sh  
#		Creates or updates (if it already exists) the collection of scripts
#
#	install_atk_scripts.sh reset   
#		Deletes the existing copy (if it exists) and reclones it
#
#	install_atk_scripts.sh uninstall
#		Uninstalls, that is deletes the directory and the respective entry in .bashrc
#
# The scripts creates (if it doesn't exist) the directory ~/progs
# clones https://github.com/merement/atk-scripts.git into there
# and adds $HOME/progs/atk-scripts to PATH, if it's not already there

comm_line="# Path to atk-scripts executables (do not modify these two lines)"
path_line="export PATH=\$PATH:\$HOME/progs/atk-scripts"

function create_clone() {
	git clone https://github.com/merement/atk-scripts.git
}

function update_clone() {
	cd atk-scripts
	git pull origin
}

if [ "$1" == "uninstall" ]; then 
	if [ -d "~/progs/atk-scripts/.git" ]; then
		rm -rf atk-scripts
	fi

	sed -i "/$comm_line/d" ~/.bashrc
	# for PATH we need to deal with slashes
	sed -i "/${path_line//\//\\/}/d" ~/.bashrc

	echo "atk-scripts were uninstalled and PATH was modified. This will be in effect in the next terminal session."
else
	cd ~
	mkdir -p progs

	cd progs
	if [ -d "atk-scripts/.git" ]; then
		# the project is there already
		if [ "$1" == "reset" ]; then 
			rm -rf atk-scripts && mkdir atk-scripts
			create_clone
			echo "Old installation is removed and scripts are reinstalled"
		else
			update_clone
			echo "Scripts are updated"
		fi
	else
		create_clone
		echo "Scripts are installed"
	fi

	# adds ~/progs/atk-scripts to PATH if it's not there
	# from http://stackoverflow.com/questions/1396066/detect-if-users-path-has-a-specific-directory-in-it
	if [[ ! ":$PATH:" == *":$HOME/progs/atk-scripts:"* ]]; then
		# check that it's not in .bashrc (the terminal was not restarted and .bashrc was possible mingled)
		if ! grep -Fxq "${comm_line}" ~/.bashrc; then
			# add an empty line to the end
			# http://unix.stackexchange.com/questions/31947/how-to-add-a-newline-to-the-end-of-a-file
			sed -i -e '$a\' ~/.bashrc
			printf "\n${comm_line}\n" >> ~/.bashrc
	#	fi
	#	if ! grep -Fxq "${path_line}" ~/.bashrc; then
			echo "${path_line}" >> ~/.bashrc
			echo "The script collection is added to PATH. The terminal must be restarted for changes to be in effect"
		fi
	fi
fi