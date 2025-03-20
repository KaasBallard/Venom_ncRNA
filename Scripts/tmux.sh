#!/bin/bash

# This script is used to start or attach to a tmux session for the Venom_ncRNA project

# Project directory
project_dir="$HOME/Dropbox/CastoeLabFolder/projects/Venom_Gene_Regulation/Venom_ncRNA"

# Set the session name
session="Venom_ncRNA"

# Check if the session exists, and attach if it does
tmux has-session -t "$session" 2>/dev/null

# $? returns 0 if the session exists, so create it if it doesn't
if [ $? != 0 ]; then
	# Change the directory to the current project directory
	cd "$project_dir" || exit

	# Create a new session
	tmux new-session -d -s "$session" -n main  # The dash -d flag is used to prevent the tmux session from immediately attaching to the terminal

	# Create a new window named 'htop'
	tmux new-window -d -n htop

	# Create a new window named 'git' for running git commands
	tmux new-window -d -n git
fi

# Attach to the session
tmux attach-session -t "$session":main