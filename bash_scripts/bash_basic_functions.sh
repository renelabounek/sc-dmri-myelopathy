#!/bin/bash
# This script contains some bash functions widely used for analysis
# Jan Valosek, fMRI laboratory Olomouc, Czech Republic, 2018-2020
#
#########################################################################
# Function for checking existency of directories or files
# USAGE:
#   check_input.sh <argument> "DIR_NAME(s) or FILE_NAME(s)"
# ARGUMENTS:
#   d - check if directory exists
#   dc - check if directory exists, if not create it
#   f - check if file exists
#   b - check if binary exists (e.g. FSL, ANTS, SCT, Matlab,...)
#   e - check if script/function exist and is executable
#########################################################################
#
# Copyright 2016-2020 Rene Labounek (1,2,3,4), Jan Valosek (1,2) and Petr Hlustik (1,2)
#
# 1 - University Hospital Olomouc, Olomouc, Czech Republic
# 2 - Palacky University Olomouc, Olomouc, Czech Republic
# 3 - University Hospital Brno, Brno, Czech Republic 
# 4 - University of Minnesota, Minneapolis, US
#
# This file is part of sc-dmri-myelopathy available at: https://github.com/renelabounek/sc-dmri-myelopathy
#
# Please, cite sc-dmri-myelopathy as:
# Labounek R, Valosek J, Horak T, Svatkova A, Bednarik P, Vojtisek L, Horakova M, Nestrasil I,
# Lenglet C, Cohen-Adad J, Bednarik J and Hlustik P. HARDI-ZOOMit protocol improves specificity
# to microstructural changes in presymptomatic myelopathy. Scientific Reports [Revised; Under review]
#
# sc-dmri-myelopathy is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or any later version.
#
# sc-dmri-myelopathy is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with sc-dmri-myelopathy.  If not, see <https://www.gnu.org/licenses/>.
#
check_input()
{

  flag="$1"
  shift

  if [[ $flag =~ "h" ]] || [[ $flag == "" ]]; then
    echo -e "Help for function for checking existency of directories or files.\nUSAGE:\n\tcheck_input <argument> DIR_NAME(s) or FILE_NAME(s) or SCRIPT_NAME(s)"
    return
  fi

  while [[ "$1" != "" ]]; do
    case $flag in
      # Check if directory exists
      d)
          if [[ ! -d "$1" ]]; then
            echo "Directory $1 does not exist."
            trap_exit
          fi
          ;;
      # Check if directory exists, if not create it
      dc)
          if [[ ! -d "$1" ]]; then
            echo "Creating $1 directory."
            mkdir $1
            if [[ ! -d "$1" ]]; then
              echo "Directory $1 could not be created."
              trap_exit
            fi
          fi
          ;;
      # Check if file exists
      f)
          if [[ ! -f "$1" ]]; then
            echo "File $1 does not exist."
            trap_exit
          fi
          ;;
      # Check if binary exists
      b)
          if [[ "$(which $1)" == "" ]]; then
            echo "Path to $1 is wrong or $1 is not installed."
            trap_exit
          fi
          ;;
      # Check if script/function exists and is executable, if not, set permissions (to 770)
      e)
          if [[ ! -x "$1" ]]; then
              if [[ ! -f "$1" ]]; then
                  echo "Function/Script $1 does not exist."
                  trap_exit
              else
                  echo "Setting permissions for $1 script/function."
                  chmod 770 $1
                  if [[ ! -x "$1" ]]; then
                      echo "Fuction/Script $1 is not executable and permissions can not be set."
                      trap_exit
                  fi
              fi
          fi
          ;;
    esac
    shift
  done
}

#########################################################################
# Function for execution of certain command and prints its output to log
# USAGE:
#   exe "command" [option]
# OPTIONS:
#   [no option] - verbose mode ON (print command' output to both terminal and to log)
#   v - turn OFF verbose (only see which command is executed but all its output is redirected to log)
#   t - turn OFF time mark

# IMPORTANT: This funtion collaborates and deppends on show function

# Jan Valosek, fMRI laboratory, Olomouc, 2019-2020. Inspired by Pavel Hok
# VER=26-02-2020
#########################################################################
exe()
{
  # Check if dependecy on show function is fullfiled
  if [[ $(which show) == " " ]]; then
      echo -e "Function show was not found. Add path to it to \$PATH variable.\nExiting..."
      trap_exit
  fi

  if [[ $1 == "h" ]] || [[ $1 == "-h" ]] || [[ $1 == "--help" ]] || [[ $1 == "" ]]; then
    echo -e "Help for function for execution of certain command and prints its output to log.\nUSAGE:\n\texe \"command_1\" [option]"
    return
  fi

  if [[ $2 =~ t ]]; then
    show "Executing: $1" $2       # Print informative message which command is executed to terminal without time mark
  else
    show "Executing: $1"          # Print informative message which command is executed to terminal with time mark
  fi
  $1 | show --stdin $2            # Run command and Redirect output of the command to log (v option) or to terminal (without v option)

  # Check if command finished with error, if so, call function show with error argument (see help of function show)
  lasterr=${PIPESTATUS[0]}                # get last error
  if [[ $lasterr != 0 ]]; then
      show "Sub-process returned error value $lasterr" e
  fi
}

#########################################################################
# Function for printing message or command's output in various ways
# USAGE:
#    show "message" [option]
# OPTIONS:
#       [no option] - verbose mode on (print output to both terminal and to log with time stamp)
#       v - turn OFF verbose (print message only to log with time stamp)
#       e - error (print message in red color to both terminal and log and exit)
#       w - warning (print message in blue color to both terminal and log but stil continue)
#       y - print standard message in yellow color
#       g - print standard message in green color
#       t - turn OFF time mark
# EXAMPLE USAGE:
#       show "Hello world"      - print "Hello world" to both terminal and log with time stamp
#       show "Hello world" t    - print "Hello world" to both terminal and log without time stamp
#       show "Hello world" g    - print "Hello world" in green color to both terminal and log with time stamp
#       show "Hello world" gt   - print "Hello world" in green color to both terminal and log without time stamp
#       show "Hello world" v    - print "Hello world" only to log (verbose off) with time stamp
#       show "Hello world" vt   - print "Hello world" only to log (verbose off) without time stamp
#       show "Hello world" e    - print "Hello world" as error in red color to both terminal and log and exit main script

# IMPORTANT: This funtion collaborates with exe function

# Jan Valosek, fMRI laboratory, Olomouc, 2019-2020. Inspired by Pavel Hok
# VER=26-02-2020
#########################################################################
# Colors definition
red=$(tput setaf 1)
yellow=$(tput setaf 3)
blue=$(tput setaf 4)
green=$(tput setaf 6)
normal=$(tput sgr0)

show()
{

  if [[ $1 == "h" ]] || [[ $1 == "-h" ]] || [[ $1 == "--help" ]] || [[ $1 == "" ]]; then
    echo -e "Help for printing message or command's output in various ways.\nUSAGE:\n\tshow \"message\" [option]"
    return
  fi

  # Switch off time mark
  if [[ $2 =~ "t" ]]; then
      print_date=""
  else
      print_date="$(date "+%T_%Y-%m-%d"): "
  fi


  # ERROR - print message in red color to both terminal and log with time stamp and exit
  if [[ $2 =~ "e" ]]; then

          echo -e "${red}${print_date}ERROR: $1\nExiting...${normal}" 1>&2
          trap_exit

  # WARNING - print message in red color to both terminal and log  with time stamp but stil continue
  elif [[ $2 =~ "w" ]]; then

          #echo -e "${blue}${print_date}WARNING: $1\nContinuing...${normal}" 1>&2
          echo -e "${blue}${print_date}WARNING: $1${normal}" 1>&2

  # Standard message in yellow color to both terminal and log with time stamp
  elif [[ $2 =~ "y" ]]; then

          echo -e "${yellow}${print_date}${1}${normal}" 1>&2

  # Standard message in green color to both terminal and log with time stamp
  elif [[ $2 =~ "g" ]]; then

          echo -e "${green}${print_date}${1}${normal}" 1>&2

  # NO OPTION - output from command or message to screen and to log (= verbose mode)
  elif [[ ! ( $2 =~ "v" ) ]]; then

          # message
          if [[ "$1" != "--stdin" ]]; then
                  echo -e "${print_date}${1}"
          # Command
          else
              while read line; do
                      echo -e "$line"
              done
          fi

  # TURN OFF verbose - output from command or message only to log
  elif [[ $2 =~ "v" ]]; then

      if [[ $LOGPATH == "" ]]; then
          echo -e "Path to log is empty. Fix it (e.g. by setting LOGPATH variable as a global in your script (export \$LOGPATH=\"log.txt\")).\nExiting..."
          trap_exit
      fi

          # message
          if [[ "$1" != "--stdin" ]]; then
                  echo -e "${print_date}${1}" >> $LOGPATH
          # Command
          else
                  while read line; do
                          echo -e "$line" >> $LOGPATH
                  done
          fi

  fi
}

#########################################################################
# Function for running matlab from command line without opening MATLAB desktop
# Explanation:
#   https://www.mathworks.com/help/matlab/ref/matlablinux.html
# USAGE:
#   run_matlab "command_1,...,command_X,exit"
# e.g.
#   run_matlab "addpath('$MATLAB_SCRIPT'),NODDI_analysis('$subdir','$SUB'),exit"

# Jan Valosek, Rene Labounek, fMRI laboratory, Olomouc, 2018-2020.
# VER=23-12-2019
#########################################################################
run_matlab()
{
  if [[ $(which matlab) == "" ]]; then
      echo "Matlab was not found."
      trap_exit
  fi

  if [[ $1 == "h" ]] || [[ $1 == "-h" ]] || [[ $1 == "--help" ]] || [[ $1 == "" ]]; then
    echo -e "Help for function for running matlab from command line without opening MATLAB desktop.\nUSAGE:\n\trun_matlab \"command_1,...,command_X,exit\""
    echo -e "EXAMPLE:\n\trun_matlab \"addpath('\$MATLAB_SCRIPT'),NODDI_analysis('\$subdir','\$SUB'),exit\""
    return
  fi

  COMMAND=$1

  if [[ $HOSTNAME == emperor ]];then
     /usr/local/MATLAB/R2014b/bin/matlab -nosplash -nodisplay -nodesktop -r "$COMMAND"           #JV 11.5.18 (R2016B does not contain Image Processing Toolbox)
  elif [[ $HOSTNAME == king ]];then
    /usr/local/MATLAB/R2016b/bin/matlab -nosplash -nodisplay -nodesktop -r "$COMMAND"           #JV 27.5.2019 (R2018B does not contain Image Processing Toolbox)
  else
    matlab -nosplash -nodisplay -nodesktop -r  "$COMMAND"
  fi
}

#########################################################################
# Function telling to the parent process that something happened (e.g. file/binary/script does not exist)
# Function finds parent process ID (ppid) and set user defined signal (SIGUSR1)
# In the parent script has to be following line for catching SIGUSR1:
#     trap "echo Exitting>&2;exit" SIGUSR1
# Inspiration:
#     https://stackoverflow.com/questions/16215196/bash-exit-parent-script-from-child-script
#########################################################################
trap_exit()
{
    OS_name=$(uname -a | awk '{print $1}')  # get OS name (Darwin or Linux)

    if [[ $OS_name == "Linux" ]]; then
        sleep 1;kill -SIGUSR1 `ps --pid $$ -o ppid=`;exit   # Find parent process ID (ppid) and set user defined signal (SIGUSR1)
    elif [[ $OS_name == "Darwin" ]]; then
        sleep 1;kill -SIGUSR1 `ps $$ -o ppid=`;exit         # Find parent process ID (ppid) and set user defined signal (SIGUSR1)
    fi
}
