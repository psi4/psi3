#!/bin/sh

# This script guess's the architecture of the host.
unames=`uname -s`
if [ "$unames" = AIX ]; then
  test -z "$vendor" && vendor=ibm
  test -z "$os" && os=aix`uname -v`.`uname -r`
  machine_type=`uname -m | sed "s/........\(..\)../\1/"`
# 590's
  if [ $machine_type = 70 ]; then
    test -z "$arch" && arch=power2
# 397's
  elif [ $machine_type = 94 ]; then
    test -z "$arch" && arch=power2
# 595's
  elif [ $machine_type = 89 ]; then
    test -z "$arch" && arch=power2
# 3CT's
  elif [ $machine_type = 59 ]; then
    test -z "$arch" && arch=power2
# 530's
  elif [ $machine_type = 10 ]; then
    test -z "$arch" && arch=power
# 540's
  elif [ $machine_type = 14 ]; then
    test -z "$arch" && arch=power
# 580's
  elif [ $machine_type = 66 ]; then
    test -z "$arch" && arch=power
# 550's
  elif [ $machine_type = 1C ]; then
    test -z "$arch" && arch=power
# default
  else
    test -z "$arch" && arch=power
  fi
elif [ "$unames" = "ULTRIX-32" ]; then
  test -z "$vendor" && vendor=dec
  test -z "$os" && os=ultrix`uname -r`
  test -z "$arch" && arch=mips
elif [ "$unames" = "HP-UX" ]; then
  test -z "$vendor" && vendor=hp
  test -z "$os" && os=hpux`uname -r`
elif [ "$unames" = "IRIX" ]; then
  unamem=`uname -m`
  test -z "$vendor" && vendor=sgi
  test -z "$os" && os=irix`uname -r`
  if [ "$unamem" = IP19 ]; then
    test -z "$arch" && arch=mips2
  elif [ "$unamem" = IP20 ]; then
    test -z "$arch" && arch=mips2
  elif [ "$unamem" = IP22 ]; then
    test -z "$arch" && arch=mips2
  elif [ "$unamem" = IP17 ]; then
    test -z "$arch" && arch=mips2
  else
    test -z "$arch" && arch=mips
  fi
elif [ "$unames" = "IRIX64" ]; then
  unamem=`uname -m`
  test -z "$vendor" && vendor=sgi
  test -z "$os" && os=irix`uname -r`
  if [ "$unamem" = IP21 ]; then
    test -z "$arch" && arch=mips4
  fi
elif [ "$unames" = "Linux" ]; then
  test -z "$vendor" && vendor=generic
  test -z "$os" && os=linux`uname -r`
  test -z "$arch" && arch=`uname -m`
  # 
  # linux patchlevels change faster than my underwear
  #
  ost=`echo $os | tr . -`
  major=`echo $ost | sed 's/\([^-]*\)-\([^-]*\)-\([^-]*\)-*.*/\1/'`
  minor=`echo $ost | sed 's/\([^-]*\)-\([^-]*\)-\([^-]*\)-*.*/\2/'`
  patch=`echo $ost | sed 's/\([^-]*\)-\([^-]*\)-\([^-]*\)-*.*/\3/'`
  os=$major.$minor
else
  test -z "$os" && os=$unames
fi

if [ -z "$arch" ]; then arch=unknown; fi
if [ -z "$vendor" ]; then vendor=unknown; fi
if [ -z "$os" ]; then os=unknown; fi

echo $arch-$vendor-$os

