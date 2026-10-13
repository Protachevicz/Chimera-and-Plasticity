#!/bin/bash
for i in {0..3}; do
  icx -fast -o z$i.z plastic_potentiation_matrix.c -lm
  echo
done

for i in {0..3}; do
   echo $i | ./z$i.z &
   echo
   sleep 1
   rm z$i.z
done
