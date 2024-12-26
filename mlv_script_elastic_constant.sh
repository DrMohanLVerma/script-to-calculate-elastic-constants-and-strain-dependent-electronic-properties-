#!/bin/bash
#======================================================#
# Script: For Elastic Constant Calculation             #
# Author: Mohan L Verma, Computational Nanomaterial    #  
# Research lab, Department of Applied Physics,         #
# Shri Shanakaracharya Technical Campus-Junwani        # 
# Bhilai(Chhattisgarh)  INDIA                          #
# Feb 24   ver: 2.0    year: 2023                      #
# it is assumed that siesta.exe binary file is linked  #
# bin directory after parallel compilation of siesta   #
#------------------------------------------------------#
# run this script using cammand                        #
#  sh mlv_script_elastic.sh                            #
# give feedback in                                     #
#                  www.drmlv.in/siesta  or             # 
#                  drmohanlv@gmail.com                 #
#===================================================== #
 
 
mkdir elastic
cd    elastic


mkdir Bands
mkdir DOS
 

 
 

mkdir cont   # read the comment at the end of this script.

for i in `seq -w 10.344557142 0.057469762 12.643347618`  
do


cp -r cont $i
cd $i
cp ../../*.psf .

 


cat > Si331.fdf <<EOF

SystemName  Si331
SystemLabel      Si331

NumberOfAtoms    18

NumberOfSpecies  1
%block ChemicalSpeciesLabel
    1   14  Si
%endblock ChemicalSpeciesLabel

LatticeConstant 1.000 Ang
 
%block LatticeVectors
        $i         0.0000000000    0.0000000000
  -5.7469761900    9.9540547510    0.0000000000
   0.0000000000    0.0000000000    6.3326220000
%endblock LatticeVectors

AtomicCoordinatesFormat Fractional
%block AtomicCoordinatesAndAtomicSpecies
     0.222222223     0.111111110     0.563009500    1
     0.111111110     0.222222223     0.436990500    1
     0.888888890     0.777777777     0.563009500    1
     0.777777777     0.888888890     0.436990500    1
     0.555555557     0.777777777     0.563009500    1
     0.444444443     0.888888890     0.436990500    1
     0.222222223     0.777777777     0.563009500    1
     0.111111110     0.888888890     0.436990500    1
     0.888888890     0.444444443     0.563009500    1
     0.777777777     0.555555557     0.436990500    1
     0.555555557     0.444444443     0.563009500    1
     0.444444443     0.555555557     0.436990500    1
     0.222222223     0.444444443     0.563009500    1
     0.111111110     0.555555557     0.436990500    1
     0.888888890     0.111111110     0.563009500    1
     0.777777777     0.222222223     0.436990500    1
     0.555555557     0.111111110     0.563009500    1
     0.444444443     0.222222223     0.436990500    1
%endblock AtomicCoordinatesAndAtomicSpecies


==================================================
==================================================
# K-points

%block kgrid_Monkhorst_Pack
5   0   0   0.0
0   5   0   0.0
0   0   1   0.0
%endblock kgrid_Monkhorst_Pack

 
 
#%block GeometryConstraints
#position from  1 to  180
#%endblock GeometryConstraints

PAO.BasisSize     DZP
PAO.EnergyShift   0.03 eV
MD.TypeOfRun      CG
MaxSCFIterations  300
SCF.MustConverge   false
MD.NumCGsteps     00
MD.MaxForceTol    0.005  eV/Ang
MeshCutoff        250 Ry
DM.MixingWeight   0.02
DM.NumberPulay   3
WriteCoorXmol   .true.
WriteMullikenPop    1
XC.functional       GGA
XC.authors          PBE
SolutionMethod  diagon
ElectronicTemperature  50 meV
SaveRho        .true.



WriteEigenvalues                yes
 
%block ProjectedDensityOfStates
-30.00 20.00 0.200 1000 eV
%endblock ProjectedDensityOfStates


WriteBands    true 
BandLinesScale ReciprocalLatticeVectors

%block BandLines
1   0.5000000000     0.5000000000     0.0000000000      K
40  0.0000000000     0.0000000000     0.0000000000     GM
50  0.5000000000     0.5000000000    -0.5000000000     M
40  0.5000000000     0.5000000000     0.0000000000     K
%endblock Bandlines



UseSaveData     true
DM.UseSaveDM    true
MD.UseSaveXV    true
MD.UseSaveCG    true
EOF
 

mpirun -np 14 siesta.exe *.fdf | tee  result.out 
Eig2DOS -f -s 0.2 -n 2000 -m -30.0 -M 20.0 *.EIG>dos_$i.dat 

cp dos_$i.dat  ./../DOS

gnubands -F *.bands>bands_$i.dat
cp bands_$i.dat ./../Bands

etot=`grep 'Total =' result.out | cut -d = -f 2`
echo $i '      ' $etot >> ../EvsA.dat 

cp EvsA.dat ../

etot=`grep 'Fermi =' result.out | cut -d = -f 2`
echo $i '      ' $etot >> ../XvsEf.dat 

cp XvsEf.dat ../

stress=`awk 'c{c--;if(!c) print $3}/siesta: Stress tensor \(static\)/{c=2}' ./*.out`

echo $i '   '$stress >> ../strain_stress.dat

cd ..
rm -rf cont 
mkdir cont

cp   ./$i/*.DM  cont  # copy these files for continuation of the next step.

 

done
cp strain_stress.dat ../

cp EvsA.dat ../

cd ..
xmgrace EvsA.dat &

 
cd..


 

