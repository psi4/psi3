Summary: PSI Quantum Chemical Program Suite
Name: psi
Version: 3.0
Release: 1
Copyright: PSITech, Inc.
Group: Applications/Scientific
Source: ftp://sirius.chem.vt.edu/psi3/psi3.tgz
URL: http://vergil.chemistry.gatech.edu/psi/
Vendor: PSITech, Inc.
Packager: T. Daniel Crawford <crawdad@vt.edu>

%description: 
The PSI3 quantum chemistry package is designed to compute properties of
small molecules using high-level ab initio techniques.

%prep
%setup

%build
./configure --exec_prefix=/usr/local/psi
make

%install
make install

%files
/usr/local/psi/bin/ccdensity
/usr/local/psi/bin/ccenergy
/usr/local/psi/bin/ccsort
/usr/local/psi/bin/cphf
/usr/local/psi/bin/extrema
/usr/local/psi/bin/intdif
/usr/local/psi/bin/normco
/usr/local/psi/bin/rgeom
/usr/local/psi/bin/cctriples
/usr/local/psi/bin/cscf
/usr/local/psi/bin/geom
/usr/local/psi/bin/makepk
/usr/local/psi/bin/oeprop
/usr/local/psi/bin/cceom
/usr/local/psi/bin/cints
/usr/local/psi/bin/detcas
/usr/local/psi/bin/input
/usr/local/psi/bin/mp2
/usr/local/psi/bin/optking
/usr/local/psi/bin/transqt
/usr/local/psi/bin/cchbar
/usr/local/psi/bin/contour
/usr/local/psi/bin/detcasman
/usr/local/psi/bin/mp2r12
/usr/local/psi/bin/psi
/usr/local/psi/bin/ugeom
/usr/local/psi/bin/cclambda
/usr/local/psi/bin/contour3d
/usr/local/psi/bin/detci
/usr/local/psi/bin/muder
/usr/local/psi/bin/psiclean
/usr/local/psi/share/pbasis.dat
/usr/local/psi/share/psi.dat
/usr/local/psi/share/tmpdisks.dat
/usr/local/psi/man/man1/cints.1
/usr/local/psi/man/man1/detcas.1
/usr/local/psi/man/man1/mp2r12.1
/usr/local/psi/man/man1/psi.1
/usr/local/psi/man/man1/ugeom.1
/usr/local/psi/man/man1/contour.1
/usr/local/psi/man/man1/detcasman.1
/usr/local/psi/man/man1/intdif.1
/usr/local/psi/man/man1/muder.1
/usr/local/psi/man/man1/psiclean.1
/usr/local/psi/man/man1/contour3d.1
/usr/local/psi/man/man1/detci.1
/usr/local/psi/man/man1/local.1
/usr/local/psi/man/man1/normco.1
/usr/local/psi/man/man1/rgeom.1
/usr/local/psi/man/man1/cphf.1
/usr/local/psi/man/man1/geom.1
/usr/local/psi/man/man1/makepk.1
/usr/local/psi/man/man1/oeprop.1
/usr/local/psi/man/man1/cscf.1
/usr/local/psi/man/man1/input.1
/usr/local/psi/man/man1/mp2.1
/usr/local/psi/man/man1/optking.1
/usr/local/psi/man/man1/transqt.1
%doc /usr/local/psi/doc/html/cints.html
%doc /usr/local/psi/doc/html/CINTS/home.html
%doc /usr/local/psi/doc/html/CINTS/commandline.html
%doc /usr/local/psi/doc/html/CINTS/description.html
%doc /usr/local/psi/doc/html/CINTS/references.html
%doc /usr/local/psi/doc/html/CINTS/compile.html
%doc /usr/local/psi/doc/html/CINTS/examples.html
%doc /usr/local/psi/doc/html/CINTS/keywords.html
%doc /usr/local/psi/doc/html/CINTS/toc.html
%doc /usr/local/psi/doc/html/userman/contents.png
%doc /usr/local/psi/doc/html/userman/img1.png
%doc /usr/local/psi/doc/html/userman/img2.png
%doc /usr/local/psi/doc/html/userman/img3.png
%doc /usr/local/psi/doc/html/userman/img4.png
%doc /usr/local/psi/doc/html/userman/img5.png
%doc /usr/local/psi/doc/html/userman/img6.png
%doc /usr/local/psi/doc/html/userman/img7.png
%doc /usr/local/psi/doc/html/userman/img8.png
%doc /usr/local/psi/doc/html/userman/next_g.png
%doc /usr/local/psi/doc/html/userman/prev.png
%doc /usr/local/psi/doc/html/userman/next.png
%doc /usr/local/psi/doc/html/userman/up_g.png
%doc /usr/local/psi/doc/html/userman/prev_g.png
%doc /usr/local/psi/doc/html/userman/up.png
%doc /usr/local/psi/doc/html/userman/index.html
%doc /usr/local/psi/doc/html/userman/node1.html
%doc /usr/local/psi/doc/html/userman/node2.html
%doc /usr/local/psi/doc/html/userman/node3.html
%doc /usr/local/psi/doc/html/userman/node4.html
%doc /usr/local/psi/doc/html/userman/node5.html
%doc /usr/local/psi/doc/html/userman/node6.html
%doc /usr/local/psi/doc/html/userman/node7.html
%doc /usr/local/psi/doc/html/userman/node8.html
%doc /usr/local/psi/doc/html/userman/node9.html
%doc /usr/local/psi/doc/html/userman/userman.html
%doc /usr/local/psi/doc/html/userman/userman.css
%doc /usr/local/psi/doc/ps/userman.ps
