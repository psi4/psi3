Summary: PSI Quantum Chemical Program Suite
Name: psi
Version: 3
Release: 2
Copyright: Virginia Tech., Georgia Tech., and others
Group: Applications/Scientific
Source: ftp://sirius.chem.vt.edu/psi3/psi-3.2.tgz
URL: http://www.psicode.org/
Vendor: Virginia Tech.
Packager: T. Daniel Crawford <crawdad@vt.edu>

%description: 
The PSI3 quantum chemistry package is designed to compute properties of
small molecules using high-level ab initio techniques.

%prep
%setup

%build
./configure
make

%install
make install

%files
/usr/local/bin/ccdensity
/usr/local/bin/ccenergy
/usr/local/bin/cceom
/usr/local/bin/cchbar
/usr/local/bin/cclambda
/usr/local/bin/ccresponse
/usr/local/bin/ccsort
/usr/local/bin/cctriples
/usr/local/bin/cints
/usr/local/bin/cis
/usr/local/bin/clag
/usr/local/bin/cphf
/usr/local/bin/cscf
/usr/local/bin/detcas
/usr/local/bin/detcasman
/usr/local/bin/detci
/usr/local/bin/extrema
/usr/local/bin/geom
/usr/local/bin/input
/usr/local/bin/localize
/usr/local/bin/mocube
/usr/local/bin/mp2
/usr/local/bin/mp2r12
/usr/local/bin/oeprop
/usr/local/bin/optking
/usr/local/bin/psi2molden
/usr/local/bin/psi3
/usr/local/bin/psiclean
/usr/local/bin/response
/usr/local/bin/stable
/usr/local/bin/tocprint
/usr/local/bin/transqt
/usr/local/share/pbasis.dat
/usr/local/share/psi.dat
/usr/local/man/man1/ccenergy.1
/usr/local/man/man1/cceom.1
/usr/local/man/man1/cints.1
/usr/local/man/man1/clag.1
/usr/local/man/man1/cscf.1
/usr/local/man/man1/cphf.1
/usr/local/man/man1/detcas.1
/usr/local/man/man1/detcasman.1
/usr/local/man/man1/detci.1
/usr/local/man/man1/geom.1
/usr/local/man/man1/input.1
/usr/local/man/man1/mocube.1
/usr/local/man/man1/mp2.1
/usr/local/man/man1/mp2r12.1
/usr/local/man/man1/oeprop.1
/usr/local/man/man1/optking.1
/usr/local/man/man1/psi3.1
/usr/local/man/man1/psiclean.1
/usr/local/man/man1/tocprint.1
/usr/local/man/man1/transqt.1
%doc /usr/local/doc/html/cints.html
%doc /usr/local/doc/html/CINTS/home.html
%doc /usr/local/doc/html/CINTS/commandline.html
%doc /usr/local/doc/html/CINTS/description.html
%doc /usr/local/doc/html/CINTS/references.html
%doc /usr/local/doc/html/CINTS/compile.html
%doc /usr/local/doc/html/CINTS/examples.html
%doc /usr/local/doc/html/CINTS/keywords.html
%doc /usr/local/doc/html/CINTS/toc.html
%doc /usr/local/doc/html/userman/contents.png
%doc /usr/local/doc/html/userman/img1.png
%doc /usr/local/doc/html/userman/img2.png
%doc /usr/local/doc/html/userman/img3.png
%doc /usr/local/doc/html/userman/img4.png
%doc /usr/local/doc/html/userman/img5.png
%doc /usr/local/doc/html/userman/img6.png
%doc /usr/local/doc/html/userman/img7.png
%doc /usr/local/doc/html/userman/img8.png
%doc /usr/local/doc/html/userman/next_g.png
%doc /usr/local/doc/html/userman/prev.png
%doc /usr/local/doc/html/userman/next.png
%doc /usr/local/doc/html/userman/up_g.png
%doc /usr/local/doc/html/userman/prev_g.png
%doc /usr/local/doc/html/userman/up.png
%doc /usr/local/doc/html/userman/index.html
%doc /usr/local/doc/html/userman/node1.html
%doc /usr/local/doc/html/userman/node2.html
%doc /usr/local/doc/html/userman/node3.html
%doc /usr/local/doc/html/userman/node4.html
%doc /usr/local/doc/html/userman/node5.html
%doc /usr/local/doc/html/userman/node6.html
%doc /usr/local/doc/html/userman/node7.html
%doc /usr/local/doc/html/userman/node8.html
%doc /usr/local/doc/html/userman/node9.html
%doc /usr/local/doc/html/userman/userman.html
%doc /usr/local/doc/html/userman/userman.css
%doc /usr/local/doc/ps/userman.ps
