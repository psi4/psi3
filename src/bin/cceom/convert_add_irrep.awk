BEGIN {
}
{

if (index($1,"buf4_init") || index($1,"file2_init")) {
$3 = "irrep, "$3
}
print $0

}

