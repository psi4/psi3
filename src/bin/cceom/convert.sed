{
s/oe_file/file2/
s/oe_copy/file2_copy/
s/buf_/buf4_/
s/, 0, outfile//
s/struct oe_dpdfile/dpdfile2/g
s/struct dpdbuf/dpdbuf4/g
s/struct dpdfile/dpdbuf4/g
s/swap12/buf4_sort/
s/swap13/buf4_sort/
s/swap14/buf4_sort/
s/swap23/buf4_sort/
s/swap24/buf4_sort/
s/swap34/buf4_sort/
s/swapbk/buf4_sort/
s/dpd_axpy/dpd_buf4_axpy/
s/dpd_copy/dpd_buf4_copy/
s/dpd_scm/dpd_buf4_scm/
s/dpd_dot/dpd_buf4_dot/
s/contract222/contract444/
s/contract221/contract424/
s/contract212/contract244/
s/contract122/contract442/
s/contract121/contract422/
s/contract111/contract222/
s/, irrep//
}

