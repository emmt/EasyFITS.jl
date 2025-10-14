module CFITSIO

using CEnum: CEnum, @cenum

using CFITSIO_jll: libcfitsio

struct Status
    code::Cint
end

fits_open_file(A, B, C, D) = ffopentest(CFITSIO_SONAME, A, B, C, D)

const off_t = Clong

function fits_parse_input_url(url, urltype, infile, outfile, extspec, rowfilter, binspec, colspec, status)
    @ccall libcfitsio.ffiurl(url::Cstring, urltype::Cstring, infile::Cstring, outfile::Cstring, extspec::Cstring, rowfilter::Cstring, binspec::Cstring, colspec::Cstring, status::Ptr{Status})::Status
end

function fits_parse_input_filename(url, urltype, infile, outfile, extspec, rowfilter, binspec, colspec, pixfilter, status)
    @ccall libcfitsio.ffifile(url::Cstring, urltype::Cstring, infile::Cstring, outfile::Cstring, extspec::Cstring, rowfilter::Cstring, binspec::Cstring, colspec::Cstring, pixfilter::Cstring, status::Ptr{Status})::Status
end

function fits_parse_rootname(url, rootname, status)
    @ccall libcfitsio.ffrtnm(url::Cstring, rootname::Cstring, status::Ptr{Status})::Status
end

function fits_file_exists(infile, exists, status)
    @ccall libcfitsio.ffexist(infile::Cstring, exists::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_parse_extspec(extspec, extnum, extname, extvers, hdutype, colname, rowexpress, status)
    @ccall libcfitsio.ffexts(extspec::Cstring, extnum::Ptr{Cint}, extname::Cstring, extvers::Ptr{Cint}, hdutype::Ptr{Cint}, colname::Cstring, rowexpress::Cstring, status::Ptr{Status})::Status
end

function fits_parse_extnum(url, extension_num, status)
    @ccall libcfitsio.ffextn(url::Cstring, extension_num::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_parse_binspec(binspec, imagetype, haxis, colname, minin, maxin, binsizein, minname, maxname, binname, weight, wtname, recip, status)
    @ccall libcfitsio.ffbins(binspec::Cstring, imagetype::Ptr{Cint}, haxis::Ptr{Cint}, colname::Ptr{NTuple{71, Cchar}}, minin::Ptr{Cdouble}, maxin::Ptr{Cdouble}, binsizein::Ptr{Cdouble}, minname::Ptr{NTuple{71, Cchar}}, maxname::Ptr{NTuple{71, Cchar}}, binname::Ptr{NTuple{71, Cchar}}, weight::Ptr{Cdouble}, wtname::Cstring, recip::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_parse_binrange(binspec, colname, minin, maxin, binsizein, minname, maxname, binname, status)
    @ccall libcfitsio.ffbinr(binspec::Ptr{Cstring}, colname::Cstring, minin::Ptr{Cdouble}, maxin::Ptr{Cdouble}, binsizein::Ptr{Cdouble}, minname::Cstring, maxname::Cstring, binname::Cstring, status::Ptr{Status})::Status
end

function fits_parse_range(rowlist, maxrows, maxranges, numranges, minrow, maxrow, status)
    @ccall libcfitsio.ffrwrg(rowlist::Cstring, maxrows::Clonglong, maxranges::Cint, numranges::Ptr{Cint}, minrow::Ptr{Clong}, maxrow::Ptr{Clong}, status::Ptr{Status})::Status
end

function fits_parse_rangell(rowlist, maxrows, maxranges, numranges, minrow, maxrow, status)
    @ccall libcfitsio.ffrwrgll(rowlist::Cstring, maxrows::Clonglong, maxranges::Cint, numranges::Ptr{Cint}, minrow::Ptr{Clonglong}, maxrow::Ptr{Clonglong}, status::Ptr{Status})::Status
end

struct tcolumn
    ttype::NTuple{70, Cchar}
    tbcol::Clonglong
    tdatatype::Cint
    trepeat::Clonglong
    tscale::Cdouble
    tzero::Cdouble
    tnull::Clonglong
    strnull::NTuple{20, Cchar}
    tform::NTuple{10, Cchar}
    twidth::Clong
end

struct FITSfile
    filehandle::Cint
    driver::Cint
    open_count::Cint
    filename::Cstring
    validcode::Cint
    only_one::Cint
    noextsyntax::Cint
    filesize::Clonglong
    logfilesize::Clonglong
    lasthdu::Cint
    bytepos::Clonglong
    io_pos::Clonglong
    curbuf::Cint
    curhdu::Cint
    hdutype::Cint
    writemode::Cint
    maxhdu::Cint
    MAXHDU::Cint
    headstart::Ptr{Clonglong}
    headend::Clonglong
    ENDpos::Clonglong
    nextkey::Clonglong
    datastart::Clonglong
    imgdim::Cint
    imgnaxis::NTuple{99, Clonglong}
    tfield::Cint
    startcol::Cint
    origrows::Clonglong
    numrows::Clonglong
    rowlength::Clonglong
    tableptr::Ptr{tcolumn}
    heapstart::Clonglong
    heapsize::Clonglong
    request_compress_type::Cint
    request_tilesize::NTuple{6, Clong}
    request_quantize_level::Cfloat
    request_quantize_method::Cint
    request_dither_seed::Cint
    request_lossy_int_compress::Cint
    request_huge_hdu::Cint
    request_hcomp_scale::Cfloat
    request_hcomp_smooth::Cint
    compress_type::Cint
    tilesize::NTuple{6, Clong}
    quantize_level::Cfloat
    quantize_method::Cint
    dither_seed::Cint
    compressimg::Cint
    zcmptype::NTuple{12, Cchar}
    zbitpix::Cint
    zndim::Cint
    znaxis::NTuple{6, Clong}
    maxtilelen::Clong
    maxelem::Clong
    cn_compressed::Cint
    cn_uncompressed::Cint
    cn_gzip_data::Cint
    cn_zscale::Cint
    cn_zzero::Cint
    cn_zblank::Cint
    zscale::Cdouble
    zzero::Cdouble
    cn_bscale::Cdouble
    cn_bzero::Cdouble
    cn_actual_bzero::Cdouble
    zblank::Cint
    rice_blocksize::Cint
    rice_bytepix::Cint
    hcomp_scale::Cfloat
    hcomp_smooth::Cint
    tilerow::Ptr{Cint}
    tiledatasize::Ptr{Clong}
    tiletype::Ptr{Cint}
    tiledata::Ptr{Ptr{Cvoid}}
    tilenullarray::Ptr{Ptr{Cvoid}}
    tileanynull::Ptr{Cint}
    iobuffer::Cstring
    bufrecnum::NTuple{40, Clong}
    dirty::NTuple{40, Cint}
    ageindex::NTuple{40, Cint}
end

struct fitsfile
    HDUposition::Cint
    Fptr::Ptr{FITSfile}
end

function fits_open_memfile(fptr, name, mode, buffptr, buffsize, deltasize, mem_realloc, status)
    @ccall libcfitsio.ffomem(fptr::Ptr{Ptr{fitsfile}}, name::Cstring, mode::Cint, buffptr::Ptr{Ptr{Cvoid}}, buffsize::Ptr{Csize_t}, deltasize::Csize_t, mem_realloc::Ptr{Cvoid}, status::Ptr{Status})::Status
end

function ffopentest(soname, fptr, filename, iomode, status)
    @ccall libcfitsio.ffopentest(soname::Cint, fptr::Ptr{Ptr{fitsfile}}, filename::Cstring, iomode::Cint, status::Ptr{Status})::Status
end

function fits_open_data(fptr, filename, iomode, status)
    @ccall libcfitsio.ffdopn(fptr::Ptr{Ptr{fitsfile}}, filename::Cstring, iomode::Cint, status::Ptr{Status})::Status
end

function fits_open_extlist(fptr, filename, iomode, extlist, hdutype, status)
    @ccall libcfitsio.ffeopn(fptr::Ptr{Ptr{fitsfile}}, filename::Cstring, iomode::Cint, extlist::Cstring, hdutype::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_open_table(fptr, filename, iomode, status)
    @ccall libcfitsio.fftopn(fptr::Ptr{Ptr{fitsfile}}, filename::Cstring, iomode::Cint, status::Ptr{Status})::Status
end

function fits_open_image(fptr, filename, iomode, status)
    @ccall libcfitsio.ffiopn(fptr::Ptr{Ptr{fitsfile}}, filename::Cstring, iomode::Cint, status::Ptr{Status})::Status
end

function fits_open_diskfile(fptr, filename, iomode, status)
    @ccall libcfitsio.ffdkopn(fptr::Ptr{Ptr{fitsfile}}, filename::Cstring, iomode::Cint, status::Ptr{Status})::Status
end

function fits_reopen_file(openfptr, newfptr, status)
    @ccall libcfitsio.ffreopen(openfptr::Ptr{fitsfile}, newfptr::Ptr{Ptr{fitsfile}}, status::Ptr{Status})::Status
end

function fits_create_file(fptr, filename, status)
    @ccall libcfitsio.ffinit(fptr::Ptr{Ptr{fitsfile}}, filename::Cstring, status::Ptr{Status})::Status
end

function fits_create_diskfile(fptr, filename, status)
    @ccall libcfitsio.ffdkinit(fptr::Ptr{Ptr{fitsfile}}, filename::Cstring, status::Ptr{Status})::Status
end

function fits_create_memfile(fptr, buffptr, buffsize, deltasize, mem_realloc, status)
    @ccall libcfitsio.ffimem(fptr::Ptr{Ptr{fitsfile}}, buffptr::Ptr{Ptr{Cvoid}}, buffsize::Ptr{Csize_t}, deltasize::Csize_t, mem_realloc::Ptr{Cvoid}, status::Ptr{Status})::Status
end

function fits_create_template(fptr, filename, tempname, status)
    @ccall libcfitsio.fftplt(fptr::Ptr{Ptr{fitsfile}}, filename::Cstring, tempname::Cstring, status::Ptr{Status})::Status
end

function fits_flush_file(fptr, status)
    @ccall libcfitsio.ffflus(fptr::Ptr{fitsfile}, status::Ptr{Status})::Status
end

function fits_flush_buffer(fptr, clearbuf, status)
    @ccall libcfitsio.ffflsh(fptr::Ptr{fitsfile}, clearbuf::Cint, status::Ptr{Status})::Status
end

function fits_close_file(fptr, status)
    @ccall libcfitsio.ffclos(fptr::Ptr{fitsfile}, status::Ptr{Status})::Status
end

function fits_delete_file(fptr, status)
    @ccall libcfitsio.ffdelt(fptr::Ptr{fitsfile}, status::Ptr{Status})::Status
end

function fits_file_name(fptr, filename, status)
    @ccall libcfitsio.ffflnm(fptr::Ptr{fitsfile}, filename::Cstring, status::Ptr{Status})::Status
end

function fits_file_mode(fptr, filemode, status)
    @ccall libcfitsio.ffflmd(fptr::Ptr{fitsfile}, filemode::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_url_type(fptr, urlType, status)
    @ccall libcfitsio.ffurlt(fptr::Ptr{fitsfile}, urlType::Cstring, status::Ptr{Status})::Status
end

function fits_get_version(version)
    @ccall libcfitsio.ffvers(version::Ptr{Cfloat})::Cfloat
end

function fits_uppercase(string)
    @ccall libcfitsio.ffupch(string::Cstring)::Cvoid
end

function fits_get_errstatus(status, errtext)
    @ccall libcfitsio.ffgerr(status::Status, errtext::Ptr{UInt8})::Cvoid
end

function fits_write_errmsg(err_message)
    @ccall libcfitsio.ffpmsg(err_message::Cstring)::Cvoid
end

function fits_write_errmark()
    @ccall libcfitsio.ffpmrk()::Cvoid
end

function fits_read_errmsg(err_message)
    @ccall libcfitsio.ffgmsg(err_message::Cstring)::Cint
end

function fits_clear_errmsg()
    @ccall libcfitsio.ffcmsg()::Cvoid
end

function fits_clear_errmark()
    @ccall libcfitsio.ffcmrk()::Cvoid
end

function fits_report_error(stream, status)
    @ccall libcfitsio.ffrprt(stream::Ptr{Libc.FILE}, status::Cint)::Cvoid
end

function fits_compare_str(templt, colname, casesen, match, exact)
    @ccall libcfitsio.ffcmps(templt::Cstring, colname::Cstring, casesen::Cint, match::Ptr{Cint}, exact::Ptr{Cint})::Cvoid
end

function fits_test_keyword(keyword, status)
    @ccall libcfitsio.fftkey(keyword::Cstring, status::Ptr{Status})::Status
end

function fits_test_record(card, status)
    @ccall libcfitsio.fftrec(card::Cstring, status::Ptr{Status})::Status
end

function fits_null_check(fptr, status)
    @ccall libcfitsio.ffnchk(fptr::Ptr{fitsfile}, status::Ptr{Status})::Status
end

function fits_make_keyn(keyroot, value, keyname, status)
    @ccall libcfitsio.ffkeyn(keyroot::Cstring, value::Cint, keyname::Cstring, status::Ptr{Status})::Status
end

function fits_make_nkey(value, keyroot, keyname, status)
    @ccall libcfitsio.ffnkey(value::Cint, keyroot::Cstring, keyname::Cstring, status::Ptr{Status})::Status
end

function fits_make_key(keyname, keyval, comm, card, status)
    @ccall libcfitsio.ffmkky(keyname::Cstring, keyval::Cstring, comm::Cstring, card::Cstring, status::Ptr{Status})::Status
end

function fits_get_keyclass(card)
    @ccall libcfitsio.ffgkcl(card::Cstring)::Cint
end

function fits_get_keytype(cval, dtype, status)
    @ccall libcfitsio.ffdtyp(cval::Cstring, dtype::Cstring, status::Ptr{Status})::Status
end

function fits_get_inttype(cval, datatype, negative, status)
    @ccall libcfitsio.ffinttyp(cval::Cstring, datatype::Ptr{Cint}, negative::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_parse_value(card, value, comm, status)
    @ccall libcfitsio.ffpsvc(card::Cstring, value::Cstring, comm::Cstring, status::Ptr{Status})::Status
end

function fits_get_keyname(card, name, length, status)
    @ccall libcfitsio.ffgknm(card::Cstring, name::Cstring, length::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_parse_template(tmplt, card, hdtype, status)
    @ccall libcfitsio.ffgthd(tmplt::Cstring, card::Cstring, hdtype::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_ascii_tform(tform, datacode, width, decim, status)
    @ccall libcfitsio.ffasfm(tform::Cstring, datacode::Ptr{Cint}, width::Ptr{Clong}, decim::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_binary_tform(tform, datacode, repeat, width, status)
    @ccall libcfitsio.ffbnfm(tform::Cstring, datacode::Ptr{Cint}, repeat::Ptr{Clong}, width::Ptr{Clong}, status::Ptr{Status})::Status
end

function fits_binary_tformll(tform, datacode, repeat, width, status)
    @ccall libcfitsio.ffbnfmll(tform::Cstring, datacode::Ptr{Cint}, repeat::Ptr{Clonglong}, width::Ptr{Clong}, status::Ptr{Status})::Status
end

function fits_get_tbcol(tfields, tform, space, rowlen, tbcol, status)
    @ccall libcfitsio.ffgabc(tfields::Cint, tform::Ptr{Cstring}, space::Cint, rowlen::Ptr{Clong}, tbcol::Ptr{Clong}, status::Ptr{Status})::Status
end

function fits_get_rowsize(fptr, nrows, status)
    @ccall libcfitsio.ffgrsz(fptr::Ptr{fitsfile}, nrows::Ptr{Clong}, status::Ptr{Status})::Status
end

function fits_get_col_display_width(fptr, colnum, width, status)
    @ccall libcfitsio.ffgcdw(fptr::Ptr{fitsfile}, colnum::Cint, width::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_write_record(fptr, card, status)
    @ccall libcfitsio.ffprec(fptr::Ptr{fitsfile}, card::Cstring, status::Ptr{Status})::Status
end

function fits_write_key(fptr, datatype, keyname, value, comm, status)
    @ccall libcfitsio.ffpky(fptr::Ptr{fitsfile}, datatype::Cint, keyname::Cstring, value::Ptr{Cvoid}, comm::Cstring, status::Ptr{Status})::Status
end

function fits_write_key_unit(fptr, keyname, unit, status)
    @ccall libcfitsio.ffpunt(fptr::Ptr{fitsfile}, keyname::Cstring, unit::Cstring, status::Ptr{Status})::Status
end

function fits_write_comment(fptr, comm, status)
    @ccall libcfitsio.ffpcom(fptr::Ptr{fitsfile}, comm::Cstring, status::Ptr{Status})::Status
end

function fits_write_history(fptr, history, status)
    @ccall libcfitsio.ffphis(fptr::Ptr{fitsfile}, history::Cstring, status::Ptr{Status})::Status
end

function fits_write_date(fptr, status)
    @ccall libcfitsio.ffpdat(fptr::Ptr{fitsfile}, status::Ptr{Status})::Status
end

function fits_get_system_time(timestr, timeref, status)
    @ccall libcfitsio.ffgstm(timestr::Cstring, timeref::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_get_system_date(day, month, year, status)
    @ccall libcfitsio.ffgsdt(day::Ptr{Cint}, month::Ptr{Cint}, year::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_date2str(year, month, day, datestr, status)
    @ccall libcfitsio.ffdt2s(year::Cint, month::Cint, day::Cint, datestr::Cstring, status::Ptr{Status})::Status
end

function fits_time2str(year, month, day, hour, minute, second, decimals, datestr, status)
    @ccall libcfitsio.fftm2s(year::Cint, month::Cint, day::Cint, hour::Cint, minute::Cint, second::Cdouble, decimals::Cint, datestr::Cstring, status::Ptr{Status})::Status
end

function fits_str2date(datestr, year, month, day, status)
    @ccall libcfitsio.ffs2dt(datestr::Cstring, year::Ptr{Cint}, month::Ptr{Cint}, day::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_str2time(datestr, year, month, day, hour, minute, second, status)
    @ccall libcfitsio.ffs2tm(datestr::Cstring, year::Ptr{Cint}, month::Ptr{Cint}, day::Ptr{Cint}, hour::Ptr{Cint}, minute::Ptr{Cint}, second::Ptr{Cdouble}, status::Ptr{Status})::Status
end

function fits_write_key_longstr(fptr, keyname, value, comm, status)
    @ccall libcfitsio.ffpkls(fptr::Ptr{fitsfile}, keyname::Cstring, value::Cstring, comm::Cstring, status::Ptr{Status})::Status
end

function fits_write_key_longwarn(fptr, status)
    @ccall libcfitsio.ffplsw(fptr::Ptr{fitsfile}, status::Ptr{Status})::Status
end

function fits_write_key_null(fptr, keyname, comm, status)
    @ccall libcfitsio.ffpkyu(fptr::Ptr{fitsfile}, keyname::Cstring, comm::Cstring, status::Ptr{Status})::Status
end

function fits_write_key_str(fptr, keyname, value, comm, status)
    @ccall libcfitsio.ffpkys(fptr::Ptr{fitsfile}, keyname::Cstring, value::Cstring, comm::Cstring, status::Ptr{Status})::Status
end

function fits_write_key_log(fptr, keyname, value, comm, status)
    @ccall libcfitsio.ffpkyl(fptr::Ptr{fitsfile}, keyname::Cstring, value::Cint, comm::Cstring, status::Ptr{Status})::Status
end

function fits_write_key_lng(fptr, keyname, value, comm, status)
    @ccall libcfitsio.ffpkyj(fptr::Ptr{fitsfile}, keyname::Cstring, value::Clonglong, comm::Cstring, status::Ptr{Status})::Status
end

function fits_write_key_ulng(fptr, keyname, value, comm, status)
    @ccall libcfitsio.ffpkyuj(fptr::Ptr{fitsfile}, keyname::Cstring, value::Culonglong, comm::Cstring, status::Ptr{Status})::Status
end

function fits_write_key_fixflt(fptr, keyname, value, decim, comm, status)
    @ccall libcfitsio.ffpkyf(fptr::Ptr{fitsfile}, keyname::Cstring, value::Cfloat, decim::Cint, comm::Cstring, status::Ptr{Status})::Status
end

function fits_write_key_flt(fptr, keyname, value, decim, comm, status)
    @ccall libcfitsio.ffpkye(fptr::Ptr{fitsfile}, keyname::Cstring, value::Cfloat, decim::Cint, comm::Cstring, status::Ptr{Status})::Status
end

function fits_write_key_fixdbl(fptr, keyname, value, decim, comm, status)
    @ccall libcfitsio.ffpkyg(fptr::Ptr{fitsfile}, keyname::Cstring, value::Cdouble, decim::Cint, comm::Cstring, status::Ptr{Status})::Status
end

function fits_write_key_dbl(fptr, keyname, value, decim, comm, status)
    @ccall libcfitsio.ffpkyd(fptr::Ptr{fitsfile}, keyname::Cstring, value::Cdouble, decim::Cint, comm::Cstring, status::Ptr{Status})::Status
end

function fits_write_key_fixcmp(fptr, keyname, value, decim, comm, status)
    @ccall libcfitsio.ffpkfc(fptr::Ptr{fitsfile}, keyname::Cstring, value::Ptr{Cfloat}, decim::Cint, comm::Cstring, status::Ptr{Status})::Status
end

function fits_write_key_cmp(fptr, keyname, value, decim, comm, status)
    @ccall libcfitsio.ffpkyc(fptr::Ptr{fitsfile}, keyname::Cstring, value::Ptr{Cfloat}, decim::Cint, comm::Cstring, status::Ptr{Status})::Status
end

function fits_write_key_fixdblcmp(fptr, keyname, value, decim, comm, status)
    @ccall libcfitsio.ffpkfm(fptr::Ptr{fitsfile}, keyname::Cstring, value::Ptr{Cdouble}, decim::Cint, comm::Cstring, status::Ptr{Status})::Status
end

function fits_write_key_dblcmp(fptr, keyname, value, decim, comm, status)
    @ccall libcfitsio.ffpkym(fptr::Ptr{fitsfile}, keyname::Cstring, value::Ptr{Cdouble}, decim::Cint, comm::Cstring, status::Ptr{Status})::Status
end

function fits_write_key_triple(fptr, keyname, intval, frac, comm, status)
    @ccall libcfitsio.ffpkyt(fptr::Ptr{fitsfile}, keyname::Cstring, intval::Clong, frac::Cdouble, comm::Cstring, status::Ptr{Status})::Status
end

function fits_write_tdim(fptr, colnum, naxis, naxes, status)
    @ccall libcfitsio.ffptdm(fptr::Ptr{fitsfile}, colnum::Cint, naxis::Cint, naxes::Ptr{Clong}, status::Ptr{Status})::Status
end

function fits_write_tdimll(fptr, colnum, naxis, naxes, status)
    @ccall libcfitsio.ffptdmll(fptr::Ptr{fitsfile}, colnum::Cint, naxis::Cint, naxes::Ptr{Clonglong}, status::Ptr{Status})::Status
end

function fits_write_keys_str(fptr, keyroot, nstart, nkey, value, comm, status)
    @ccall libcfitsio.ffpkns(fptr::Ptr{fitsfile}, keyroot::Cstring, nstart::Cint, nkey::Cint, value::Ptr{Cstring}, comm::Ptr{Cstring}, status::Ptr{Status})::Status
end

function fits_write_keys_log(fptr, keyroot, nstart, nkey, value, comm, status)
    @ccall libcfitsio.ffpknl(fptr::Ptr{fitsfile}, keyroot::Cstring, nstart::Cint, nkey::Cint, value::Ptr{Cint}, comm::Ptr{Cstring}, status::Ptr{Status})::Status
end

function fits_write_keys_lng(fptr, keyroot, nstart, nkey, value, comm, status)
    @ccall libcfitsio.ffpknj(fptr::Ptr{fitsfile}, keyroot::Cstring, nstart::Cint, nkey::Cint, value::Ptr{Clong}, comm::Ptr{Cstring}, status::Ptr{Status})::Status
end

function fits_write_keys_fixflt(fptr, keyroot, nstart, nkey, value, decim, comm, status)
    @ccall libcfitsio.ffpknf(fptr::Ptr{fitsfile}, keyroot::Cstring, nstart::Cint, nkey::Cint, value::Ptr{Cfloat}, decim::Cint, comm::Ptr{Cstring}, status::Ptr{Status})::Status
end

function fits_write_keys_flt(fptr, keyroot, nstart, nkey, value, decim, comm, status)
    @ccall libcfitsio.ffpkne(fptr::Ptr{fitsfile}, keyroot::Cstring, nstart::Cint, nkey::Cint, value::Ptr{Cfloat}, decim::Cint, comm::Ptr{Cstring}, status::Ptr{Status})::Status
end

function fits_write_keys_fixdbl(fptr, keyroot, nstart, nkey, value, decim, comm, status)
    @ccall libcfitsio.ffpkng(fptr::Ptr{fitsfile}, keyroot::Cstring, nstart::Cint, nkey::Cint, value::Ptr{Cdouble}, decim::Cint, comm::Ptr{Cstring}, status::Ptr{Status})::Status
end

function fits_write_keys_dbl(fptr, keyroot, nstart, nkey, value, decim, comm, status)
    @ccall libcfitsio.ffpknd(fptr::Ptr{fitsfile}, keyroot::Cstring, nstart::Cint, nkey::Cint, value::Ptr{Cdouble}, decim::Cint, comm::Ptr{Cstring}, status::Ptr{Status})::Status
end

function fits_copy_key(infptr, outfptr, incol, outcol, rootname, status)
    @ccall libcfitsio.ffcpky(infptr::Ptr{fitsfile}, outfptr::Ptr{fitsfile}, incol::Cint, outcol::Cint, rootname::Cstring, status::Ptr{Status})::Status
end

function fits_write_imghdr(fptr, bitpix, naxis, naxes, status)
    @ccall libcfitsio.ffphps(fptr::Ptr{fitsfile}, bitpix::Cint, naxis::Cint, naxes::Ptr{Clong}, status::Ptr{Status})::Status
end

function fits_write_imghdrll(fptr, bitpix, naxis, naxes, status)
    @ccall libcfitsio.ffphpsll(fptr::Ptr{fitsfile}, bitpix::Cint, naxis::Cint, naxes::Ptr{Clonglong}, status::Ptr{Status})::Status
end

function fits_write_grphdr(fptr, simple, bitpix, naxis, naxes, pcount, gcount, extend, status)
    @ccall libcfitsio.ffphpr(fptr::Ptr{fitsfile}, simple::Cint, bitpix::Cint, naxis::Cint, naxes::Ptr{Clong}, pcount::Clonglong, gcount::Clonglong, extend::Cint, status::Ptr{Status})::Status
end

function fits_write_grphdrll(fptr, simple, bitpix, naxis, naxes, pcount, gcount, extend, status)
    @ccall libcfitsio.ffphprll(fptr::Ptr{fitsfile}, simple::Cint, bitpix::Cint, naxis::Cint, naxes::Ptr{Clonglong}, pcount::Clonglong, gcount::Clonglong, extend::Cint, status::Ptr{Status})::Status
end

function fits_write_atblhdr(fptr, naxis1, naxis2, tfields, ttype, tbcol, tform, tunit, extname, status)
    @ccall libcfitsio.ffphtb(fptr::Ptr{fitsfile}, naxis1::Clonglong, naxis2::Clonglong, tfields::Cint, ttype::Ptr{Cstring}, tbcol::Ptr{Clong}, tform::Ptr{Cstring}, tunit::Ptr{Cstring}, extname::Cstring, status::Ptr{Status})::Status
end

function fits_write_btblhdr(fptr, naxis2, tfields, ttype, tform, tunit, extname, pcount, status)
    @ccall libcfitsio.ffphbn(fptr::Ptr{fitsfile}, naxis2::Clonglong, tfields::Cint, ttype::Ptr{Cstring}, tform::Ptr{Cstring}, tunit::Ptr{Cstring}, extname::Cstring, pcount::Clonglong, status::Ptr{Status})::Status
end

function fits_write_exthdr(fptr, xtension, bitpix, naxis, naxes, pcount, gcount, status)
    @ccall libcfitsio.ffphext(fptr::Ptr{fitsfile}, xtension::Cstring, bitpix::Cint, naxis::Cint, naxes::Ptr{Clong}, pcount::Clonglong, gcount::Clonglong, status::Ptr{Status})::Status
end

function fits_write_key_template(fptr, filename, status)
    @ccall libcfitsio.ffpktp(fptr::Ptr{fitsfile}, filename::Cstring, status::Ptr{Status})::Status
end

function fits_get_hdrspace(fptr, nexist, nmore, status)
    @ccall libcfitsio.ffghsp(fptr::Ptr{fitsfile}, nexist::Ptr{Cint}, nmore::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_get_hdrpos(fptr, nexist, position, status)
    @ccall libcfitsio.ffghps(fptr::Ptr{fitsfile}, nexist::Ptr{Cint}, position::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_movabs_key(fptr, nrec, status)
    @ccall libcfitsio.ffmaky(fptr::Ptr{fitsfile}, nrec::Cint, status::Ptr{Status})::Status
end

function fits_movrel_key(fptr, nrec, status)
    @ccall libcfitsio.ffmrky(fptr::Ptr{fitsfile}, nrec::Cint, status::Ptr{Status})::Status
end

function fits_find_nextkey(fptr, inclist, ninc, exclist, nexc, card, status)
    @ccall libcfitsio.ffgnxk(fptr::Ptr{fitsfile}, inclist::Ptr{Cstring}, ninc::Cint, exclist::Ptr{Cstring}, nexc::Cint, card::Cstring, status::Ptr{Status})::Status
end

function fits_read_record(fptr, nrec, card, status)
    @ccall libcfitsio.ffgrec(fptr::Ptr{fitsfile}, nrec::Cint, card::Cstring, status::Ptr{Status})::Status
end

function fits_read_card(fptr, keyname, card, status)
    @ccall libcfitsio.ffgcrd(fptr::Ptr{fitsfile}, keyname::Cstring, card::Cstring, status::Ptr{Status})::Status
end

function fits_read_str(fptr, string, card, status)
    @ccall libcfitsio.ffgstr(fptr::Ptr{fitsfile}, string::Cstring, card::Cstring, status::Ptr{Status})::Status
end

function fits_read_key_unit(fptr, keyname, unit, status)
    @ccall libcfitsio.ffgunt(fptr::Ptr{fitsfile}, keyname::Cstring, unit::Cstring, status::Ptr{Status})::Status
end

function fits_read_keyn(fptr, nkey, keyname, keyval, comm, status)
    @ccall libcfitsio.ffgkyn(fptr::Ptr{fitsfile}, nkey::Cint, keyname::Cstring, keyval::Cstring, comm::Cstring, status::Ptr{Status})::Status
end

function fits_read_key(fptr, datatype, keyname, value, comm, status)
    @ccall libcfitsio.ffgky(fptr::Ptr{fitsfile}, datatype::Cint, keyname::Cstring, value::Ptr{Cvoid}, comm::Cstring, status::Ptr{Status})::Status
end

function fits_read_keyword(fptr, keyname, keyval, comm, status)
    @ccall libcfitsio.ffgkey(fptr::Ptr{fitsfile}, keyname::Cstring, keyval::Cstring, comm::Cstring, status::Ptr{Status})::Status
end

function fits_read_key_str(fptr, keyname, value, comm, status)
    @ccall libcfitsio.ffgkys(fptr::Ptr{fitsfile}, keyname::Cstring, value::Cstring, comm::Cstring, status::Ptr{Status})::Status
end

function fits_read_key_log(fptr, keyname, value, comm, status)
    @ccall libcfitsio.ffgkyl(fptr::Ptr{fitsfile}, keyname::Cstring, value::Ptr{Cint}, comm::Cstring, status::Ptr{Status})::Status
end

function fits_read_key_lng(fptr, keyname, value, comm, status)
    @ccall libcfitsio.ffgkyj(fptr::Ptr{fitsfile}, keyname::Cstring, value::Ptr{Clong}, comm::Cstring, status::Ptr{Status})::Status
end

function fits_read_key_lnglng(fptr, keyname, value, comm, status)
    @ccall libcfitsio.ffgkyjj(fptr::Ptr{fitsfile}, keyname::Cstring, value::Ptr{Clonglong}, comm::Cstring, status::Ptr{Status})::Status
end

function fits_read_key_ulnglng(fptr, keyname, value, comm, status)
    @ccall libcfitsio.ffgkyujj(fptr::Ptr{fitsfile}, keyname::Cstring, value::Ptr{Culonglong}, comm::Cstring, status::Ptr{Status})::Status
end

function fits_read_key_flt(fptr, keyname, value, comm, status)
    @ccall libcfitsio.ffgkye(fptr::Ptr{fitsfile}, keyname::Cstring, value::Ptr{Cfloat}, comm::Cstring, status::Ptr{Status})::Status
end

function fits_read_key_dbl(fptr, keyname, value, comm, status)
    @ccall libcfitsio.ffgkyd(fptr::Ptr{fitsfile}, keyname::Cstring, value::Ptr{Cdouble}, comm::Cstring, status::Ptr{Status})::Status
end

function fits_read_key_cmp(fptr, keyname, value, comm, status)
    @ccall libcfitsio.ffgkyc(fptr::Ptr{fitsfile}, keyname::Cstring, value::Ptr{Cfloat}, comm::Cstring, status::Ptr{Status})::Status
end

function fits_read_key_dblcmp(fptr, keyname, value, comm, status)
    @ccall libcfitsio.ffgkym(fptr::Ptr{fitsfile}, keyname::Cstring, value::Ptr{Cdouble}, comm::Cstring, status::Ptr{Status})::Status
end

function fits_read_key_triple(fptr, keyname, ivalue, dvalue, comm, status)
    @ccall libcfitsio.ffgkyt(fptr::Ptr{fitsfile}, keyname::Cstring, ivalue::Ptr{Clong}, dvalue::Ptr{Cdouble}, comm::Cstring, status::Ptr{Status})::Status
end

function fits_get_key_strlen(fptr, keyname, length, status)
    @ccall libcfitsio.ffgksl(fptr::Ptr{fitsfile}, keyname::Cstring, length::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_get_key_com_strlen(fptr, keyname, length, comlength, status)
    @ccall libcfitsio.ffgkcsl(fptr::Ptr{fitsfile}, keyname::Cstring, length::Ptr{Cint}, comlength::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_key_longstr(fptr, keyname, value, comm, status)
    @ccall libcfitsio.ffgkls(fptr::Ptr{fitsfile}, keyname::Cstring, value::Ptr{Cstring}, comm::Cstring, status::Ptr{Status})::Status
end

function fits_read_string_key(fptr, keyname, firstchar, maxchar, value, valuelen, comm, status)
    @ccall libcfitsio.ffgsky(fptr::Ptr{fitsfile}, keyname::Cstring, firstchar::Cint, maxchar::Cint, value::Cstring, valuelen::Ptr{Cint}, comm::Cstring, status::Ptr{Status})::Status
end

function fits_read_string_key_com(fptr, keyname, firstchar, maxchar, maxcomchar, value, valuelen, comm, comlen, status)
    @ccall libcfitsio.ffgskyc(fptr::Ptr{fitsfile}, keyname::Cstring, firstchar::Cint, maxchar::Cint, maxcomchar::Cint, value::Cstring, valuelen::Ptr{Cint}, comm::Cstring, comlen::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_free_memory(value, status)
    @ccall libcfitsio.fffree(value::Ptr{Cvoid}, status::Ptr{Status})::Status
end

function fits_read_tdim(fptr, colnum, maxdim, naxis, naxes, status)
    @ccall libcfitsio.ffgtdm(fptr::Ptr{fitsfile}, colnum::Cint, maxdim::Cint, naxis::Ptr{Cint}, naxes::Ptr{Clong}, status::Ptr{Status})::Status
end

function fits_read_tdimll(fptr, colnum, maxdim, naxis, naxes, status)
    @ccall libcfitsio.ffgtdmll(fptr::Ptr{fitsfile}, colnum::Cint, maxdim::Cint, naxis::Ptr{Cint}, naxes::Ptr{Clonglong}, status::Ptr{Status})::Status
end

function fits_decode_tdim(fptr, tdimstr, colnum, maxdim, naxis, naxes, status)
    @ccall libcfitsio.ffdtdm(fptr::Ptr{fitsfile}, tdimstr::Cstring, colnum::Cint, maxdim::Cint, naxis::Ptr{Cint}, naxes::Ptr{Clong}, status::Ptr{Status})::Status
end

function fits_decode_tdimll(fptr, tdimstr, colnum, maxdim, naxis, naxes, status)
    @ccall libcfitsio.ffdtdmll(fptr::Ptr{fitsfile}, tdimstr::Cstring, colnum::Cint, maxdim::Cint, naxis::Ptr{Cint}, naxes::Ptr{Clonglong}, status::Ptr{Status})::Status
end

function fits_read_keys_str(fptr, keyname, nstart, nmax, value, nfound, status)
    @ccall libcfitsio.ffgkns(fptr::Ptr{fitsfile}, keyname::Cstring, nstart::Cint, nmax::Cint, value::Ptr{Cstring}, nfound::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_keys_log(fptr, keyname, nstart, nmax, value, nfound, status)
    @ccall libcfitsio.ffgknl(fptr::Ptr{fitsfile}, keyname::Cstring, nstart::Cint, nmax::Cint, value::Ptr{Cint}, nfound::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_keys_lng(fptr, keyname, nstart, nmax, value, nfound, status)
    @ccall libcfitsio.ffgknj(fptr::Ptr{fitsfile}, keyname::Cstring, nstart::Cint, nmax::Cint, value::Ptr{Clong}, nfound::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_keys_lnglng(fptr, keyname, nstart, nmax, value, nfound, status)
    @ccall libcfitsio.ffgknjj(fptr::Ptr{fitsfile}, keyname::Cstring, nstart::Cint, nmax::Cint, value::Ptr{Clonglong}, nfound::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_keys_flt(fptr, keyname, nstart, nmax, value, nfound, status)
    @ccall libcfitsio.ffgkne(fptr::Ptr{fitsfile}, keyname::Cstring, nstart::Cint, nmax::Cint, value::Ptr{Cfloat}, nfound::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_keys_dbl(fptr, keyname, nstart, nmax, value, nfound, status)
    @ccall libcfitsio.ffgknd(fptr::Ptr{fitsfile}, keyname::Cstring, nstart::Cint, nmax::Cint, value::Ptr{Cdouble}, nfound::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_imghdr(fptr, maxdim, simple, bitpix, naxis, naxes, pcount, gcount, extend, status)
    @ccall libcfitsio.ffghpr(fptr::Ptr{fitsfile}, maxdim::Cint, simple::Ptr{Cint}, bitpix::Ptr{Cint}, naxis::Ptr{Cint}, naxes::Ptr{Clong}, pcount::Ptr{Clong}, gcount::Ptr{Clong}, extend::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_imghdrll(fptr, maxdim, simple, bitpix, naxis, naxes, pcount, gcount, extend, status)
    @ccall libcfitsio.ffghprll(fptr::Ptr{fitsfile}, maxdim::Cint, simple::Ptr{Cint}, bitpix::Ptr{Cint}, naxis::Ptr{Cint}, naxes::Ptr{Clonglong}, pcount::Ptr{Clong}, gcount::Ptr{Clong}, extend::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_atblhdr(fptr, maxfield, naxis1, naxis2, tfields, ttype, tbcol, tform, tunit, extname, status)
    @ccall libcfitsio.ffghtb(fptr::Ptr{fitsfile}, maxfield::Cint, naxis1::Ptr{Clong}, naxis2::Ptr{Clong}, tfields::Ptr{Cint}, ttype::Ptr{Cstring}, tbcol::Ptr{Clong}, tform::Ptr{Cstring}, tunit::Ptr{Cstring}, extname::Cstring, status::Ptr{Status})::Status
end

function fits_read_btblhdr(fptr, maxfield, naxis2, tfields, ttype, tform, tunit, extname, pcount, status)
    @ccall libcfitsio.ffghbn(fptr::Ptr{fitsfile}, maxfield::Cint, naxis2::Ptr{Clong}, tfields::Ptr{Cint}, ttype::Ptr{Cstring}, tform::Ptr{Cstring}, tunit::Ptr{Cstring}, extname::Cstring, pcount::Ptr{Clong}, status::Ptr{Status})::Status
end

function fits_read_atblhdrll(fptr, maxfield, naxis1, naxis2, tfields, ttype, tbcol, tform, tunit, extname, status)
    @ccall libcfitsio.ffghtbll(fptr::Ptr{fitsfile}, maxfield::Cint, naxis1::Ptr{Clonglong}, naxis2::Ptr{Clonglong}, tfields::Ptr{Cint}, ttype::Ptr{Cstring}, tbcol::Ptr{Clonglong}, tform::Ptr{Cstring}, tunit::Ptr{Cstring}, extname::Cstring, status::Ptr{Status})::Status
end

function fits_read_btblhdrll(fptr, maxfield, naxis2, tfields, ttype, tform, tunit, extname, pcount, status)
    @ccall libcfitsio.ffghbnll(fptr::Ptr{fitsfile}, maxfield::Cint, naxis2::Ptr{Clonglong}, tfields::Ptr{Cint}, ttype::Ptr{Cstring}, tform::Ptr{Cstring}, tunit::Ptr{Cstring}, extname::Cstring, pcount::Ptr{Clonglong}, status::Ptr{Status})::Status
end

function fits_hdr2str(fptr, exclude_comm, exclist, nexc, header, nkeys, status)
    @ccall libcfitsio.ffhdr2str(fptr::Ptr{fitsfile}, exclude_comm::Cint, exclist::Ptr{Cstring}, nexc::Cint, header::Ptr{Cstring}, nkeys::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_convert_hdr2str(fptr, exclude_comm, exclist, nexc, header, nkeys, status)
    @ccall libcfitsio.ffcnvthdr2str(fptr::Ptr{fitsfile}, exclude_comm::Cint, exclist::Ptr{Cstring}, nexc::Cint, header::Ptr{Cstring}, nkeys::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_update_card(fptr, keyname, card, status)
    @ccall libcfitsio.ffucrd(fptr::Ptr{fitsfile}, keyname::Cstring, card::Cstring, status::Ptr{Status})::Status
end

function fits_update_key(fptr, datatype, keyname, value, comm, status)
    @ccall libcfitsio.ffuky(fptr::Ptr{fitsfile}, datatype::Cint, keyname::Cstring, value::Ptr{Cvoid}, comm::Cstring, status::Ptr{Status})::Status
end

function fits_update_key_null(fptr, keyname, comm, status)
    @ccall libcfitsio.ffukyu(fptr::Ptr{fitsfile}, keyname::Cstring, comm::Cstring, status::Ptr{Status})::Status
end

function fits_update_key_str(fptr, keyname, value, comm, status)
    @ccall libcfitsio.ffukys(fptr::Ptr{fitsfile}, keyname::Cstring, value::Cstring, comm::Cstring, status::Ptr{Status})::Status
end

function fits_update_key_longstr(fptr, keyname, value, comm, status)
    @ccall libcfitsio.ffukls(fptr::Ptr{fitsfile}, keyname::Cstring, value::Cstring, comm::Cstring, status::Ptr{Status})::Status
end

function fits_update_key_log(fptr, keyname, value, comm, status)
    @ccall libcfitsio.ffukyl(fptr::Ptr{fitsfile}, keyname::Cstring, value::Cint, comm::Cstring, status::Ptr{Status})::Status
end

function fits_update_key_lng(fptr, keyname, value, comm, status)
    @ccall libcfitsio.ffukyj(fptr::Ptr{fitsfile}, keyname::Cstring, value::Clonglong, comm::Cstring, status::Ptr{Status})::Status
end

function fits_update_key_ulng(fptr, keyname, value, comm, status)
    @ccall libcfitsio.ffukyuj(fptr::Ptr{fitsfile}, keyname::Cstring, value::Culonglong, comm::Cstring, status::Ptr{Status})::Status
end

function fits_update_key_fixflt(fptr, keyname, value, decim, comm, status)
    @ccall libcfitsio.ffukyf(fptr::Ptr{fitsfile}, keyname::Cstring, value::Cfloat, decim::Cint, comm::Cstring, status::Ptr{Status})::Status
end

function fits_update_key_flt(fptr, keyname, value, decim, comm, status)
    @ccall libcfitsio.ffukye(fptr::Ptr{fitsfile}, keyname::Cstring, value::Cfloat, decim::Cint, comm::Cstring, status::Ptr{Status})::Status
end

function fits_update_key_fixdbl(fptr, keyname, value, decim, comm, status)
    @ccall libcfitsio.ffukyg(fptr::Ptr{fitsfile}, keyname::Cstring, value::Cdouble, decim::Cint, comm::Cstring, status::Ptr{Status})::Status
end

function fits_update_key_dbl(fptr, keyname, value, decim, comm, status)
    @ccall libcfitsio.ffukyd(fptr::Ptr{fitsfile}, keyname::Cstring, value::Cdouble, decim::Cint, comm::Cstring, status::Ptr{Status})::Status
end

function fits_update_key_fixcmp(fptr, keyname, value, decim, comm, status)
    @ccall libcfitsio.ffukfc(fptr::Ptr{fitsfile}, keyname::Cstring, value::Ptr{Cfloat}, decim::Cint, comm::Cstring, status::Ptr{Status})::Status
end

function fits_update_key_cmp(fptr, keyname, value, decim, comm, status)
    @ccall libcfitsio.ffukyc(fptr::Ptr{fitsfile}, keyname::Cstring, value::Ptr{Cfloat}, decim::Cint, comm::Cstring, status::Ptr{Status})::Status
end

function fits_update_key_fixdblcmp(fptr, keyname, value, decim, comm, status)
    @ccall libcfitsio.ffukfm(fptr::Ptr{fitsfile}, keyname::Cstring, value::Ptr{Cdouble}, decim::Cint, comm::Cstring, status::Ptr{Status})::Status
end

function fits_update_key_dblcmp(fptr, keyname, value, decim, comm, status)
    @ccall libcfitsio.ffukym(fptr::Ptr{fitsfile}, keyname::Cstring, value::Ptr{Cdouble}, decim::Cint, comm::Cstring, status::Ptr{Status})::Status
end

function fits_modify_record(fptr, nkey, card, status)
    @ccall libcfitsio.ffmrec(fptr::Ptr{fitsfile}, nkey::Cint, card::Cstring, status::Ptr{Status})::Status
end

function fits_modify_card(fptr, keyname, card, status)
    @ccall libcfitsio.ffmcrd(fptr::Ptr{fitsfile}, keyname::Cstring, card::Cstring, status::Ptr{Status})::Status
end

function fits_modify_name(fptr, oldname, newname, status)
    @ccall libcfitsio.ffmnam(fptr::Ptr{fitsfile}, oldname::Cstring, newname::Cstring, status::Ptr{Status})::Status
end

function fits_modify_comment(fptr, keyname, comm, status)
    @ccall libcfitsio.ffmcom(fptr::Ptr{fitsfile}, keyname::Cstring, comm::Cstring, status::Ptr{Status})::Status
end

function fits_modify_key_null(fptr, keyname, comm, status)
    @ccall libcfitsio.ffmkyu(fptr::Ptr{fitsfile}, keyname::Cstring, comm::Cstring, status::Ptr{Status})::Status
end

function fits_modify_key_str(fptr, keyname, value, comm, status)
    @ccall libcfitsio.ffmkys(fptr::Ptr{fitsfile}, keyname::Cstring, value::Cstring, comm::Cstring, status::Ptr{Status})::Status
end

function fits_modify_key_longstr(fptr, keyname, value, comm, status)
    @ccall libcfitsio.ffmkls(fptr::Ptr{fitsfile}, keyname::Cstring, value::Cstring, comm::Cstring, status::Ptr{Status})::Status
end

function fits_modify_key_log(fptr, keyname, value, comm, status)
    @ccall libcfitsio.ffmkyl(fptr::Ptr{fitsfile}, keyname::Cstring, value::Cint, comm::Cstring, status::Ptr{Status})::Status
end

function fits_modify_key_lng(fptr, keyname, value, comm, status)
    @ccall libcfitsio.ffmkyj(fptr::Ptr{fitsfile}, keyname::Cstring, value::Clonglong, comm::Cstring, status::Ptr{Status})::Status
end

function fits_modify_key_ulng(fptr, keyname, value, comm, status)
    @ccall libcfitsio.ffmkyuj(fptr::Ptr{fitsfile}, keyname::Cstring, value::Culonglong, comm::Cstring, status::Ptr{Status})::Status
end

function fits_modify_key_fixflt(fptr, keyname, value, decim, comm, status)
    @ccall libcfitsio.ffmkyf(fptr::Ptr{fitsfile}, keyname::Cstring, value::Cfloat, decim::Cint, comm::Cstring, status::Ptr{Status})::Status
end

function fits_modify_key_flt(fptr, keyname, value, decim, comm, status)
    @ccall libcfitsio.ffmkye(fptr::Ptr{fitsfile}, keyname::Cstring, value::Cfloat, decim::Cint, comm::Cstring, status::Ptr{Status})::Status
end

function fits_modify_key_fixdbl(fptr, keyname, value, decim, comm, status)
    @ccall libcfitsio.ffmkyg(fptr::Ptr{fitsfile}, keyname::Cstring, value::Cdouble, decim::Cint, comm::Cstring, status::Ptr{Status})::Status
end

function fits_modify_key_dbl(fptr, keyname, value, decim, comm, status)
    @ccall libcfitsio.ffmkyd(fptr::Ptr{fitsfile}, keyname::Cstring, value::Cdouble, decim::Cint, comm::Cstring, status::Ptr{Status})::Status
end

function fits_modify_key_fixcmp(fptr, keyname, value, decim, comm, status)
    @ccall libcfitsio.ffmkfc(fptr::Ptr{fitsfile}, keyname::Cstring, value::Ptr{Cfloat}, decim::Cint, comm::Cstring, status::Ptr{Status})::Status
end

function fits_modify_key_cmp(fptr, keyname, value, decim, comm, status)
    @ccall libcfitsio.ffmkyc(fptr::Ptr{fitsfile}, keyname::Cstring, value::Ptr{Cfloat}, decim::Cint, comm::Cstring, status::Ptr{Status})::Status
end

function fits_modify_key_fixdblcmp(fptr, keyname, value, decim, comm, status)
    @ccall libcfitsio.ffmkfm(fptr::Ptr{fitsfile}, keyname::Cstring, value::Ptr{Cdouble}, decim::Cint, comm::Cstring, status::Ptr{Status})::Status
end

function fits_modify_key_dblcmp(fptr, keyname, value, decim, comm, status)
    @ccall libcfitsio.ffmkym(fptr::Ptr{fitsfile}, keyname::Cstring, value::Ptr{Cdouble}, decim::Cint, comm::Cstring, status::Ptr{Status})::Status
end

function fits_insert_record(fptr, nkey, card, status)
    @ccall libcfitsio.ffirec(fptr::Ptr{fitsfile}, nkey::Cint, card::Cstring, status::Ptr{Status})::Status
end

function fits_insert_card(fptr, card, status)
    @ccall libcfitsio.ffikey(fptr::Ptr{fitsfile}, card::Cstring, status::Ptr{Status})::Status
end

function fits_insert_key_null(fptr, keyname, comm, status)
    @ccall libcfitsio.ffikyu(fptr::Ptr{fitsfile}, keyname::Cstring, comm::Cstring, status::Ptr{Status})::Status
end

function fits_insert_key_str(fptr, keyname, value, comm, status)
    @ccall libcfitsio.ffikys(fptr::Ptr{fitsfile}, keyname::Cstring, value::Cstring, comm::Cstring, status::Ptr{Status})::Status
end

function fits_insert_key_longstr(fptr, keyname, value, comm, status)
    @ccall libcfitsio.ffikls(fptr::Ptr{fitsfile}, keyname::Cstring, value::Cstring, comm::Cstring, status::Ptr{Status})::Status
end

function fits_insert_key_log(fptr, keyname, value, comm, status)
    @ccall libcfitsio.ffikyl(fptr::Ptr{fitsfile}, keyname::Cstring, value::Cint, comm::Cstring, status::Ptr{Status})::Status
end

function fits_insert_key_lng(fptr, keyname, value, comm, status)
    @ccall libcfitsio.ffikyj(fptr::Ptr{fitsfile}, keyname::Cstring, value::Clonglong, comm::Cstring, status::Ptr{Status})::Status
end

function fits_insert_key_fixflt(fptr, keyname, value, decim, comm, status)
    @ccall libcfitsio.ffikyf(fptr::Ptr{fitsfile}, keyname::Cstring, value::Cfloat, decim::Cint, comm::Cstring, status::Ptr{Status})::Status
end

function fits_insert_key_flt(fptr, keyname, value, decim, comm, status)
    @ccall libcfitsio.ffikye(fptr::Ptr{fitsfile}, keyname::Cstring, value::Cfloat, decim::Cint, comm::Cstring, status::Ptr{Status})::Status
end

function fits_insert_key_fixdbl(fptr, keyname, value, decim, comm, status)
    @ccall libcfitsio.ffikyg(fptr::Ptr{fitsfile}, keyname::Cstring, value::Cdouble, decim::Cint, comm::Cstring, status::Ptr{Status})::Status
end

function fits_insert_key_dbl(fptr, keyname, value, decim, comm, status)
    @ccall libcfitsio.ffikyd(fptr::Ptr{fitsfile}, keyname::Cstring, value::Cdouble, decim::Cint, comm::Cstring, status::Ptr{Status})::Status
end

function fits_insert_key_fixcmp(fptr, keyname, value, decim, comm, status)
    @ccall libcfitsio.ffikfc(fptr::Ptr{fitsfile}, keyname::Cstring, value::Ptr{Cfloat}, decim::Cint, comm::Cstring, status::Ptr{Status})::Status
end

function fits_insert_key_cmp(fptr, keyname, value, decim, comm, status)
    @ccall libcfitsio.ffikyc(fptr::Ptr{fitsfile}, keyname::Cstring, value::Ptr{Cfloat}, decim::Cint, comm::Cstring, status::Ptr{Status})::Status
end

function fits_insert_key_fixdblcmp(fptr, keyname, value, decim, comm, status)
    @ccall libcfitsio.ffikfm(fptr::Ptr{fitsfile}, keyname::Cstring, value::Ptr{Cdouble}, decim::Cint, comm::Cstring, status::Ptr{Status})::Status
end

function fits_insert_key_dblcmp(fptr, keyname, value, decim, comm, status)
    @ccall libcfitsio.ffikym(fptr::Ptr{fitsfile}, keyname::Cstring, value::Ptr{Cdouble}, decim::Cint, comm::Cstring, status::Ptr{Status})::Status
end

function fits_delete_key(fptr, keyname, status)
    @ccall libcfitsio.ffdkey(fptr::Ptr{fitsfile}, keyname::Cstring, status::Ptr{Status})::Status
end

function fits_delete_str(fptr, string, status)
    @ccall libcfitsio.ffdstr(fptr::Ptr{fitsfile}, string::Cstring, status::Ptr{Status})::Status
end

function fits_delete_record(fptr, keypos, status)
    @ccall libcfitsio.ffdrec(fptr::Ptr{fitsfile}, keypos::Cint, status::Ptr{Status})::Status
end

function fits_get_hdu_num(fptr, chdunum)
    @ccall libcfitsio.ffghdn(fptr::Ptr{fitsfile}, chdunum::Ptr{Cint})::Cint
end

function fits_get_hdu_type(fptr, exttype, status)
    @ccall libcfitsio.ffghdt(fptr::Ptr{fitsfile}, exttype::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_get_hduaddr(fptr, headstart, datastart, dataend, status)
    @ccall libcfitsio.ffghad(fptr::Ptr{fitsfile}, headstart::Ptr{Clong}, datastart::Ptr{Clong}, dataend::Ptr{Clong}, status::Ptr{Status})::Status
end

function fits_get_hduaddrll(fptr, headstart, datastart, dataend, status)
    @ccall libcfitsio.ffghadll(fptr::Ptr{fitsfile}, headstart::Ptr{Clonglong}, datastart::Ptr{Clonglong}, dataend::Ptr{Clonglong}, status::Ptr{Status})::Status
end

function fits_get_hduoff(fptr, headstart, datastart, dataend, status)
    @ccall libcfitsio.ffghof(fptr::Ptr{fitsfile}, headstart::Ptr{off_t}, datastart::Ptr{off_t}, dataend::Ptr{off_t}, status::Ptr{Status})::Status
end

function fits_get_img_param(fptr, maxaxis, imgtype, naxis, naxes, status)
    @ccall libcfitsio.ffgipr(fptr::Ptr{fitsfile}, maxaxis::Cint, imgtype::Ptr{Cint}, naxis::Ptr{Cint}, naxes::Ptr{Clong}, status::Ptr{Status})::Status
end

function fits_get_img_paramll(fptr, maxaxis, imgtype, naxis, naxes, status)
    @ccall libcfitsio.ffgiprll(fptr::Ptr{fitsfile}, maxaxis::Cint, imgtype::Ptr{Cint}, naxis::Ptr{Cint}, naxes::Ptr{Clonglong}, status::Ptr{Status})::Status
end

function fits_get_img_type(fptr, imgtype, status)
    @ccall libcfitsio.ffgidt(fptr::Ptr{fitsfile}, imgtype::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_get_img_equivtype(fptr, imgtype, status)
    @ccall libcfitsio.ffgiet(fptr::Ptr{fitsfile}, imgtype::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_get_img_dim(fptr, naxis, status)
    @ccall libcfitsio.ffgidm(fptr::Ptr{fitsfile}, naxis::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_get_img_size(fptr, nlen, naxes, status)
    @ccall libcfitsio.ffgisz(fptr::Ptr{fitsfile}, nlen::Cint, naxes::Ptr{Clong}, status::Ptr{Status})::Status
end

function fits_get_img_sizell(fptr, nlen, naxes, status)
    @ccall libcfitsio.ffgiszll(fptr::Ptr{fitsfile}, nlen::Cint, naxes::Ptr{Clonglong}, status::Ptr{Status})::Status
end

function fits_movabs_hdu(fptr, hdunum, exttype, status)
    @ccall libcfitsio.ffmahd(fptr::Ptr{fitsfile}, hdunum::Cint, exttype::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_movrel_hdu(fptr, hdumov, exttype, status)
    @ccall libcfitsio.ffmrhd(fptr::Ptr{fitsfile}, hdumov::Cint, exttype::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_movnam_hdu(fptr, exttype, hduname, hduvers, status)
    @ccall libcfitsio.ffmnhd(fptr::Ptr{fitsfile}, exttype::Cint, hduname::Cstring, hduvers::Cint, status::Ptr{Status})::Status
end

function fits_get_num_hdus(fptr, nhdu, status)
    @ccall libcfitsio.ffthdu(fptr::Ptr{fitsfile}, nhdu::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_create_img(fptr, bitpix, naxis, naxes, status)
    @ccall libcfitsio.ffcrim(fptr::Ptr{fitsfile}, bitpix::Cint, naxis::Cint, naxes::Ptr{Clong}, status::Ptr{Status})::Status
end

function fits_create_imgll(fptr, bitpix, naxis, naxes, status)
    @ccall libcfitsio.ffcrimll(fptr::Ptr{fitsfile}, bitpix::Cint, naxis::Cint, naxes::Ptr{Clonglong}, status::Ptr{Status})::Status
end

function fits_create_tbl(fptr, tbltype, naxis2, tfields, ttype, tform, tunit, extname, status)
    @ccall libcfitsio.ffcrtb(fptr::Ptr{fitsfile}, tbltype::Cint, naxis2::Clonglong, tfields::Cint, ttype::Ptr{Cstring}, tform::Ptr{Cstring}, tunit::Ptr{Cstring}, extname::Cstring, status::Ptr{Status})::Status
end

function fits_create_hdu(fptr, status)
    @ccall libcfitsio.ffcrhd(fptr::Ptr{fitsfile}, status::Ptr{Status})::Status
end

function fits_insert_img(fptr, bitpix, naxis, naxes, status)
    @ccall libcfitsio.ffiimg(fptr::Ptr{fitsfile}, bitpix::Cint, naxis::Cint, naxes::Ptr{Clong}, status::Ptr{Status})::Status
end

function fits_insert_imgll(fptr, bitpix, naxis, naxes, status)
    @ccall libcfitsio.ffiimgll(fptr::Ptr{fitsfile}, bitpix::Cint, naxis::Cint, naxes::Ptr{Clonglong}, status::Ptr{Status})::Status
end

function fits_insert_atbl(fptr, naxis1, naxis2, tfields, ttype, tbcol, tform, tunit, extname, status)
    @ccall libcfitsio.ffitab(fptr::Ptr{fitsfile}, naxis1::Clonglong, naxis2::Clonglong, tfields::Cint, ttype::Ptr{Cstring}, tbcol::Ptr{Clong}, tform::Ptr{Cstring}, tunit::Ptr{Cstring}, extname::Cstring, status::Ptr{Status})::Status
end

function fits_insert_btbl(fptr, naxis2, tfields, ttype, tform, tunit, extname, pcount, status)
    @ccall libcfitsio.ffibin(fptr::Ptr{fitsfile}, naxis2::Clonglong, tfields::Cint, ttype::Ptr{Cstring}, tform::Ptr{Cstring}, tunit::Ptr{Cstring}, extname::Cstring, pcount::Clonglong, status::Ptr{Status})::Status
end

function fits_resize_img(fptr, bitpix, naxis, naxes, status)
    @ccall libcfitsio.ffrsim(fptr::Ptr{fitsfile}, bitpix::Cint, naxis::Cint, naxes::Ptr{Clong}, status::Ptr{Status})::Status
end

function fits_resize_imgll(fptr, bitpix, naxis, naxes, status)
    @ccall libcfitsio.ffrsimll(fptr::Ptr{fitsfile}, bitpix::Cint, naxis::Cint, naxes::Ptr{Clonglong}, status::Ptr{Status})::Status
end

function fits_delete_hdu(fptr, hdutype, status)
    @ccall libcfitsio.ffdhdu(fptr::Ptr{fitsfile}, hdutype::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_copy_hdu(infptr, outfptr, morekeys, status)
    @ccall libcfitsio.ffcopy(infptr::Ptr{fitsfile}, outfptr::Ptr{fitsfile}, morekeys::Cint, status::Ptr{Status})::Status
end

function fits_copy_file(infptr, outfptr, prev, cur, follow, status)
    @ccall libcfitsio.ffcpfl(infptr::Ptr{fitsfile}, outfptr::Ptr{fitsfile}, prev::Cint, cur::Cint, follow::Cint, status::Ptr{Status})::Status
end

function fits_copy_header(infptr, outfptr, status)
    @ccall libcfitsio.ffcphd(infptr::Ptr{fitsfile}, outfptr::Ptr{fitsfile}, status::Ptr{Status})::Status
end

function fits_copy_hdutab(infptr, outfptr, firstrow, nrows, status)
    @ccall libcfitsio.ffcpht(infptr::Ptr{fitsfile}, outfptr::Ptr{fitsfile}, firstrow::Clonglong, nrows::Clonglong, status::Ptr{Status})::Status
end

function fits_copy_data(infptr, outfptr, status)
    @ccall libcfitsio.ffcpdt(infptr::Ptr{fitsfile}, outfptr::Ptr{fitsfile}, status::Ptr{Status})::Status
end

function fits_write_hdu(fptr, outstream, status)
    @ccall libcfitsio.ffwrhdu(fptr::Ptr{fitsfile}, outstream::Ptr{Libc.FILE}, status::Ptr{Status})::Status
end

function fits_set_hdustruc(fptr, status)
    @ccall libcfitsio.ffrdef(fptr::Ptr{fitsfile}, status::Ptr{Status})::Status
end

function fits_set_hdrsize(fptr, morekeys, status)
    @ccall libcfitsio.ffhdef(fptr::Ptr{fitsfile}, morekeys::Cint, status::Ptr{Status})::Status
end

function fits_write_theap(fptr, theap, status)
    @ccall libcfitsio.ffpthp(fptr::Ptr{fitsfile}, theap::Clong, status::Ptr{Status})::Status
end

function fits_encode_chksum(sum, complm, ascii)
    @ccall libcfitsio.ffesum(sum::Culong, complm::Cint, ascii::Cstring)::Cvoid
end

function fits_decode_chksum(ascii, complm, sum)
    @ccall libcfitsio.ffdsum(ascii::Cstring, complm::Cint, sum::Ptr{Culong})::Culong
end

function fits_write_chksum(fptr, status)
    @ccall libcfitsio.ffpcks(fptr::Ptr{fitsfile}, status::Ptr{Status})::Status
end

function fits_update_chksum(fptr, status)
    @ccall libcfitsio.ffupck(fptr::Ptr{fitsfile}, status::Ptr{Status})::Status
end

function fits_verify_chksum(fptr, datastatus, hdustatus, status)
    @ccall libcfitsio.ffvcks(fptr::Ptr{fitsfile}, datastatus::Ptr{Cint}, hdustatus::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_get_chksum(fptr, datasum, hdusum, status)
    @ccall libcfitsio.ffgcks(fptr::Ptr{fitsfile}, datasum::Ptr{Culong}, hdusum::Ptr{Culong}, status::Ptr{Status})::Status
end

function fits_set_bscale(fptr, scale, zeroval, status)
    @ccall libcfitsio.ffpscl(fptr::Ptr{fitsfile}, scale::Cdouble, zeroval::Cdouble, status::Ptr{Status})::Status
end

function fits_set_tscale(fptr, colnum, scale, zeroval, status)
    @ccall libcfitsio.fftscl(fptr::Ptr{fitsfile}, colnum::Cint, scale::Cdouble, zeroval::Cdouble, status::Ptr{Status})::Status
end

function fits_set_imgnull(fptr, nulvalue, status)
    @ccall libcfitsio.ffpnul(fptr::Ptr{fitsfile}, nulvalue::Clonglong, status::Ptr{Status})::Status
end

function fits_set_btblnull(fptr, colnum, nulvalue, status)
    @ccall libcfitsio.fftnul(fptr::Ptr{fitsfile}, colnum::Cint, nulvalue::Clonglong, status::Ptr{Status})::Status
end

function fits_set_atblnull(fptr, colnum, nulstring, status)
    @ccall libcfitsio.ffsnul(fptr::Ptr{fitsfile}, colnum::Cint, nulstring::Cstring, status::Ptr{Status})::Status
end

function fits_get_colnum(fptr, casesen, templt, colnum, status)
    @ccall libcfitsio.ffgcno(fptr::Ptr{fitsfile}, casesen::Cint, templt::Cstring, colnum::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_get_colname(fptr, casesen, templt, colname, colnum, status)
    @ccall libcfitsio.ffgcnn(fptr::Ptr{fitsfile}, casesen::Cint, templt::Cstring, colname::Cstring, colnum::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_get_coltype(fptr, colnum, typecode, repeat, width, status)
    @ccall libcfitsio.ffgtcl(fptr::Ptr{fitsfile}, colnum::Cint, typecode::Ptr{Cint}, repeat::Ptr{Clong}, width::Ptr{Clong}, status::Ptr{Status})::Status
end

function fits_get_coltypell(fptr, colnum, typecode, repeat, width, status)
    @ccall libcfitsio.ffgtclll(fptr::Ptr{fitsfile}, colnum::Cint, typecode::Ptr{Cint}, repeat::Ptr{Clonglong}, width::Ptr{Clonglong}, status::Ptr{Status})::Status
end

function fits_get_eqcoltype(fptr, colnum, typecode, repeat, width, status)
    @ccall libcfitsio.ffeqty(fptr::Ptr{fitsfile}, colnum::Cint, typecode::Ptr{Cint}, repeat::Ptr{Clong}, width::Ptr{Clong}, status::Ptr{Status})::Status
end

function fits_get_eqcoltypell(fptr, colnum, typecode, repeat, width, status)
    @ccall libcfitsio.ffeqtyll(fptr::Ptr{fitsfile}, colnum::Cint, typecode::Ptr{Cint}, repeat::Ptr{Clonglong}, width::Ptr{Clonglong}, status::Ptr{Status})::Status
end

function fits_get_num_rows(fptr, nrows, status)
    @ccall libcfitsio.ffgnrw(fptr::Ptr{fitsfile}, nrows::Ptr{Clong}, status::Ptr{Status})::Status
end

function fits_get_num_rowsll(fptr, nrows, status)
    @ccall libcfitsio.ffgnrwll(fptr::Ptr{fitsfile}, nrows::Ptr{Clonglong}, status::Ptr{Status})::Status
end

function fits_get_num_cols(fptr, ncols, status)
    @ccall libcfitsio.ffgncl(fptr::Ptr{fitsfile}, ncols::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_get_acolparms(fptr, colnum, ttype, tbcol, tunit, tform, tscal, tzero, tnull, tdisp, status)
    @ccall libcfitsio.ffgacl(fptr::Ptr{fitsfile}, colnum::Cint, ttype::Cstring, tbcol::Ptr{Clong}, tunit::Cstring, tform::Cstring, tscal::Ptr{Cdouble}, tzero::Ptr{Cdouble}, tnull::Cstring, tdisp::Cstring, status::Ptr{Status})::Status
end

function fits_get_bcolparms(fptr, colnum, ttype, tunit, dtype, repeat, tscal, tzero, tnull, tdisp, status)
    @ccall libcfitsio.ffgbcl(fptr::Ptr{fitsfile}, colnum::Cint, ttype::Cstring, tunit::Cstring, dtype::Cstring, repeat::Ptr{Clong}, tscal::Ptr{Cdouble}, tzero::Ptr{Cdouble}, tnull::Ptr{Clong}, tdisp::Cstring, status::Ptr{Status})::Status
end

function fits_get_bcolparmsll(fptr, colnum, ttype, tunit, dtype, repeat, tscal, tzero, tnull, tdisp, status)
    @ccall libcfitsio.ffgbclll(fptr::Ptr{fitsfile}, colnum::Cint, ttype::Cstring, tunit::Cstring, dtype::Cstring, repeat::Ptr{Clonglong}, tscal::Ptr{Cdouble}, tzero::Ptr{Cdouble}, tnull::Ptr{Clonglong}, tdisp::Cstring, status::Ptr{Status})::Status
end

struct iteratorCol
    fptr::Ptr{fitsfile}
    colnum::Cint
    colname::NTuple{70, Cchar}
    datatype::Cint
    iotype::Cint
    array::Ptr{Cvoid}
    repeat::Clong
    tlmin::Clong
    tlmax::Clong
    tunit::NTuple{70, Cchar}
    tdisp::NTuple{70, Cchar}
end

function fits_iterate_data(ncols, data, offset, nPerLoop, workFn, userPointer, status)
    @ccall libcfitsio.ffiter(ncols::Cint, data::Ptr{iteratorCol}, offset::Clong, nPerLoop::Clong, workFn::Ptr{Cvoid}, userPointer::Ptr{Cvoid}, status::Ptr{Status})::Status
end

function fits_read_grppar_byt(fptr, group, firstelem, nelem, array, status)
    @ccall libcfitsio.ffggpb(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clong, nelem::Clong, array::Ptr{Cuchar}, status::Ptr{Status})::Status
end

function fits_read_grppar_sbyt(fptr, group, firstelem, nelem, array, status)
    @ccall libcfitsio.ffggpsb(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clong, nelem::Clong, array::Ptr{Int8}, status::Ptr{Status})::Status
end

function fits_read_grppar_usht(fptr, group, firstelem, nelem, array, status)
    @ccall libcfitsio.ffggpui(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clong, nelem::Clong, array::Ptr{Cushort}, status::Ptr{Status})::Status
end

function fits_read_grppar_ulng(fptr, group, firstelem, nelem, array, status)
    @ccall libcfitsio.ffggpuj(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clong, nelem::Clong, array::Ptr{Culong}, status::Ptr{Status})::Status
end

function fits_read_grppar_ulnglng(fptr, group, firstelem, nelem, array, status)
    @ccall libcfitsio.ffggpujj(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clong, nelem::Clong, array::Ptr{Culonglong}, status::Ptr{Status})::Status
end

function fits_read_grppar_sht(fptr, group, firstelem, nelem, array, status)
    @ccall libcfitsio.ffggpi(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clong, nelem::Clong, array::Ptr{Cshort}, status::Ptr{Status})::Status
end

function fits_read_grppar_lng(fptr, group, firstelem, nelem, array, status)
    @ccall libcfitsio.ffggpj(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clong, nelem::Clong, array::Ptr{Clong}, status::Ptr{Status})::Status
end

function fits_read_grppar_lnglng(fptr, group, firstelem, nelem, array, status)
    @ccall libcfitsio.ffggpjj(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clong, nelem::Clong, array::Ptr{Clonglong}, status::Ptr{Status})::Status
end

function fits_read_grppar_int(fptr, group, firstelem, nelem, array, status)
    @ccall libcfitsio.ffggpk(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clong, nelem::Clong, array::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_grppar_uint(fptr, group, firstelem, nelem, array, status)
    @ccall libcfitsio.ffggpuk(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clong, nelem::Clong, array::Ptr{Cuint}, status::Ptr{Status})::Status
end

function fits_read_grppar_flt(fptr, group, firstelem, nelem, array, status)
    @ccall libcfitsio.ffggpe(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clong, nelem::Clong, array::Ptr{Cfloat}, status::Ptr{Status})::Status
end

function fits_read_grppar_dbl(fptr, group, firstelem, nelem, array, status)
    @ccall libcfitsio.ffggpd(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clong, nelem::Clong, array::Ptr{Cdouble}, status::Ptr{Status})::Status
end

function fits_read_pix(fptr, datatype, firstpix, nelem, nulval, array, anynul, status)
    @ccall libcfitsio.ffgpxv(fptr::Ptr{fitsfile}, datatype::Cint, firstpix::Ptr{Clong}, nelem::Clonglong, nulval::Ptr{Cvoid}, array::Ptr{Cvoid}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_pixll(fptr, datatype, firstpix, nelem, nulval, array, anynul, status)
    @ccall libcfitsio.ffgpxvll(fptr::Ptr{fitsfile}, datatype::Cint, firstpix::Ptr{Clonglong}, nelem::Clonglong, nulval::Ptr{Cvoid}, array::Ptr{Cvoid}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_pixnull(fptr, datatype, firstpix, nelem, array, nullarray, anynul, status)
    @ccall libcfitsio.ffgpxf(fptr::Ptr{fitsfile}, datatype::Cint, firstpix::Ptr{Clong}, nelem::Clonglong, array::Ptr{Cvoid}, nullarray::Cstring, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_pixnullll(fptr, datatype, firstpix, nelem, array, nullarray, anynul, status)
    @ccall libcfitsio.ffgpxfll(fptr::Ptr{fitsfile}, datatype::Cint, firstpix::Ptr{Clonglong}, nelem::Clonglong, array::Ptr{Cvoid}, nullarray::Cstring, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_img(fptr, datatype, firstelem, nelem, nulval, array, anynul, status)
    @ccall libcfitsio.ffgpv(fptr::Ptr{fitsfile}, datatype::Cint, firstelem::Clonglong, nelem::Clonglong, nulval::Ptr{Cvoid}, array::Ptr{Cvoid}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_imgnull(fptr, datatype, firstelem, nelem, array, nullarray, anynul, status)
    @ccall libcfitsio.ffgpf(fptr::Ptr{fitsfile}, datatype::Cint, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cvoid}, nullarray::Cstring, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_img_byt(fptr, group, firstelem, nelem, nulval, array, anynul, status)
    @ccall libcfitsio.ffgpvb(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, nulval::Cuchar, array::Ptr{Cuchar}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_img_sbyt(fptr, group, firstelem, nelem, nulval, array, anynul, status)
    @ccall libcfitsio.ffgpvsb(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, nulval::Int8, array::Ptr{Int8}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_img_usht(fptr, group, firstelem, nelem, nulval, array, anynul, status)
    @ccall libcfitsio.ffgpvui(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, nulval::Cushort, array::Ptr{Cushort}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_img_ulng(fptr, group, firstelem, nelem, nulval, array, anynul, status)
    @ccall libcfitsio.ffgpvuj(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, nulval::Culong, array::Ptr{Culong}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_img_sht(fptr, group, firstelem, nelem, nulval, array, anynul, status)
    @ccall libcfitsio.ffgpvi(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, nulval::Cshort, array::Ptr{Cshort}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_img_lng(fptr, group, firstelem, nelem, nulval, array, anynul, status)
    @ccall libcfitsio.ffgpvj(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, nulval::Clong, array::Ptr{Clong}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_img_ulnglng(fptr, group, firstelem, nelem, nulval, array, anynul, status)
    @ccall libcfitsio.ffgpvujj(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, nulval::Culonglong, array::Ptr{Culonglong}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_img_lnglng(fptr, group, firstelem, nelem, nulval, array, anynul, status)
    @ccall libcfitsio.ffgpvjj(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, nulval::Clonglong, array::Ptr{Clonglong}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_img_uint(fptr, group, firstelem, nelem, nulval, array, anynul, status)
    @ccall libcfitsio.ffgpvuk(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, nulval::Cuint, array::Ptr{Cuint}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_img_int(fptr, group, firstelem, nelem, nulval, array, anynul, status)
    @ccall libcfitsio.ffgpvk(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, nulval::Cint, array::Ptr{Cint}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_img_flt(fptr, group, firstelem, nelem, nulval, array, anynul, status)
    @ccall libcfitsio.ffgpve(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, nulval::Cfloat, array::Ptr{Cfloat}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_img_dbl(fptr, group, firstelem, nelem, nulval, array, anynul, status)
    @ccall libcfitsio.ffgpvd(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, nulval::Cdouble, array::Ptr{Cdouble}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_imgnull_byt(fptr, group, firstelem, nelem, array, nularray, anynul, status)
    @ccall libcfitsio.ffgpfb(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cuchar}, nularray::Cstring, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_imgnull_sbyt(fptr, group, firstelem, nelem, array, nularray, anynul, status)
    @ccall libcfitsio.ffgpfsb(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Int8}, nularray::Cstring, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_imgnull_usht(fptr, group, firstelem, nelem, array, nularray, anynul, status)
    @ccall libcfitsio.ffgpfui(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cushort}, nularray::Cstring, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_imgnull_ulng(fptr, group, firstelem, nelem, array, nularray, anynul, status)
    @ccall libcfitsio.ffgpfuj(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Culong}, nularray::Cstring, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_imgnull_sht(fptr, group, firstelem, nelem, array, nularray, anynul, status)
    @ccall libcfitsio.ffgpfi(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cshort}, nularray::Cstring, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_imgnull_lng(fptr, group, firstelem, nelem, array, nularray, anynul, status)
    @ccall libcfitsio.ffgpfj(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Clong}, nularray::Cstring, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_imgnull_ulnglng(fptr, group, firstelem, nelem, array, nularray, anynul, status)
    @ccall libcfitsio.ffgpfujj(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Culonglong}, nularray::Cstring, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_imgnull_lnglng(fptr, group, firstelem, nelem, array, nularray, anynul, status)
    @ccall libcfitsio.ffgpfjj(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Clonglong}, nularray::Cstring, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_imgnull_uint(fptr, group, firstelem, nelem, array, nularray, anynul, status)
    @ccall libcfitsio.ffgpfuk(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cuint}, nularray::Cstring, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_imgnull_int(fptr, group, firstelem, nelem, array, nularray, anynul, status)
    @ccall libcfitsio.ffgpfk(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cint}, nularray::Cstring, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_imgnull_flt(fptr, group, firstelem, nelem, array, nularray, anynul, status)
    @ccall libcfitsio.ffgpfe(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cfloat}, nularray::Cstring, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_imgnull_dbl(fptr, group, firstelem, nelem, array, nularray, anynul, status)
    @ccall libcfitsio.ffgpfd(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cdouble}, nularray::Cstring, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_2d_byt(fptr, group, nulval, ncols, naxis1, naxis2, array, anynul, status)
    @ccall libcfitsio.ffg2db(fptr::Ptr{fitsfile}, group::Clong, nulval::Cuchar, ncols::Clonglong, naxis1::Clonglong, naxis2::Clonglong, array::Ptr{Cuchar}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_2d_sbyt(fptr, group, nulval, ncols, naxis1, naxis2, array, anynul, status)
    @ccall libcfitsio.ffg2dsb(fptr::Ptr{fitsfile}, group::Clong, nulval::Int8, ncols::Clonglong, naxis1::Clonglong, naxis2::Clonglong, array::Ptr{Int8}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_2d_usht(fptr, group, nulval, ncols, naxis1, naxis2, array, anynul, status)
    @ccall libcfitsio.ffg2dui(fptr::Ptr{fitsfile}, group::Clong, nulval::Cushort, ncols::Clonglong, naxis1::Clonglong, naxis2::Clonglong, array::Ptr{Cushort}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_2d_ulng(fptr, group, nulval, ncols, naxis1, naxis2, array, anynul, status)
    @ccall libcfitsio.ffg2duj(fptr::Ptr{fitsfile}, group::Clong, nulval::Culong, ncols::Clonglong, naxis1::Clonglong, naxis2::Clonglong, array::Ptr{Culong}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_2d_sht(fptr, group, nulval, ncols, naxis1, naxis2, array, anynul, status)
    @ccall libcfitsio.ffg2di(fptr::Ptr{fitsfile}, group::Clong, nulval::Cshort, ncols::Clonglong, naxis1::Clonglong, naxis2::Clonglong, array::Ptr{Cshort}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_2d_lng(fptr, group, nulval, ncols, naxis1, naxis2, array, anynul, status)
    @ccall libcfitsio.ffg2dj(fptr::Ptr{fitsfile}, group::Clong, nulval::Clong, ncols::Clonglong, naxis1::Clonglong, naxis2::Clonglong, array::Ptr{Clong}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_2d_ulnglng(fptr, group, nulval, ncols, naxis1, naxis2, array, anynul, status)
    @ccall libcfitsio.ffg2dujj(fptr::Ptr{fitsfile}, group::Clong, nulval::Culonglong, ncols::Clonglong, naxis1::Clonglong, naxis2::Clonglong, array::Ptr{Culonglong}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_2d_lnglng(fptr, group, nulval, ncols, naxis1, naxis2, array, anynul, status)
    @ccall libcfitsio.ffg2djj(fptr::Ptr{fitsfile}, group::Clong, nulval::Clonglong, ncols::Clonglong, naxis1::Clonglong, naxis2::Clonglong, array::Ptr{Clonglong}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_2d_uint(fptr, group, nulval, ncols, naxis1, naxis2, array, anynul, status)
    @ccall libcfitsio.ffg2duk(fptr::Ptr{fitsfile}, group::Clong, nulval::Cuint, ncols::Clonglong, naxis1::Clonglong, naxis2::Clonglong, array::Ptr{Cuint}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_2d_int(fptr, group, nulval, ncols, naxis1, naxis2, array, anynul, status)
    @ccall libcfitsio.ffg2dk(fptr::Ptr{fitsfile}, group::Clong, nulval::Cint, ncols::Clonglong, naxis1::Clonglong, naxis2::Clonglong, array::Ptr{Cint}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_2d_flt(fptr, group, nulval, ncols, naxis1, naxis2, array, anynul, status)
    @ccall libcfitsio.ffg2de(fptr::Ptr{fitsfile}, group::Clong, nulval::Cfloat, ncols::Clonglong, naxis1::Clonglong, naxis2::Clonglong, array::Ptr{Cfloat}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_2d_dbl(fptr, group, nulval, ncols, naxis1, naxis2, array, anynul, status)
    @ccall libcfitsio.ffg2dd(fptr::Ptr{fitsfile}, group::Clong, nulval::Cdouble, ncols::Clonglong, naxis1::Clonglong, naxis2::Clonglong, array::Ptr{Cdouble}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_3d_byt(fptr, group, nulval, ncols, nrows, naxis1, naxis2, naxis3, array, anynul, status)
    @ccall libcfitsio.ffg3db(fptr::Ptr{fitsfile}, group::Clong, nulval::Cuchar, ncols::Clonglong, nrows::Clonglong, naxis1::Clonglong, naxis2::Clonglong, naxis3::Clonglong, array::Ptr{Cuchar}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_3d_sbyt(fptr, group, nulval, ncols, nrows, naxis1, naxis2, naxis3, array, anynul, status)
    @ccall libcfitsio.ffg3dsb(fptr::Ptr{fitsfile}, group::Clong, nulval::Int8, ncols::Clonglong, nrows::Clonglong, naxis1::Clonglong, naxis2::Clonglong, naxis3::Clonglong, array::Ptr{Int8}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_3d_usht(fptr, group, nulval, ncols, nrows, naxis1, naxis2, naxis3, array, anynul, status)
    @ccall libcfitsio.ffg3dui(fptr::Ptr{fitsfile}, group::Clong, nulval::Cushort, ncols::Clonglong, nrows::Clonglong, naxis1::Clonglong, naxis2::Clonglong, naxis3::Clonglong, array::Ptr{Cushort}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_3d_ulng(fptr, group, nulval, ncols, nrows, naxis1, naxis2, naxis3, array, anynul, status)
    @ccall libcfitsio.ffg3duj(fptr::Ptr{fitsfile}, group::Clong, nulval::Culong, ncols::Clonglong, nrows::Clonglong, naxis1::Clonglong, naxis2::Clonglong, naxis3::Clonglong, array::Ptr{Culong}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_3d_sht(fptr, group, nulval, ncols, nrows, naxis1, naxis2, naxis3, array, anynul, status)
    @ccall libcfitsio.ffg3di(fptr::Ptr{fitsfile}, group::Clong, nulval::Cshort, ncols::Clonglong, nrows::Clonglong, naxis1::Clonglong, naxis2::Clonglong, naxis3::Clonglong, array::Ptr{Cshort}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_3d_lng(fptr, group, nulval, ncols, nrows, naxis1, naxis2, naxis3, array, anynul, status)
    @ccall libcfitsio.ffg3dj(fptr::Ptr{fitsfile}, group::Clong, nulval::Clong, ncols::Clonglong, nrows::Clonglong, naxis1::Clonglong, naxis2::Clonglong, naxis3::Clonglong, array::Ptr{Clong}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_3d_ulnglng(fptr, group, nulval, ncols, nrows, naxis1, naxis2, naxis3, array, anynul, status)
    @ccall libcfitsio.ffg3dujj(fptr::Ptr{fitsfile}, group::Clong, nulval::Culonglong, ncols::Clonglong, nrows::Clonglong, naxis1::Clonglong, naxis2::Clonglong, naxis3::Clonglong, array::Ptr{Culonglong}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_3d_lnglng(fptr, group, nulval, ncols, nrows, naxis1, naxis2, naxis3, array, anynul, status)
    @ccall libcfitsio.ffg3djj(fptr::Ptr{fitsfile}, group::Clong, nulval::Clonglong, ncols::Clonglong, nrows::Clonglong, naxis1::Clonglong, naxis2::Clonglong, naxis3::Clonglong, array::Ptr{Clonglong}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_3d_uint(fptr, group, nulval, ncols, nrows, naxis1, naxis2, naxis3, array, anynul, status)
    @ccall libcfitsio.ffg3duk(fptr::Ptr{fitsfile}, group::Clong, nulval::Cuint, ncols::Clonglong, nrows::Clonglong, naxis1::Clonglong, naxis2::Clonglong, naxis3::Clonglong, array::Ptr{Cuint}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_3d_int(fptr, group, nulval, ncols, nrows, naxis1, naxis2, naxis3, array, anynul, status)
    @ccall libcfitsio.ffg3dk(fptr::Ptr{fitsfile}, group::Clong, nulval::Cint, ncols::Clonglong, nrows::Clonglong, naxis1::Clonglong, naxis2::Clonglong, naxis3::Clonglong, array::Ptr{Cint}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_3d_flt(fptr, group, nulval, ncols, nrows, naxis1, naxis2, naxis3, array, anynul, status)
    @ccall libcfitsio.ffg3de(fptr::Ptr{fitsfile}, group::Clong, nulval::Cfloat, ncols::Clonglong, nrows::Clonglong, naxis1::Clonglong, naxis2::Clonglong, naxis3::Clonglong, array::Ptr{Cfloat}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_3d_dbl(fptr, group, nulval, ncols, nrows, naxis1, naxis2, naxis3, array, anynul, status)
    @ccall libcfitsio.ffg3dd(fptr::Ptr{fitsfile}, group::Clong, nulval::Cdouble, ncols::Clonglong, nrows::Clonglong, naxis1::Clonglong, naxis2::Clonglong, naxis3::Clonglong, array::Ptr{Cdouble}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_subset(fptr, datatype, blc, trc, inc, nulval, array, anynul, status)
    @ccall libcfitsio.ffgsv(fptr::Ptr{fitsfile}, datatype::Cint, blc::Ptr{Clong}, trc::Ptr{Clong}, inc::Ptr{Clong}, nulval::Ptr{Cvoid}, array::Ptr{Cvoid}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_subset_byt(fptr, colnum, naxis, naxes, blc, trc, inc, nulval, array, anynul, status)
    @ccall libcfitsio.ffgsvb(fptr::Ptr{fitsfile}, colnum::Cint, naxis::Cint, naxes::Ptr{Clong}, blc::Ptr{Clong}, trc::Ptr{Clong}, inc::Ptr{Clong}, nulval::Cuchar, array::Ptr{Cuchar}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_subset_sbyt(fptr, colnum, naxis, naxes, blc, trc, inc, nulval, array, anynul, status)
    @ccall libcfitsio.ffgsvsb(fptr::Ptr{fitsfile}, colnum::Cint, naxis::Cint, naxes::Ptr{Clong}, blc::Ptr{Clong}, trc::Ptr{Clong}, inc::Ptr{Clong}, nulval::Int8, array::Ptr{Int8}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_subset_usht(fptr, colnum, naxis, naxes, blc, trc, inc, nulval, array, anynul, status)
    @ccall libcfitsio.ffgsvui(fptr::Ptr{fitsfile}, colnum::Cint, naxis::Cint, naxes::Ptr{Clong}, blc::Ptr{Clong}, trc::Ptr{Clong}, inc::Ptr{Clong}, nulval::Cushort, array::Ptr{Cushort}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_subset_ulng(fptr, colnum, naxis, naxes, blc, trc, inc, nulval, array, anynul, status)
    @ccall libcfitsio.ffgsvuj(fptr::Ptr{fitsfile}, colnum::Cint, naxis::Cint, naxes::Ptr{Clong}, blc::Ptr{Clong}, trc::Ptr{Clong}, inc::Ptr{Clong}, nulval::Culong, array::Ptr{Culong}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_subset_sht(fptr, colnum, naxis, naxes, blc, trc, inc, nulval, array, anynul, status)
    @ccall libcfitsio.ffgsvi(fptr::Ptr{fitsfile}, colnum::Cint, naxis::Cint, naxes::Ptr{Clong}, blc::Ptr{Clong}, trc::Ptr{Clong}, inc::Ptr{Clong}, nulval::Cshort, array::Ptr{Cshort}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_subset_lng(fptr, colnum, naxis, naxes, blc, trc, inc, nulval, array, anynul, status)
    @ccall libcfitsio.ffgsvj(fptr::Ptr{fitsfile}, colnum::Cint, naxis::Cint, naxes::Ptr{Clong}, blc::Ptr{Clong}, trc::Ptr{Clong}, inc::Ptr{Clong}, nulval::Clong, array::Ptr{Clong}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_subset_ulnglng(fptr, colnum, naxis, naxes, blc, trc, inc, nulval, array, anynul, status)
    @ccall libcfitsio.ffgsvujj(fptr::Ptr{fitsfile}, colnum::Cint, naxis::Cint, naxes::Ptr{Clong}, blc::Ptr{Clong}, trc::Ptr{Clong}, inc::Ptr{Clong}, nulval::Culonglong, array::Ptr{Culonglong}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_subset_lnglng(fptr, colnum, naxis, naxes, blc, trc, inc, nulval, array, anynul, status)
    @ccall libcfitsio.ffgsvjj(fptr::Ptr{fitsfile}, colnum::Cint, naxis::Cint, naxes::Ptr{Clong}, blc::Ptr{Clong}, trc::Ptr{Clong}, inc::Ptr{Clong}, nulval::Clonglong, array::Ptr{Clonglong}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_subset_uint(fptr, colnum, naxis, naxes, blc, trc, inc, nulval, array, anynul, status)
    @ccall libcfitsio.ffgsvuk(fptr::Ptr{fitsfile}, colnum::Cint, naxis::Cint, naxes::Ptr{Clong}, blc::Ptr{Clong}, trc::Ptr{Clong}, inc::Ptr{Clong}, nulval::Cuint, array::Ptr{Cuint}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_subset_int(fptr, colnum, naxis, naxes, blc, trc, inc, nulval, array, anynul, status)
    @ccall libcfitsio.ffgsvk(fptr::Ptr{fitsfile}, colnum::Cint, naxis::Cint, naxes::Ptr{Clong}, blc::Ptr{Clong}, trc::Ptr{Clong}, inc::Ptr{Clong}, nulval::Cint, array::Ptr{Cint}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_subset_flt(fptr, colnum, naxis, naxes, blc, trc, inc, nulval, array, anynul, status)
    @ccall libcfitsio.ffgsve(fptr::Ptr{fitsfile}, colnum::Cint, naxis::Cint, naxes::Ptr{Clong}, blc::Ptr{Clong}, trc::Ptr{Clong}, inc::Ptr{Clong}, nulval::Cfloat, array::Ptr{Cfloat}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_subset_dbl(fptr, colnum, naxis, naxes, blc, trc, inc, nulval, array, anynul, status)
    @ccall libcfitsio.ffgsvd(fptr::Ptr{fitsfile}, colnum::Cint, naxis::Cint, naxes::Ptr{Clong}, blc::Ptr{Clong}, trc::Ptr{Clong}, inc::Ptr{Clong}, nulval::Cdouble, array::Ptr{Cdouble}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_subsetnull_byt(fptr, colnum, naxis, naxes, blc, trc, inc, array, flagval, anynul, status)
    @ccall libcfitsio.ffgsfb(fptr::Ptr{fitsfile}, colnum::Cint, naxis::Cint, naxes::Ptr{Clong}, blc::Ptr{Clong}, trc::Ptr{Clong}, inc::Ptr{Clong}, array::Ptr{Cuchar}, flagval::Cstring, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_subsetnull_sbyt(fptr, colnum, naxis, naxes, blc, trc, inc, array, flagval, anynul, status)
    @ccall libcfitsio.ffgsfsb(fptr::Ptr{fitsfile}, colnum::Cint, naxis::Cint, naxes::Ptr{Clong}, blc::Ptr{Clong}, trc::Ptr{Clong}, inc::Ptr{Clong}, array::Ptr{Int8}, flagval::Cstring, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_subsetnull_usht(fptr, colnum, naxis, naxes, blc, trc, inc, array, flagval, anynul, status)
    @ccall libcfitsio.ffgsfui(fptr::Ptr{fitsfile}, colnum::Cint, naxis::Cint, naxes::Ptr{Clong}, blc::Ptr{Clong}, trc::Ptr{Clong}, inc::Ptr{Clong}, array::Ptr{Cushort}, flagval::Cstring, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_subsetnull_ulng(fptr, colnum, naxis, naxes, blc, trc, inc, array, flagval, anynul, status)
    @ccall libcfitsio.ffgsfuj(fptr::Ptr{fitsfile}, colnum::Cint, naxis::Cint, naxes::Ptr{Clong}, blc::Ptr{Clong}, trc::Ptr{Clong}, inc::Ptr{Clong}, array::Ptr{Culong}, flagval::Cstring, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_subsetnull_sht(fptr, colnum, naxis, naxes, blc, trc, inc, array, flagval, anynul, status)
    @ccall libcfitsio.ffgsfi(fptr::Ptr{fitsfile}, colnum::Cint, naxis::Cint, naxes::Ptr{Clong}, blc::Ptr{Clong}, trc::Ptr{Clong}, inc::Ptr{Clong}, array::Ptr{Cshort}, flagval::Cstring, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_subsetnull_lng(fptr, colnum, naxis, naxes, blc, trc, inc, array, flagval, anynul, status)
    @ccall libcfitsio.ffgsfj(fptr::Ptr{fitsfile}, colnum::Cint, naxis::Cint, naxes::Ptr{Clong}, blc::Ptr{Clong}, trc::Ptr{Clong}, inc::Ptr{Clong}, array::Ptr{Clong}, flagval::Cstring, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_subsetnull_ulnglng(fptr, colnum, naxis, naxes, blc, trc, inc, array, flagval, anynul, status)
    @ccall libcfitsio.ffgsfujj(fptr::Ptr{fitsfile}, colnum::Cint, naxis::Cint, naxes::Ptr{Clong}, blc::Ptr{Clong}, trc::Ptr{Clong}, inc::Ptr{Clong}, array::Ptr{Culonglong}, flagval::Cstring, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_subsetnull_lnglng(fptr, colnum, naxis, naxes, blc, trc, inc, array, flagval, anynul, status)
    @ccall libcfitsio.ffgsfjj(fptr::Ptr{fitsfile}, colnum::Cint, naxis::Cint, naxes::Ptr{Clong}, blc::Ptr{Clong}, trc::Ptr{Clong}, inc::Ptr{Clong}, array::Ptr{Clonglong}, flagval::Cstring, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_subsetnull_uint(fptr, colnum, naxis, naxes, blc, trc, inc, array, flagval, anynul, status)
    @ccall libcfitsio.ffgsfuk(fptr::Ptr{fitsfile}, colnum::Cint, naxis::Cint, naxes::Ptr{Clong}, blc::Ptr{Clong}, trc::Ptr{Clong}, inc::Ptr{Clong}, array::Ptr{Cuint}, flagval::Cstring, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_subsetnull_int(fptr, colnum, naxis, naxes, blc, trc, inc, array, flagval, anynul, status)
    @ccall libcfitsio.ffgsfk(fptr::Ptr{fitsfile}, colnum::Cint, naxis::Cint, naxes::Ptr{Clong}, blc::Ptr{Clong}, trc::Ptr{Clong}, inc::Ptr{Clong}, array::Ptr{Cint}, flagval::Cstring, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_subsetnull_flt(fptr, colnum, naxis, naxes, blc, trc, inc, array, flagval, anynul, status)
    @ccall libcfitsio.ffgsfe(fptr::Ptr{fitsfile}, colnum::Cint, naxis::Cint, naxes::Ptr{Clong}, blc::Ptr{Clong}, trc::Ptr{Clong}, inc::Ptr{Clong}, array::Ptr{Cfloat}, flagval::Cstring, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_subsetnull_dbl(fptr, colnum, naxis, naxes, blc, trc, inc, array, flagval, anynul, status)
    @ccall libcfitsio.ffgsfd(fptr::Ptr{fitsfile}, colnum::Cint, naxis::Cint, naxes::Ptr{Clong}, blc::Ptr{Clong}, trc::Ptr{Clong}, inc::Ptr{Clong}, array::Ptr{Cdouble}, flagval::Cstring, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_copy_image_section(infptr, outfile, imagesection, status)
    @ccall libcfitsio.fits_copy_image_section(infptr::Ptr{fitsfile}, outfile::Ptr{fitsfile}, imagesection::Cstring, status::Ptr{Status})::Status
end

function fits_comp_img(infptr, outfptr, compress_type, tilesize, parm1, parm2, status)
    @ccall libcfitsio.fits_comp_img(infptr::Ptr{fitsfile}, outfptr::Ptr{fitsfile}, compress_type::Cint, tilesize::Ptr{Clong}, parm1::Cint, parm2::Cint, status::Ptr{Status})::Status
end

function fits_decomp_img(infptr, outfptr, status)
    @ccall libcfitsio.fits_decomp_img(infptr::Ptr{fitsfile}, outfptr::Ptr{fitsfile}, status::Ptr{Status})::Status
end

function fits_read_col(fptr, datatype, colnum, firstrow, firstelem, nelem, nulval, array, anynul, status)
    @ccall libcfitsio.ffgcv(fptr::Ptr{fitsfile}, datatype::Cint, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, nulval::Ptr{Cvoid}, array::Ptr{Cvoid}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_cols(fptr, ncols, datatype, colnum, firstrow, nrows, nulval, array, anynul, status)
    @ccall libcfitsio.ffgcvn(fptr::Ptr{fitsfile}, ncols::Cint, datatype::Ptr{Cint}, colnum::Ptr{Cint}, firstrow::Clonglong, nrows::Clonglong, nulval::Ptr{Ptr{Cvoid}}, array::Ptr{Ptr{Cvoid}}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_colnull(fptr, datatype, colnum, firstrow, firstelem, nelem, array, nullarray, anynul, status)
    @ccall libcfitsio.ffgcf(fptr::Ptr{fitsfile}, datatype::Cint, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cvoid}, nullarray::Cstring, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_col_str(fptr, colnum, firstrow, firstelem, nelem, nulval, array, anynul, status)
    @ccall libcfitsio.ffgcvs(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, nulval::Cstring, array::Ptr{Cstring}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_col_log(fptr, colnum, firstrow, firstelem, nelem, nulval, array, anynul, status)
    @ccall libcfitsio.ffgcvl(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, nulval::Cchar, array::Cstring, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_col_byt(fptr, colnum, firstrow, firstelem, nelem, nulval, array, anynul, status)
    @ccall libcfitsio.ffgcvb(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, nulval::Cuchar, array::Ptr{Cuchar}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_col_sbyt(fptr, colnum, firstrow, firstelem, nelem, nulval, array, anynul, status)
    @ccall libcfitsio.ffgcvsb(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, nulval::Int8, array::Ptr{Int8}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_col_usht(fptr, colnum, firstrow, firstelem, nelem, nulval, array, anynul, status)
    @ccall libcfitsio.ffgcvui(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, nulval::Cushort, array::Ptr{Cushort}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_col_ulng(fptr, colnum, firstrow, firstelem, nelem, nulval, array, anynul, status)
    @ccall libcfitsio.ffgcvuj(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, nulval::Culong, array::Ptr{Culong}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_col_sht(fptr, colnum, firstrow, firstelem, nelem, nulval, array, anynul, status)
    @ccall libcfitsio.ffgcvi(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, nulval::Cshort, array::Ptr{Cshort}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_col_lng(fptr, colnum, firstrow, firstelem, nelem, nulval, array, anynul, status)
    @ccall libcfitsio.ffgcvj(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, nulval::Clong, array::Ptr{Clong}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_col_ulnglng(fptr, colnum, firstrow, firstelem, nelem, nulval, array, anynul, status)
    @ccall libcfitsio.ffgcvujj(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, nulval::Culonglong, array::Ptr{Culonglong}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_col_lnglng(fptr, colnum, firstrow, firstelem, nelem, nulval, array, anynul, status)
    @ccall libcfitsio.ffgcvjj(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, nulval::Clonglong, array::Ptr{Clonglong}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_col_uint(fptr, colnum, firstrow, firstelem, nelem, nulval, array, anynul, status)
    @ccall libcfitsio.ffgcvuk(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, nulval::Cuint, array::Ptr{Cuint}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_col_int(fptr, colnum, firstrow, firstelem, nelem, nulval, array, anynul, status)
    @ccall libcfitsio.ffgcvk(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, nulval::Cint, array::Ptr{Cint}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_col_flt(fptr, colnum, firstrow, firstelem, nelem, nulval, array, anynul, status)
    @ccall libcfitsio.ffgcve(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, nulval::Cfloat, array::Ptr{Cfloat}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_col_dbl(fptr, colnum, firstrow, firstelem, nelem, nulval, array, anynul, status)
    @ccall libcfitsio.ffgcvd(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, nulval::Cdouble, array::Ptr{Cdouble}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_col_cmp(fptr, colnum, firstrow, firstelem, nelem, nulval, array, anynul, status)
    @ccall libcfitsio.ffgcvc(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, nulval::Cfloat, array::Ptr{Cfloat}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_col_dblcmp(fptr, colnum, firstrow, firstelem, nelem, nulval, array, anynul, status)
    @ccall libcfitsio.ffgcvm(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, nulval::Cdouble, array::Ptr{Cdouble}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_col_bit(fptr, colnum, firstrow, firstbit, nbits, larray, status)
    @ccall libcfitsio.ffgcx(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstbit::Clonglong, nbits::Clonglong, larray::Cstring, status::Ptr{Status})::Status
end

function fits_read_col_bit_usht(fptr, colnum, firstrow, nrows, firstbit, nbits, array, status)
    @ccall libcfitsio.ffgcxui(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, nrows::Clonglong, firstbit::Clong, nbits::Cint, array::Ptr{Cushort}, status::Ptr{Status})::Status
end

function fits_read_col_bit_uint(fptr, colnum, firstrow, nrows, firstbit, nbits, array, status)
    @ccall libcfitsio.ffgcxuk(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, nrows::Clonglong, firstbit::Clong, nbits::Cint, array::Ptr{Cuint}, status::Ptr{Status})::Status
end

function fits_read_colnull_str(fptr, colnum, firstrow, firstelem, nelem, array, nularray, anynul, status)
    @ccall libcfitsio.ffgcfs(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cstring}, nularray::Cstring, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_colnull_log(fptr, colnum, firstrow, firstelem, nelem, array, nularray, anynul, status)
    @ccall libcfitsio.ffgcfl(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Cstring, nularray::Cstring, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_colnull_byt(fptr, colnum, firstrow, firstelem, nelem, array, nularray, anynul, status)
    @ccall libcfitsio.ffgcfb(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cuchar}, nularray::Cstring, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_colnull_sbyt(fptr, colnum, firstrow, firstelem, nelem, array, nularray, anynul, status)
    @ccall libcfitsio.ffgcfsb(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Int8}, nularray::Cstring, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_colnull_usht(fptr, colnum, firstrow, firstelem, nelem, array, nularray, anynul, status)
    @ccall libcfitsio.ffgcfui(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cushort}, nularray::Cstring, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_colnull_ulng(fptr, colnum, firstrow, firstelem, nelem, array, nularray, anynul, status)
    @ccall libcfitsio.ffgcfuj(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Culong}, nularray::Cstring, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_colnull_sht(fptr, colnum, firstrow, firstelem, nelem, array, nularray, anynul, status)
    @ccall libcfitsio.ffgcfi(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cshort}, nularray::Cstring, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_colnull_lng(fptr, colnum, firstrow, firstelem, nelem, array, nularray, anynul, status)
    @ccall libcfitsio.ffgcfj(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Clong}, nularray::Cstring, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_colnull_ulnglng(fptr, colnum, firstrow, firstelem, nelem, array, nularray, anynul, status)
    @ccall libcfitsio.ffgcfujj(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Culonglong}, nularray::Cstring, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_colnull_lnglng(fptr, colnum, firstrow, firstelem, nelem, array, nularray, anynul, status)
    @ccall libcfitsio.ffgcfjj(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Clonglong}, nularray::Cstring, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_colnull_uint(fptr, colnum, firstrow, firstelem, nelem, array, nularray, anynul, status)
    @ccall libcfitsio.ffgcfuk(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cuint}, nularray::Cstring, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_colnull_int(fptr, colnum, firstrow, firstelem, nelem, array, nularray, anynul, status)
    @ccall libcfitsio.ffgcfk(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cint}, nularray::Cstring, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_colnull_flt(fptr, colnum, firstrow, firstelem, nelem, array, nularray, anynul, status)
    @ccall libcfitsio.ffgcfe(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cfloat}, nularray::Cstring, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_colnull_dbl(fptr, colnum, firstrow, firstelem, nelem, array, nularray, anynul, status)
    @ccall libcfitsio.ffgcfd(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cdouble}, nularray::Cstring, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_colnull_cmp(fptr, colnum, firstrow, firstelem, nelem, array, nularray, anynul, status)
    @ccall libcfitsio.ffgcfc(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cfloat}, nularray::Cstring, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_colnull_dblcmp(fptr, colnum, firstrow, firstelem, nelem, array, nularray, anynul, status)
    @ccall libcfitsio.ffgcfm(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cdouble}, nularray::Cstring, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_read_descript(fptr, colnum, rownum, length, heapaddr, status)
    @ccall libcfitsio.ffgdes(fptr::Ptr{fitsfile}, colnum::Cint, rownum::Clonglong, length::Ptr{Clong}, heapaddr::Ptr{Clong}, status::Ptr{Status})::Status
end

function fits_read_descriptll(fptr, colnum, rownum, length, heapaddr, status)
    @ccall libcfitsio.ffgdesll(fptr::Ptr{fitsfile}, colnum::Cint, rownum::Clonglong, length::Ptr{Clonglong}, heapaddr::Ptr{Clonglong}, status::Ptr{Status})::Status
end

function fits_read_descripts(fptr, colnum, firstrow, nrows, length, heapaddr, status)
    @ccall libcfitsio.ffgdess(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, nrows::Clonglong, length::Ptr{Clong}, heapaddr::Ptr{Clong}, status::Ptr{Status})::Status
end

function fits_read_descriptsll(fptr, colnum, firstrow, nrows, length, heapaddr, status)
    @ccall libcfitsio.ffgdessll(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, nrows::Clonglong, length::Ptr{Clonglong}, heapaddr::Ptr{Clonglong}, status::Ptr{Status})::Status
end

function fits_read_tblbytes(fptr, firstrow, firstchar, nchars, values, status)
    @ccall libcfitsio.ffgtbb(fptr::Ptr{fitsfile}, firstrow::Clonglong, firstchar::Clonglong, nchars::Clonglong, values::Ptr{Cuchar}, status::Ptr{Status})::Status
end

function fits_write_grppar_byt(fptr, group, firstelem, nelem, array, status)
    @ccall libcfitsio.ffpgpb(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clong, nelem::Clong, array::Ptr{Cuchar}, status::Ptr{Status})::Status
end

function fits_write_grppar_sbyt(fptr, group, firstelem, nelem, array, status)
    @ccall libcfitsio.ffpgpsb(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clong, nelem::Clong, array::Ptr{Int8}, status::Ptr{Status})::Status
end

function fits_write_grppar_usht(fptr, group, firstelem, nelem, array, status)
    @ccall libcfitsio.ffpgpui(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clong, nelem::Clong, array::Ptr{Cushort}, status::Ptr{Status})::Status
end

function fits_write_grppar_ulng(fptr, group, firstelem, nelem, array, status)
    @ccall libcfitsio.ffpgpuj(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clong, nelem::Clong, array::Ptr{Culong}, status::Ptr{Status})::Status
end

function fits_write_grppar_sht(fptr, group, firstelem, nelem, array, status)
    @ccall libcfitsio.ffpgpi(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clong, nelem::Clong, array::Ptr{Cshort}, status::Ptr{Status})::Status
end

function fits_write_grppar_lng(fptr, group, firstelem, nelem, array, status)
    @ccall libcfitsio.ffpgpj(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clong, nelem::Clong, array::Ptr{Clong}, status::Ptr{Status})::Status
end

function fits_write_grppar_ulnglng(fptr, group, firstelem, nelem, array, status)
    @ccall libcfitsio.ffpgpujj(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clong, nelem::Clong, array::Ptr{Culonglong}, status::Ptr{Status})::Status
end

function fits_write_grppar_lnglng(fptr, group, firstelem, nelem, array, status)
    @ccall libcfitsio.ffpgpjj(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clong, nelem::Clong, array::Ptr{Clonglong}, status::Ptr{Status})::Status
end

function fits_write_grppar_uint(fptr, group, firstelem, nelem, array, status)
    @ccall libcfitsio.ffpgpuk(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clong, nelem::Clong, array::Ptr{Cuint}, status::Ptr{Status})::Status
end

function fits_write_grppar_int(fptr, group, firstelem, nelem, array, status)
    @ccall libcfitsio.ffpgpk(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clong, nelem::Clong, array::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_write_grppar_flt(fptr, group, firstelem, nelem, array, status)
    @ccall libcfitsio.ffpgpe(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clong, nelem::Clong, array::Ptr{Cfloat}, status::Ptr{Status})::Status
end

function fits_write_grppar_dbl(fptr, group, firstelem, nelem, array, status)
    @ccall libcfitsio.ffpgpd(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clong, nelem::Clong, array::Ptr{Cdouble}, status::Ptr{Status})::Status
end

function fits_write_pix(fptr, datatype, firstpix, nelem, array, status)
    @ccall libcfitsio.ffppx(fptr::Ptr{fitsfile}, datatype::Cint, firstpix::Ptr{Clong}, nelem::Clonglong, array::Ptr{Cvoid}, status::Ptr{Status})::Status
end

function fits_write_pixll(fptr, datatype, firstpix, nelem, array, status)
    @ccall libcfitsio.ffppxll(fptr::Ptr{fitsfile}, datatype::Cint, firstpix::Ptr{Clonglong}, nelem::Clonglong, array::Ptr{Cvoid}, status::Ptr{Status})::Status
end

function fits_write_pixnull(fptr, datatype, firstpix, nelem, array, nulval, status)
    @ccall libcfitsio.ffppxn(fptr::Ptr{fitsfile}, datatype::Cint, firstpix::Ptr{Clong}, nelem::Clonglong, array::Ptr{Cvoid}, nulval::Ptr{Cvoid}, status::Ptr{Status})::Status
end

function fits_write_pixnullll(fptr, datatype, firstpix, nelem, array, nulval, status)
    @ccall libcfitsio.ffppxnll(fptr::Ptr{fitsfile}, datatype::Cint, firstpix::Ptr{Clonglong}, nelem::Clonglong, array::Ptr{Cvoid}, nulval::Ptr{Cvoid}, status::Ptr{Status})::Status
end

function fits_write_img(fptr, datatype, firstelem, nelem, array, status)
    @ccall libcfitsio.ffppr(fptr::Ptr{fitsfile}, datatype::Cint, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cvoid}, status::Ptr{Status})::Status
end

function fits_write_img_byt(fptr, group, firstelem, nelem, array, status)
    @ccall libcfitsio.ffpprb(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cuchar}, status::Ptr{Status})::Status
end

function fits_write_img_sbyt(fptr, group, firstelem, nelem, array, status)
    @ccall libcfitsio.ffpprsb(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Int8}, status::Ptr{Status})::Status
end

function fits_write_img_usht(fptr, group, firstelem, nelem, array, status)
    @ccall libcfitsio.ffpprui(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cushort}, status::Ptr{Status})::Status
end

function fits_write_img_ulng(fptr, group, firstelem, nelem, array, status)
    @ccall libcfitsio.ffppruj(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Culong}, status::Ptr{Status})::Status
end

function fits_write_img_sht(fptr, group, firstelem, nelem, array, status)
    @ccall libcfitsio.ffppri(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cshort}, status::Ptr{Status})::Status
end

function fits_write_img_lng(fptr, group, firstelem, nelem, array, status)
    @ccall libcfitsio.ffpprj(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Clong}, status::Ptr{Status})::Status
end

function fits_write_img_ulnglng(fptr, group, firstelem, nelem, array, status)
    @ccall libcfitsio.ffpprujj(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Culonglong}, status::Ptr{Status})::Status
end

function fits_write_img_lnglng(fptr, group, firstelem, nelem, array, status)
    @ccall libcfitsio.ffpprjj(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Clonglong}, status::Ptr{Status})::Status
end

function fits_write_img_uint(fptr, group, firstelem, nelem, array, status)
    @ccall libcfitsio.ffppruk(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cuint}, status::Ptr{Status})::Status
end

function fits_write_img_int(fptr, group, firstelem, nelem, array, status)
    @ccall libcfitsio.ffpprk(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_write_img_flt(fptr, group, firstelem, nelem, array, status)
    @ccall libcfitsio.ffppre(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cfloat}, status::Ptr{Status})::Status
end

function fits_write_img_dbl(fptr, group, firstelem, nelem, array, status)
    @ccall libcfitsio.ffpprd(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cdouble}, status::Ptr{Status})::Status
end

function fits_write_imgnull(fptr, datatype, firstelem, nelem, array, nulval, status)
    @ccall libcfitsio.ffppn(fptr::Ptr{fitsfile}, datatype::Cint, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cvoid}, nulval::Ptr{Cvoid}, status::Ptr{Status})::Status
end

function fits_write_imgnull_byt(fptr, group, firstelem, nelem, array, nulval, status)
    @ccall libcfitsio.ffppnb(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cuchar}, nulval::Cuchar, status::Ptr{Status})::Status
end

function fits_write_imgnull_sbyt(fptr, group, firstelem, nelem, array, nulval, status)
    @ccall libcfitsio.ffppnsb(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Int8}, nulval::Int8, status::Ptr{Status})::Status
end

function fits_write_imgnull_usht(fptr, group, firstelem, nelem, array, nulval, status)
    @ccall libcfitsio.ffppnui(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cushort}, nulval::Cushort, status::Ptr{Status})::Status
end

function fits_write_imgnull_ulng(fptr, group, firstelem, nelem, array, nulval, status)
    @ccall libcfitsio.ffppnuj(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Culong}, nulval::Culong, status::Ptr{Status})::Status
end

function fits_write_imgnull_sht(fptr, group, firstelem, nelem, array, nulval, status)
    @ccall libcfitsio.ffppni(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cshort}, nulval::Cshort, status::Ptr{Status})::Status
end

function fits_write_imgnull_lng(fptr, group, firstelem, nelem, array, nulval, status)
    @ccall libcfitsio.ffppnj(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Clong}, nulval::Clong, status::Ptr{Status})::Status
end

function fits_write_imgnull_ulnglng(fptr, group, firstelem, nelem, array, nulval, status)
    @ccall libcfitsio.ffppnujj(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Culonglong}, nulval::Culonglong, status::Ptr{Status})::Status
end

function fits_write_imgnull_lnglng(fptr, group, firstelem, nelem, array, nulval, status)
    @ccall libcfitsio.ffppnjj(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Clonglong}, nulval::Clonglong, status::Ptr{Status})::Status
end

function fits_write_imgnull_uint(fptr, group, firstelem, nelem, array, nulval, status)
    @ccall libcfitsio.ffppnuk(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cuint}, nulval::Cuint, status::Ptr{Status})::Status
end

function fits_write_imgnull_int(fptr, group, firstelem, nelem, array, nulval, status)
    @ccall libcfitsio.ffppnk(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cint}, nulval::Cint, status::Ptr{Status})::Status
end

function fits_write_imgnull_flt(fptr, group, firstelem, nelem, array, nulval, status)
    @ccall libcfitsio.ffppne(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cfloat}, nulval::Cfloat, status::Ptr{Status})::Status
end

function fits_write_imgnull_dbl(fptr, group, firstelem, nelem, array, nulval, status)
    @ccall libcfitsio.ffppnd(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cdouble}, nulval::Cdouble, status::Ptr{Status})::Status
end

function fits_write_img_null(fptr, group, firstelem, nelem, status)
    @ccall libcfitsio.ffppru(fptr::Ptr{fitsfile}, group::Clong, firstelem::Clonglong, nelem::Clonglong, status::Ptr{Status})::Status
end

function fits_write_null_img(fptr, firstelem, nelem, status)
    @ccall libcfitsio.ffpprn(fptr::Ptr{fitsfile}, firstelem::Clonglong, nelem::Clonglong, status::Ptr{Status})::Status
end

function fits_write_2d_byt(fptr, group, ncols, naxis1, naxis2, array, status)
    @ccall libcfitsio.ffp2db(fptr::Ptr{fitsfile}, group::Clong, ncols::Clonglong, naxis1::Clonglong, naxis2::Clonglong, array::Ptr{Cuchar}, status::Ptr{Status})::Status
end

function fits_write_2d_sbyt(fptr, group, ncols, naxis1, naxis2, array, status)
    @ccall libcfitsio.ffp2dsb(fptr::Ptr{fitsfile}, group::Clong, ncols::Clonglong, naxis1::Clonglong, naxis2::Clonglong, array::Ptr{Int8}, status::Ptr{Status})::Status
end

function fits_write_2d_usht(fptr, group, ncols, naxis1, naxis2, array, status)
    @ccall libcfitsio.ffp2dui(fptr::Ptr{fitsfile}, group::Clong, ncols::Clonglong, naxis1::Clonglong, naxis2::Clonglong, array::Ptr{Cushort}, status::Ptr{Status})::Status
end

function fits_write_2d_ulng(fptr, group, ncols, naxis1, naxis2, array, status)
    @ccall libcfitsio.ffp2duj(fptr::Ptr{fitsfile}, group::Clong, ncols::Clonglong, naxis1::Clonglong, naxis2::Clonglong, array::Ptr{Culong}, status::Ptr{Status})::Status
end

function fits_write_2d_sht(fptr, group, ncols, naxis1, naxis2, array, status)
    @ccall libcfitsio.ffp2di(fptr::Ptr{fitsfile}, group::Clong, ncols::Clonglong, naxis1::Clonglong, naxis2::Clonglong, array::Ptr{Cshort}, status::Ptr{Status})::Status
end

function fits_write_2d_lng(fptr, group, ncols, naxis1, naxis2, array, status)
    @ccall libcfitsio.ffp2dj(fptr::Ptr{fitsfile}, group::Clong, ncols::Clonglong, naxis1::Clonglong, naxis2::Clonglong, array::Ptr{Clong}, status::Ptr{Status})::Status
end

function fits_write_2d_ulnglng(fptr, group, ncols, naxis1, naxis2, array, status)
    @ccall libcfitsio.ffp2dujj(fptr::Ptr{fitsfile}, group::Clong, ncols::Clonglong, naxis1::Clonglong, naxis2::Clonglong, array::Ptr{Culonglong}, status::Ptr{Status})::Status
end

function fits_write_2d_lnglng(fptr, group, ncols, naxis1, naxis2, array, status)
    @ccall libcfitsio.ffp2djj(fptr::Ptr{fitsfile}, group::Clong, ncols::Clonglong, naxis1::Clonglong, naxis2::Clonglong, array::Ptr{Clonglong}, status::Ptr{Status})::Status
end

function fits_write_2d_uint(fptr, group, ncols, naxis1, naxis2, array, status)
    @ccall libcfitsio.ffp2duk(fptr::Ptr{fitsfile}, group::Clong, ncols::Clonglong, naxis1::Clonglong, naxis2::Clonglong, array::Ptr{Cuint}, status::Ptr{Status})::Status
end

function fits_write_2d_int(fptr, group, ncols, naxis1, naxis2, array, status)
    @ccall libcfitsio.ffp2dk(fptr::Ptr{fitsfile}, group::Clong, ncols::Clonglong, naxis1::Clonglong, naxis2::Clonglong, array::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_write_2d_flt(fptr, group, ncols, naxis1, naxis2, array, status)
    @ccall libcfitsio.ffp2de(fptr::Ptr{fitsfile}, group::Clong, ncols::Clonglong, naxis1::Clonglong, naxis2::Clonglong, array::Ptr{Cfloat}, status::Ptr{Status})::Status
end

function fits_write_2d_dbl(fptr, group, ncols, naxis1, naxis2, array, status)
    @ccall libcfitsio.ffp2dd(fptr::Ptr{fitsfile}, group::Clong, ncols::Clonglong, naxis1::Clonglong, naxis2::Clonglong, array::Ptr{Cdouble}, status::Ptr{Status})::Status
end

function fits_write_3d_byt(fptr, group, ncols, nrows, naxis1, naxis2, naxis3, array, status)
    @ccall libcfitsio.ffp3db(fptr::Ptr{fitsfile}, group::Clong, ncols::Clonglong, nrows::Clonglong, naxis1::Clonglong, naxis2::Clonglong, naxis3::Clonglong, array::Ptr{Cuchar}, status::Ptr{Status})::Status
end

function fits_write_3d_sbyt(fptr, group, ncols, nrows, naxis1, naxis2, naxis3, array, status)
    @ccall libcfitsio.ffp3dsb(fptr::Ptr{fitsfile}, group::Clong, ncols::Clonglong, nrows::Clonglong, naxis1::Clonglong, naxis2::Clonglong, naxis3::Clonglong, array::Ptr{Int8}, status::Ptr{Status})::Status
end

function fits_write_3d_usht(fptr, group, ncols, nrows, naxis1, naxis2, naxis3, array, status)
    @ccall libcfitsio.ffp3dui(fptr::Ptr{fitsfile}, group::Clong, ncols::Clonglong, nrows::Clonglong, naxis1::Clonglong, naxis2::Clonglong, naxis3::Clonglong, array::Ptr{Cushort}, status::Ptr{Status})::Status
end

function fits_write_3d_ulng(fptr, group, ncols, nrows, naxis1, naxis2, naxis3, array, status)
    @ccall libcfitsio.ffp3duj(fptr::Ptr{fitsfile}, group::Clong, ncols::Clonglong, nrows::Clonglong, naxis1::Clonglong, naxis2::Clonglong, naxis3::Clonglong, array::Ptr{Culong}, status::Ptr{Status})::Status
end

function fits_write_3d_sht(fptr, group, ncols, nrows, naxis1, naxis2, naxis3, array, status)
    @ccall libcfitsio.ffp3di(fptr::Ptr{fitsfile}, group::Clong, ncols::Clonglong, nrows::Clonglong, naxis1::Clonglong, naxis2::Clonglong, naxis3::Clonglong, array::Ptr{Cshort}, status::Ptr{Status})::Status
end

function fits_write_3d_lng(fptr, group, ncols, nrows, naxis1, naxis2, naxis3, array, status)
    @ccall libcfitsio.ffp3dj(fptr::Ptr{fitsfile}, group::Clong, ncols::Clonglong, nrows::Clonglong, naxis1::Clonglong, naxis2::Clonglong, naxis3::Clonglong, array::Ptr{Clong}, status::Ptr{Status})::Status
end

function fits_write_3d_ulnglng(fptr, group, ncols, nrows, naxis1, naxis2, naxis3, array, status)
    @ccall libcfitsio.ffp3dujj(fptr::Ptr{fitsfile}, group::Clong, ncols::Clonglong, nrows::Clonglong, naxis1::Clonglong, naxis2::Clonglong, naxis3::Clonglong, array::Ptr{Culonglong}, status::Ptr{Status})::Status
end

function fits_write_3d_lnglng(fptr, group, ncols, nrows, naxis1, naxis2, naxis3, array, status)
    @ccall libcfitsio.ffp3djj(fptr::Ptr{fitsfile}, group::Clong, ncols::Clonglong, nrows::Clonglong, naxis1::Clonglong, naxis2::Clonglong, naxis3::Clonglong, array::Ptr{Clonglong}, status::Ptr{Status})::Status
end

function fits_write_3d_uint(fptr, group, ncols, nrows, naxis1, naxis2, naxis3, array, status)
    @ccall libcfitsio.ffp3duk(fptr::Ptr{fitsfile}, group::Clong, ncols::Clonglong, nrows::Clonglong, naxis1::Clonglong, naxis2::Clonglong, naxis3::Clonglong, array::Ptr{Cuint}, status::Ptr{Status})::Status
end

function fits_write_3d_int(fptr, group, ncols, nrows, naxis1, naxis2, naxis3, array, status)
    @ccall libcfitsio.ffp3dk(fptr::Ptr{fitsfile}, group::Clong, ncols::Clonglong, nrows::Clonglong, naxis1::Clonglong, naxis2::Clonglong, naxis3::Clonglong, array::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_write_3d_flt(fptr, group, ncols, nrows, naxis1, naxis2, naxis3, array, status)
    @ccall libcfitsio.ffp3de(fptr::Ptr{fitsfile}, group::Clong, ncols::Clonglong, nrows::Clonglong, naxis1::Clonglong, naxis2::Clonglong, naxis3::Clonglong, array::Ptr{Cfloat}, status::Ptr{Status})::Status
end

function fits_write_3d_dbl(fptr, group, ncols, nrows, naxis1, naxis2, naxis3, array, status)
    @ccall libcfitsio.ffp3dd(fptr::Ptr{fitsfile}, group::Clong, ncols::Clonglong, nrows::Clonglong, naxis1::Clonglong, naxis2::Clonglong, naxis3::Clonglong, array::Ptr{Cdouble}, status::Ptr{Status})::Status
end

function fits_write_subset(fptr, datatype, fpixel, lpixel, array, status)
    @ccall libcfitsio.ffpss(fptr::Ptr{fitsfile}, datatype::Cint, fpixel::Ptr{Clong}, lpixel::Ptr{Clong}, array::Ptr{Cvoid}, status::Ptr{Status})::Status
end

function fits_write_subset_byt(fptr, group, naxis, naxes, fpixel, lpixel, array, status)
    @ccall libcfitsio.ffpssb(fptr::Ptr{fitsfile}, group::Clong, naxis::Clong, naxes::Ptr{Clong}, fpixel::Ptr{Clong}, lpixel::Ptr{Clong}, array::Ptr{Cuchar}, status::Ptr{Status})::Status
end

function fits_write_subset_sbyt(fptr, group, naxis, naxes, fpixel, lpixel, array, status)
    @ccall libcfitsio.ffpsssb(fptr::Ptr{fitsfile}, group::Clong, naxis::Clong, naxes::Ptr{Clong}, fpixel::Ptr{Clong}, lpixel::Ptr{Clong}, array::Ptr{Int8}, status::Ptr{Status})::Status
end

function fits_write_subset_usht(fptr, group, naxis, naxes, fpixel, lpixel, array, status)
    @ccall libcfitsio.ffpssui(fptr::Ptr{fitsfile}, group::Clong, naxis::Clong, naxes::Ptr{Clong}, fpixel::Ptr{Clong}, lpixel::Ptr{Clong}, array::Ptr{Cushort}, status::Ptr{Status})::Status
end

function fits_write_subset_ulng(fptr, group, naxis, naxes, fpixel, lpixel, array, status)
    @ccall libcfitsio.ffpssuj(fptr::Ptr{fitsfile}, group::Clong, naxis::Clong, naxes::Ptr{Clong}, fpixel::Ptr{Clong}, lpixel::Ptr{Clong}, array::Ptr{Culong}, status::Ptr{Status})::Status
end

function fits_write_subset_sht(fptr, group, naxis, naxes, fpixel, lpixel, array, status)
    @ccall libcfitsio.ffpssi(fptr::Ptr{fitsfile}, group::Clong, naxis::Clong, naxes::Ptr{Clong}, fpixel::Ptr{Clong}, lpixel::Ptr{Clong}, array::Ptr{Cshort}, status::Ptr{Status})::Status
end

function fits_write_subset_lng(fptr, group, naxis, naxes, fpixel, lpixel, array, status)
    @ccall libcfitsio.ffpssj(fptr::Ptr{fitsfile}, group::Clong, naxis::Clong, naxes::Ptr{Clong}, fpixel::Ptr{Clong}, lpixel::Ptr{Clong}, array::Ptr{Clong}, status::Ptr{Status})::Status
end

function fits_write_subset_ulnglng(fptr, group, naxis, naxes, fpixel, lpixel, array, status)
    @ccall libcfitsio.ffpssujj(fptr::Ptr{fitsfile}, group::Clong, naxis::Clong, naxes::Ptr{Clong}, fpixel::Ptr{Clong}, lpixel::Ptr{Clong}, array::Ptr{Culonglong}, status::Ptr{Status})::Status
end

function fits_write_subset_lnglng(fptr, group, naxis, naxes, fpixel, lpixel, array, status)
    @ccall libcfitsio.ffpssjj(fptr::Ptr{fitsfile}, group::Clong, naxis::Clong, naxes::Ptr{Clong}, fpixel::Ptr{Clong}, lpixel::Ptr{Clong}, array::Ptr{Clonglong}, status::Ptr{Status})::Status
end

function fits_write_subset_uint(fptr, group, naxis, naxes, fpixel, lpixel, array, status)
    @ccall libcfitsio.ffpssuk(fptr::Ptr{fitsfile}, group::Clong, naxis::Clong, naxes::Ptr{Clong}, fpixel::Ptr{Clong}, lpixel::Ptr{Clong}, array::Ptr{Cuint}, status::Ptr{Status})::Status
end

function fits_write_subset_int(fptr, group, naxis, naxes, fpixel, lpixel, array, status)
    @ccall libcfitsio.ffpssk(fptr::Ptr{fitsfile}, group::Clong, naxis::Clong, naxes::Ptr{Clong}, fpixel::Ptr{Clong}, lpixel::Ptr{Clong}, array::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_write_subset_flt(fptr, group, naxis, naxes, fpixel, lpixel, array, status)
    @ccall libcfitsio.ffpsse(fptr::Ptr{fitsfile}, group::Clong, naxis::Clong, naxes::Ptr{Clong}, fpixel::Ptr{Clong}, lpixel::Ptr{Clong}, array::Ptr{Cfloat}, status::Ptr{Status})::Status
end

function fits_write_subset_dbl(fptr, group, naxis, naxes, fpixel, lpixel, array, status)
    @ccall libcfitsio.ffpssd(fptr::Ptr{fitsfile}, group::Clong, naxis::Clong, naxes::Ptr{Clong}, fpixel::Ptr{Clong}, lpixel::Ptr{Clong}, array::Ptr{Cdouble}, status::Ptr{Status})::Status
end

function fits_write_col(fptr, datatype, colnum, firstrow, firstelem, nelem, array, status)
    @ccall libcfitsio.ffpcl(fptr::Ptr{fitsfile}, datatype::Cint, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cvoid}, status::Ptr{Status})::Status
end

function fits_write_cols(fptr, ncols, datatype, colnum, firstrow, nrows, array, nulval, status)
    @ccall libcfitsio.ffpcln(fptr::Ptr{fitsfile}, ncols::Cint, datatype::Ptr{Cint}, colnum::Ptr{Cint}, firstrow::Clonglong, nrows::Clonglong, array::Ptr{Ptr{Cvoid}}, nulval::Ptr{Ptr{Cvoid}}, status::Ptr{Status})::Status
end

function fits_write_col_str(fptr, colnum, firstrow, firstelem, nelem, array, status)
    @ccall libcfitsio.ffpcls(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cstring}, status::Ptr{Status})::Status
end

function fits_write_col_log(fptr, colnum, firstrow, firstelem, nelem, array, status)
    @ccall libcfitsio.ffpcll(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Cstring, status::Ptr{Status})::Status
end

function fits_write_col_byt(fptr, colnum, firstrow, firstelem, nelem, array, status)
    @ccall libcfitsio.ffpclb(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cuchar}, status::Ptr{Status})::Status
end

function fits_write_col_sbyt(fptr, colnum, firstrow, firstelem, nelem, array, status)
    @ccall libcfitsio.ffpclsb(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Int8}, status::Ptr{Status})::Status
end

function fits_write_col_usht(fptr, colnum, firstrow, firstelem, nelem, array, status)
    @ccall libcfitsio.ffpclui(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cushort}, status::Ptr{Status})::Status
end

function fits_write_col_ulng(fptr, colnum, firstrow, firstelem, nelem, array, status)
    @ccall libcfitsio.ffpcluj(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Culong}, status::Ptr{Status})::Status
end

function fits_write_col_sht(fptr, colnum, firstrow, firstelem, nelem, array, status)
    @ccall libcfitsio.ffpcli(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cshort}, status::Ptr{Status})::Status
end

function fits_write_col_lng(fptr, colnum, firstrow, firstelem, nelem, array, status)
    @ccall libcfitsio.ffpclj(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Clong}, status::Ptr{Status})::Status
end

function fits_write_col_ulnglng(fptr, colnum, firstrow, firstelem, nelem, array, status)
    @ccall libcfitsio.ffpclujj(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Culonglong}, status::Ptr{Status})::Status
end

function fits_write_col_lnglng(fptr, colnum, firstrow, firstelem, nelem, array, status)
    @ccall libcfitsio.ffpcljj(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Clonglong}, status::Ptr{Status})::Status
end

function fits_write_col_uint(fptr, colnum, firstrow, firstelem, nelem, array, status)
    @ccall libcfitsio.ffpcluk(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cuint}, status::Ptr{Status})::Status
end

function fits_write_col_int(fptr, colnum, firstrow, firstelem, nelem, array, status)
    @ccall libcfitsio.ffpclk(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_write_col_flt(fptr, colnum, firstrow, firstelem, nelem, array, status)
    @ccall libcfitsio.ffpcle(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cfloat}, status::Ptr{Status})::Status
end

function fits_write_col_dbl(fptr, colnum, firstrow, firstelem, nelem, array, status)
    @ccall libcfitsio.ffpcld(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cdouble}, status::Ptr{Status})::Status
end

function fits_write_col_cmp(fptr, colnum, firstrow, firstelem, nelem, array, status)
    @ccall libcfitsio.ffpclc(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cfloat}, status::Ptr{Status})::Status
end

function fits_write_col_dblcmp(fptr, colnum, firstrow, firstelem, nelem, array, status)
    @ccall libcfitsio.ffpclm(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cdouble}, status::Ptr{Status})::Status
end

function fits_write_col_null(fptr, colnum, firstrow, firstelem, nelem, status)
    @ccall libcfitsio.ffpclu(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, status::Ptr{Status})::Status
end

function fits_write_col_bit(fptr, colnum, frow, fbit, nbit, larray, status)
    @ccall libcfitsio.ffpclx(fptr::Ptr{fitsfile}, colnum::Cint, frow::Clonglong, fbit::Clong, nbit::Clong, larray::Cstring, status::Ptr{Status})::Status
end

function fits_write_nullrows(fptr, firstrow, nrows, status)
    @ccall libcfitsio.ffprwu(fptr::Ptr{fitsfile}, firstrow::Clonglong, nrows::Clonglong, status::Ptr{Status})::Status
end

function fits_write_colnull(fptr, datatype, colnum, firstrow, firstelem, nelem, array, nulval, status)
    @ccall libcfitsio.ffpcn(fptr::Ptr{fitsfile}, datatype::Cint, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cvoid}, nulval::Ptr{Cvoid}, status::Ptr{Status})::Status
end

function fits_write_colnull_str(fptr, colnum, firstrow, firstelem, nelem, array, nulvalue, status)
    @ccall libcfitsio.ffpcns(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cstring}, nulvalue::Cstring, status::Ptr{Status})::Status
end

function fits_write_colnull_log(fptr, colnum, firstrow, firstelem, nelem, array, nulvalue, status)
    @ccall libcfitsio.ffpcnl(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Cstring, nulvalue::Cchar, status::Ptr{Status})::Status
end

function fits_write_colnull_byt(fptr, colnum, firstrow, firstelem, nelem, array, nulvalue, status)
    @ccall libcfitsio.ffpcnb(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cuchar}, nulvalue::Cuchar, status::Ptr{Status})::Status
end

function fits_write_colnull_sbyt(fptr, colnum, firstrow, firstelem, nelem, array, nulvalue, status)
    @ccall libcfitsio.ffpcnsb(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Int8}, nulvalue::Int8, status::Ptr{Status})::Status
end

function fits_write_colnull_usht(fptr, colnum, firstrow, firstelem, nelem, array, nulvalue, status)
    @ccall libcfitsio.ffpcnui(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cushort}, nulvalue::Cushort, status::Ptr{Status})::Status
end

function fits_write_colnull_ulng(fptr, colnum, firstrow, firstelem, nelem, array, nulvalue, status)
    @ccall libcfitsio.ffpcnuj(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Culong}, nulvalue::Culong, status::Ptr{Status})::Status
end

function fits_write_colnull_sht(fptr, colnum, firstrow, firstelem, nelem, array, nulvalue, status)
    @ccall libcfitsio.ffpcni(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cshort}, nulvalue::Cshort, status::Ptr{Status})::Status
end

function fits_write_colnull_lng(fptr, colnum, firstrow, firstelem, nelem, array, nulvalue, status)
    @ccall libcfitsio.ffpcnj(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Clong}, nulvalue::Clong, status::Ptr{Status})::Status
end

function fits_write_colnull_ulnglng(fptr, colnum, firstrow, firstelem, nelem, array, nulvalue, status)
    @ccall libcfitsio.ffpcnujj(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Culonglong}, nulvalue::Culonglong, status::Ptr{Status})::Status
end

function fits_write_colnull_lnglng(fptr, colnum, firstrow, firstelem, nelem, array, nulvalue, status)
    @ccall libcfitsio.ffpcnjj(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Clonglong}, nulvalue::Clonglong, status::Ptr{Status})::Status
end

function fits_write_colnull_uint(fptr, colnum, firstrow, firstelem, nelem, array, nulvalue, status)
    @ccall libcfitsio.ffpcnuk(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cuint}, nulvalue::Cuint, status::Ptr{Status})::Status
end

function fits_write_colnull_int(fptr, colnum, firstrow, firstelem, nelem, array, nulvalue, status)
    @ccall libcfitsio.ffpcnk(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cint}, nulvalue::Cint, status::Ptr{Status})::Status
end

function fits_write_colnull_flt(fptr, colnum, firstrow, firstelem, nelem, array, nulvalue, status)
    @ccall libcfitsio.ffpcne(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cfloat}, nulvalue::Cfloat, status::Ptr{Status})::Status
end

function fits_write_colnull_dbl(fptr, colnum, firstrow, firstelem, nelem, array, nulvalue, status)
    @ccall libcfitsio.ffpcnd(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Ptr{Cdouble}, nulvalue::Cdouble, status::Ptr{Status})::Status
end

function fits_write_ext(fptr, offset, nelem, array, status)
    @ccall libcfitsio.ffpextn(fptr::Ptr{fitsfile}, offset::Clonglong, nelem::Clonglong, array::Ptr{Cvoid}, status::Ptr{Status})::Status
end

function fits_read_ext(fptr, offset, nelem, array, status)
    @ccall libcfitsio.ffgextn(fptr::Ptr{fitsfile}, offset::Clonglong, nelem::Clonglong, array::Ptr{Cvoid}, status::Ptr{Status})::Status
end

function fits_write_descript(fptr, colnum, rownum, length, heapaddr, status)
    @ccall libcfitsio.ffpdes(fptr::Ptr{fitsfile}, colnum::Cint, rownum::Clonglong, length::Clonglong, heapaddr::Clonglong, status::Ptr{Status})::Status
end

function fits_compress_heap(fptr, status)
    @ccall libcfitsio.ffcmph(fptr::Ptr{fitsfile}, status::Ptr{Status})::Status
end

function fits_test_heap(fptr, heapsize, unused, overlap, valid, status)
    @ccall libcfitsio.fftheap(fptr::Ptr{fitsfile}, heapsize::Ptr{Clonglong}, unused::Ptr{Clonglong}, overlap::Ptr{Clonglong}, valid::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_write_tblbytes(fptr, firstrow, firstchar, nchars, values, status)
    @ccall libcfitsio.ffptbb(fptr::Ptr{fitsfile}, firstrow::Clonglong, firstchar::Clonglong, nchars::Clonglong, values::Ptr{Cuchar}, status::Ptr{Status})::Status
end

function fits_insert_rows(fptr, firstrow, nrows, status)
    @ccall libcfitsio.ffirow(fptr::Ptr{fitsfile}, firstrow::Clonglong, nrows::Clonglong, status::Ptr{Status})::Status
end

function fits_delete_rows(fptr, firstrow, nrows, status)
    @ccall libcfitsio.ffdrow(fptr::Ptr{fitsfile}, firstrow::Clonglong, nrows::Clonglong, status::Ptr{Status})::Status
end

function fits_delete_rowrange(fptr, ranges, status)
    @ccall libcfitsio.ffdrrg(fptr::Ptr{fitsfile}, ranges::Cstring, status::Ptr{Status})::Status
end

function fits_delete_rowlist(fptr, rownum, nrows, status)
    @ccall libcfitsio.ffdrws(fptr::Ptr{fitsfile}, rownum::Ptr{Clong}, nrows::Clong, status::Ptr{Status})::Status
end

function fits_delete_rowlistll(fptr, rownum, nrows, status)
    @ccall libcfitsio.ffdrwsll(fptr::Ptr{fitsfile}, rownum::Ptr{Clonglong}, nrows::Clonglong, status::Ptr{Status})::Status
end

function fits_insert_col(fptr, numcol, ttype, tform, status)
    @ccall libcfitsio.fficol(fptr::Ptr{fitsfile}, numcol::Cint, ttype::Cstring, tform::Cstring, status::Ptr{Status})::Status
end

function fits_insert_cols(fptr, firstcol, ncols, ttype, tform, status)
    @ccall libcfitsio.fficls(fptr::Ptr{fitsfile}, firstcol::Cint, ncols::Cint, ttype::Ptr{Cstring}, tform::Ptr{Cstring}, status::Ptr{Status})::Status
end

function fits_delete_col(fptr, numcol, status)
    @ccall libcfitsio.ffdcol(fptr::Ptr{fitsfile}, numcol::Cint, status::Ptr{Status})::Status
end

function fits_copy_col(infptr, outfptr, incol, outcol, create_col, status)
    @ccall libcfitsio.ffcpcl(infptr::Ptr{fitsfile}, outfptr::Ptr{fitsfile}, incol::Cint, outcol::Cint, create_col::Cint, status::Ptr{Status})::Status
end

function fits_copy_cols(infptr, outfptr, incol, outcol, ncols, create_col, status)
    @ccall libcfitsio.ffccls(infptr::Ptr{fitsfile}, outfptr::Ptr{fitsfile}, incol::Cint, outcol::Cint, ncols::Cint, create_col::Cint, status::Ptr{Status})::Status
end

function fits_copy_rows(infptr, outfptr, firstrow, nrows, status)
    @ccall libcfitsio.ffcprw(infptr::Ptr{fitsfile}, outfptr::Ptr{fitsfile}, firstrow::Clonglong, nrows::Clonglong, status::Ptr{Status})::Status
end

function fits_copy_selrows(infptr, outfptr, firstrow, nrows, row_status, status)
    @ccall libcfitsio.ffcpsr(infptr::Ptr{fitsfile}, outfptr::Ptr{fitsfile}, firstrow::Clonglong, nrows::Clonglong, row_status::Cstring, status::Ptr{Status})::Status
end

function fits_modify_vector_len(fptr, colnum, newveclen, status)
    @ccall libcfitsio.ffmvec(fptr::Ptr{fitsfile}, colnum::Cint, newveclen::Clonglong, status::Ptr{Status})::Status
end

function fits_read_img_coord(fptr, xrval, yrval, xrpix, yrpix, xinc, yinc, rot, type, status)
    @ccall libcfitsio.ffgics(fptr::Ptr{fitsfile}, xrval::Ptr{Cdouble}, yrval::Ptr{Cdouble}, xrpix::Ptr{Cdouble}, yrpix::Ptr{Cdouble}, xinc::Ptr{Cdouble}, yinc::Ptr{Cdouble}, rot::Ptr{Cdouble}, type::Cstring, status::Ptr{Status})::Status
end

function fits_read_img_coord_version(fptr, version, xrval, yrval, xrpix, yrpix, xinc, yinc, rot, type, status)
    @ccall libcfitsio.ffgicsa(fptr::Ptr{fitsfile}, version::Cchar, xrval::Ptr{Cdouble}, yrval::Ptr{Cdouble}, xrpix::Ptr{Cdouble}, yrpix::Ptr{Cdouble}, xinc::Ptr{Cdouble}, yinc::Ptr{Cdouble}, rot::Ptr{Cdouble}, type::Cstring, status::Ptr{Status})::Status
end

function fits_read_tbl_coord(fptr, xcol, ycol, xrval, yrval, xrpix, yrpix, xinc, yinc, rot, type, status)
    @ccall libcfitsio.ffgtcs(fptr::Ptr{fitsfile}, xcol::Cint, ycol::Cint, xrval::Ptr{Cdouble}, yrval::Ptr{Cdouble}, xrpix::Ptr{Cdouble}, yrpix::Ptr{Cdouble}, xinc::Ptr{Cdouble}, yinc::Ptr{Cdouble}, rot::Ptr{Cdouble}, type::Cstring, status::Ptr{Status})::Status
end

function fits_pix_to_world(xpix, ypix, xref, yref, xrefpix, yrefpix, xinc, yinc, rot, type, xpos, ypos, status)
    @ccall libcfitsio.ffwldp(xpix::Cdouble, ypix::Cdouble, xref::Cdouble, yref::Cdouble, xrefpix::Cdouble, yrefpix::Cdouble, xinc::Cdouble, yinc::Cdouble, rot::Cdouble, type::Cstring, xpos::Ptr{Cdouble}, ypos::Ptr{Cdouble}, status::Ptr{Status})::Status
end

function fits_world_to_pix(xpos, ypos, xref, yref, xrefpix, yrefpix, xinc, yinc, rot, type, xpix, ypix, status)
    @ccall libcfitsio.ffxypx(xpos::Cdouble, ypos::Cdouble, xref::Cdouble, yref::Cdouble, xrefpix::Cdouble, yrefpix::Cdouble, xinc::Cdouble, yinc::Cdouble, rot::Cdouble, type::Cstring, xpix::Ptr{Cdouble}, ypix::Ptr{Cdouble}, status::Ptr{Status})::Status
end

function fits_get_image_wcs_keys(fptr, header, status)
    @ccall libcfitsio.ffgiwcs(fptr::Ptr{fitsfile}, header::Ptr{Cstring}, status::Ptr{Status})::Status
end

function fits_get_table_wcs_keys(fptr, xcol, ycol, header, status)
    @ccall libcfitsio.ffgtwcs(fptr::Ptr{fitsfile}, xcol::Cint, ycol::Cint, header::Ptr{Cstring}, status::Ptr{Status})::Status
end

function fits_find_rows(infptr, expr, firstrow, nrows, n_good_rows, row_status, status)
    @ccall libcfitsio.fffrow(infptr::Ptr{fitsfile}, expr::Cstring, firstrow::Clong, nrows::Clong, n_good_rows::Ptr{Clong}, row_status::Cstring, status::Ptr{Status})::Status
end

function fits_find_first_row(fptr, expr, rownum, status)
    @ccall libcfitsio.ffffrw(fptr::Ptr{fitsfile}, expr::Cstring, rownum::Ptr{Clong}, status::Ptr{Status})::Status
end

function fits_find_rows_cmp(fptr, expr, timeCol, parCol, valCol, ntimes, times, time_status, status)
    @ccall libcfitsio.fffrwc(fptr::Ptr{fitsfile}, expr::Cstring, timeCol::Cstring, parCol::Cstring, valCol::Cstring, ntimes::Clong, times::Ptr{Cdouble}, time_status::Cstring, status::Ptr{Status})::Status
end

function fits_select_rows(infptr, outfptr, expr, status)
    @ccall libcfitsio.ffsrow(infptr::Ptr{fitsfile}, outfptr::Ptr{fitsfile}, expr::Cstring, status::Ptr{Status})::Status
end

function fits_calc_rows(fptr, datatype, expr, firstrow, nelements, nulval, array, anynul, status)
    @ccall libcfitsio.ffcrow(fptr::Ptr{fitsfile}, datatype::Cint, expr::Cstring, firstrow::Clong, nelements::Clong, nulval::Ptr{Cvoid}, array::Ptr{Cvoid}, anynul::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_calculator(infptr, expr, outfptr, parName, parInfo, status)
    @ccall libcfitsio.ffcalc(infptr::Ptr{fitsfile}, expr::Cstring, outfptr::Ptr{fitsfile}, parName::Cstring, parInfo::Cstring, status::Ptr{Status})::Status
end

function fits_calculator_rng(infptr, expr, outfptr, parName, parInfo, nRngs, start, _end, status)
    @ccall libcfitsio.ffcalc_rng(infptr::Ptr{fitsfile}, expr::Cstring, outfptr::Ptr{fitsfile}, parName::Cstring, parInfo::Cstring, nRngs::Cint, start::Ptr{Clong}, _end::Ptr{Clong}, status::Ptr{Status})::Status
end

function fits_test_expr(fptr, expr, maxdim, datatype, nelem, naxis, naxes, status)
    @ccall libcfitsio.fftexp(fptr::Ptr{fitsfile}, expr::Cstring, maxdim::Cint, datatype::Ptr{Cint}, nelem::Ptr{Clong}, naxis::Ptr{Cint}, naxes::Ptr{Clong}, status::Ptr{Status})::Status
end

function fits_create_group(fptr, grpname, grouptype, status)
    @ccall libcfitsio.ffgtcr(fptr::Ptr{fitsfile}, grpname::Cstring, grouptype::Cint, status::Ptr{Status})::Status
end

function fits_insert_group(fptr, grpname, grouptype, status)
    @ccall libcfitsio.ffgtis(fptr::Ptr{fitsfile}, grpname::Cstring, grouptype::Cint, status::Ptr{Status})::Status
end

function fits_change_group(gfptr, grouptype, status)
    @ccall libcfitsio.ffgtch(gfptr::Ptr{fitsfile}, grouptype::Cint, status::Ptr{Status})::Status
end

function fits_remove_group(gfptr, rmopt, status)
    @ccall libcfitsio.ffgtrm(gfptr::Ptr{fitsfile}, rmopt::Cint, status::Ptr{Status})::Status
end

function fits_copy_group(infptr, outfptr, cpopt, status)
    @ccall libcfitsio.ffgtcp(infptr::Ptr{fitsfile}, outfptr::Ptr{fitsfile}, cpopt::Cint, status::Ptr{Status})::Status
end

function fits_merge_groups(infptr, outfptr, mgopt, status)
    @ccall libcfitsio.ffgtmg(infptr::Ptr{fitsfile}, outfptr::Ptr{fitsfile}, mgopt::Cint, status::Ptr{Status})::Status
end

function fits_compact_group(gfptr, cmopt, status)
    @ccall libcfitsio.ffgtcm(gfptr::Ptr{fitsfile}, cmopt::Cint, status::Ptr{Status})::Status
end

function fits_verify_group(gfptr, firstfailed, status)
    @ccall libcfitsio.ffgtvf(gfptr::Ptr{fitsfile}, firstfailed::Ptr{Clong}, status::Ptr{Status})::Status
end

function fits_open_group(mfptr, group, gfptr, status)
    @ccall libcfitsio.ffgtop(mfptr::Ptr{fitsfile}, group::Cint, gfptr::Ptr{Ptr{fitsfile}}, status::Ptr{Status})::Status
end

function fits_add_group_member(gfptr, mfptr, hdupos, status)
    @ccall libcfitsio.ffgtam(gfptr::Ptr{fitsfile}, mfptr::Ptr{fitsfile}, hdupos::Cint, status::Ptr{Status})::Status
end

function fits_get_num_members(gfptr, nmembers, status)
    @ccall libcfitsio.ffgtnm(gfptr::Ptr{fitsfile}, nmembers::Ptr{Clong}, status::Ptr{Status})::Status
end

function fits_get_num_groups(mfptr, nmembers, status)
    @ccall libcfitsio.ffgmng(mfptr::Ptr{fitsfile}, nmembers::Ptr{Clong}, status::Ptr{Status})::Status
end

function fits_open_member(gfptr, member, mfptr, status)
    @ccall libcfitsio.ffgmop(gfptr::Ptr{fitsfile}, member::Clong, mfptr::Ptr{Ptr{fitsfile}}, status::Ptr{Status})::Status
end

function fits_copy_member(gfptr, mfptr, member, cpopt, status)
    @ccall libcfitsio.ffgmcp(gfptr::Ptr{fitsfile}, mfptr::Ptr{fitsfile}, member::Clong, cpopt::Cint, status::Ptr{Status})::Status
end

function fits_transfer_member(infptr, outfptr, member, tfopt, status)
    @ccall libcfitsio.ffgmtf(infptr::Ptr{fitsfile}, outfptr::Ptr{fitsfile}, member::Clong, tfopt::Cint, status::Ptr{Status})::Status
end

function fits_remove_member(fptr, member, rmopt, status)
    @ccall libcfitsio.ffgmrm(fptr::Ptr{fitsfile}, member::Clong, rmopt::Cint, status::Ptr{Status})::Status
end

function fits_init_https()
    @ccall libcfitsio.ffihtps()::Cint
end

function fits_cleanup_https()
    @ccall libcfitsio.ffchtps()::Cint
end

function fits_verbose_https(flag)
    @ccall libcfitsio.ffvhtps(flag::Cint)::Cvoid
end

function fits_show_download_progress(flag)
    @ccall libcfitsio.ffshdwn(flag::Cint)::Cvoid
end

function fits_get_timeout()
    @ccall libcfitsio.ffgtmo()::Cint
end

function fits_set_timeout(sec, status)
    @ccall libcfitsio.ffstmo(sec::Cint, status::Ptr{Status})::Status
end

struct wtbarr
    i::Cint
    m::Cint
    kind::Cint
    extnam::NTuple{72, Cchar}
    extver::Cint
    extlev::Cint
    ttype::NTuple{72, Cchar}
    row::Clong
    ndim::Cint
    dimlen::Ptr{Cint}
    arrayp::Ptr{Ptr{Cdouble}}
end

function fits_read_wcstab(fptr, nwtb, wtb, status)
    @ccall libcfitsio.fits_read_wcstab(fptr::Ptr{fitsfile}, nwtb::Cint, wtb::Ptr{wtbarr}, status::Ptr{Status})::Status
end

function CFITS2Unit(fptr)
    @ccall libcfitsio.CFITS2Unit(fptr::Ptr{fitsfile})::Cint
end

function CUnit2FITS(unit)
    @ccall libcfitsio.CUnit2FITS(unit::Cint)::Ptr{fitsfile}
end

function fits_get_token(ptr, delimiter, token, isanumber)
    @ccall libcfitsio.fits_get_token(ptr::Ptr{Cstring}, delimiter::Cstring, token::Cstring, isanumber::Ptr{Cint})::Cint
end

function fits_get_token2(ptr, delimiter, token, isanumber, status)
    @ccall libcfitsio.fits_get_token2(ptr::Ptr{Cstring}, delimiter::Cstring, token::Ptr{Cstring}, isanumber::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_split_names(list)
    @ccall libcfitsio.fits_split_names(list::Cstring)::Cstring
end

function ffifile2(url, urltype, infile, outfile, extspec, rowfilter, binspec, colspec, pixfilter, compspec, status)
    @ccall libcfitsio.ffifile2(url::Cstring, urltype::Cstring, infile::Cstring, outfile::Cstring, extspec::Cstring, rowfilter::Cstring, binspec::Cstring, colspec::Cstring, pixfilter::Cstring, compspec::Cstring, status::Ptr{Status})::Status
end

function fits_copy_cell2image(fptr, newptr, colname, rownum, status)
    @ccall libcfitsio.fits_copy_cell2image(fptr::Ptr{fitsfile}, newptr::Ptr{fitsfile}, colname::Cstring, rownum::Clong, status::Ptr{Status})::Status
end

function fits_copy_image2cell(fptr, newptr, colname, rownum, copykeyflag, status)
    @ccall libcfitsio.fits_copy_image2cell(fptr::Ptr{fitsfile}, newptr::Ptr{fitsfile}, colname::Cstring, rownum::Clong, copykeyflag::Cint, status::Ptr{Status})::Status
end

function fits_copy_pixlist2image(infptr, outfptr, firstkey, naxis, colnum, status)
    @ccall libcfitsio.fits_copy_pixlist2image(infptr::Ptr{fitsfile}, outfptr::Ptr{fitsfile}, firstkey::Cint, naxis::Cint, colnum::Ptr{Cint}, status::Ptr{Status})::Status
end

function ffimport_file(filename, contents, status)
    @ccall libcfitsio.ffimport_file(filename::Cstring, contents::Ptr{Cstring}, status::Ptr{Status})::Status
end

function fits_init_cfitsio()
    @ccall libcfitsio.fits_init_cfitsio()::Cint
end

function ffopen(fptr, filename, iomode, status)
    @ccall libcfitsio.ffopen(fptr::Ptr{Ptr{fitsfile}}, filename::Cstring, iomode::Cint, status::Ptr{Status})::Status
end

function fits_delete_iraf_file(filename, status)
    @ccall libcfitsio.fits_delete_iraf_file(filename::Cstring, status::Ptr{Status})::Status
end

function fits_translate_keyword(inrec, outrec, patterns, npat, n_value, n_offset, n_range, pat_num, i, j, m, n, status)
    @ccall libcfitsio.fits_translate_keyword(inrec::Cstring, outrec::Cstring, patterns::Ptr{NTuple{2, Cstring}}, npat::Cint, n_value::Cint, n_offset::Cint, n_range::Cint, pat_num::Ptr{Cint}, i::Ptr{Cint}, j::Ptr{Cint}, m::Ptr{Cint}, n::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_translate_keywords(infptr, outfptr, firstkey, patterns, npat, n_value, n_offset, n_range, status)
    @ccall libcfitsio.fits_translate_keywords(infptr::Ptr{fitsfile}, outfptr::Ptr{fitsfile}, firstkey::Cint, patterns::Ptr{NTuple{2, Cstring}}, npat::Cint, n_value::Cint, n_offset::Cint, n_range::Cint, status::Ptr{Status})::Status
end

function fits_get_section_range(ptr, secmin, secmax, incre, status)
    @ccall libcfitsio.fits_get_section_range(ptr::Ptr{Cstring}, secmin::Ptr{Clong}, secmax::Ptr{Clong}, incre::Ptr{Clong}, status::Ptr{Status})::Status
end

function ffmbyt(fptr, bytpos, ignore_err, status)
    @ccall libcfitsio.ffmbyt(fptr::Ptr{fitsfile}, bytpos::Clonglong, ignore_err::Cint, status::Ptr{Status})::Status
end

function ffverifydate(year, month, day, status)
    @ccall libcfitsio.ffverifydate(year::Cint, month::Cint, day::Cint, status::Ptr{Status})::Status
end

function ffpknjj(fptr, keyroot, nstart, nkey, value, comm, status)
    @ccall libcfitsio.ffpknjj(fptr::Ptr{fitsfile}, keyroot::Cstring, nstart::Cint, nkey::Cint, value::Ptr{Clonglong}, comm::Ptr{Cstring}, status::Ptr{Status})::Status
end

function fffkls(value, status)
    @ccall libcfitsio.fffkls(value::Cstring, status::Ptr{Status})::Status
end

function ffh2st(fptr, header, status)
    @ccall libcfitsio.ffh2st(fptr::Ptr{fitsfile}, header::Ptr{Cstring}, status::Ptr{Status})::Status
end

function ffchfl(fptr, status)
    @ccall libcfitsio.ffchfl(fptr::Ptr{fitsfile}, status::Ptr{Status})::Status
end

function ffcdfl(fptr, status)
    @ccall libcfitsio.ffcdfl(fptr::Ptr{fitsfile}, status::Ptr{Status})::Status
end

function ffrhdu(fptr, hdutype, status)
    @ccall libcfitsio.ffrhdu(fptr::Ptr{fitsfile}, hdutype::Ptr{Cint}, status::Ptr{Status})::Status
end

function ffcsum(fptr, nrec, sum, status)
    @ccall libcfitsio.ffcsum(fptr::Ptr{fitsfile}, nrec::Clong, sum::Ptr{Culong}, status::Ptr{Status})::Status
end

function ffgcl(fptr, colnum, firstrow, firstelem, nelem, array, status)
    @ccall libcfitsio.ffgcl(fptr::Ptr{fitsfile}, colnum::Cint, firstrow::Clonglong, firstelem::Clonglong, nelem::Clonglong, array::Cstring, status::Ptr{Status})::Status
end

function fits_iter_set_by_name(col, fptr, colname, datatype, iotype)
    @ccall libcfitsio.fits_iter_set_by_name(col::Ptr{iteratorCol}, fptr::Ptr{fitsfile}, colname::Cstring, datatype::Cint, iotype::Cint)::Cint
end

function fits_iter_set_by_num(col, fptr, colnum, datatype, iotype)
    @ccall libcfitsio.fits_iter_set_by_num(col::Ptr{iteratorCol}, fptr::Ptr{fitsfile}, colnum::Cint, datatype::Cint, iotype::Cint)::Cint
end

function fits_iter_set_file(col, fptr)
    @ccall libcfitsio.fits_iter_set_file(col::Ptr{iteratorCol}, fptr::Ptr{fitsfile})::Cint
end

function fits_iter_set_colname(col, colname)
    @ccall libcfitsio.fits_iter_set_colname(col::Ptr{iteratorCol}, colname::Cstring)::Cint
end

function fits_iter_set_colnum(col, colnum)
    @ccall libcfitsio.fits_iter_set_colnum(col::Ptr{iteratorCol}, colnum::Cint)::Cint
end

function fits_iter_set_datatype(col, datatype)
    @ccall libcfitsio.fits_iter_set_datatype(col::Ptr{iteratorCol}, datatype::Cint)::Cint
end

function fits_iter_set_iotype(col, iotype)
    @ccall libcfitsio.fits_iter_set_iotype(col::Ptr{iteratorCol}, iotype::Cint)::Cint
end

function fits_iter_get_file(col)
    @ccall libcfitsio.fits_iter_get_file(col::Ptr{iteratorCol})::Ptr{fitsfile}
end

function fits_iter_get_colname(col)
    @ccall libcfitsio.fits_iter_get_colname(col::Ptr{iteratorCol})::Cstring
end

function fits_iter_get_colnum(col)
    @ccall libcfitsio.fits_iter_get_colnum(col::Ptr{iteratorCol})::Cint
end

function fits_iter_get_datatype(col)
    @ccall libcfitsio.fits_iter_get_datatype(col::Ptr{iteratorCol})::Cint
end

function fits_iter_get_iotype(col)
    @ccall libcfitsio.fits_iter_get_iotype(col::Ptr{iteratorCol})::Cint
end

function fits_iter_get_array(col)
    @ccall libcfitsio.fits_iter_get_array(col::Ptr{iteratorCol})::Ptr{Cvoid}
end

function fits_iter_get_tlmin(col)
    @ccall libcfitsio.fits_iter_get_tlmin(col::Ptr{iteratorCol})::Clong
end

function fits_iter_get_tlmax(col)
    @ccall libcfitsio.fits_iter_get_tlmax(col::Ptr{iteratorCol})::Clong
end

function fits_iter_get_repeat(col)
    @ccall libcfitsio.fits_iter_get_repeat(col::Ptr{iteratorCol})::Clong
end

function fits_iter_get_tunit(col)
    @ccall libcfitsio.fits_iter_get_tunit(col::Ptr{iteratorCol})::Cstring
end

function fits_iter_get_tdisp(col)
    @ccall libcfitsio.fits_iter_get_tdisp(col::Ptr{iteratorCol})::Cstring
end

function ffhist(fptr, outfile, imagetype, naxis, colname, minin, maxin, binsizein, minname, maxname, binname, weightin, wtcol, recip, rowselect, status)
    @ccall libcfitsio.ffhist(fptr::Ptr{Ptr{fitsfile}}, outfile::Cstring, imagetype::Cint, naxis::Cint, colname::Ptr{NTuple{71, Cchar}}, minin::Ptr{Cdouble}, maxin::Ptr{Cdouble}, binsizein::Ptr{Cdouble}, minname::Ptr{NTuple{71, Cchar}}, maxname::Ptr{NTuple{71, Cchar}}, binname::Ptr{NTuple{71, Cchar}}, weightin::Cdouble, wtcol::Ptr{Cchar}, recip::Cint, rowselect::Cstring, status::Ptr{Status})::Status
end

function ffhist2(fptr, outfile, imagetype, naxis, colname, minin, maxin, binsizein, minname, maxname, binname, weightin, wtcol, recip, rowselect, status)
    @ccall libcfitsio.ffhist2(fptr::Ptr{Ptr{fitsfile}}, outfile::Cstring, imagetype::Cint, naxis::Cint, colname::Ptr{NTuple{71, Cchar}}, minin::Ptr{Cdouble}, maxin::Ptr{Cdouble}, binsizein::Ptr{Cdouble}, minname::Ptr{NTuple{71, Cchar}}, maxname::Ptr{NTuple{71, Cchar}}, binname::Ptr{NTuple{71, Cchar}}, weightin::Cdouble, wtcol::Ptr{Cchar}, recip::Cint, rowselect::Cstring, status::Ptr{Status})::Status
end

function ffhist3(fptr, outfile, imagetype, naxis, colname, minin, maxin, binsizein, minname, maxname, binname, weightin, wtcol, recip, selectrow, status)
    @ccall libcfitsio.ffhist3(fptr::Ptr{fitsfile}, outfile::Cstring, imagetype::Cint, naxis::Cint, colname::Ptr{NTuple{71, Cchar}}, minin::Ptr{Cdouble}, maxin::Ptr{Cdouble}, binsizein::Ptr{Cdouble}, minname::Ptr{NTuple{71, Cchar}}, maxname::Ptr{NTuple{71, Cchar}}, binname::Ptr{NTuple{71, Cchar}}, weightin::Cdouble, wtcol::Ptr{Cchar}, recip::Cint, selectrow::Cstring, status::Ptr{Cint})::Ptr{fitsfile}
end

function fits_select_image_section(fptr, outfile, imagesection, status)
    @ccall libcfitsio.fits_select_image_section(fptr::Ptr{Ptr{fitsfile}}, outfile::Cstring, imagesection::Cstring, status::Ptr{Status})::Status
end

function fits_calc_binning(fptr, naxis, colname, minin, maxin, binsizein, minname, maxname, binname, colnum, haxes, amin, amax, binsize, status)
    @ccall libcfitsio.fits_calc_binning(fptr::Ptr{fitsfile}, naxis::Cint, colname::Ptr{NTuple{71, Cchar}}, minin::Ptr{Cdouble}, maxin::Ptr{Cdouble}, binsizein::Ptr{Cdouble}, minname::Ptr{NTuple{71, Cchar}}, maxname::Ptr{NTuple{71, Cchar}}, binname::Ptr{NTuple{71, Cchar}}, colnum::Ptr{Cint}, haxes::Ptr{Clong}, amin::Ptr{Cfloat}, amax::Ptr{Cfloat}, binsize::Ptr{Cfloat}, status::Ptr{Status})::Status
end

function fits_calc_binningd(fptr, naxis, colname, minin, maxin, binsizein, minname, maxname, binname, colnum, haxes, amin, amax, binsize, status)
    @ccall libcfitsio.fits_calc_binningd(fptr::Ptr{fitsfile}, naxis::Cint, colname::Ptr{NTuple{71, Cchar}}, minin::Ptr{Cdouble}, maxin::Ptr{Cdouble}, binsizein::Ptr{Cdouble}, minname::Ptr{NTuple{71, Cchar}}, maxname::Ptr{NTuple{71, Cchar}}, binname::Ptr{NTuple{71, Cchar}}, colnum::Ptr{Cint}, haxes::Ptr{Clong}, amin::Ptr{Cdouble}, amax::Ptr{Cdouble}, binsize::Ptr{Cdouble}, status::Ptr{Status})::Status
end

function fits_write_keys_histo(fptr, histptr, naxis, colnum, status)
    @ccall libcfitsio.fits_write_keys_histo(fptr::Ptr{fitsfile}, histptr::Ptr{fitsfile}, naxis::Cint, colnum::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_rebin_wcs(fptr, naxis, amin, binsize, status)
    @ccall libcfitsio.fits_rebin_wcs(fptr::Ptr{fitsfile}, naxis::Cint, amin::Ptr{Cfloat}, binsize::Ptr{Cfloat}, status::Ptr{Status})::Status
end

function fits_rebin_wcsd(fptr, naxis, amin, binsize, status)
    @ccall libcfitsio.fits_rebin_wcsd(fptr::Ptr{fitsfile}, naxis::Cint, amin::Ptr{Cdouble}, binsize::Ptr{Cdouble}, status::Ptr{Status})::Status
end

function fits_make_hist(fptr, histptr, bitpix, naxis, naxes, colnum, amin, amax, binsize, weight, wtcolnum, recip, selectrow, status)
    @ccall libcfitsio.fits_make_hist(fptr::Ptr{fitsfile}, histptr::Ptr{fitsfile}, bitpix::Cint, naxis::Cint, naxes::Ptr{Clong}, colnum::Ptr{Cint}, amin::Ptr{Cfloat}, amax::Ptr{Cfloat}, binsize::Ptr{Cfloat}, weight::Cfloat, wtcolnum::Cint, recip::Cint, selectrow::Cstring, status::Ptr{Status})::Status
end

function fits_make_histd(fptr, histptr, bitpix, naxis, naxes, colnum, amin, amax, binsize, weight, wtcolnum, recip, selectrow, status)
    @ccall libcfitsio.fits_make_histd(fptr::Ptr{fitsfile}, histptr::Ptr{fitsfile}, bitpix::Cint, naxis::Cint, naxes::Ptr{Clong}, colnum::Ptr{Cint}, amin::Ptr{Cdouble}, amax::Ptr{Cdouble}, binsize::Ptr{Cdouble}, weight::Cdouble, wtcolnum::Cint, recip::Cint, selectrow::Cstring, status::Ptr{Status})::Status
end

struct PixelFilter
    count::Cint
    path::Ptr{Cstring}
    tag::Ptr{Cstring}
    ifptr::Ptr{Ptr{fitsfile}}
    expression::Cstring
    bitpix::Cint
    blank::Clong
    ofptr::Ptr{fitsfile}
    keyword::NTuple{75, Cchar}
    comment::NTuple{73, Cchar}
end

function fits_pixel_filter(filter, status)
    @ccall libcfitsio.fits_pixel_filter(filter::Ptr{PixelFilter}, status::Ptr{Status})::Status
end

function fits_execute_template(ff, ngp_template, status)
    @ccall libcfitsio.fits_execute_template(ff::Ptr{fitsfile}, ngp_template::Cstring, status::Ptr{Status})::Status
end

function fits_img_stats_short(array, nx, ny, nullcheck, nullvalue, ngoodpix, minvalue, maxvalue, mean, sigma, noise1, noise2, noise3, noise5, status)
    @ccall libcfitsio.fits_img_stats_short(array::Ptr{Cshort}, nx::Clong, ny::Clong, nullcheck::Cint, nullvalue::Cshort, ngoodpix::Ptr{Clong}, minvalue::Ptr{Cshort}, maxvalue::Ptr{Cshort}, mean::Ptr{Cdouble}, sigma::Ptr{Cdouble}, noise1::Ptr{Cdouble}, noise2::Ptr{Cdouble}, noise3::Ptr{Cdouble}, noise5::Ptr{Cdouble}, status::Ptr{Status})::Status
end

function fits_img_stats_int(array, nx, ny, nullcheck, nullvalue, ngoodpix, minvalue, maxvalue, mean, sigma, noise1, noise2, noise3, noise5, status)
    @ccall libcfitsio.fits_img_stats_int(array::Ptr{Cint}, nx::Clong, ny::Clong, nullcheck::Cint, nullvalue::Cint, ngoodpix::Ptr{Clong}, minvalue::Ptr{Cint}, maxvalue::Ptr{Cint}, mean::Ptr{Cdouble}, sigma::Ptr{Cdouble}, noise1::Ptr{Cdouble}, noise2::Ptr{Cdouble}, noise3::Ptr{Cdouble}, noise5::Ptr{Cdouble}, status::Ptr{Status})::Status
end

function fits_img_stats_float(array, nx, ny, nullcheck, nullvalue, ngoodpix, minvalue, maxvalue, mean, sigma, noise1, noise2, noise3, noise5, status)
    @ccall libcfitsio.fits_img_stats_float(array::Ptr{Cfloat}, nx::Clong, ny::Clong, nullcheck::Cint, nullvalue::Cfloat, ngoodpix::Ptr{Clong}, minvalue::Ptr{Cfloat}, maxvalue::Ptr{Cfloat}, mean::Ptr{Cdouble}, sigma::Ptr{Cdouble}, noise1::Ptr{Cdouble}, noise2::Ptr{Cdouble}, noise3::Ptr{Cdouble}, noise5::Ptr{Cdouble}, status::Ptr{Status})::Status
end

function fits_set_compression_type(fptr, ctype, status)
    @ccall libcfitsio.fits_set_compression_type(fptr::Ptr{fitsfile}, ctype::Cint, status::Ptr{Status})::Status
end

function fits_set_tile_dim(fptr, ndim, dims, status)
    @ccall libcfitsio.fits_set_tile_dim(fptr::Ptr{fitsfile}, ndim::Cint, dims::Ptr{Clong}, status::Ptr{Status})::Status
end

function fits_set_noise_bits(fptr, noisebits, status)
    @ccall libcfitsio.fits_set_noise_bits(fptr::Ptr{fitsfile}, noisebits::Cint, status::Ptr{Status})::Status
end

function fits_set_quantize_level(fptr, qlevel, status)
    @ccall libcfitsio.fits_set_quantize_level(fptr::Ptr{fitsfile}, qlevel::Cfloat, status::Ptr{Status})::Status
end

function fits_set_hcomp_scale(fptr, scale, status)
    @ccall libcfitsio.fits_set_hcomp_scale(fptr::Ptr{fitsfile}, scale::Cfloat, status::Ptr{Status})::Status
end

function fits_set_hcomp_smooth(fptr, smooth, status)
    @ccall libcfitsio.fits_set_hcomp_smooth(fptr::Ptr{fitsfile}, smooth::Cint, status::Ptr{Status})::Status
end

function fits_set_quantize_method(fptr, method, status)
    @ccall libcfitsio.fits_set_quantize_method(fptr::Ptr{fitsfile}, method::Cint, status::Ptr{Status})::Status
end

function fits_set_quantize_dither(fptr, dither, status)
    @ccall libcfitsio.fits_set_quantize_dither(fptr::Ptr{fitsfile}, dither::Cint, status::Ptr{Status})::Status
end

function fits_set_dither_seed(fptr, seed, status)
    @ccall libcfitsio.fits_set_dither_seed(fptr::Ptr{fitsfile}, seed::Cint, status::Ptr{Status})::Status
end

function fits_set_dither_offset(fptr, offset, status)
    @ccall libcfitsio.fits_set_dither_offset(fptr::Ptr{fitsfile}, offset::Cint, status::Ptr{Status})::Status
end

function fits_set_lossy_int(fptr, lossy_int, status)
    @ccall libcfitsio.fits_set_lossy_int(fptr::Ptr{fitsfile}, lossy_int::Cint, status::Ptr{Status})::Status
end

function fits_set_huge_hdu(fptr, huge, status)
    @ccall libcfitsio.fits_set_huge_hdu(fptr::Ptr{fitsfile}, huge::Cint, status::Ptr{Status})::Status
end

function fits_set_compression_pref(infptr, outfptr, status)
    @ccall libcfitsio.fits_set_compression_pref(infptr::Ptr{fitsfile}, outfptr::Ptr{fitsfile}, status::Ptr{Status})::Status
end

function fits_get_compression_type(fptr, ctype, status)
    @ccall libcfitsio.fits_get_compression_type(fptr::Ptr{fitsfile}, ctype::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_get_tile_dim(fptr, ndim, dims, status)
    @ccall libcfitsio.fits_get_tile_dim(fptr::Ptr{fitsfile}, ndim::Cint, dims::Ptr{Clong}, status::Ptr{Status})::Status
end

function fits_get_quantize_level(fptr, qlevel, status)
    @ccall libcfitsio.fits_get_quantize_level(fptr::Ptr{fitsfile}, qlevel::Ptr{Cfloat}, status::Ptr{Status})::Status
end

function fits_get_noise_bits(fptr, noisebits, status)
    @ccall libcfitsio.fits_get_noise_bits(fptr::Ptr{fitsfile}, noisebits::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_get_hcomp_scale(fptr, scale, status)
    @ccall libcfitsio.fits_get_hcomp_scale(fptr::Ptr{fitsfile}, scale::Ptr{Cfloat}, status::Ptr{Status})::Status
end

function fits_get_hcomp_smooth(fptr, smooth, status)
    @ccall libcfitsio.fits_get_hcomp_smooth(fptr::Ptr{fitsfile}, smooth::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_get_dither_seed(fptr, seed, status)
    @ccall libcfitsio.fits_get_dither_seed(fptr::Ptr{fitsfile}, seed::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_img_compress(infptr, outfptr, status)
    @ccall libcfitsio.fits_img_compress(infptr::Ptr{fitsfile}, outfptr::Ptr{fitsfile}, status::Ptr{Status})::Status
end

function fits_is_compressed_image(fptr, status)
    @ccall libcfitsio.fits_is_compressed_image(fptr::Ptr{fitsfile}, status::Ptr{Status})::Status
end

function fits_is_reentrant()
    @ccall libcfitsio.fits_is_reentrant()::Cint
end

function fits_img_decompress_header(infptr, outfptr, status)
    @ccall libcfitsio.fits_img_decompress_header(infptr::Ptr{fitsfile}, outfptr::Ptr{fitsfile}, status::Ptr{Status})::Status
end

function fits_img_decompress(infptr, outfptr, status)
    @ccall libcfitsio.fits_img_decompress(infptr::Ptr{fitsfile}, outfptr::Ptr{fitsfile}, status::Ptr{Status})::Status
end

function fits_hcompress(a, nx, ny, scale, output, nbytes, status)
    @ccall libcfitsio.fits_hcompress(a::Ptr{Cint}, nx::Cint, ny::Cint, scale::Cint, output::Cstring, nbytes::Ptr{Clong}, status::Ptr{Status})::Status
end

function fits_hcompress64(a, nx, ny, scale, output, nbytes, status)
    @ccall libcfitsio.fits_hcompress64(a::Ptr{Clonglong}, nx::Cint, ny::Cint, scale::Cint, output::Cstring, nbytes::Ptr{Clong}, status::Ptr{Status})::Status
end

function fits_hdecompress(input, smooth, a, nx, ny, scale, status)
    @ccall libcfitsio.fits_hdecompress(input::Ptr{Cuchar}, smooth::Cint, a::Ptr{Cint}, nx::Ptr{Cint}, ny::Ptr{Cint}, scale::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_hdecompress64(input, smooth, a, nx, ny, scale, status)
    @ccall libcfitsio.fits_hdecompress64(input::Ptr{Cuchar}, smooth::Cint, a::Ptr{Clonglong}, nx::Ptr{Cint}, ny::Ptr{Cint}, scale::Ptr{Cint}, status::Ptr{Status})::Status
end

function fits_compress_table(infptr, outfptr, status)
    @ccall libcfitsio.fits_compress_table(infptr::Ptr{fitsfile}, outfptr::Ptr{fitsfile}, status::Ptr{Status})::Status
end

function fits_uncompress_table(infptr, outfptr, status)
    @ccall libcfitsio.fits_uncompress_table(infptr::Ptr{fitsfile}, outfptr::Ptr{fitsfile}, status::Ptr{Status})::Status
end

# Skipping MacroDefinition: CFITSIO_VERSION 4.6.2

const CFITSIO_MICRO = 2

const CFITSIO_MINOR = 6

const CFITSIO_MAJOR = 4

const CFITSIO_SONAME = 10

const OFF_T = off_t

const USE_LL_SUFFIX = 1

const LONGLONG_MAX = typemax(Clonglong)

const LONGLONG_MIN = typemin(Clonglong)

# const ffcpimg = fits_copy_image_section

# const fits_compress_img = fits_comp_img

# const fits_decompress_img = fits_decomp_img

# const fits_write_nulrows = ffprwu

const NIOBUF = 40

const IOBUFLEN = 2880

const FLEN_FILENAME = 1025

const FLEN_KEYWORD = 75

const FLEN_CARD = 81

const FLEN_VALUE = 71

const FLEN_COMMENT = 73

const FLEN_ERRMSG = 81

const FLEN_STATUS = 31

const TBIT = 1

const TBYTE = 11

const TSBYTE = 12

const TLOGICAL = 14

const TSTRING = 16

const TUSHORT = 20

const TSHORT = 21

const TUINT = 30

const TINT = 31

const TULONG = 40

const TLONG = 41

const TINT32BIT = 41

const TFLOAT = 42

const TULONGLONG = 80

const TLONGLONG = 81

const TDOUBLE = 82

const TCOMPLEX = 83

const TDBLCOMPLEX = 163

const TYP_STRUC_KEY = 10

const TYP_CMPRS_KEY = 20

const TYP_SCAL_KEY = 30

const TYP_NULL_KEY = 40

const TYP_DIM_KEY = 50

const TYP_RANG_KEY = 60

const TYP_UNIT_KEY = 70

const TYP_DISP_KEY = 80

const TYP_HDUID_KEY = 90

const TYP_CKSUM_KEY = 100

const TYP_WCS_KEY = 110

const TYP_REFSYS_KEY = 120

const TYP_COMM_KEY = 130

const TYP_CONT_KEY = 140

const TYP_USER_KEY = 150

const INT32BIT = Cint

const BYTE_IMG = 8

const SHORT_IMG = 16

const LONG_IMG = 32

const LONGLONG_IMG = 64

const FLOAT_IMG = -32

const DOUBLE_IMG = -64

const SBYTE_IMG = 10

const USHORT_IMG = 20

const ULONG_IMG = 40

const ULONGLONG_IMG = 80

const IMAGE_HDU = 0

const ASCII_TBL = 1

const BINARY_TBL = 2

const ANY_HDU = -1

const READONLY = 0

const READWRITE = 1

const FLOATNULLVALUE = -(9.11912f-36)

const DOUBLENULLVALUE = -9.1191291391491e-36

const NO_DITHER = -1

const SUBTRACTIVE_DITHER_1 = 1

const SUBTRACTIVE_DITHER_2 = 2

const MAX_COMPRESS_DIM = 6

const RICE_1 = 11

const GZIP_1 = 21

const GZIP_2 = 22

const PLIO_1 = 31

const HCOMPRESS_1 = 41

const BZIP2_1 = 51

const NOCOMPRESS = -1

const TRUE = 1

const FALSE = 0

const CASESEN = 1

const CASEINSEN = 0

const GT_ID_ALL_URI = 0

const GT_ID_REF = 1

const GT_ID_POS = 2

const GT_ID_ALL = 3

const GT_ID_REF_URI = 11

const GT_ID_POS_URI = 12

const OPT_RM_GPT = 0

const OPT_RM_ENTRY = 1

const OPT_RM_MBR = 2

const OPT_RM_ALL = 3

const OPT_GCP_GPT = 0

const OPT_GCP_MBR = 1

const OPT_GCP_ALL = 2

const OPT_MCP_ADD = 0

const OPT_MCP_NADD = 1

const OPT_MCP_REPL = 2

const OPT_MCP_MOV = 3

const OPT_MRG_COPY = 0

const OPT_MRG_MOV = 1

const OPT_CMT_MBR = 1

const OPT_CMT_MBR_DEL = 11

const VALIDSTRUC = 555

const InputCol = 0

const InputOutputCol = 1

const OutputCol = 2

const TemporaryCol = 3

const CREATE_DISK_FILE = -106

const OPEN_DISK_FILE = -105

const SKIP_TABLE = -104

const SKIP_IMAGE = -103

const SKIP_NULL_PRIMARY = -102

const USE_MEM_BUFF = -101

const OVERFLOW_ERR = -11

const PREPEND_PRIMARY = -9

const SAME_FILE = 101

const TOO_MANY_FILES = 103

const FILE_NOT_OPENED = 104

const FILE_NOT_CREATED = 105

const WRITE_ERROR = 106

const END_OF_FILE = 107

const READ_ERROR = 108

const FILE_NOT_CLOSED = 110

const ARRAY_TOO_BIG = 111

const READONLY_FILE = 112

const MEMORY_ALLOCATION = 113

const BAD_FILEPTR = 114

const NULL_INPUT_PTR = 115

const SEEK_ERROR = 116

const BAD_NETTIMEOUT = 117

const BAD_URL_PREFIX = 121

const TOO_MANY_DRIVERS = 122

const DRIVER_INIT_FAILED = 123

const NO_MATCHING_DRIVER = 124

const URL_PARSE_ERROR = 125

const RANGE_PARSE_ERROR = 126

const SHARED_ERRBASE = 150

const SHARED_BADARG = SHARED_ERRBASE + 1

const SHARED_NULPTR = SHARED_ERRBASE + 2

const SHARED_TABFULL = SHARED_ERRBASE + 3

const SHARED_NOTINIT = SHARED_ERRBASE + 4

const SHARED_IPCERR = SHARED_ERRBASE + 5

const SHARED_NOMEM = SHARED_ERRBASE + 6

const SHARED_AGAIN = SHARED_ERRBASE + 7

const SHARED_NOFILE = SHARED_ERRBASE + 8

const SHARED_NORESIZE = SHARED_ERRBASE + 9

const HEADER_NOT_EMPTY = 201

const KEY_NO_EXIST = 202

const KEY_OUT_BOUNDS = 203

const VALUE_UNDEFINED = 204

const NO_QUOTE = 205

const BAD_INDEX_KEY = 206

const BAD_KEYCHAR = 207

const BAD_ORDER = 208

const NOT_POS_INT = 209

const NO_END = 210

const BAD_BITPIX = 211

const BAD_NAXIS = 212

const BAD_NAXES = 213

const BAD_PCOUNT = 214

const BAD_GCOUNT = 215

const BAD_TFIELDS = 216

const NEG_WIDTH = 217

const NEG_ROWS = 218

const COL_NOT_FOUND = 219

const BAD_SIMPLE = 220

const NO_SIMPLE = 221

const NO_BITPIX = 222

const NO_NAXIS = 223

const NO_NAXES = 224

const NO_XTENSION = 225

const NOT_ATABLE = 226

const NOT_BTABLE = 227

const NO_PCOUNT = 228

const NO_GCOUNT = 229

const NO_TFIELDS = 230

const NO_TBCOL = 231

const NO_TFORM = 232

const NOT_IMAGE = 233

const BAD_TBCOL = 234

const NOT_TABLE = 235

const COL_TOO_WIDE = 236

const COL_NOT_UNIQUE = 237

const BAD_ROW_WIDTH = 241

const UNKNOWN_EXT = 251

const UNKNOWN_REC = 252

const END_JUNK = 253

const BAD_HEADER_FILL = 254

const BAD_DATA_FILL = 255

const BAD_TFORM = 261

const BAD_TFORM_DTYPE = 262

const BAD_TDIM = 263

const BAD_HEAP_PTR = 264

const BAD_HDU_NUM = 301

const BAD_COL_NUM = 302

const NEG_FILE_POS = 304

const NEG_BYTES = 306

const BAD_ROW_NUM = 307

const BAD_ELEM_NUM = 308

const NOT_ASCII_COL = 309

const NOT_LOGICAL_COL = 310

const BAD_ATABLE_FORMAT = 311

const BAD_BTABLE_FORMAT = 312

const NO_NULL = 314

const NOT_VARI_LEN = 317

const BAD_DIMEN = 320

const BAD_PIX_NUM = 321

const ZERO_SCALE = 322

const NEG_AXIS = 323

const NOT_GROUP_TABLE = 340

const HDU_ALREADY_MEMBER = 341

const MEMBER_NOT_FOUND = 342

const GROUP_NOT_FOUND = 343

const BAD_GROUP_ID = 344

const TOO_MANY_HDUS_TRACKED = 345

const HDU_ALREADY_TRACKED = 346

const BAD_OPTION = 347

const IDENTICAL_POINTERS = 348

const BAD_GROUP_ATTACH = 349

const BAD_GROUP_DETACH = 350

const BAD_I2C = 401

const BAD_F2C = 402

const BAD_INTKEY = 403

const BAD_LOGICALKEY = 404

const BAD_FLOATKEY = 405

const BAD_DOUBLEKEY = 406

const BAD_C2I = 407

const BAD_C2F = 408

const BAD_C2D = 409

const BAD_DATATYPE = 410

const BAD_DECIM = 411

const NUM_OVERFLOW = 412

const DATA_COMPRESSION_ERR = 413

const DATA_DECOMPRESSION_ERR = 414

const NO_COMPRESSED_TILE = 415

const BAD_DATE = 420

const PARSE_SYNTAX_ERR = 431

const PARSE_BAD_TYPE = 432

const PARSE_LRG_VECTOR = 433

const PARSE_NO_OUTPUT = 434

const PARSE_BAD_COL = 435

const PARSE_BAD_OUTPUT = 436

const ANGLE_TOO_BIG = 501

const BAD_WCS_VAL = 502

const WCS_ERROR = 503

const BAD_WCS_PROJ = 504

const NO_WCS_KEY = 505

const APPROX_WCS_KEY = 506

const NO_CLOSE_ERROR = 999

const NGP_ERRBASE = 360

const NGP_OK = 0

const NGP_NO_MEMORY = NGP_ERRBASE + 0

const NGP_READ_ERR = NGP_ERRBASE + 1

const NGP_NUL_PTR = NGP_ERRBASE + 2

const NGP_EMPTY_CURLINE = NGP_ERRBASE + 3

const NGP_UNREAD_QUEUE_FULL = NGP_ERRBASE + 4

const NGP_INC_NESTING = NGP_ERRBASE + 5

const NGP_ERR_FOPEN = NGP_ERRBASE + 6

const NGP_EOF = NGP_ERRBASE + 7

const NGP_BAD_ARG = NGP_ERRBASE + 8

const NGP_TOKEN_NOT_EXPECT = NGP_ERRBASE + 9

# Retrieve version number from the library.
function fits_get_version()
    version = round(Int, CFITSIO.fits_get_version(Ref{Cfloat}())*10_000.0)
    major = div(version, 10_000)
    minor, micro = divrem(version - 10_000*major, 100)
    return VersionNumber(major, minor, micro)
end

function __init__()
    # Version number according to header files used to generate this code.
    gen_version = VersionNumber(CFITSIO_MAJOR,CFITSIO_MINOR,CFITSIO_MICRO)

    # Version number of the library provided by the artifact.
    lib_version = fits_get_version()

    # Check for compatibility.
    (lib_version.major == gen_version.major && lib_version.minor ≥ gen_version.minor) || @warn """
`CFITSIO_jll` library has version $(lib_version) while the headers used to generate the code
of `CFITSIO.jl` have version $(gen_version). You should regenerate code in `CFITSIO.jl`
following instructions in `$(@__DIR__)/README.md`.
"""

end

end # module
