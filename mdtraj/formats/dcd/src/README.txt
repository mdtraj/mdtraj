dcdplugin.c
=========

VMD DCD Plugin with a few changes

- Update to latest version of dcdplugin.c from VMD
- Removed static keyword from several functions
- Added "dcd_nsets" and "dcd_rewind" from the original mdtraj copy of dcdplugin.c
- nsets is returned (via function parameter) from open_dcd_read
- Omitted VMD "initialization stuff"
- void * changed to dcdhandle *
- some parameters marked const in open_dcd_write
- with_unitcell is passed as an argument to open_dcd_write, rather than depending on
  an environment variable
