HTSLIB_static_LDFLAGS = -Wl,-O2 -Wl,--sort-common -Wl,--as-needed -Wl,-z,relro -Wl,-z,now -Wl,--disable-new-dtags -Wl,--gc-sections -Wl,-rpath,/eep/softwares/miniconda/lib -Wl,-rpath-link,/eep/softwares/miniconda/lib -L/eep/softwares/miniconda/lib
HTSLIB_static_LIBS = -lpthread -lz -lm -lbz2 -llzma -ldeflate
