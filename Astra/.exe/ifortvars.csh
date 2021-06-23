if !($?PATH) then
    setenv PATH /opt/intel/fc/9.1.039/bin
else
    set TEST = `echo $PATH | grep /opt/intel/fc/9.1.039/bin`
    if ($TEST == '') setenv PATH /opt/intel/fc/9.1.039/bin:$PATH
    unset TEST
endif

if !($?LD_LIBRARY_PATH) then
    setenv LD_LIBRARY_PATH /opt/intel/fc/9.1.039/lib
else
    set TEST = `echo $LD_LIBRARY_PATH | grep /opt/intel/fc/9.1.039/lib`
    if ($TEST == '') setenv LD_LIBRARY_PATH /opt/intel/fc/9.1.039/lib:$LD_LIBRARY_PATH
    unset TEST
endif

if !($?MANPATH) then
    setenv MANPATH /opt/intel/fc/9.1.039/man:`man -w`
else
    set TEST = `echo $MANPATH | grep /opt/intel/fc/9.1.039/man`
    if ($TEST == '') setenv MANPATH /opt/intel/fc/9.1.039/man:$MANPATH
    unset TEST
endif

if !($?INTEL_LICENSE_FILE) then
    setenv INTEL_LICENSE_FILE /opt/intel/fc/9.1.039/licenses
else
    set TEST = `echo $INTEL_LICENSE_FILE | grep /opt/intel/fc/9.1.039/licenses`
    if ($TEST == '') setenv INTEL_LICENSE_FILE /opt/intel/fc/9.1.039/licenses:$INTEL_LICENSE_FILE
    unset TEST
endif


