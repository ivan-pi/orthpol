set(SRC
    r1mach.f90
    d1mach.f90
    scheb.f
    schri.f
    sgauss.f
    sgchri.f
    skern.f
    sknum.f
    slancz.f
    slob.f
    smccheb.f
    smcdis.f
    sqgp.f
    sradau.f
    srecur.f
    ssti.f
    dcheb.f
    dchri.f
    dgauss.f
    dgchri.f
    dkern.f
    dknum.f
    dlancz.f
    dlob.f
    dmcheb.f
    dmcdis.f
    dqgp.f
    dradau.f
    drecur.f
    dsti.f
    nu0jac.f
    nu0lag.f
    nu0her.f
    orthpol.f90
    )

add_library(orthpol ${SRC})
target_compile_options(orthpol PRIVATE -Wno-implicit-interface -Wno-compare-reals)
add_library(p::orthpol ALIAS orthpol)