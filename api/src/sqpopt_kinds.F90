    module sqpopt_kinds

    use, intrinsic :: iso_fortran_env, only: real32, real64, real128

    implicit none

    private

#ifdef REAL32
    integer,parameter,public :: sqpopt_module_wp = real32   !! Real working precision [4 bytes]
#elif REAL64
    integer,parameter,public :: sqpopt_module_wp = real64   !! Real working precision [8 bytes]
#elif REAL128
    integer,parameter,public :: sqpopt_module_wp = real128  !! Real working precision [16 bytes]
#else
    integer,parameter,public :: sqpopt_module_wp = real64   !! Real working precision if not specified [8 bytes]
#endif

    end module sqpopt_kinds
