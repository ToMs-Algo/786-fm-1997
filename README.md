# Algorithm 786: multiple-precision complex arithmetic and functions

> David M. Smith. 1998.  
> Algorithm 786: multiple-precision complex arithmetic and functions.  
> ACM Trans. Math. Softw. 24, 4 (Dec. 1998), 359–367.  
> https://doi.org/10.1145/293686.293687

## FM Package

```for
!     FM 1.1                  David M. Smith               5-19-97

!  The FM routines in this package perform floating-point
!  multiple-precision arithmetic, and the IM routines perform
!  integer multiple-precision arithmetic.
```

### LIST OF ROUTINES

```for

!  These are the FM routines that are designed to be called by
!  the user.  All are subroutines except logical function FMCOMP.
!  MA, MB, MC refer to FM format numbers.

!  In each case it is permissible to use the same array more than
!  once in the calling sequence.  The statement MA = MA*MA can
!  be written CALL FMMPY(MA,MA,MA).

!  For each of these routines there is also a version available for
!  which the argument list is the same but all FM numbers are in packed
!  format.  The routines using packed numbers have the same names except
!  'FM' is replaced by 'FP' at the start of each name.


!  FMABS(MA,MB)         MB = ABS(MA)

!  FMACOS(MA,MB)        MB = ACOS(MA)

!  FMADD(MA,MB,MC)      MC = MA + MB

!  FMADDI(MA,IVAL)      MA = MA + IVAL   Increment an FM number by a one
!                                        word integer.  Note this call
!                                        does not have an "MB" result
!                                        like FMDIVI and FMMPYI.

!  FMASIN(MA,MB)        MB = ASIN(MA)

!  FMATAN(MA,MB)        MB = ATAN(MA)

!  FMATN2(MA,MB,MC)     MC = ATAN2(MA,MB)

!  FMBIG(MA)            MA = Biggest FM number less than overflow.

!  FMCHSH(MA,MB,MC)     MB = COSH(MA),  MC = SINH(MA).  Faster than
!                            making two separate calls.

!  FMCOMP(MA,LREL,MB)        Logical comparison of MA and MB.
!                            LREL is a CHARACTER*2 value identifying
!                            which comparison is made.
!                            Example:  IF (FMCOMP(MA,'GE',MB)) ...

!  FMCONS                    Set several saved constants that depend
!                            on MBASE, the base being used.  FMCONS
!                            should be called immediately after
!                            changing MBASE.

!  FMCOS(MA,MB)         MB = COS(MA)

!  FMCOSH(MA,MB)        MB = COSH(MA)

!  FMCSSN(MA,MB,MC)     MB = COS(MA),  MC = SIN(MA).  Faster than
!                            making two separate calls.

!  FMDIG(NSTACK,KST)         Find a set of precisions to use during
!                            Newton iteration for finding a simple
!                            root starting with about double
!                            precision accuracy.

!  FMDIM(MA,MB,MC)      MC = DIM(MA,MB)

!  FMDIV(MA,MB,MC)      MC = MA/MB

!  FMDIVI(MA,IVAL,MB)   MB = MA/IVAL   IVAL is a one word integer.

!  FMDP2M(X,MA)         MA = X    Convert from double precision to FM.

!  FMDPM(X,MA)          MA = X    Convert from double precision to FM.
!                                 Much faster than FMDP2M, but MA agrees
!                                 with X only to D.P. accuracy.  See
!                                 the comments in the two routines.

!  FMEQ(MA,MB)          MB = MA   Both have precision NDIG.
!                                 This is the version to use for
!                                 standard  B = A  statements.

!  FMEQU(MA,MB,NA,NB)   MB = MA   Version for changing precision.
!                                 MA has NA digits (i.e., MA was
!                                 computed using NDIG = NA), and MB
!                                 will be defined having NB digits.
!                                 MB is zero-padded if NB.GT.NA
!                                 MB is rounded if NB.LT.NA

!  FMEXP(MA,MB)         MB = EXP(MA)

!  FMFORM(FORM,MA,STRING)    MA is converted to a character string
!                               using format FORM and returned in
!                               STRING.  FORM can represent I, F,
!                               E, or 1PE formats.  Example:
!                               CALL FMFORM('F60.40',MA,STRING)

!  FMFPRT(FORM,MA)           Print MA on unit KW using FORM format.

!  FMI2M(IVAL,MA)       MA = IVAL   Convert from one word integer
!                                   to FM.

!  FMINP(LINE,MA,LA,LB) MA = LINE   Input conversion.
!                                   Convert LINE(LA) through LINE(LB)
!                                   from characters to FM.

!  FMINT(MA,MB)         MB = INT(MA)    Integer part of MA.

!  FMIPWR(MA,IVAL,MB)   MB = MA**IVAL   Raise an FM number to a one
!                                       word integer power.

!  FMLG10(MA,MB)        MB = LOG10(MA)

!  FMLN(MA,MB)          MB = LOG(MA)

!  FMLNI(IVAL,MA)       MA = LOG(IVAL)   Natural log of a one word
!                                        integer.

!  FMM2DP(MA,X)         X  = MA     Convert from FM to double precision.

!  FMM2I(MA,IVAL)       IVAL = MA   Convert from FM to integer.

!  FMM2SP(MA,X)         X  = MA     Convert from FM to single precision.

!  FMMAX(MA,MB,MC)      MC = MAX(MA,MB)

!  FMMIN(MA,MB,MC)      MC = MIN(MA,MB)

!  FMMOD(MA,MB,MC)      MC = MA mod MB

!  FMMPY(MA,MB,MC)      MC = MA*MB

!  FMMPYI(MA,IVAL,MB)   MB = MA*IVAL    Multiply by a one word integer.

!  FMNINT(MA,MB)        MB = NINT(MA)   Nearest FM integer.

!  FMOUT(MA,LINE,LB)    LINE = MA   Convert from FM to character.
!                                   LINE is a character array of
!                                   length LB.

!  FMPI(MA)             MA = pi

!  FMPRNT(MA)                Print MA on unit KW using current format.

!  FMPWR(MA,MB,MC)      MC = MA**MB

!  FMREAD(KREAD,MA)     MA   is returned after reading one (possibly
!                            multi-line) FM number on unit KREAD.  This
!                            routine reads numbers written by FMWRIT.

!  FMRPWR(MA,K,J,MB)    MB = MA**(K/J)  Rational power.  Faster than
!                            FMPWR for functions like the cube root.

!  FMSET(NPREC)              Set default values and machine-dependent
!                            variables to give at least NPREC base 10
!                            digits plus three base 10 guard digits.
!                            Must be called to initialize FM package.

!  FMSIGN(MA,MB,MC)     MC = SIGN(MA,MB)   Sign transfer.

!  FMSIN(MA,MB)         MB = SIN(MA)

!  FMSINH(MA,MB)        MB = SINH(MA)

!  FMSP2M(X,MA)         MA = X   Convert from single precision to FM.

!  FMSQR(MA,MB)         MB = MA*MA   Faster than FMMPY.

!  FMSQRT(MA,MB)        MB = SQRT(MA)

!  FMST2M(STRING,MA)    MA = STRING
!                            Convert from character string to FM.
!                            Often more convenient than FMINP, which
!                            converts an array of CHARACTER*1 values.
!                            Example:   CALL FMST2M('123.4',MA).

!  FMSUB(MA,MB,MC)      MC = MA - MB

!  FMTAN(MA,MB)         MB = TAN(MA)

!  FMTANH(MA,MB)        MB = TANH(MA)

!  FMULP(MA,MB)         MB = One Unit in the Last Place of MA.

!  FMWRIT(KWRITE,MA)         Write MA on unit KWRITE.
!                            Multi-line numbers will have '&' as the
!                            last nonblank character on all but the last
!                            line.  These numbers can then be read
!                            easily using FMREAD.


!  These are the integer routines that are designed to be called by
!  the user.  All are subroutines except logical function IMCOMP.
!  MA, MB, MC refer to IM format numbers.  In each case the version
!  of the routine to handle packed IM numbers has the same name,
!  with 'IM' replaced by 'IP'.

!  IMABS(MA,MB)         MB = ABS(MA)

!  IMADD(MA,MB,MC)      MC = MA + MB

!  IMBIG(MA)            MA = Biggest IM number less than overflow.

!  IMCOMP(MA,LREL,MB)        Logical comparison of MA and MB.
!                            LREL is a CHARACTER*2 value identifying
!                            which comparison is made.
!                            Example:  IF (IMCOMP(MA,'GE',MB)) ...

!  IMDIM(MA,MB,MC)      MC = DIM(MA,MB)

!  IMDIV(MA,MB,MC)      MC = int(MA/MB)
!                            Use IMDIVR if the remainder is also needed.

!  IMDIVI(MA,IVAL,MB)   MB = int(MA/IVAL)
!                            IVAL is a one word integer.  Use IMDVIR
!                            to get the remainder also.

!  IMDIVR(MA,MB,MC,MD)  MC = int(MA/MB),   MD = MA mod MB
!                            When both the quotient and remainder are
!                            needed, this routine is twice as fast as
!                            calling both IMDIV and IMMOD.

!  IMDVIR(MA,IVAL,MB,IREM)   MB = int(MA/IVAL),   IREM = MA mod IVAL
!                            IVAL and IREM are one word integers.

!  IMEQ(MA,MB)          MB = MA

!  IMFM2I(MAFM,MB)      MB = MAFM  Convert from real (FM) format
!                                  to integer (IM) format.

!  IMFORM(FORM,MA,STRING)    MA is converted to a character string
!                               using format FORM and returned in
!                               STRING.  FORM can represent I, F,
!                               E, or 1PE formats.  Example:
!                               CALL IMFORM('I70',MA,STRING)

!  IMFPRT(FORM,MA)           Print MA on unit KW using FORM format.

!  IMGCD(MA,MB,MC)      MC = greatest common divisor of MA and MB.

!  IMI2FM(MA,MBFM)    MBFM = MA  Convert from integer (IM) format
!                                to real (FM) format.

!  IMI2M(IVAL,MA)       MA = IVAL   Convert from one word integer
!                                   to IM.

!  IMINP(LINE,MA,LA,LB) MA = LINE   Input conversion.
!                                   Convert LINE(LA) through LINE(LB)
!                                   from characters to IM.

!  IMM2DP(MA,X)         X  = MA     Convert from IM to double precision.

!  IMM2I(MA,IVAL)       IVAL = MA   Convert from IM to one word integer.

!  IMMAX(MA,MB,MC)      MC = MAX(MA,MB)

!  IMMIN(MA,MB,MC)      MC = MIN(MA,MB)

!  IMMOD(MA,MB,MC)      MC = MA mod MB

!  IMMPY(MA,MB,MC)      MC = MA*MB

!  IMMPYI(MA,IVAL,MB)   MB = MA*IVAL    Multiply by a one word integer.

!  IMMPYM(MA,MB,MC,MD)  MD = MA*MB mod MC
!                            Slightly faster than calling IMMPY and
!                            IMMOD separately, and it works for cases
!                            where IMMPY would return OVERFLOW.

!  IMOUT(MA,LINE,LB)    LINE = MA   Convert from IM to character.
!                                   LINE is a character array of
!                                   length LB.

!  IMPMOD(MA,MB,MC,MD)       MD = MA**MB mod MC

!  IMPRNT(MA)                Print MA on unit KW.

!  IMPWR(MA,MB,MC)      MC = MA**MB

!  IMREAD(KREAD,MA)     MA   is returned after reading one (possibly
!                            multi-line) IM number on unit KREAD.  This
!                            routine reads numbers written by IMWRIT.

!  IMSIGN(MA,MB,MC)     MC = SIGN(MA,MB)   Sign transfer.

!  IMSQR(MA,MB)         MB = MA*MA   Faster than IMMPY.

!  IMST2M(STRING,MA)    MA = STRING
!                            Convert from character string to IM.
!                            Often more convenient than IMINP, which
!                            converts an array of CHARACTER*1 values.
!                            Example:   CALL IMST2M('12345678901',MA).

!  IMSUB(MA,MB,MC)      MC = MA - MB

!  IMWRIT(KWRITE,MA)         Write MA on unit KWRITE.
!                            Multi-line numbers will have '&' as the
!                            last nonblank character on all but the last
!                            line.  These numbers can then be read
!                            easily using IMREAD.

!  Many of the IM routines call FM routines, but none of the FM
!  routines call IM routines, so the IM routines can be omitted
!  if none are called explicitly from a program.
```
