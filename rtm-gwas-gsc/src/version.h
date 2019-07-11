#ifndef VERSION_H
#define VERSION_H


#ifndef RTM_GWAS_VERSION
#define RTM_GWAS_VERSION    unknown
#endif


#define STRINGIZING1(x)             #x
#define STRINGIZING2(x)             STRINGIZING1(x)
#define RTM_GWAS_VERSION_STRING     STRINGIZING2(RTM_GWAS_VERSION)


#endif // VERSION_H
