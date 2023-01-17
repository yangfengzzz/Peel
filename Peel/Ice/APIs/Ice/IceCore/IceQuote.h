///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Cool tip from Flipcode.
 *	\file		IceQuote.h
 *	\author		Alberto Garc�a-Baquero Vega
 *	\date		January, 11, 2001
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Include Guard
#ifndef ICEQUOTE_H
#define ICEQUOTE_H

#define _QUOTE(x) #x
#define QUOTE(x) _QUOTE(x)
#define __FILE__LINE__ __FILE__ "(" QUOTE(__LINE__) ") : "

#define __NOTE(x) message(x)
#define __FILE_LINE message(__FILE__LINE__)

#define __TODO(x)                                                 \
    message(__FILE__LINE__                                        \
            "\n"                                                  \
            " ------------------------------------------------\n" \
            "|  TODO :   " #x                                     \
            "\n"                                                  \
            " -------------------------------------------------\n")

#define __FIXME(x)                                                \
    message(__FILE__LINE__                                        \
            "\n"                                                  \
            " ------------------------------------------------\n" \
            "|  FIXME :  " #x                                     \
            "\n"                                                  \
            " -------------------------------------------------\n")

#define __todo(x) message(__FILE__LINE__ " TODO :   " #x "\n")
#define __fixme(x) message(__FILE__LINE__ " FIXME:   " #x "\n")

#endif  // ICEQUOTE_H
