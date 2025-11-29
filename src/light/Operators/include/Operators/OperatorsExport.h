#pragma once

#if defined(_WIN32) || defined(__CYGWIN__)

#ifdef OPERATORS_BUILD
#define OPERATORS_API __declspec(dllexport)
#else
#define OPERATORS_API __declspec(dllimport)
#endif

#else

#define OPERATORS_API __attribute__((visibility("default")))

#endif