#pragma once

#if defined(_WIN32) || defined(__CYGWIN__)

#ifdef CIRCUITS_BUILD
#define CIRCUITS_API __declspec(dllexport)
#else
#define CIRCUITS_API __declspec(dllimport)
#endif

#else

#define CIRCUITS_API __attribute__((visibility("default")))

#endif