#pragma once

#if defined(_WIN32) || defined(__CYGWIN__)

#ifdef MEASUREMENT_BUILD
#define MEASUREMENTS_API __declspec(dllexport)
#else
#define MEASUREMENTS_API __declspec(dllimport)
#endif

#else

#define MEASUREMENTS_API __attribute__((visibility("default")))

#endif