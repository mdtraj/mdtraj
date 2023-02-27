#include "center.h"
#include "msvccompat.h"
#ifdef __NO_INTRINSICS
#include "center_generic.h"
#elif defined(__ARM_NEON)
#include "center_arm.h"
#elif defined(__SSE2__) 
#include "center_sse.h"
#else
#include "center_generic.h"
#endif