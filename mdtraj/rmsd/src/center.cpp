#include "center.h"
#include "msvccompat.h"
#ifdef __ARM_NEON
#include "center_arm.h"
#else
#include "center_sse.h"
#endif