// Common definitions

#ifndef _EYE_DEFS_
#define _EYE_DEFS_

#include "defs.h"
#include "mdefs.h"

#include "cstrng.h"

/*#build_stop*/

namespace eye {

using namespace lwml;

// local exceptions

DEF_EX_CLASS(ex_base, ex_epl);

DEF_EX_TYPE(ex_epl, ex_not_find, "incorrect called function follow, but we didn't find linb yet");

//dfsdf

// Интерфейс для логгера событий, возникающих при обработке

}; // namespace eye

#endif // _EYE_DEFS_
