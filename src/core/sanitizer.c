#include "sanitizer.h"

#if __has_feature(address_sanitizer) || defined(__SANITIZE_ADDRESS__)

// NOLINTBEGIN(bugprone-reserved-identifier)

const char *__lsan_default_options(void)
{
    return "print_suppressions=0";
}

const char *__lsan_default_suppressions(void)
{
    return "leak:libmpi.so\n"
           "leak:libpmix.so\n";
}

// NOLINTEND(bugprone-reserved-identifier)

#endif
