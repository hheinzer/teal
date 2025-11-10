const char *__asan_default_options(void)  // NOLINT(bugprone-reserved-identifier)
{
    return "abort_on_error=1:"
           "disable_coredump=0:"
           "unmap_shadow_on_exit=1";
}
