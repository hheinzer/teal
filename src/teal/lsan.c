const char *__lsan_default_options(void)  // NOLINT(bugprone-reserved-identifier)
{
    return "print_suppressions=0";
}

const char *__lsan_default_suppressions(void)  // NOLINT(bugprone-reserved-identifier)
{
    return "leak:libmpi.so\n"
           "leak:libpmix.so\n";
}
