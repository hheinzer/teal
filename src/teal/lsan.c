const char *__lsan_default_options(void)
{
    return "print_suppressions=0";
}

const char *__lsan_default_suppressions(void)
{
    return "leak:libpmix.so\n"
           "leak:libmpi.so\n";
}
