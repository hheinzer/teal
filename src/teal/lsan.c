// NOLINTBEGIN(bugprone-reserved-identifier)
__attribute((unused)) const char *__lsan_default_options(void)
{
    return "print_suppressions=0";
}
__attribute((unused)) const char *__lsan_default_suppressions(void)
{
    return "leak:libmpi.so\n"
           "leak:libpmix.so\n";
}
// NOLINTEND(bugprone-reserved-identifier)
