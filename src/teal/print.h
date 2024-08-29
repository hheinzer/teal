#pragma once

[[gnu::format(printf, 2, 3)]] void print_key(const char *key, const char *format, ...);

void print_size(double size_bytes);
