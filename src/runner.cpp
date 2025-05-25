#include <gmp.h>

int main() {
    mpz_t n;
    int flag;

    mpz_init(n);
    mpz_set_ui(n, 100);
    mpz_out_str(stdout, 10, n);

    return 0;
}
