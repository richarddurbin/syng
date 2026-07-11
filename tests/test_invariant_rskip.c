#include <check.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <setjmp.h>

/* Import the error handling function from rskip.c */
extern void rsError(int type, const char *funcName, const char *format, ...);

static sigjmp_buf jump_buffer;
static volatile sig_atomic_t got_signal = 0;

static void signal_handler(int sig) {
    got_signal = 1;
    siglongjmp(jump_buffer, 1);
}

START_TEST(test_buffer_overflow_protection)
{
    /* Invariant: Buffer reads/writes never exceed declared length */
    char overflow_2x[512];
    char overflow_10x[2560];
    
    memset(overflow_2x, 'A', sizeof(overflow_2x) - 1);
    overflow_2x[sizeof(overflow_2x) - 1] = '\0';
    
    memset(overflow_10x, 'B', sizeof(overflow_10x) - 1);
    overflow_10x[sizeof(overflow_10x) - 1] = '\0';
    
    const char *payloads[] = {
        overflow_10x,                    /* 10x buffer size - exploit case */
        overflow_2x,                     /* 2x buffer size - boundary case */
        "validFunc"                      /* Valid short input */
    };
    int num_payloads = sizeof(payloads) / sizeof(payloads[0]);

    struct sigaction sa, old_sa;
    sa.sa_handler = signal_handler;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = 0;
    sigaction(SIGSEGV, &sa, &old_sa);
    sigaction(SIGBUS, &sa, &old_sa);

    for (int i = 0; i < num_payloads; i++) {
        got_signal = 0;
        if (sigsetjmp(jump_buffer, 1) == 0) {
            /* Call with oversized funcName - should not crash or corrupt stack */
            rsError(1, payloads[i], "test message %s", "arg");
        }
        ck_assert_msg(got_signal == 0, 
            "Buffer overflow detected with payload %d (len=%zu)", 
            i, strlen(payloads[i]));
    }

    sigaction(SIGSEGV, &old_sa, NULL);
    sigaction(SIGBUS, &old_sa, NULL);
}
END_TEST

Suite *security_suite(void)
{
    Suite *s;
    TCase *tc_core;

    s = suite_create("Security");
    tc_core = tcase_create("Core");

    tcase_add_test(tc_core, test_buffer_overflow_protection);
    suite_add_tcase(s, tc_core);

    return s;
}

int main(void)
{
    int number_failed;
    Suite *s;
    SRunner *sr;

    s = security_suite();
    sr = srunner_create(s);

    srunner_run_all(sr, CK_NORMAL);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);

    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}